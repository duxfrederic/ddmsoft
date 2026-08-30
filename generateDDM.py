#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A small GUI program to interface a custom DDM setup.

@author: Frédéric Dux, biosoft intern@IPC with Jerome Crassous
"""

import  numpy                as      np
from    scipy.signal.windows import  tukey
from    os.path              import  exists, basename, dirname, join, splitext
from    os                   import  makedirs, remove
import  cv2
from    joblib               import  Parallel, delayed
from    multiprocessing      import  cpu_count

from    utilities            import  ddm_matrices, RadialAverager


def _as_gray(frame):
    """Return a video frame as a two-dimensional float-compatible array."""
    frame = np.asarray(frame)
    if frame.ndim == 3:
        if frame.shape[-1] == 1:
            frame = frame[..., 0]
        else:
            frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
    if frame.ndim != 2:
        raise ValueError(f"Expected a grayscale or colour image, got shape {frame.shape}")
    return frame


def _video_frame_count(filename):
    capture = cv2.VideoCapture(filename)
    if not capture.isOpened():
        capture.release()
        raise IOError(f"Could not open video: {filename}")
    count = int(capture.get(cv2.CAP_PROP_FRAME_COUNT))
    capture.release()
    if count <= 0:
        raise ValueError(f"Could not determine the frame count for: {filename}")
    return count


def _video_frames(filename):
    capture = cv2.VideoCapture(filename)
    if not capture.isOpened():
        capture.release()
        raise IOError(f"Could not open video: {filename}")
    try:
        while True:
            success, frame = capture.read()
            if not success:
                break
            yield _as_gray(frame)
    finally:
        capture.release()


def partition_frame_counts(nframes, npartitions):
    """Split a frame count into non-empty partitions without dropping frames."""
    nframes = int(nframes)
    npartitions = int(npartitions)
    if nframes < 1:
        raise ValueError("A video must contain at least one frame")
    if not 1 <= npartitions <= nframes:
        raise ValueError("The number of partitions must be between 1 and the frame count")
    base, remainder = divmod(nframes, npartitions)
    return [base + (index < remainder) for index in range(npartitions)]

def tukey_twoD(width, alpha):
    """2D tukey lowpass window with a circular support
    """
    base = np.zeros((width, width))
    tuk = tukey(width, alpha)
    tuk = tuk[int(len(tuk)/2)-1:]  # Second half of tukey window
    x = np.linspace(-width/2, width/2, width)
    y = np.linspace(-width/2, width/2, width)
    for x_index in range(0, width):
        for y_index in range(0, width):
            # Only plot tukey value with in circle of radius width
            if int(np.sqrt(x[x_index]**2 + y[y_index]**2)) <= width/2:
                base[x_index, y_index] = tuk[int(np.sqrt(x[x_index]**2
                     + y[y_index]**2))]
                    # Based on k**2 find tukey window value and place in matrix
    return base

class timeDependantDDM():
    
    def __init__(self, freq, pixelsize, maxCouples, ptPerDecade, Npartitions):
        self.freq        = freq
        self.pixelsize   = pixelsize
        self.maxCouples  = maxCouples
        self.ptPerDecade = ptPerDecade
        self.Npartitions = Npartitions
        self.completed   = False
        self.fftdone     = False
        self.loaded      = False
        self.partitions  = {}
        self.ddmmatrices = {}
        self.Nperpart    = 1
    
    def loadVideo(self, path):
        Nframes = _video_frame_count(path)
        counts = partition_frame_counts(Nframes, self.Npartitions)
        self.partitions.clear()
        frame_offset = 0
        frames = _video_frames(path)
        for count in counts:
            stack = FFTStack(self.freq, self.pixelsize, self.maxCouples,
                             self.ptPerDecade, t0=frame_offset, alone=False)
            stack.loadVideoFromGenerator(frames, path, count)
            self.partitions[frame_offset] = stack
            frame_offset += count
        self.Nperpart = counts[0]
        self.loaded = True
        # loaded the videos, useless to attempt doing it in parallel though.
        # (io speed limited, serial nature of the generator)
        
        
    def fftAllStacks(self):
        for stack in self.partitions.values():
            stack.fftVideo()
        self.fftdone = True
        # also not doing that in parallel, as the fft routines are already
        # heavily optimized and will use 100% of the cpu
        
        
    def __ddmOneStack(self, i):
        return self.partitions[i].stackToDDM()
    
    
    def ddmAllStacks(self):
        offsets = list(self.partitions)
        listofDics = Parallel(n_jobs=cpu_count())(
            delayed(self.__ddmOneStack)(offset) for offset in offsets
        )
        del self.partitions
        for offset, matrix in zip(offsets, listofDics):
            self.ddmmatrices[offset] = matrix
        self.completed = True
    
    def getStatus(self):
        return self.loaded, self.fftdone, self.completed


class FFTStack():
    
    def __init__(self, freq, pixelsize, maxCouples, ptPerDecade, Nangle=1, 
                 t0=0, alone=True, debug=False, windowing=False,
                 progress_callback=None):
        self.freq        = freq
        self.pixelsize   = pixelsize
        self.maxCouples  = maxCouples
        self.ptPerDecade = ptPerDecade
        self.t0          = t0
        self.progress    = 0
        self.completed   = False
        self.fftdone     = False
        self.averaged    = {}
        self.Nbimages    = 0
        self.alone       = alone
        self.Nangle      = Nangle
        # time average params:
        self.interval    = 1
        # average fft:
        self.averageFFT  = None
        self.debug       = debug
        self.windowing   = windowing
        self.progress_callback = progress_callback
        
    def __len__(self):
        return self.Nbimages
    
    def __getitem__(self, t):
        """returns the image at time t"""
        if t<0: t= len(self)+t
        if t > len(self): t = t - self.t0
        return self.data[t,:,:]
    
    def loadVideoFromGenerator(self, vreader_generator, filename, nframes):
        self.Nbimages = int(nframes)
        if self.Nbimages < 1:
            raise ValueError("A stack must contain at least one frame")
        self.filename = filename
        frames = iter(vreader_generator)
        try:
            firstFrame = _as_gray(next(frames))
        except StopIteration as error:
            raise ValueError("The video ended before the requested frames were read") from error
        x, y       = firstFrame.shape[:2]
        self.shape = (x,y)
        self.data  = np.zeros((nframes, x, y), dtype=np.complex64)
        self.data[0, :, :] = firstFrame
        for i in range(1, nframes):
            try:
                frame = _as_gray(next(frames))
            except StopIteration as error:
                raise ValueError("The video ended before the requested frames were read") from error
            if frame.shape != self.shape:
                raise ValueError("All video frames must have the same dimensions")
            self.data[i, :, :] = frame
    
    def loadVideo(self, filename, t0=0):
        self.filename    = filename
        self.t0          = t0
        frames = list(_video_frames(filename))
        if not frames:
            raise ValueError(f"No frames found in video: {filename}")
        shape = frames[0].shape
        if any(frame.shape != shape for frame in frames):
            raise ValueError("All video frames must have the same dimensions")
        self.data        = np.asarray(frames, dtype=np.complex64)
        self.Nbimages    = self.data.shape[0]
        # get the images shape while checking that the last image does exist
        self.shape       = self.data.shape[1:]
        
    def getTotalNumberOfOperations(self):
            if self.maxCouples > 0:
                mult = self.maxCouples
            else:
                mult = self.Nbimages
            return self.Nbimages + len(logSpaced(self.Nbimages, self.ptPerDecade))*mult
        
    def fftVideo(self):
        if self.windowing:
            if self.shape[0] == self.shape[1]:
                self.window = tukey_twoD(self.shape[0], self.windowing)
            else:
                self.window = np.outer(
                    tukey(self.shape[0], self.windowing),
                    tukey(self.shape[1], self.windowing),
                )
        for t in range(self.Nbimages):
            #self.data[t,:,:] = dctn(self.data[t,:,:])
            if self.windowing:
                self.data[t,:,:] = np.fft.fft2(self.window*self.data[t,:,:])
            else:
                self.data[t,:,:] = np.fft.fft2(self.data[t,:,:])
            self.progress   += 1 
            self._report_progress()
        self.fftdone = True
        
    def getAverageFFT(self):
        if not self.fftdone:
            print('do fft first')
            return 0
        if self.averageFFT is None:
            self.averageFFT =  np.mean(self.data, axis=0)
        return self.averageFFT        
    
    def ddm(self, idts, maxNCouples=1000):
        """Perform time averaged and radial averaged DDM for given time intervals.
        Returns the DDM matrix."""
        # time averaging parameters:
        N             = self.Nbimages
        if maxNCouples == 0:
            self.increment = 1
            print("Max number of increments. (max # of couples set to 0)")
        else:
            self.increment = N//maxNCouples
        if self.increment == 0:
            self.increment = 1
        elif self.increment > 1:
            print(f"The size of the stack is {N}, but only {maxNCouples} couples of frames\
 were allowed. (More statistics can be extracted from the stack by setting\
 a higher maximal number of couples. (slower)")
        ra   = RadialAverager(self.shape, self.Nangle)
        nq = min(self.shape) // 2
        DDMs= [np.zeros((len(idts), nq)) for _ in range(self.Nangle)]

        for i, idt in enumerate(idts):
            curves = ra(self.timeAverage(idt))
            for j, curve in enumerate(curves):
                DDMs[j][i] = curve[:nq]
            if self.maxCouples > 0 :
                self.progress += self.maxCouples
            else:
                self.progress += self.Nbimages
            self._report_progress()
        return DDMs
        
    def timeAverage(self, dt):
        """Does at most maxNCouples spectreDiff on regularly spaced couples of images. 
        Separation within couple is dt."""
        #Spread initial times over the available range
        initialTimes  = np.arange(0, len(self)-dt, self.increment)
        if initialTimes.size == 0:
            raise ValueError(f"Lag time {dt} is not available for a {len(self)} frame stack")
        inverseleng   = 1./ initialTimes.size
        #perform the time average
        avgFFT        = np.zeros(self.shape)
        for t in initialTimes:
            # divide in the loop to avoid precision loss, slower but heh
            avgFFT   += spectrumDiff(self[int(t)], self[int(t+dt)]) * inverseleng
        if self.debug:
            np.save(join(dirname(self.filename), f"{basename(self.filename)}_time_averaged_diff_fft_tau={dt*1000:.01f}ms"), avgFFT)
        return avgFFT 


    def stackToDDM(self):
        idts       = logSpaced(self.Nbimages, self.ptPerDecade)
        dts        = idts/float(self.freq)
    
        DDMs       = self.ddm(idts, self.maxCouples)
        
        qs         = np.pi/(DDMs[0].shape[-1]*self.pixelsize) * np.arange(DDMs[0].shape[-1])
        # move the values to the center of each bin (see utilities.RadialAverager)
        qs         = qs + 0.5 * (qs[1]-qs[0])
        
        expname    = splitext(self.filename)[0]
        expdir     = dirname(expname)
        savedir    = join(expdir, ddm_matrices)
        expname    = basename(expname)
        makedirs(savedir, exist_ok=True)
        if self.alone:
            extra = ''
        else:
            extra = '__i='+str(self.t0)+'__'

        if len(DDMs)==1:
            np.save(join( savedir, expname + extra + '_QS'), qs)
            np.save(join( savedir, expname + extra + '_deltaTs'), dts)
            np.save(join( savedir, expname + extra + '_DDM_matrix'), DDMs[0])
        else:
            angles = np.arange(0, self.Nangle + 1) / self.Nangle * 180
            for i, DDM in enumerate(DDMs):
                angle_param = f"_{angles[i]:.01f}_"
                np.save(join( savedir, expname + extra + angle_param + '_QS'), qs)
                np.save(join( savedir, expname + extra + angle_param + '_deltaTs'), dts)
                np.save(join( savedir, expname + extra + angle_param + '_DDM_matrix'), DDMs[i])
                
        result = (DDMs, dts, qs)
        del DDMs, qs, idts, dts, self.data
        self.completed = True
        return result

    def _report_progress(self):
        if self.progress_callback is not None:
            self.progress_callback(self.progress, self.getTotalNumberOfOperations())
        
        
    def getProgress(self):
        return self.progress
    def isFFTdone(self):
        return self.fftdone
    def isCompleted(self):
        return self.completed
    
def spectrumDiff(imfourier0,imfourier1):
    """
        simply returns the power spectrum of the difference of two fourier images.
    """
    diff = imfourier1-imfourier0
    diff = np.real( diff*np.conj(diff) )
    return diff


def logSpaced(L, pointsPerDecade=20):
    """
    Generate an array of log spaced integers smaller than L.
    taken from https://github.com/MathieuLeocmach/colloids/blob/master/python/colloids/ddm.py
    """
    L = int(L)
    pointsPerDecade = int(pointsPerDecade)
    if L < 2:
        raise ValueError("At least two frames are required to calculate lag times")
    if pointsPerDecade < 1:
        raise ValueError("pointsPerDecade must be positive")
    nbdecades = np.log10(L - 1)
    npoints = max(2, int(np.ceil(np.log10(L) * pointsPerDecade)))
    return np.unique(np.logspace(
        start=0, stop=nbdecades,
        num=npoints,
        base=10, endpoint=True
    ).astype(int))
    
    
def readVideoFrame(filename, framenumber):
    """
        reads a gray image from the framenumberth image of the video contained at
        filename
    """
    if not exists(filename):
        raise IOError
    cap = cv2.VideoCapture(filename)
    cap.set(1, framenumber) # 2 is the CV_CAP_PROP_POS_FRAMES flag
    res, frame = cap.read()
    cap.release()
    if not res:
        raise IOError(f"Could not read frame {framenumber} from {filename}")
    gray = _as_gray(frame)
    return gray

def concatenateVideos(pathsToVideos, pathToConcatenated, as_grey=True):
    """Concatenate videos with OpenCV, without requiring an ffmpeg binary."""
    if not pathsToVideos:
        raise ValueError("At least one video is required")
    first = cv2.VideoCapture(pathsToVideos[0])
    if not first.isOpened():
        first.release()
        raise IOError(f"Could not open video: {pathsToVideos[0]}")
    width = int(first.get(cv2.CAP_PROP_FRAME_WIDTH))
    height = int(first.get(cv2.CAP_PROP_FRAME_HEIGHT))
    fps = first.get(cv2.CAP_PROP_FPS)
    first.release()
    if width < 1 or height < 1:
        raise ValueError("Could not determine the video dimensions")
    fps = fps if fps > 0 else 30.0
    writer = cv2.VideoWriter(
        pathToConcatenated,
        cv2.VideoWriter_fourcc(*"MJPG"),
        fps,
        (width, height),
    )
    if not writer.isOpened():
        writer.release()
        raise IOError(f"Could not create video: {pathToConcatenated}")
    try:
        for path in pathsToVideos:
            capture = cv2.VideoCapture(path)
            if not capture.isOpened():
                raise IOError(f"Could not open video: {path}")
            try:
                while True:
                    success, frame = capture.read()
                    if not success:
                        break
                    if frame.shape[1::-1] != (width, height):
                        raise ValueError("All videos must have the same dimensions")
                    writer.write(frame)
            finally:
                capture.release()
    except Exception:
        writer.release()
        if exists(pathToConcatenated):
            remove(pathToConcatenated)
        raise
    writer.release()
    
        

































