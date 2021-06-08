#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A small GUI program to interface a custom DDM setup.

@author: Frédéric Dux, biosoft intern@IPC with Jerome Crassous
"""

import  numpy              as      np
from    scipy.signal       import  tukey
from    skvideo.io         import  vread, vreader, FFmpegReader
from    os.path            import  exists, basename, dirname, join
from    os                 import  makedirs, remove
import  cv2
from    joblib             import  Parallel, delayed
from    multiprocessing    import  cpu_count
from    subprocess         import  call

from    utilities          import  ddm_matrices, RadialAverager

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
        vid_reader    = FFmpegReader(path) 
        Nframes       = vid_reader.getShape()[0]
        self.Nperpart = Nframes // self.Npartitions
        del vid_reader
        vid_reader    = vreader(path, as_grey=True)
        for i in range(self.Npartitions):
            stack     = FFTStack(self.freq, self.pixelsize, self.maxCouples, 
                                 self.ptPerDecade, t0=i, alone=False)
            stack.loadVideoFromGenerator(vid_reader, path, self.Nperpart)
            
            self.partitions[i] = stack
        del vid_reader
        self.loaded = True
        # loaded the videos, useless to attempt doing it in parallel though.
        # (io speed limited, serial nature of the generator)
        
        
    def fftAllStacks(self):
        for start in self.partitions:
            self.partitions[start].fftVideo()
        self.fftdone = True
        # also not doing that in parallel, as the fft routines are already
        # heavily optimized and will use 100% of the cpu
        
        
    def __ddmOneStack(self, i):
        return self.partitions[i].stackToDDM()
    
    
    def ddmAllStacks(self):
        listofDics =  Parallel(n_jobs=cpu_count()) ( delayed(self.__ddmOneStack) (i) for i in range(self.Npartitions) )
        del self.partitions
        for i in range(self.Npartitions):
            self.ddmmatrices[i*self.Nperpart] = listofDics[i]
        self.completed = True
    
    def getStatus(self):
        return self.loaded, self.fftdone, self.completed


class FFTStack():
    
    def __init__(self, freq, pixelsize, maxCouples, ptPerDecade, Nangle=1, 
                 t0=0, alone=True, debug=False, windowing=False):
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
        
    def __len__(self):
        return self.Nbimages
    
    def __getitem__(self, t):
        """returns the image at time t"""
        if t<0: t= len(self)+t
        if t > len(self): t = t - self.t0
        return self.data[t,:,:]
    
    def loadVideoFromGenerator(self, vreader_generator, filename, nframes):
        self.Nbimages = nframes
        self.filename = filename
        firstFrame = vreader_generator.send(None)[0,:,:,0]
        
        x, y       = firstFrame.shape[:2]
        self.shape = (x,y)
        self.data  = np.zeros((nframes, x, y), dtype=np.complex64)
        self.data[0, :, :] = firstFrame
        for i in range(1, nframes):
            self.data[i, :, :] = vreader_generator.send(None)[0,:,:,0]
    
    def loadVideo(self, filename, t0=0):
        self.filename    = filename
        self.t0          = t0
        self.data        = vread(filename, as_grey=1)[:,:,:,0].astype(np.complex64)
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
            self.window = tukey_twoD(self.shape[0], self.windowing)
        for t in range(self.Nbimages):
            #self.data[t,:,:] = dctn(self.data[t,:,:])
            if self.windowing:
                self.data[t,:,:] = np.fft.fft2(self.window*self.data[t,:,:])
            else:
                self.data[t,:,:] = np.fft.fft2(self.data[t,:,:])
            self.progress   += 1 
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
        DDMs= [np.zeros((len(idts), self.shape[0]//2)) for _ in range(self.Nangle)]

        for i, idt in enumerate(idts):
            curves = ra(self.timeAverage(idt))
            for j, curve in enumerate(curves):
                DDMs[j][i] = curve[:self.shape[0]//2]
            if self.maxCouples > 0 :
                self.progress += self.maxCouples
            else:
                self.progress += self.Nbimages
        return DDMs
        
    def timeAverage(self, dt):
        """Does at most maxNCouples spectreDiff on regularly spaced couples of images. 
        Separation within couple is dt."""
        #Spread initial times over the available range
        initialTimes  = np.arange(0, len(self)-dt, self.increment)
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
        
        expname    = self.filename.replace('.avi', '')
        expdir     = dirname(expname)
        savedir    = join(expdir, ddm_matrices)
        expname    = basename(expname)
        makedirs(savedir, exist_ok=True)
        if self.alone:
            extra = ''
        else:
            extra = '__i='+str(self.t0*self.Nbimages)+'__'

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
                
        del DDMs, qs, idts, dts, self.data
        self.completed = True
        
        
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
    nbdecades = np.log10(L)
    return np.unique(np.logspace(
        start=0, stop=nbdecades, 
        num= int(nbdecades) * pointsPerDecade, 
        base=10, endpoint=False
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
    gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
    return gray

def concatenateVideos(pathsToVideos, pathToConcatenated, as_grey=True):
    workingdir = dirname(pathsToVideos[0])
    listofvidsforffmpeg = join(workingdir, "my_list_of_videos_to_stack_3215.txt")
    with open(listofvidsforffmpeg, 'w') as f:
        for path in pathsToVideos:
            f.writelines(f"file {path}\n")
    call(['ffmpeg', '-safe', '0', '-f', 'concat', '-i', listofvidsforffmpeg, \
          '-c', 'copy', pathToConcatenated ])
    remove(listofvidsforffmpeg)
    
        




































