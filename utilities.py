#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A small GUI program to interface a custom DDM setup.

@author: Frédéric Dux, biosoft intern@IPC with Jerome Crassous
"""
import        re
import        numpy           as      np
from          os.path         import  exists, join, basename, isfile, normpath, \
                                      split, dirname
from          glob            import  glob
from          shutil          import  move


musthaves    = ['_DDM_matrix.npy', '_deltaTs.npy', '_QS.npy']
ddm_matrices = "ddm_matrices" # subdirectory in which computed matrices are stored

def water_viscosity(T):
    """
       https://en.wikipedia.org/wiki/Viscosity#Water
       not sourced, but it works
       (also: converted to the natural number basis)
    """
    return 2.414e-5 * np.exp(570.6 / (T - 140.))


def effectiveTemperature(thermistance_measured):
    """
        fit from the getTemperatureLA115.py
    """
    ### +2. : correction from the swelling curves, make the transition
    ### happen at 32°C (pnipam). Temperature to be re-calibrated with a better
    ### probed-radiator contact eventually.
    ### also, from uncertainties: offset 111 was +- 7, to be reduced.
    return 0.624 * thermistance_measured + 111 + 1.5

def extractCrudef(ddm, A=None, B=None):
    if (A is None) and (B is None):
        aplusb = (ddm[-1,:]+ddm[-2,:]+ddm[-3,:])/3
        B      = ddm[0, :]
        A      = aplusb-B
    f = 1 - (ddm-B)/(A[np.newaxis,:])
    return f

def loadDirectory(path):
    params = {}
    if exists(path):
        files    = glob(join(path, '*.txt'))
        videos   = glob(join(path, '*.avi'))
        files.sort()
        videos.sort()
        if len(videos)==0:
            # case no video found
            print(f'No video found in directory {path}')
            return {}
        if len(files)==0:
            # case no config file
            print(f'No config file found in directory {path}')
            for vid in videos:
                params[vid] =  {'pixelsize':'?', 'temperature':'?', 'framerate':'?'}
            return params
        elif len(files) == len(videos) and len(videos)>1:
            print(f'one config file per video found in directory {path}')
            for file, vid in zip(files, videos):
                if not file.replace('.txt','')==vid.replace('.avi',''):
                    # case wrong config file
                    return 0
                params[vid] = readParams(file)
            return params
        elif len(files) == 1:
            print(f'one common config file found for {len(videos)} videos in directory {path}')
            for vid in videos:
                params[vid] = readParams(files[0])
            return params
        else:
            print(len(files), len(videos))
    else:
        # case directory does not exist
        return {}
    
def loadAnalyzedVideos(path):
    path = join(path, ddm_matrices)
    if exists(path):
        files = glob(join(path, '*'+musthaves[0]))
        files.sort()
        exps  = [file.replace(musthaves[0], '') for file in files]
        computeds = {}
        
        for exp, file in zip(exps, files):
            ok = True
            for musthave in musthaves:
                if not exists(exp+musthave):
                    ok = False
            if ok:
                data = [np.load(exp+musthave) for musthave in musthaves]
                computeds[file] = data
        displays  = {basename(computed):computed for computed in computeds.keys()}
        return computeds, displays
    else:
        return {}, {}
    
    
def readParams(path, existingparams=None):
    if existingparams:
        convertedkey = path.replace(musthaves[0], '.avi').replace(ddm_matrices, '')
        convertedkey = normpath(convertedkey)
        if convertedkey in existingparams:
            return existingparams[convertedkey]
    
    generic_name = 'acquisition_parameters.txt'
    if isfile(path):
        
        # go back to the main directory:
        pathdir, taildir       = split(path)
        pathdirdir, taildirdir = split(pathdir)
        if ddm_matrices in taildir or ddm_matrices in taildirdir:
            path = normpath(path.replace(ddm_matrices, ''))
        # and proceed
        path = path.replace('.avi', '.txt')
        for musthave in musthaves:
            path = path.replace(musthave, '.txt')
        if not exists(path):
            path = path.replace(basename(path), generic_name)
    else:
        path = join(path, generic_name)
        
    if not exists(path):
        path += '.txt' # sometimes windows hides the extension of the file
    if not exists(path):
        try:
            path = glob(join(dirname(path), '*.txt'))[0]
        except:
            print("no appropriate config file found")
            pass # then we're screwed, allow the traceback from below to unroll
        
        
    with open(path, 'r') as f:
        text   = f.readlines()
        params = {}
        for element in text:
            key, value  = element.split(':')
            key         = key.strip()
            value       = value.replace('\n', '')
            value       = value.strip()
            params[key] = value
    return params


def saveMatrixCSV(path, ddm_dts_qs):
    base         = basename(path)
    musthavescsv = [must.replace('.npy', '.csv') for must in musthaves]
    for item, ext in zip(ddm_dts_qs, musthavescsv):
        savepath = path.replace(base, base+ext)
        np.savetxt(savepath, item, fmt='%.6e', delimiter='\t')
     

def saveAutocorrelationCSV(path, ddm_dts_qs, A, B, qsest):
    ddm, dts, qs = ddm_dts_qs
    qindexmin = np.argmin(np.abs(np.min(qsest)-qs))
    qindexmax = np.argmin(np.abs(np.max(qsest)-qs))
    frefined  = extractCrudef(ddm[:, qindexmin:qindexmax+1], A, B)
    np.savetxt(path + '_autocorrelationmatrix.csv', frefined, fmt='%.6e', delimiter='\t')
    np.savetxt(path + '_qs.csv', qs[qindexmin:qindexmax+1], fmt='%.6e', delimiter='\t')
    np.savetxt(path + '_dts.csv', dts, fmt='%.6e', delimiter='\t')


def saveFitTextFile(path, qs, A, B, params, paramnames, viscosity=None,\
                    temperature=None):
    
    if temperature and viscosity:
        D        = params[0]
        kb       = 1.381e-23
        RHeff    = kb*temperature/(6*np.pi*viscosity*(D))*1e9
    else:
        RHeff = None
    path = normpath(path)
    if not (path.endswith('.csv') or path.endswith('.txt')):
        path += '.txt'
    delimiter = '\t'
    header = ["q [m^-1]", "A", "B"] + paramnames
    if not RHeff is  None:
        header += ["R_H eff (nm)"]
    header = delimiter.join(header) + "\n"
    with open(path, 'w') as savefile:
        savefile.write(header)
        for i in range(len(qs)):
            qab       = [qs[i], A[i], B[i]]
            paramline = [paramset[i] for paramset in params] 
            towrite   = qab + paramline
            if not RHeff is None:
                towrite += [RHeff[i]]
            towrite   = [f"{p:.3e}" for p in towrite]
            savefile.write(delimiter.join(towrite)+"\n")
            
            

def saveSingleCONTINfit(path, sol, video, q):
    with open(path, 'w') as f:
        f.write(f"Matrix:\t{video}\n")
        f.write(f"Wavenumber [1/m]:\t{q:.03e}\n")
        f.write(f"alpha with smallest residual:\t{sol.chosen_alpha:.03e}\n")
        for i in range(len(sol.alphas)):
            alpha  = sol.alphas[i]
            res    = sol.alpha_residuals[i]
            dist   = sol.alpha_g[i]
            noise  = sol.alpha_noise[i]
            ampl   = sol.alpha_amplitude[i]
            ddmfit = sol.alpha_ddmfit[i] 
            f.write(3*"###############################"+"\n")
            f.write(f"alpha: {alpha:.03e} with residuals {res:.05e}\n")
            f.write(f"amplitude:\t{sol.amplitude:.03e}\n")
            f.write(f"noise:\t{sol.noise:.03e}\n\n")
            f.write("Distribution of decay rates\n")
            if sol.sizes is None:
                f.write("Gamma [m^2/s] \t intensity []\n")
                for gamma, p in zip(sol.gamma_range, dist):
                    f.write(f"{gamma:.03e}\t{p:.03e}\n")
            else:
                f.write("Gamma [m^2/s] \t size [nm] \t intensity []\n")
                for gamma, size, p in zip(sol.gamma_range, sol.sizes, dist):
                    f.write(f"{gamma:.03e}\t{1e9*size:.01e}\t{p:.03e}\n")
                    
            f.write("\ntau q^2 [s/m^2]\texp data\t CONTIN fit\n")
            intensity = ampl * (1-ddmfit) + noise
            for tau, exp, fit in zip(sol.tau, sol.ddmdata, intensity):
                f.write(f"{tau:.03e}\t{exp:.03e}\t{fit:.03e}\n")

def saveCONTINfit(path, CONTINsolutions, CONTINwindow, video, q):
    if bool(CONTINwindow.Element('continsaveall').Get()):
        for video in CONTINsolutions:
            for q in CONTINsolutions[video]:
                path_mod = path + f"_{basename(video.replace('.npy', ''))}_q={q:.03e}.csv"
                sol = CONTINsolutions[video][q]
                saveSingleCONTINfit(path_mod, sol, video, q)
    else:       
        sol = CONTINsolutions[video][q]
        saveSingleCONTINfit(path, sol, video, q)
            
class RadialAverager(object):
    """
    Radial average of a 2D array centred on (0,0), like the result of fft2d.
    
    adapted from https://github.com/MathieuLeocmach/colloids/blob/master/python/colloids/ddm.py
    
    """
    def __init__(self, shape, N=1):
        """
            A RadialAverager instance can process only arrays of a given shape, fixed at instanciation.
            N is the number of portions the q-wheel must be divided into.
        """
        assert len(shape) == 2
        self.shape       = shape
        self.N           = N
        # matrix of distances (we never do fftshift)
        self.dists       = np.sqrt(np.fft.fftfreq(shape[0])[:,None]**2 +  np.fft.fftfreq(shape[1])[None,:]**2)

        self.radbins     = np.arange(max(shape)//2+1)/float(max(shape))
        if N > 1:
            # matrix of arguments
            self.args    = np.arctan(np.fft.fftfreq(shape[1])[None,:] / np.fft.fftfreq(shape[0])[:,None]) + np.pi/2 
            # angular division of the wheel:
            self.argbins = np.arange(-0.5, N) / N * (np.pi) 
            self.args[self.args>(N-0.5)/N*np.pi] -= np.pi
            self.argbins[0] -= 1e-10
            
            # number of pixels at each distance
            self.hd          = []
            self.whos        = []
            for i in range(self.N):
                who = np.where( ( self.args > self.argbins[i]) * (self.args <= self.argbins[i+1]) * (self.dists <= 0.5) + np.isnan(self.args) ) 
                self.whos.append(who)
                self.hd.append( np.histogram(self.dists[who], self.radbins)[0] ) 
        else:
            self.hd = [np.histogram(self.dists, self.radbins)[0]]
    
    def __call__(self, im):
        """Perform and return the radial average of the specrum 'im'"""
        assert im.shape == self.dists.shape
        if self.N > 1:
            avgs = []
            for who, hd in zip(self.whos, self.hd):
                hw = np.histogram(self.dists[who], self.radbins, weights=im[who])[0]
                avgs.append(hw/hd)
            return avgs
        else:
            hw = np.histogram(self.dists, self.radbins, weights=im)[0]
            return [hw/self.hd[0]]
        
class RadialAverager_test(object):
    def __init__(self, shape, N=1, centred_theta=True):
        """
        Initialise the averager class according to the shape of the dft array to be averaged over.
        Slice the image into 2N equally large parts, before and after Pi. centred_theta centres
        slices on the angle, otherwise start at 0 degrees.
        """
        assert len(shape) == 2
        self.shape = shape
        self.N = N
        self.centred_theta = centred_theta
        self.angles = np.pi / N * np.arange(N)
        self.centre = tuple(int(length / 2) for length in shape)
        self.make_masks()

    def __call__(self, spectrum):
        """
        Call the class as a function after initalisation. Multiplies each pixel with either 0 or 1,
        according to slice, then sum each slice and normalise intensity through the number of pixels
        of that slice.
        """
        shifted_spectrum = np.fft.fftshift(spectrum)
        shifted_spectrum[0, 0] = 0
        angle_vectors = np.sum(self.masks * shifted_spectrum, axis=(2, 3)) / self.counts
        return angle_vectors

    def make_masks(self):
        """
        Makes an N x max(shape) array of shape[0] x shape[1] boolean arrays (masks),
        used o switch the pixels on or off.
        """
        self.masks = np.zeros((self.N, self.shape[0], self.shape[0], self.shape[1]), dtype=bool)
        x, y       = np.ogrid[:self.shape[0], :self.shape[1]]
        cx, cy     = self.centre
        dx         = x - cx
        dy         = y - cy
        r2         = dx**2 + dy**2
        for i, angle in enumerate(self.angles):
            start_angle, stop_angle = angle, angle+np.pi/self.N
            if self.centred_theta:
                start_angle -= np.pi / (2 * self.N)
                stop_angle  -= np.pi / (2 * self.N)
            theta = np.arctan2(dx, dy) - start_angle
            # wrap angles between 0 and pi
            theta %= 2*np.pi
            angle_mask = ((theta < (stop_angle - start_angle - 1e-2)) & (theta > 0)) \
                         | ((theta < (stop_angle - start_angle + np.pi - 1e-2)) & (theta > np.pi))
            for r in np.arange(self.shape[0]):
                radius_mask = (r2 >= np.square(r)) & (r2 <= np.square(r+1))
                self.masks[i, r]   = radius_mask * angle_mask
        self.counts = np.sum(self.masks, axis=(2, 3))


def DtoRH(D, temp, dD=None, effectivetemp=False, viscosity=None):
    kb           = 1.381e-23
    if temp < 150:
        temp        += 273.15
    if effectivetemp:
        temperature  = effectiveTemperature(temp)
    else:
        temperature  = temp
    if not viscosity:
        viscosity    = water_viscosity(temperature)
    RH           = kb*temperature/(6*np.pi*viscosity*D)*1e9
    if dD:
        dRH = kb*temperature/(6*np.pi*viscosity*D**2) * dD * 1e9
        return RH, dRH
    return RH
            
def renameTimestamp(directory):
    ts = "[0-1][0-9][0-9][0-9]2[0-9][0-9][0-9][0-2][0-9][0-6][0-9][0-9][0-9]"
    finder = re.compile(ts)
    files = glob(join(directory, '*'))
    i = 1
    known_patterns = {}
    for file in files:
        match = finder.search(file)
        if match:
            pattern = file[match.start():match.end()]
            if not pattern in known_patterns:
                known_patterns[pattern] = i, [file]
            else:
                known_patterns[pattern][1].append(file)
            i += 1
    
    for foundpattern in known_patterns:
        i, filestomove = known_patterns[foundpattern]
        for filetomove in filestomove:
            move(filetomove, filetomove.replace(foundpattern, f"{i:02d}"))
            
            
