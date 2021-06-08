#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A small GUI program to interface a custom DDM setup.

@author: Frédéric Dux, biosoft intern@IPC with Jerome Crassous
"""

import numpy              as      np

from   scipy.optimize     import  minimize
from   os.path            import  join, dirname
from   math               import  factorial
from   utilities          import  musthaves


mycolormap = 'plasma' 


def mergeDDM(computedData, mode='merge', title=''):
    params       = {}
    keys  = list(computedData.keys())
    keys  = [key for key in keys if (not "Merged" in key)]
    framerateold = -10
    for vid in keys:
        param        = {}
        framerate    = 1/(computedData[vid][1][1]-computedData[vid][1][0])
        param['framerate'] = framerate
        print(framerate)
        if framerate == framerateold and mode == 'merge':
            print("Videos to be merged with the same framerate: error?")
        framerateold = framerate
        params[vid]  = param
    keys.sort(key=lambda vid: -float(params[vid]['framerate']))
    first = keys[0]
    ddm, dts, qs = computedData[first]
    framerateold = params[first]['framerate']
    if mode == 'average':
        i = 1
        while i < len(keys):
            N = 1
            while i < len(keys) and params[keys[i]]['framerate'] == framerateold :
                if not "Averaged(" in keys[i]:
                    ddm += computedData[keys[i]][0] 
                    N   += 1
                    framerateold = params[keys[i]]['framerate']
                i   += 1
            ddm = ddm / N
            tosave = [ddm, dts, qs]
            for musthave, arr in zip(musthaves, tosave):
                np.save(join(dirname(vid),f"Averaged({framerateold})_{title}"+musthave), arr) 
            i += 1
    else:
        for i in range(1,len(keys)):
            ddmnew, dtsnew, qsnew = computedData[keys[i]]
            ddm, dts = merge(ddm, ddmnew, dts, dtsnew)
        tosave = [ddm, dts, qs]
        for musthave, arr in zip(musthaves, tosave):
            np.save(join(dirname(vid),f"Merged_{title}"+musthave), arr) 


def merge(ddmf, ddms, dtsf, dtss):
    """
        adapted from http://perso.ens-lyon.fr/thomas.gibaud/ddm
        ( https://aapt.scitation.org/doi/10.1119/1.4939516 )
    """ 
    # Find the closest time at the fast freq to the smallest time at the small freq
    boundary     = np.argmin(np.abs(dtsf - dtss[0]))
    # Rescale the value at the slow freq according to the value at t=boundary for the fast freq
    ddms        *= ddmf[boundary] / ddms[0]
    # find the first third of their overlap
    overlap0     = (len(ddmf)-1 - boundary)
    overlap1     = np.argmin(np.abs(dtss - dtsf[boundary+overlap0]))
    #interpolate on this first third the DDM at 4Hz on the times at 400Hz
    interpolated = np.transpose([
        np.interp(
            dtsf[boundary:boundary+overlap0],
            dtss[:overlap1], 
            v)
        for v in ddms[:overlap1].T])
    #do a smooth transition on this first third
    x          = ((dtsf[boundary:boundary+overlap0]-dtsf[boundary])/(dtsf[boundary+overlap0]-dtsf[boundary]))[:,None]
    transition = (1-x) * ddmf[boundary:boundary+overlap0] + x * interpolated
    # Merge f Hz, transition to s Hz
    dts = np.concatenate([dtsf[:boundary+overlap0], dtss[overlap1:]])
    ddm = np.concatenate([ddmf[:boundary], transition, ddms[overlap1:]], axis=0)
    #"""
    return ddm, dts     




def getTemperature(D, viscosity, known_radius):
    kb  = 1.381e-23
    eta = viscosity
    T   = D * 6 * np.pi * eta * known_radius / kb
    return T

def getRadius(D, viscosity, temperature):
    kb  = 1.381e-23
    eta = viscosity
    T   = temperature
    r   = kb * T / ( 6 * np.pi * eta * D)
    return r
    

def fToDDM(f, A, B):
    return A * (1 - f) + B


def single_exponential(params, QS, DTS):
    D       = params[2]
    Gamma   = D*QS**2
    f       = np.exp(-DTS*Gamma)
    return f


def cumulant_exponential(params, QS, DTS):
    numorder  = len(params) - 2
    Gamma     = params[2]*QS**2
    cumulants = []
    for order in range(0, numorder-1):
        cumu = params[ 3+order ]
        if (order+1) % 2 == 0:
            cumu = np.abs(cumu)
        cumulants.append(cumu)
    mean_decay = np.exp(-DTS*Gamma)
    lineardev  = 1
    for i, cumulant in enumerate(cumulants):
        exponent   = i + 2
        lineardev += (-cumulant*DTS*QS**2)**exponent / factorial(exponent)
    mean_decay    *= lineardev
    return mean_decay


def stretch_exponential(params, QS, DTS):
    Gamma   = params[2]*QS**2
    beta    = params[3]
    decay   = np.exp(-(DTS*Gamma)**beta)
    return decay


def dbl_exponential_stretch(params, QS, DTS):
    Gamma1 = params[2]*QS**2
    Gamma2 = params[3]*QS**2
    beta2  = params[4]
    alpha  = params[5]
    decay  = alpha * np.exp(-DTS*Gamma1) + (1-alpha) * np.exp(- (DTS*Gamma2)**beta2 )
    return decay


def exponential_with_flow(params, QS, DTS):
    Gamma  = params[2]*QS**2
    # v_flow is the effective, projected speed of the particles defining the flow
    # aka, once we obtain v_flow, must divide it by cos(theta) 
    v_flow = params[3]
    decay  = np.exp(-DTS*Gamma) * np.cos(QS * v_flow * DTS)
    return decay

def stretch_exponential_with_flow(params, QS, DTS):
    Gamma  = params[2]*QS**2
    # v_flow is the effective, projected speed of the particles defining the flow
    # aka, once we obtain v_flow, must divide it by cos(theta) 
    v_flow = params[3]
    beta   = params[4]
    decay  = np.exp(-(DTS*Gamma)**beta) * np.cos(QS * v_flow * DTS)
    return decay

def dbl_exponential_stretch_with_flow(params, QS, DTS):
    Gamma1  = params[2]*QS**2
    Gamma2  = params[3]*QS**2
    beta2   = params[4]
    alpha   = params[5]
    v_flow  = params[6]
    # v_flow is the effective, projected speed of the particles defining the flow
    # aka, once we obtain v_flow, must divide it by cos(theta) 
    decay  = (alpha * np.exp(-DTS*Gamma1) + (1-alpha) * np.exp(- (DTS*Gamma2)**beta2 )) * np.cos(QS * v_flow * DTS)
    return decay


def ddm_penalty(params, ddm, QS, DTS, model, fixed, ini):
    params       = [p if f == False else i  for p,f,i in zip(params, fixed, ini)]
    A, B         = params[0], params[1]
    base_penalty = ( ddm-fToDDM(model(params, QS, DTS), A, B) )**2 / (np.median(DTS) + DTS)
    return (base_penalty).flatten()

def ddm_penalty_sum(params, ddm, QS, DTS, model, fixed, ini):
    return ddm_penalty(params, ddm, QS, DTS, model, fixed, ini).sum()

def ddm_penalty_dbl_exponential(params, ddm, QS, DTS, model, fixed, ini):
    params   = [p if f == False else i  for p,f,i in zip(params, fixed, ini)]
    B        = params[1]
    negB     = B if B < 0 else 0
    alpha    = params[5]
    notgood  = np.abs(alpha) if alpha > 1 or alpha < 0 else 0
    beta2    = params[4]
    notgood2 = np.abs(beta2) if beta2 > 1 or beta2 < 0 else 0
    base     = ddm_penalty(params, ddm, QS, DTS, model, fixed, ini).sum()
    reg      = 1e16*np.abs(notgood) + 1e16*np.abs(notgood2) + 1e7*(negB**2)
    return base+reg

def est_A_B(ddm, a='', b=''):
    aplusb      = (ddm[-1,:]+ddm[-2,:]+ddm[-3,:] ) / 3
    if b == '':
        b       = ddm[0, :]
    else:
       b = np.ones(ddm.shape[1])*b
    if a == '':
        a       = aplusb-b
    else:
        a = np.ones(ddm.shape[1])*a
    return a, b

def fitOneDDMmatrix(ddm_dts_qs, model, ini, fixed, qmin=0, 
                    qmax=None, dtmin=0, dtmax=None):
    """
        returns A_fit, B_fit, modelparams_fit, f_analytical, opt_object
    """
    ddm, dts, qs   = ddm_dts_qs
    ddmopt         = ddm[dtmin:dtmax:,qmin:qmax]
    
    dtsopt         = dts[dtmin:dtmax]
    qsopt          = qs[qmin:qmax]
    QSopt, DTSopt  = np.meshgrid(qsopt, dts)

    a, b = ini[-2:]
    if a == '':
        fixed[-2] = False
    if b == '':
        fixed[-1] = False
    a, b = est_A_B(ddm, a, b)
    # now rotate: a, b are to the front in the optimization. (whereas they
    # were put at the end in the front end)
    ini   = list(np.roll(ini, 2))
    fixed = list(np.roll(fixed, 2))
    # an opt object stating "False" to trick the first optimization in the loop
    # to use the default ini values
    class optdecoy():
        def __init__(self):
            self.success = False
    opt = optdecoy()
    # save the given ini value:
    ini_arg = ini.copy()
    
    ## now we are model dependant!
    
    if model == 'single_exponential':
        A, B, D, fs = [], [], [], []
        for i, (q, ai, bi) in enumerate(zip(qsopt, a, b)):
            ddmrow = ddmopt[:, i]
            if opt.success == True:
                ini_arg[:2] = ai, bi
                ini_arg = [float(e) for e in ini_arg]
                ini = opt.x
            else:
                ini_arg[:2] = ai, bi
                ini_arg = [float(e) for e in ini_arg]
                ini =  ini_arg
            opt    = minimize(ddm_penalty_sum, ini, method='Nelder-Mead',
                           args=(ddmrow, q, dtsopt, single_exponential, fixed, ini_arg))
            opt.x = [e if f == False else i for e,f,i in zip(opt.x, fixed, ini_arg)]
            A.append(opt.x[0])
            B.append(opt.x[1])
            D.append(opt.x[2])
            f = single_exponential(opt.x, q, dts)
            fs.append(f)
        f = np.vstack(fs).T
        A, B, D = np.array(A), np.array(B), np.array(D)
        return A, B, [D], f, opt
       
    
    elif model.startswith('cumulant_'):
        # desired order
        order = int(model.replace('cumulant_', ''))
        A, B, cumulants, fs = [], [], [[] for _ in range(order)], []
        for i, (q, ai, bi) in enumerate(zip(qsopt, a, b)):
            if opt.success == True:
                ini_arg[:2] = ai, bi
                ini_arg = [float(e) for e in ini_arg]
                ini = opt.x
            else:
                ini_arg[:2] = ai, bi
                ini_arg = [float(e) for e in ini_arg]
                ini =  ini_arg
            ddmrow = ddmopt[:, i]
            opt    = minimize(ddm_penalty_sum, ini, method='Nelder-Mead', tol=1e-9,
                           args=(ddmrow, q, dtsopt, cumulant_exponential, fixed, ini_arg))
            opt.x = [e if f == False else i for e,f,i in zip(opt.x, fixed, ini_arg)]
            A.append(opt.x[0])
            B.append(opt.x[1])
            [cumulants[i].append(opt.x[i+2]) for i in range(order)]
            cumulants = [list(np.abs(cumulants[i])) if (i+1)%2 == 0 else cumulants[i] for i in range(order)]
            f = cumulant_exponential(opt.x, q, dts)
            fs.append(f)
        f = np.vstack(fs).T
        A, B = np.array(A), np.array(B)
        cumulants = [np.array(cumulants[i]) for i in range(order)]
        return A, B, cumulants, f, opt

        
    elif model == 'stretch':
        A, B, D, beta, fs = [], [], [], [], []
        for i, (q, ai, bi) in enumerate(zip(qsopt, a, b)):
            if opt.success == True:
                ini_arg[:2] = ai, bi
                ini_arg = [float(e) for e in ini_arg]
                ini = opt.x
            else:
                ini_arg[:2] = ai, bi
                ini_arg = [float(e) for e in ini_arg]
                ini =  ini_arg
            ddmrow = ddmopt[:, i]
            opt    = minimize(ddm_penalty_sum, ini, method='Nelder-Mead',
                           args=(ddmrow, q, dtsopt, stretch_exponential, fixed, ini_arg))
            opt.x = [e if f == False else i for e,f,i in zip(opt.x, fixed, ini_arg)]
            A.append(opt.x[0])
            B.append(opt.x[1])
            D.append(opt.x[2])
            beta.append(opt.x[3])
            f = stretch_exponential(opt.x, q, dts)
            fs.append(f)
        f = np.vstack(fs).T
        A, B, D, beta = np.array(A), np.array(B), np.array(D), np.array(beta)
        return A, B, [D, beta], f, opt
        

    elif model == "dblexp_2ndstretched":
        """
        "dblexp_2ndstretched": ["Diffusion coefficient 1", "Diffusion coefficient 2",\
                                         "stretch coefficient (2)", "weighting parameter"]
        """
        A, B, D1, D2, beta2, alpha, fs = [], [], [], [], [], [], []
        for i, (q, ai, bi) in enumerate(zip(qsopt, a, b)):
            if opt.success == True:
                ini_arg[:2] = ai, bi
                ini_arg = [float(e) for e in ini_arg]
                ini = opt.x
            else:
                ini_arg[:2] = ai, bi
                ini_arg = [float(e) for e in ini_arg]
                ini =  ini_arg
            ddmrow = ddmopt[:, i]
            opt    = minimize(ddm_penalty_dbl_exponential, ini, method='Nelder-Mead',
                           args=(ddmrow, q, dtsopt, dbl_exponential_stretch, fixed, ini_arg))
            opt.x = [e if f == False else i for e,f,i in zip(opt.x, fixed, ini_arg)]
            A.append(opt.x[0])
            B.append(opt.x[1])
            D1.append(opt.x[2])
            D2.append(opt.x[3])
            beta2.append(opt.x[4])
            alpha.append(opt.x[5])
            f = dbl_exponential_stretch(opt.x, q, dts)
            fs.append(f)
        f = np.vstack(fs).T
        A, B, D1, D2 = np.array(A), np.array(B), np.array(D1), np.array(D2)
        beta2, alpha = np.array(beta2), np.array(alpha)
        return A, B, [D1, D2, beta2, alpha], f, opt
            
    elif model == 'expcos':
        A, B, D, v_flow, fs = [], [], [], [], []
        for i, (q, ai, bi) in enumerate(zip(qsopt, a, b)):
            if opt.success == True:
                ini_arg[:2] = ai, bi
                ini_arg = [float(e) for e in ini_arg]
                ini = opt.x
            else:
                ini_arg[:2] = ai, bi
                ini_arg = [float(e) for e in ini_arg]
                ini =  ini_arg
            ddmrow = ddmopt[:, i]
            opt    = minimize(ddm_penalty_sum, ini, method='Nelder-Mead',
                           args=(ddmrow, q, dtsopt, exponential_with_flow, fixed, ini_arg))
            opt.x = [e if f == False else i for e,f,i in zip(opt.x, fixed, ini_arg)]
            A.append(opt.x[0])
            B.append(opt.x[1])
            D.append(opt.x[2])
            v_flow.append(opt.x[3])
            f = exponential_with_flow(opt.x, q, dts)
            fs.append(f)
        f = np.vstack(fs).T
        A, B, D, v_flow = np.array(A), np.array(B), np.array(D), np.abs(v_flow)
        return A, B, [D, v_flow], f, opt
    
    elif model == 'expcosstretch':
        A, B, D, v_flow, beta, fs = [], [], [], [], [], []
        for i, (q, ai, bi) in enumerate(zip(qsopt, a, b)):
            if opt.success == True:
                ini_arg[:2] = ai, bi
                ini_arg = [float(e) for e in ini_arg]
                ini = opt.x
            else:
                ini_arg[:2] = ai, bi
                ini_arg = [float(e) for e in ini_arg]
                ini =  ini_arg
            ddmrow = ddmopt[:, i]
            opt    = minimize(ddm_penalty_sum, ini, method='Nelder-Mead',
                           args=(ddmrow, q, dtsopt, stretch_exponential_with_flow, fixed, ini_arg))
            opt.x = [e if f == False else i for e,f,i in zip(opt.x, fixed, ini_arg)]
            A.append(opt.x[0])
            B.append(opt.x[1])
            D.append(opt.x[2])
            v_flow.append(opt.x[3])
            beta.append(opt.x[4])
            f = stretch_exponential_with_flow(opt.x, q, dts)
            fs.append(f)
        f = np.vstack(fs).T
        A, B, D, v_flow, beta = np.array(A), np.array(B), np.array(D), np.abs(v_flow), np.array(beta)
        return A, B, [D, v_flow, beta], f, opt
    
    elif model == "dblexpcosstretch":
        """
        "dblexpcosstretch": ["Diffusion coefficient 1", "Diffusion coefficient 2",\
                               "stretch coefficient (2)", "weighting parameter", "effective speed"]
        """
        A, B, D1, D2, beta2, alpha, vs, fs = [], [], [], [], [], [], [], []
        for i, (q, ai, bi) in enumerate(zip(qsopt, a, b)):
            if opt.success == True:
                ini_arg[:2] = ai, bi
                ini_arg = [float(e) for e in ini_arg]
                ini = opt.x
            else:
                ini_arg[:2] = ai, bi
                ini_arg = [float(e) for e in ini_arg]
                ini =  ini_arg
            ddmrow = ddmopt[:, i]
            opt    = minimize(ddm_penalty_dbl_exponential, ini, method='Nelder-Mead',
                           args=(ddmrow, q, dtsopt, dbl_exponential_stretch_with_flow, fixed, ini_arg))
            opt.x = [e if f == False else i for e,f,i in zip(opt.x, fixed, ini_arg)]
            A.append(opt.x[0])
            B.append(opt.x[1])
            D1.append(opt.x[2])
            D2.append(opt.x[3])
            beta2.append(opt.x[4])
            alpha.append(opt.x[5])
            vs.append(opt.x[6])
            f = dbl_exponential_stretch_with_flow(opt.x, q, dts)
            fs.append(f)
        f = np.vstack(fs).T
        A, B, D1, D2 = np.array(A), np.array(B), np.array(D1), np.array(D2)
        beta2, alpha, vs = np.array(beta2), np.array(alpha), np.array(vs)
        return A, B, [D1, D2, beta2, alpha, vs], f, opt
    
    elif model == 'None':
        f      = 1 - (ddm-b)/(a[np.newaxis,:])
        return a, b, [None], f, None
    
    else:
        return "Not implemented"
