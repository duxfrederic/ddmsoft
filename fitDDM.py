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
    keys = [key for key in computedData if "Merged" not in str(key)]
    if not keys:
        raise ValueError("No matrices were supplied")

    def frame_rate(key):
        dts = np.asarray(computedData[key][1])
        if dts.size < 2 or dts[1] <= dts[0]:
            raise ValueError(f"Invalid time grid for matrix {key}")
        return 1 / (dts[1] - dts[0])

    for key in keys:
        ddm, dts, qs = computedData[key]
        if np.asarray(ddm).ndim != 2 or len(dts) != len(ddm) or len(qs) != ddm.shape[1]:
            raise ValueError(f"Incompatible DDM matrix: {key}")

    if mode == 'average':
        groups = []
        for key in keys:
            rate = frame_rate(key)
            group = next((group for group in groups if np.isclose(group[0], rate)), None)
            if group is None:
                group = [rate, []]
                groups.append(group)
            group[1].append(key)

        output_paths = []
        for rate, group_keys in groups:
            first_ddm, dts, qs = computedData[group_keys[0]]
            if any(
                not np.array_equal(dts, computedData[key][1]) or
                not np.array_equal(qs, computedData[key][2]) or
                np.shape(first_ddm) != np.shape(computedData[key][0])
                for key in group_keys[1:]
            ):
                raise ValueError("Matrices being averaged must have matching grids")
            ddm = np.mean(
                [np.asarray(computedData[key][0], dtype=float) for key in group_keys],
                axis=0,
            )
            output_paths.append(_save_merged(
                dirname(group_keys[0]), f"Averaged({rate:g})_{title}", ddm, dts, qs
            ))
        return output_paths

    keys.sort(key=frame_rate, reverse=True)
    ddm, dts, qs = (np.array(value, copy=True) for value in computedData[keys[0]])
    for key in keys[1:]:
        ddmnew, dtsnew, qsnew = computedData[key]
        if not np.array_equal(qs, qsnew):
            raise ValueError("Matrices being merged must use the same q grid")
        ddm, dts = merge(ddm, ddmnew, dts, dtsnew)
    return _save_merged(dirname(keys[0]), f"Merged_{title}", ddm, dts, qs)


def _save_merged(directory, prefix, ddm, dts, qs):
    paths = []
    for musthave, arr in zip(musthaves, (ddm, dts, qs)):
        path = join(directory, prefix + musthave)
        np.save(path, arr)
        paths.append(path)
    return paths


def merge(ddmf, ddms, dtsf, dtss):
    """
        Taken from http://perso.ens-lyon.fr/thomas.gibaud/ddm
        ( https://aapt.scitation.org/doi/10.1119/1.4939516 )
    """ 
    ddmf = np.asarray(ddmf, dtype=float)
    ddms = np.asarray(ddms, dtype=float)
    dtsf = np.asarray(dtsf, dtype=float)
    dtss = np.asarray(dtss, dtype=float)
    if ddmf.ndim != 2 or ddms.ndim != 2:
        raise ValueError("DDM matrices must be two-dimensional")
    if ddmf.shape[1] != ddms.shape[1]:
        raise ValueError("DDM matrices must have the same number of q values")
    if len(dtsf) != len(ddmf) or len(dtss) != len(ddms):
        raise ValueError("Each DDM matrix must match its time grid")
    if len(dtsf) < 2 or len(dtss) < 2 or np.any(np.diff(dtsf) <= 0) or np.any(np.diff(dtss) <= 0):
        raise ValueError("Time grids must be strictly increasing")
    if dtss[0] > dtsf[-1] or dtsf[0] > dtss[-1]:
        raise ValueError("The DDM time grids do not overlap")

    boundary = int(np.argmin(np.abs(dtsf - dtss[0])))
    overlap_end = min(dtsf[-1], dtss[-1])
    fast_end = int(np.searchsorted(dtsf, overlap_end, side='right') - 1)
    if fast_end < boundary:
        raise ValueError("The DDM time grids do not overlap")
    overlap_times = dtsf[boundary:fast_end + 1]

    if np.any(ddms[0] == 0):
        raise ValueError("Cannot merge a matrix with zero initial amplitude")
    scale = ddmf[boundary] / ddms[0]
    slow_scaled = ddms * scale
    interpolated = np.column_stack([
        np.interp(overlap_times, dtss, slow_scaled[:, column])
        for column in range(ddms.shape[1])
    ])
    if len(overlap_times) == 1:
        transition = interpolated
    else:
        blend = np.linspace(0, 1, len(overlap_times))[:, None]
        transition = (1 - blend) * ddmf[boundary:fast_end + 1] + blend * interpolated

    slow_tail = np.searchsorted(dtss, overlap_times[-1], side='right')
    dts = np.concatenate([dtsf[:boundary], overlap_times, dtss[slow_tail:]])
    ddm = np.concatenate([ddmf[:boundary], transition, slow_scaled[slow_tail:]], axis=0)
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
    tauq      = DTS * QS**2
    Gamma     = params[2]  # Gamma mean, rescaled.
    cumulants = params[3:]  # mu2, mu3, mu4 ...
    #  mean decay (first moment)
    meandecay = np.exp(- Gamma * tauq)
    # add more
    deviation = 1.
    for i, cumulant in enumerate(cumulants):
        order = i + 2
        deviation += (-1)**order * cumulant * tauq**order / factorial(order)
    return deviation * meandecay


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
    if len(ddm) < 1:
        raise ValueError("A DDM matrix must contain at least one time point")
    aplusb      = np.mean(ddm[-min(3, len(ddm)):], axis=0)
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
    ini            = list(ini)
    fixed          = list(fixed)
    ddmopt         = ddm[dtmin:dtmax, qmin:qmax]
    
    dtsopt         = dts[dtmin:dtmax]
    qsopt          = qs[qmin:qmax]
    if ddmopt.size == 0 or len(qsopt) == 0 or len(dtsopt) == 0:
        raise ValueError("The selected q and time ranges must contain data")
    QSopt, DTSopt  = np.meshgrid(qsopt, dtsopt)

    a, b = ini[-2:]
    if a == '':
        fixed[-2] = False
    if b == '':
        fixed[-1] = False
    a, b = est_A_B(ddmopt, a, b)
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
            for cumulant_index in range(order):
                if (cumulant_index + 1) % 2 == 0:
                    opt.x[cumulant_index + 2] = abs(opt.x[cumulant_index + 2])
                cumulants[cumulant_index].append(opt.x[cumulant_index + 2])
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
