#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 17:07:08 2021

@author: fred
"""

import   numpy              as      np
from     scipy.optimize     import  minimize

def regularizedSquareSum(g, y, A, alpha, Omega, weights):
    # unpack everything we put in "g":
    amplitude, noise = g[:2]
    scale = (amplitude + noise)
    g = g[2:]
    
    # the extremes of g should be pushed towards 0.
    reg1  = 0.1 * scale * ( np.abs(g[0]) + np.abs(g[-1]) )
    if np.isnan(reg1):
        reg1 = 0
    
    # g is positive anyways
    g     = np.abs(g)
    
    # the result as is: (@ is the matrix multiplication)
    yfit  = amplitude * (1 - A @ g) + noise
    
    # sum the of the weighted square residuals:
    var   = np.nansum( weights * (y-yfit)**2 )
    
    # regularization: curvature; take the 2nd derivative of g
    r     = np.diff(g, 2)
    reg2  = scale * alpha**2 * np.nansum( r**2 ) 
    # if an additional regularization or prior is encoded in the matrix Omega,
    # uncomment the following line and comment out the above one.
    # reg2  = alpha**2 * np.nansum( (r - Omega@g)**2 )
        
    # also keep sum(g) = 1:
    # reg3 = 0.1 * scale * np.abs(np.sum(g)-1)
    # if np.isnan(reg3):
        # reg3 = 0
    # add everything to the score to minimize:
    
    return var + reg2 + reg1



def CONTIN(tau, ddmdata, gammaRange, g0=None, weights=None, alpha=0.1, 
           maxiter=10, tol=1e-3):
    # build the matrix that convolutes g to f.
    # s is the decay rate, tau is (time * q²).
    S, T    = np.meshgrid(gammaRange, tau)
    A       = np.exp(-T * S)
    # so A@g is indeed a sum of exponential decays with rates distributed 
    # according to g. 
    

    # equal weights unless specified:
    if weights is None: 
        weights = np.ones_like(ddmdata)
        weights = 1 / np.sqrt(np.arange(1,len(ddmdata)+1))
        
    
    # in the literature, the regularization is often encoded in a matrix
    # Omega that acts on the distribution g 
    # (see e.g. Scotti, Liu et al 2015: journal of phys. chem.  142, 234905)
    # however, here the curvature will simply be included by taking
    # finite differences twice. Thus I'll set it to zero.
    # see as well the end of the regularizedSquareSum function.
    Omega     = np.zeros((len(gammaRange)-2, len(gammaRange)))

    
    # estimate the amplitude as usual in DDM. 
    aplusb    = (ddmdata[-1]+ddmdata[-2])/2
    b         = ddmdata[0]
    amplitude = aplusb - b
    noise     = b
    
    # we'll store some stuff along the way. 
    alpha_residuals = []
    alpha_g         = []
    alpha_ddmfit    = []
    alpha_amplitude = []
    alpha_noise     = []
    # make sure alpha is iterable:
    try:
        iter(alpha)
    except:
        alpha   = [alpha]
    alphas      = list(alpha)
    # Nested loops here for the sake of robustness: 
    # scan over a range of the weight of the regularization alpha;
    # and guide the optimizer by feeding it its previous estimation as initial
    # guess everal times. 
    for alpha in alphas:
        yield alpha
        # build some stupid initial guess for the distribution of 
        # decay rates, g:
        if g0 is None:
            g     = np.ones_like(gammaRange)
            g[0]  = 0; g[-1] = 0
        else:
            g     = g0
        # normalize it:
        # g = ( np.sum(A@g)  / np.sum(ddmdata) ) * g
        g = g / np.sum(g)
        # the amplitude and noise will be fed to the solver in the same array:
        # (because no way around it and all variables are to be tuned)
        g = np.concatenate( ([amplitude, noise], g) )
        trackresiduals = []
        for i in range(maxiter):
            # everytime, minimize with the previous solution as a minimal guess
            sol = minimize(regularizedSquareSum, g, 
                           args=(ddmdata, A, alpha, Omega, weights),
                           method='Nelder-Mead')
            g  = np.abs(sol.x)
            
            residuals = sol.fun 
            trackresiduals.append(residuals)
            # when the residuals are stable, stop:
            if len(trackresiduals) > 1:
                diff = np.abs(np.diff(trackresiduals))[-1]
                if diff < np.abs(trackresiduals[-1]) * tol:
                    break
        # for this value of the regularization weight alpha, store the 
        # score we obtained:
        alpha_residuals.append(trackresiduals[-1])
        # and as well the distribution we got:
        alpha_g.append(g[2:])
        alpha_ddmfit.append(A @ g[2:])
        alpha_amplitude.append(g[0])
        alpha_noise.append(g[1])
    # now, using the L-curve criterion, select the "optimal" alpha:
    minimum          = np.argmin(alpha_residuals)
    # and give the user the corresponding solution
    g                = alpha_g[minimum]
    alpha            = alphas[minimum]
    ddmfit           = alpha_ddmfit[minimum]
    amplitude        = alpha_amplitude[minimum]
    noise            = alpha_noise[minimum]
    sol = CONTINsolution(g, alpha, tau, ddmfit, ddmdata, amplitude, noise, 
                         alpha_residuals, alphas, gamma_range=gammaRange,
                         alpha_g=alpha_g, alpha_ddmfit=alpha_ddmfit,
                         alpha_amplitude=alpha_amplitude, alpha_noise=alpha_noise)
    yield sol


class CONTINsolution():
    def __init__(self, g, chosen_alpha, tau, ddmfit, ddmdata, amplitude, noise, 
                 alpha_residuals, alphas, gamma_range, alpha_g, alpha_ddmfit,
                 alpha_amplitude, alpha_noise):
        self.g               = g
        self.chosen_alpha    = chosen_alpha
        self.ddmfit          = ddmfit
        self.amplitude       = amplitude
        self.noise           = noise
        self.alpha_residuals = alpha_residuals 
        self.alphas          = alphas
        self.gamma_range     = gamma_range
        self.tau             = tau
        self.ddmdata         = ddmdata
        self.alpha_g         = alpha_g
        self.alpha_ddmfit    = alpha_ddmfit
        self.alpha_amplitude = alpha_amplitude
        self.alpha_noise     = alpha_noise
        self.sizes           = None
        