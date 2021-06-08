#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A small GUI program to interface a custom DDM setup.

@author: Frédéric Dux, biosoft intern@IPC with Jerome Crassous
"""

import     numpy                     as      np
import     matplotlib.pyplot         as      plt
from       matplotlib.widgets        import  Slider
from       datetime                  import  datetime
from       utilities                 import  musthaves, extractCrudef
plt.switch_backend('TkAgg')
mycolormap = 'plasma'

def getUniqueNum(comment, title):
    title = title.replace(musthaves[0], '')
    date  = datetime.strftime(datetime.now(), '%H:%S.%f')[:-4] + ' | '
    return date + comment + ' ' + title

def plotAutocorrelation(qs, dts, f, title=None):
    if title:
        plt.figure(num=getUniqueNum('autocor grid:' , title))
    else:
        plt.figure()
    QS, DTS = np.meshgrid(qs, dts)
    plt.pcolormesh(QS/1e6,DTS,f,vmin=0,vmax=1,cmap=mycolormap)
    plt.xlabel('$q$ [μm$^{-1}$]')
    plt.ylabel('$\Delta t$ [s]')
    plt.tight_layout()
    plt.show(block=False)


def plotSomefs(qmin, qmax, datasets, title=None):
    qmin, qmax = int(qmin), int(qmax)
    if title:
        plt.figure(num=getUniqueNum('some fs:', title))
    else:
        plt.figure()
    # we would like to plot 5 qs.
    if (qmax - qmin +1) >= 5:
        increment = (qmax-qmin + 1)//5
    else:
        increment = 1
        
    # now we plot the raw stuff:
    tsraw, qsraw, ddmraw, fraw = datasets[0]
    if fraw is None:
        if len(datasets)>1:
            fraw  = extractCrudef(ddmraw[:,qmin:qmax], A=datasets[1][4], B=datasets[1][5])
        else:
            fraw  = extractCrudef(ddmraw[:,qmin:qmax])
        datasets[0][3] = fraw
    qsraw = qsraw[qmin:qmax]
    # we'll need the colors in case we also plot the fits:
    colors = []
    
    for qindex in range(0, qmax-qmin, increment):
        q = qsraw[qindex]
        label      = r"{0:.02f} μm$^{{-1}}$".format(float(q)/1e6)
        p = plt.semilogx(tsraw*q**2, fraw[:,qindex], marker='D', label=label, mfc="None",
                     ls="None", zorder=-qindex)
        colors.append(p[0].get_color())
    
    if len(datasets) > 1:
        tsest, qsest, ddmest, f_fit, A_fit, B_fit = datasets[1]
        for qindex, color in zip(range(0, qmax-qmin, increment), colors):
            q = qsest[qindex]
            plt.semilogx(tsest*q**2, f_fit[:,qindex], color=color, alpha=0.9, zorder=1+qindex)
    
    plt.legend(frameon=False)
    plt.xlabel(r"$\tau \,q^2$ [s/m$^2$]")
    plt.ylabel(r"$f(q,\tau)$")
    plt.tight_layout()
    plt.show(block=False)



def plotFitParams(qs, fitparams, names, title=None):
    if title:
        fig = plt.figure(num=getUniqueNum('fit params:', title))
    else:
        fig = plt.figure()
    Np      = len(fitparams)
    axs     = fig.subplots(Np, 1, sharex=True)
    fig.set_figheight(3+2*Np)
    if Np  == 1:
        axs = [axs]
    for name, params, ax in zip(names, fitparams, axs):
        if name == "2nd cumulant" or name == "3rd cumulant":
            to_plot = params / fitparams[0]
            name += " / D"
        else :
            to_plot = params
        if name == "effective flow speed":
            # to_plot *= 1e6
            name += " [m/s]"
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        elif "Diffusion coefficient" in name:
            name += " [m$^2$/s]"
        ax.plot(qs/1e6, to_plot, 'o', mfc="None",ls="None")
        ax.set_ylabel(name)
    ax.set_xlabel('q [μm$^{-1}$]')
    plt.show(block=False)



def loadDDM(path):
    path_dts = path.replace('_DDM_matrix', '_deltaTs')
    path_qs  = path.replace('_DDM_matrix', '_QS')
    return np.load(path_dts), np.load(path_qs), np.load(path)


def plotABandR(qs, A, B, D, temperature=None, viscosity=None, title=None):
    if title:
        fig = plt.figure(num=getUniqueNum('A, B, D or R:', title), figsize=(5,5))
    else:
        fig = plt.figure(figsize=(5,5))
    if temperature and viscosity:
        kb       = 1.381e-23
        yval     = kb*temperature/(6*np.pi*viscosity*(D))*1e9
        title    = "Amplitude, noise and hydrodynamic radii"
        ylabel   = "$R_H$ [nm]" 
    else:
        title    = "Amplitude, noise and diffusion coefficient"
        ylabel   = "$D$ [m$^2$/s]"
        yval     = D
    ax1, ax2 = fig.subplots(2,1, sharex=True)
    ax1.semilogy(qs/1e6, A, '.', label='$A(q)$ (amplitude)')
    ax1.semilogy(qs/1e6, B, '.', label='$B(q)$ (noise)')
    ax1.set_ylabel('Amplitude [~]')
    ax2.set_xlabel('$q$ [μm$^{-1}$]')
    ax2.plot(qs/1e6, yval,'.',color='crimson')
    yvalmean = np.mean(yval)
    yvalstd  = np.std(yval)
    ax2.axhline(yvalmean, color='gray', lw=2.5, zorder=-1000)
    ax2.axhline(yvalmean-yvalstd, color='gray', lw=2, ls=':', zorder=-1000)
    ax2.axhline(yvalmean+yvalstd, color='gray', lw=2, ls=':', zorder=-1000)
    ax2.set_ylabel(ylabel)
    ax1.legend(frameon=False)
    plt.tight_layout()
    plt.show(block=False)
    

def plotDDMmeasAndFit(qs, dts, ddm, qsfit=None, dtsfit=None, ddmfit=None, title=None):
    if title:
        fig = plt.figure(num=getUniqueNum('DDM grid:', title))
    else:
        fig = plt.figure()
    QS, DTS       = np.meshgrid(qs, dts)
    
    ax1, ax2      = fig.subplots(1,2, sharey=True)
    ax1.set_title('DDM matrix: measure')
    ax2.set_title('DDM matrix: fit')
    ax1.set_ylabel('$\Delta t$ [s]')
    ax1.set_xlabel('$q$ [μm$^{-1}$]')
    vmin, vmax    = 0.9 * np.percentile(ddm, 5), 1.1 * np.percentile(ddm, 95)
    ax1.pcolormesh(QS/1e6, DTS, ddm, vmin=vmin, vmax=vmax, cmap=mycolormap)
    ax2.set_ylabel('$\Delta t$ [s]')
    ax2.set_xlabel('$q$ [μm$^{-1}$]')
    if not qsfit is None:
        QSfit, DTSfit = np.meshgrid(qsfit, dtsfit)
        ax2.pcolormesh(QSfit/1e6, DTSfit, ddmfit, vmin=vmin, vmax=vmax, cmap=mycolormap)
        ax1.set_xlim((qsfit.min()/1e6,qsfit.max()/1e6))
        ax2.set_xlim((qsfit.min()/1e6,qsfit.max()/1e6))
    plt.tight_layout()
    plt.show(block=False)


def plotDDMmatrix(datasets, amplitude_and_noise=None, title=None, dtsmin=0, dtsmax=-1):
    if title:
        fig = plt.figure(num=getUniqueNum('main plot:', title))
    else:
        fig = plt.figure()
    global datasetsPush, lautoraw, lddmraw, lautoest, lauto_fit, lddmest, sfreq, t1,\
           dtsminlim, dtsmaxlim, dtsminpush, dtsmaxpush
    dtsminpush   = dtsmin
    dtsmaxpush   = dtsmax
    datasetsPush = datasets
    axs          = fig.subplots(1,2)
    fig.set_figheight(6); fig.set_figwidth(12)
    axs[1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    tsraw, qsraw, ddmraw, fraw = datasetsPush[0]
    if fraw is None:
        fraw  = extractCrudef(ddmraw)
        datasetsPush[0][3] = fraw
    
    if len(datasetsPush) > 1:
        qi    = np.argmin(np.abs(datasetsPush[1][1][0]-qsraw))
        qf    = np.argmin(np.abs(datasetsPush[1][1][-1]-qsraw))
        print("qi at plotDDM::plotDDMmatrix:", qi)
    else:
        qi    = len(qsraw)//3
        
    lautoraw, = axs[0].semilogx(qsraw[qi]**2*tsraw, fraw[:,qi], 's')
    lddmraw,  = axs[1].semilogx(tsraw, ddmraw[:,qi], 's')
    dtsminlim = axs[0].axvline(tsraw[dtsminpush]*qsraw[qi]**2, ls='--', color='0.5')
    dtsmaxlim = axs[0].axvline(tsraw[dtsmaxpush]*qsraw[qi]**2, ls='--', color='0.5')    
    axs[1].axvline(tsraw[dtsminpush], ls='--', color='0.5')
    axs[1].axvline(tsraw[dtsmaxpush], ls='--', color='0.5')  
    
    axcolor   = 'lightgoldenrodyellow'
    axfreq    = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    sfreq     = Slider(axfreq, 'Freq', 0, len(qsraw)-1, valinit=qi, valstep=1)
    t1        = plt.text(0.5, 1.2, fr'$q$: {qsraw[qi]/1e6:.01f} μm$^{{-1}}$')
    
    #### aesthetics
    axs[0].set_ylim((-0.1,1.1))
    axs[0].set_xlabel(r"$\tau q^2$ [s/m$^2$]")
    axs[0].set_ylabel(r"$f(q,\tau)$")
    axs[0].set_title(r"Correlation function $f(q, \tau)$")
    axs[1].set_xlabel(r"$\tau$ [s]")
    axs[1].set_ylabel(r"$\mathcal{D}(q,\tau)$")
    axs[1].set_title(r"$\mathcal{D}(q,\tau)=A(q)\,(1-f(q,\tau) + B(q)$")
    
    plt.subplots_adjust(left=0.1, bottom=0.25)

    if len(datasetsPush) > 1:
        tsest, qsest, ddmest, f_fit, A_fit, B_fit = datasetsPush[1]
        fest       = extractCrudef(ddmraw[:,qi:qf+1], A_fit, B_fit)
        qi1        = np.argmin(np.abs(qsraw[qi]-qsest))
        lauto_fit, = axs[0].semilogx(qsest[qi1]**2*tsest, f_fit[:,qi1])
        lautoest,  = axs[0].semilogx(qsest[qi1]**2*tsest, fest[:,qi1], 'D', mfc='None')
        lddmest,   = axs[1].semilogx(tsest, ddmest[:,qi1])
    
    def press(event):
        try:
            button = event.button
        except:
            button = 'None'
        if event.key == 'right' or button == 'down':
            if sfreq.val < len(qsraw) - 1:
                sfreq.set_val(sfreq.val + 1)
        elif event.key == 'left' or button == 'up':
            if sfreq.val > 0:
                sfreq.set_val(sfreq.val - 1)
        qi0 = int(sfreq.val)
        t1.set_text(fr'$q$: {qsraw[qi0]/1e6:.01f} μm$^{{-1}}$')
        # im.set_data(dic_of_images[datetimes[indexvalue]])
        update(sfreq.val)
        fig.canvas.draw_idle()
    
    def reset(event):
        sfreq.reset()
    def update(val):
        qi0 = int(sfreq.val)
        tsraw, qsraw, ddmraw, fraw = datasetsPush[0]
        q_value = qsraw[qi0]
        t1.set_text(fr'$q$: {q_value/1e6:.01f} μm$^{{-1}}$')
        
        lautoraw.set_xdata(qsraw[qi0]**2*tsraw)
        lautoraw.set_ydata(fraw[:,qi0])
        lddmraw.set_xdata(tsraw)
        lddmraw.set_ydata(ddmraw[:,qi0])
        
        dtsminlim.set_xdata(q_value**2*tsraw[dtsminpush])
        dtsmaxlim.set_xdata(q_value**2*tsraw[dtsmaxpush])
        
        if len(datasetsPush) > 1:
            tsest, qsest, ddmest, f_fit, A, B = datasetsPush[1]
            qi1 = np.argmin(np.abs(qsraw[qi0]-qsest))
            if qsraw[qi0] != qsest[qi1]:
                lauto_fit.set_xdata([])
                lauto_fit.set_ydata([])
                lautoest.set_xdata([])
                lautoest.set_ydata([])
                lddmest.set_xdata([])
                lddmest.set_ydata([])
            else:
                lauto_fit.set_xdata(qsest[qi1]**2*tsest)
                lauto_fit.set_ydata(f_fit[:,qi1])
                lautoest.set_xdata(qsest[qi1]**2*tsest)
                lautoest.set_ydata(fest[:,qi1])
                lddmest.set_xdata(tsest)
                lddmest.set_ydata(ddmest[:,qi1])
            
        fig.canvas.draw_idle()
        axs[1].set_ylim( [np.min(ddmraw[:,qi0]), np.max(ddmraw[:,qi0])] )
        axs[0].set_xlim((0.9*qsraw[qi0]**2*np.min(tsraw), 1.1*qsraw[qi0]**2*np.max(tsraw)))    
        t1.set_text(fr'$q$: {qsraw[qi0]/1e6:.01f} μm$^{{-1}}$')
    fig.canvas.mpl_connect('key_press_event', press)
    fig.canvas.mpl_connect('scroll_event', press)
    sfreq.on_changed(update)
    plt.show(block=False)



def plotCONTIN(sol):
    alphas    = sol.alphas
    alpha     = sol.chosen_alpha
    ini       = alphas.index(alpha)
    amplitude = sol.alpha_amplitude[ini]
    noise     = sol.alpha_noise[ini]
    f_fit     = sol.alpha_ddmfit[ini]
    g         = sol.alpha_g[ini]
    tau       = sol.tau
    ddmdata   = sol.ddmdata
    s         = sol.gamma_range
    
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(15,7))
    axcolor   = 'lightgoldenrodyellow'
    axalpha   = plt.axes([0.35, 0.9, 0.3, 0.03], facecolor=axcolor)
    salpha    = Slider(axalpha, '', 0, len(alphas)-1, valinit=ini, valstep=1)
    fitplot,  = ax1.semilogx(tau, amplitude * (1-f_fit) + noise, 
                           label=r'CONTIN')
    dataplot, = ax1.semilogx(tau, ddmdata, '.', label='data')
    ax1.set_xlabel(r"$\tau q^2$ [s/m²]")
    g         = g / np.nansum(g)
    if sol.sizes is None:
        distplot, = ax2.plot(s, g, label=r"$\Gamma$ distribution per CONTIN")
        ax2.set_xlabel(r"$\Gamma$ [m²/s]")
    else:
        distplot, = ax2.plot(sol.sizes*1e9, g, label=r"$R_h$ distribution per CONTIN")
        ax2.set_xlabel(r"$R_h$ [nm]")
    ax2.set_ylabel("Intensity")
    ax2.set_ylim((0,1.5*np.max(g)))
    
    
    ax1.legend()
    ax2.legend()
    t1 = plt.text(0.5, 1.2, fr'$\alpha$: {alpha:.02e}')
    def press(event):
        try:
            button = event.button
        except:
            button = 'None'
        if event.key == 'right' or button == 'down':
            if salpha.val < len(alphas) - 1:
                salpha.set_val(salpha.val + 1)
        elif event.key == 'left' or button == 'up':
            if salpha.val > 0:
                salpha.set_val(salpha.val - 1)
        # im.set_data(dic_of_images[datetimes[indexvalue]])
        update(salpha.val)
        fig.canvas.draw_idle()
    
    def reset(event):
        salpha.reset()
    def update(val):
        val = int(val)
        alpha = alphas[val]
        t1.set_text(fr'$\alpha$: {alpha:.02e}')
        amplitude = sol.alpha_amplitude[val]
        noise     = sol.alpha_noise[val]
        f         = sol.alpha_ddmfit[val]
        fitplot.set_ydata(amplitude * (1 - f) + noise)
        g         = sol.alpha_g[val] / np.nansum(sol.alpha_g[val])
        distplot.set_ydata(g)
        fig.canvas.draw_idle()
    fig.canvas.mpl_connect('key_press_event', press)
    fig.canvas.mpl_connect('scroll_event', press)
    salpha.on_changed(update)
    plt.tight_layout(rect=(0,0,1,0.93))
    plt.show(block=False)
    
    
    
    
    
    
    


if __name__ == "__main__":
    ######## test the main plot function
    """ (comment this line to activate)
    path = "/home/fred/Nextcloud/ddm_matrices/Merged__DDM_matrix.npy"
    path = "/home/frederic/Desktop/ddm_matrices/20°C-531fps_-01_DDM_matrix.npy"
    dts, qs, ddm = loadDDM(path)
    plotDDMmatrix([[dts, qs, ddm, ddm]])
    #"""
