#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A small GUI program to interface a custom DDM setup.

@author: Frédéric Dux, biosoft intern@IPC with Jerome Crassous
"""
from          os.path         import  basename, exists, dirname, join
from          glob            import  glob
import        numpy           as      np
import        threading
import        re
import        traceback



from          utilities       import  loadAnalyzedVideos, loadDirectory,\
                                      water_viscosity, saveFitTextFile, \
                                      saveMatrixCSV, musthaves, saveAutocorrelationCSV,\
                                      renameTimestamp, saveCONTINfit
from          generateDDM     import  FFTStack, timeDependantDDM, concatenateVideos
from          plotDDM         import  plotDDMmatrix, plotDDMmeasAndFit, \
                                      plotABandR, plotFitParams, plotSomefs,\
                                      plotCONTIN
from          fitDDM          import  fitOneDDMmatrix, getRadius, mergeDDM
from          contin          import  CONTIN
from          layoutFile      import  DEFAULTFONT, layout, DEFAULTDATADIRECTORY, sg
from          namedConcepts   import  FITMODELS, FITPARAMNAMES, FITPARAMDEFAULTS



FITPARAMGUESS = {}
for key, value in FITPARAMDEFAULTS.items():
    fixed = [False for v in value]
    FITPARAMGUESS[key] = value, fixed



##### generate a layout for the input of the guess fit params:
winparamsactive = False # by default, that second window won't be active
def getInitialGuessParamsInputWindowLayout(fitmodel):
    currentguess, fixed = FITPARAMGUESS[fitmodel]
    values = [v if not v == '' else e for v, e in zip(currentguess, FITPARAMDEFAULTS[fitmodel])]
    display_names = FITPARAMNAMES[fitmodel]
    layout = [[sg.Text("Initial guess for the fitting solver", font=DEFAULTFONT)]]
    nametitle = 'parameters'
    valuetitle = 'values'
    layout.append([sg.Text(f"{nametitle}", font=DEFAULTFONT, size=(25,1)),\
                   sg.Text(f"{valuetitle}", font=DEFAULTFONT, size=(10,1)),\
                   sg.Text("fixed?")])
    for i, (name, value, check) in enumerate(zip(display_names, values, fixed)):
        display_element = sg.Text(f"{name}", font=DEFAULTFONT, size=(25,1))
        input_element   = sg.InputText(value, font=DEFAULTFONT, key=i, size=(10,1))
        fixed_checkbox  = sg.Checkbox('', default=check, key=f"{i}_fixed")
        layout.append([display_element, input_element, fixed_checkbox])
    layout.append([ sg.Button('submit', font=DEFAULTFONT, key='submitinitialguess',\
                             enable_events=True)])
    return layout

winCONTINactive = False
CONTINsolutions = {}
def getCONTINWindowLayout(q_index):
    layout = [[sg.Text("CONTIN", font=DEFAULTFONT)]]
    layout.append([sg.Text("q index", font=DEFAULTFONT, size=(35,1)),  
                       sg.InputText(f"{q_index:.0f}", size=(9,1), justification='right', key='continqindex', enable_events=True,\
                       tooltip="CONTIN is expensive and will act on a single wavenumber",\
                       font=DEFAULTFONT)])
    layout.append([sg.Text("Smallest decay rate", font=DEFAULTFONT, size=(35,1)),  
                       sg.InputText('1e-13', size=(9,1), justification='right', key='continmindecay', enable_events=True,\
                       tooltip="should be smaller than the smallest expected decay rate in the distribution",\
                       font=DEFAULTFONT)])
    layout.append([sg.Text("Biggest decay rate", font=DEFAULTFONT, size=(35,1)),  
                       sg.InputText('1e-11', size=(9,1), justification='right', key='continmaxdecay', enable_events=True,\
                       tooltip="should be bigger than the biggest expected decay rate in the distribution",\
                       font=DEFAULTFONT)])
    layout.append([sg.Text("N decay rate", font=DEFAULTFONT, size=(35,1)),  
                       sg.InputText('10', size=(9,1), justification='right', key='continNdecay', enable_events=True,\
                       tooltip="The number of decay rates used to build the total decay",\
                       font=DEFAULTFONT)])
    layout.append([sg.Text("Smallest alpha", font=DEFAULTFONT, size=(35,1)),  
                       sg.InputText('0.01', size=(9,1), justification='right', key='continminalpha', enable_events=True,\
                       tooltip="alpha is the regularization parameter: will try each from smallest to biggest in N steps",\
                       font=DEFAULTFONT)])
    layout.append([sg.Text("Biggest alpha", font=DEFAULTFONT, size=(35,1)),  
                       sg.InputText('1', size=(9,1), justification='right', key='continmaxalpha', enable_events=True,\
                       tooltip="alpha is the regularization parameter: will try each from smallest to biggest in N steps",\
                       font=DEFAULTFONT)])
    layout.append([sg.Text("N alpha", font=DEFAULTFONT, size=(35,1)),  
                       sg.InputText('3', size=(9,1), justification='right', key='continNalpha', enable_events=True,\
                       tooltip="alpha is the regularization parameter: will try each from smallest to biggest in N steps",\
                       font=DEFAULTFONT)])
    layout.append([sg.Text("max iter", font=DEFAULTFONT, size=(35,1)),  
                       sg.InputText('10', size=(9,1), justification='right', key='continmaxiter', enable_events=True,\
                       tooltip="max number of iterations when optimizing the distribution",\
                       font=DEFAULTFONT)])
    layout.append([ sg.Button('Estimate', font=DEFAULTFONT, key='goCONTIN',\
                             enable_events=True),
                    sg.Button('Plot', font=DEFAULTFONT, key='plotCONTIN',\
                             enable_events=True),
                    sg.Button('Save', font=DEFAULTFONT, key='saveCONTIN',\
                             enable_events=True),
                    sg.Checkbox('save all the fits?', 
                                default=0, key="continsaveall")])
    return layout


def addOrRemoveSizesCONTINSolution(solution):
    # if temperature and viscosity entered, carry the sizes too:
    temperature, viscosity = readTempAndVisco()
    if temperature and viscosity:
        sizes = getRadius(solution.gamma_range, viscosity, temperature)
        solution.sizes = sizes
    else:
        solution.sizes = None

def runCONTIN(continwindow, videotofit):
    q_index  = int(float(continwindow.Element('continqindex').Get()))
    alphamin = float(continwindow.Element('continminalpha').Get())
    alphamax = float(continwindow.Element('continmaxalpha').Get())
    alphaN   = int(float(continwindow.Element('continNalpha').Get()))
    gammamin = float(continwindow.Element('continmindecay').Get())
    gammamax = float(continwindow.Element('continmaxdecay').Get())
    gammaN   = int(float(continwindow.Element('continNdecay').Get()))
    maxiter  = int(float(continwindow.Element('continmaxiter').Get()))
    
    alphas   = np.linspace(alphamin, alphamax, alphaN)
    gammas   = np.linspace(gammamin, gammamax, gammaN)
    dtmin, dtmax     = int(value['dtsminslider']), int(value['dtsmaxslider'])
    matrix, time, qs = computeddata[videotofit]
    ddmdata = matrix[dtmin:dtmax, q_index]
    q       = qs[q_index]
    tau     = time[dtmin:dtmax] * q**2
    state.Update(f"Fitting the matrix {basename(videotofit)} \n using CONTIN.")
    window.Refresh()
    run     = CONTIN(tau, ddmdata, gammas, alpha=alphas, maxiter=maxiter)
    for i, e in enumerate(run):
        pb.UpdateBar(i/alphaN*pb.MaxValue)
    state.Update("Idle")
    window.Refresh()
    e.q_index = q_index
    # if temperature and viscosity entered, carry the sizes too:
    addOrRemoveSizesCONTINSolution(e)
    return e
    
def readAllGuessParamsFromInputWindow(window, fitmodel):
    values = []
    fixed  = []
    # all except A and B:
    for i, paramname in enumerate(FITPARAMNAMES[fitmodel]):
        try:
            try:
                values.append(float(window.Element(i).Get()))
            except:
                values.append('')
            fixed.append(bool(window.Element(f"{i}_fixed").Get()))
        except:
            sg.Popup(f"Problem with reading the value of the param {paramname}")
    values[:-2] = [v if not v == '' else e for v, e in zip(values[:-2], FITPARAMDEFAULTS[fitmodel][:-2])]
    return values, fixed

def fitCurrentMatrixWithCurrentFitModel(videotofit):
    fitmethod    = value['fitmodel']
    fitmethod    = FITMODELS[fitmethod]
    state.Update(f"Fitting the matrix {basename(videotofit)} \n using {fitmethod}.")
    window.Refresh()
    qmin, qmax   = int(value['qminslider']), int(value['qmaxslider'])
    dtmin, dtmax = int(value['dtsminslider']), int(value['dtsmaxslider'])
    A, B, cumulants, f, opt = fitOneDDMmatrix(computeddata[videotofit], 
                                              model=fitmethod, 
                                              ini=FITPARAMGUESS[fitmethod][0],
                                              fixed=FITPARAMGUESS[fitmethod][1],
                                              qmin=qmin, qmax=qmax,
                                              dtmin=dtmin, dtmax=dtmax)                                           
    _, dtsraw, qsraw = computeddata[videotofit]
    ddmest = A*(1-f)+B
    dtsest = dtsraw
    qsest  = qsraw[qmin:qmax]
    fitteddata[videotofit] = [dtsest, qsest, ddmest, A, B, cumulants, f, fitmethod]
    state.Update("Idle")
    window.Refresh()


def popupCheckboxListOfItems(listOfItems, default=True, title="", usebasename=True):
    filter_ = sg.PopupGetText("Filter to apply? (portions of text the matrices must have)")
    filter_not = sg.PopupGetText("Filter to apply? (portions of the text the matrices must not have)")
    listOfItems = [e for e in listOfItems if filter_ in str(e)]
    if filter_not :
        listOfItems = [e for e in listOfItems if (not filter_not in str(e))]
    layout = []
    layout.append( [sg.Button("Submit", key="Submit", enable_events=True)] )
    for item in listOfItems:
        if usebasename:
            displayitem = basename(item)
        else:
            displayitem = str(item)
        layout.append([sg.Checkbox(displayitem, default=default, key=item,\
                                   enable_events=True)])
    choicepopup  = sg.Window(title=title, layout=layout)
    event, value = choicepopup.Read(timeout=1)
    usefulvalues = value
    while True:
        event, value = choicepopup.Read()
        if event == "Submit":
            break
        elif event == "filter":
            pass
        elif event == None:
            return []
        else:
            usefulvalues = value
    choicepopup.Close()
    return usefulvalues

def readTempAndVisco():
    try:
        temperature, viscosity = float(tempinput.Get())+273.1, viscoinput.Get()
        if viscosity == 'water':
            viscosity = water_viscosity(temperature)
        else:
            viscosity = float(viscosity)
    except:
        temperature, viscosity = None, None
    return temperature, viscosity

def renameTimestampGUI():
    pathToVideos = sg.PopupGetFolder("Please select the directory where to rename the ugly timestamps")
    if pathToVideos:
        renameTimestamp(pathToVideos)

def concatenateVideosGUI():
    pathToVideos = sg.PopupGetFolder("Please select the directory containing the videos to be concatenated")
    if not pathToVideos:
        return
    videos       = glob(join(pathToVideos, '*.avi'))
    values       = popupCheckboxListOfItems(videos, title="List of videos to be appended", \
                                            usebasename=True)
    if not values:
        return
    toappend     = [e for e in values if values[e] == True]
    pathToNew    = sg.PopupGetFile("Also, choose a location and name for the concatenated video.",\
                                   save_as=True)
    if not pathToNew:
        sg.Popup("Concatenation aborted: no video name provided.")
        return
    concatenateVideos(toappend, pathToNew)
    
def mergeDDMGUI(computeddata, event, value):
    global actuallist
    actuallist = popupCheckboxListOfItems(list(computeddata.keys()), title=f"Matrices to {event.replace('button', '')}")
    if not actuallist:
        computeddata, displaycomputeddata = loadAnalyzedVideos(path)
        try:
            videoathand = displaycomputeddata[value['computedlist']]
        except:
            event, value = window.Read(timeout=1)
            videoathand = displaycomputeddata[value['computedlist']]
        return computeddata, displaycomputeddata, videoathand
    actualdata = {e:computeddata[e] for e in actuallist if actuallist[e] == True}
    name       = sg.PopupGetText("Enter a name for the merged matrix")
    if not name:
        name = ''
    mergeDDM(actualdata, mode=event.replace('button', ''), title=name)
    computeddata, displaycomputeddata = loadAnalyzedVideos(path)
    computedlist.Update(values=[key for key in displaycomputeddata.keys()])
    try:
        videoathand = displaycomputeddata[value['computedlist']]
    except:
        event, value = window.Read(timeout=1)
        videoathand = displaycomputeddata[value['computedlist']]
    return computeddata, displaycomputeddata, videoathand    


# Show the Window to the user
window = sg.Window('DDMSoft', layout, return_keyboard_events=True, icon='logo.png')
# shortcuts to the important elements:
pb            = window.FindElement('progressbar')
infow         = window.FindElement('videodescriptions')
state         = window.FindElement('state')
recompute     = window.FindElement('recompute')
maxcouples    = window.FindElement('maxcouples')
ptperdecade   = window.FindElement('ptperdecade')
computedlist  = window.FindElement('computedlist')
fitmodel      = window.FindElement('fitmodel')
mergebutton   = window.FindElement('mergebutton')
fitbutton     = window.FindElement('fitbutton')
qminslider    = window.FindElement('qminslider')
qmaxslider    = window.FindElement('qmaxslider')
dtsminslider  = window.FindElement('dtsminslider')
dtsmaxslider  = window.FindElement('dtsmaxslider')
noiseandamp   = window.FindElement('noiseandamp')
viscoinput    = window.FindElement('viscosity')
tempinput     = window.FindElement('temperature')
nangle        = window.FindElement('Nangle')



params     = {}
### some variables related to the fitting part:
fitteddata = {}
temperature, viscosity = None, None

videoathand = None  # no video selected by default.
last_path   = DEFAULTDATADIRECTORY

# Event loop. Read buttons, make callbacks
while True:
    # Read the Window
    event, value = window.Read(timeout=100)
   
    if event in ('Quit', 'Exit', None):
        break
    try:
        ##### grab the parameters
        try:
            recomputebool = bool(recompute.Get())
            maxCouples    = int(maxcouples.Get())
            ptPerDecade   = int(ptperdecade.Get())
            Nangle        = int(nangle.Get())
        except:
            pass
        
        ##### load everything that's to be loaded
        if  event == 'Open directory' or event == "" or event == "o:79" or event == "o:32":
            path   = sg.PopupGetFolder("Please choose the directory containing your videos", \
                                       default_path=last_path)
            if not path:
                continue
            if not exists(path):
                sg.Popup("Error", f"This directory ({path}) does not exist!")
                continue
            params = loadDirectory(path)
            info  = f'{"name":<46}{"frame rate (Hz)":<18}{"pixel size (meter)":<18}\n'
            for vid in params.keys():
                p     = params[vid]
                info += f"\n{basename(vid):<46}{p['framerate']:<18}{p['pixelsize']:<18}"
            infow.Update(info)
            # also loading the already computed matrices:
            computeddata, displaycomputeddata = loadAnalyzedVideos(path)
            if len(computeddata)>0:
                computedlist.Update(values=[key for key in displaycomputeddata.keys()],
                                    set_to_index=0)
                try:
                    videoathand = displaycomputeddata[value['computedlist']]
                except:
                    event, value = window.Read(timeout=1)
                    videoathand = displaycomputeddata[value['computedlist']]
                
                qsathand     = computeddata[videoathand][2]
                dtsathand    = computeddata[videoathand][1]
                _maxqval     = len(qsathand)-1
                _dtsmaxval   = len(dtsathand)-1
                try:
                    _qmin    = qmin
                    _qmax    = qmax
                    _dtmin   = dtmin
                    _dtmax   = dtmax
                except:
                    _qmin    = 0
                    _qmax    = _maxqval
                    _dtmin   = 0
                    _dtmax   = _dtsmaxval
                qminslider.Update(range=(0,_maxqval), value=_qmin)
                qmaxslider.Update(range=(0,_maxqval), value=_qmax)
                dtsminslider.Update(range=(0, _dtsmaxval), value=_dtmin)
                dtsmaxslider.Update(range=(0, _dtsmaxval), value=_dtmax)

                last_path             = path
            
        ##### if the video parameters are changed, take that into account
        elif event == 'videodescriptions':
            news  = infow.Get().split('\n')[2:]
            news  = [new for new in news if new != '']
            _keys = []
            try:
                for vid, new in zip(list(params.keys()), news):
                    new = re.sub(' +', ' ', new).strip()
                    vid_written, framerate, pixelsize  = new.split(' ')
                    _keys.append(vid_written)
                    assert(basename(vid)==vid_written)
                    params[vid]['framerate'] = framerate
                    params[vid]['pixelsize'] = pixelsize
            except AssertionError :
                sg.Popup("You shoud not change the name of the video here, it will end badly.")
            except:
                pass
                

        ##### main course, do the crunching of the listed videos
        elif event == 'Process':
            Nvids  =   len(params)
            i      =   0
            for vid in params.keys():
                try:
                    freq           = float(params[vid]['framerate'])
                except ValueError:
                    sg.Popup(f"The value '{params[vid]['framerate']}' can not be interpreted as a frame rate.",\
                               title="Bad parameters")
                    break
                try:
                    pixelsize      = float(params[vid]['pixelsize'])
                except ValueError:
                    sg.Popup(f"The value '{params[vid]['pixelsize']}' can not be interpreted as a pixel size.",\
                               title="Bad parameters")
                    break
                if (not recomputebool) and \
                   (not basename(vid.replace('.avi', '_DDM_matrix.npy')) in displaycomputeddata.keys() ):
                    goahead = True
                elif recomputebool:
                    goahead = True
                else :
                    goahead = False
                if goahead:
                    fftstack       = FFTStack(freq, pixelsize, maxCouples, ptPerDecade, 
                                              Nangle=Nangle)
                    i += 1
                    state.Update(f"Processing video {i} out of {Nvids}: loading")
                    fftstack.loadVideo(vid)
                    print(f'length of the stack: {fftstack.Nbimages}')
                    MaxValue       = fftstack.getTotalNumberOfOperations()
                    pbmax          = pb.MaxValue
                    workThread     = threading.Thread(target=fftstack.fftVideo)
                    workThread.start()
                    state.Update(f"Processing video {i} out of {Nvids}: FFT")                
                    while not fftstack.isFFTdone():
                        pb.UpdateBar(fftstack.getProgress()/MaxValue*pbmax)
                        event, values = window.Read(timeout=50)
                        if event == 'Quit':
                            workThread.join()
                            break
                    workThread     = threading.Thread(target=fftstack.stackToDDM)
                    workThread.start()
                    state.Update(f"Processing video {i} out of {Nvids}: time averaging")
                    while not fftstack.isCompleted():
                        pb.UpdateBar(fftstack.getProgress()/MaxValue*pbmax)
                        event, values = window.Read(timeout=50)
                        if event == 'Quit':
                            workThread.join()
                            break
                    pb.UpdateBar(pbmax)
                    # at the end, update the list of computed videos:
                    computeddata, displaycomputeddata = loadAnalyzedVideos(path)
                    computedlist.Update(values=[key for key in displaycomputeddata.keys()],
                                        set_to_index=0)
                    if len(computeddata) > 0:
                        try:
                            videoathand = displaycomputeddata[value['computedlist']]
                        except:
                            event, value = window.Read(timeout=1)
                            videoathand = displaycomputeddata[value['computedlist']]
                        qsathand    = computeddata[videoathand][2]
                        dtsathand   = computeddata[videoathand][1]
                        qminslider.Update(range=(0,len(qsathand)-1))
                        qmaxslider.Update(range=(0,len(qsathand)-1))
                        dtsminslider.Update(range=(0, len(dtsathand)-1))
                        dtsmaxslider.Update(range=(0, len(dtsathand)-1))
                        state.Update('Idle')
                else:
                    print(f"Video {basename(vid)} already processed.")
                    state.Update(f"Video {basename(vid)} already processed.")
        
        ##### change the focus of the fitting routines
        elif event == 'computedlist':
            try:
                videoathand = displaycomputeddata[value['computedlist']]
            except:
                event, value = window.Read(timeout=1)
                videoathand = displaycomputeddata[value['computedlist']]
            qsathand    = computeddata[videoathand][2]
            dtsathand   = computeddata[videoathand][1]
            qminslider.Update(range=(0, len(qsathand)-1))
            qmaxslider.Update(range=(0, len(qsathand)-1))
            dtsminslider.Update(range=(0, len(dtsathand)-1))
            dtsmaxslider.Update(range=(0, len(dtsathand)-1))
            dtsmaxslider.Update(value=len(dtsathand)-1)
            if int(value['dtsmaxslider']) > len(dtsathand)-1:
                dtsmaxslider.Update(range=(0, len(dtsathand)-1), value=len(dtsathand)-1)
            else:
                dtsmaxslider.Update(range=(0, len(dtsathand)-1), value=int(value['dtsmaxslider']))
            if int(value['dtsminslider']) > len(dtsathand)-1:
                dtsminslider.Update(range=(0, len(dtsathand)-1), value=0)
            else:
                dtsminslider.Update(range=(0, len(dtsathand)-1), value=int(value['dtsminslider']))
            window.Finalize()
        
        ##### in case merging is desired, merge all the videos in the list together 
        ##### (except those that are themselves a result of a merging)
        elif event == 'mergebutton' or event == 'averagebutton':
            computeddata, displaycomputeddata, videoathand = mergeDDMGUI(computeddata, event, value)

            
        ##### plotting the main graph:
        elif event == "mainplotbutton":
            ddmraw, dtsraw, qsraw = computeddata[videoathand]
            try:
                dtsest, qsest, ddmest, A, B, _, f, plotfitmethod = fitteddata[videoathand]
                plotDDMmatrix([[dtsraw, qsraw, ddmraw, None], [dtsest, qsest, ddmest, f, A, B]],\
                              title=basename(videoathand), dtsmin=int(value['dtsminslider']),\
                              dtsmax=int(value['dtsmaxslider']))
                
            except:
                print(traceback.format_exc())
                plotDDMmatrix([[dtsraw, qsraw, ddmraw, None]], title=basename(videoathand),\
                              dtsmin=int(value['dtsminslider']), dtsmax=int(value['dtsmaxslider']))
                
        elif event == "Show the DDM matrix (image)":
            ddmraw, dtsraw, qsraw = computeddata[videoathand]
            try:
                dtsest, qsest, ddmest, _, _, _, f, plotfitmethod = fitteddata[videoathand]
                plotDDMmeasAndFit(qsraw, dtsraw, ddmraw, qsfit=qsest, dtsfit=dtsest, ddmfit=ddmest,\
                                  title=basename(videoathand))
            except:
                print(traceback.format_exc())
                plotDDMmeasAndFit(qsraw, dtsraw, ddmraw, title=basename(videoathand))
        ##### plotting A, B and (either) Rh or D
        elif event == "noiseandamp":
            _, qsest, _, A, B, cumulants, _, plotfitmethod = fitteddata[videoathand]
            temperature, viscosity = readTempAndVisco()
            plotABandR(qsest, A, B, cumulants[0], temperature=temperature, viscosity=viscosity,\
                       title=basename(videoathand))
            
        ##### fitting
        elif event == "CONTIN":
            if not winCONTINactive:
                winCONTIN = sg.Window("CONTIN", getCONTINWindowLayout(value['qminslider']))
                winCONTINactive = True
            else:
                try:
                    winCONTIN.BringToFront()
                except Exception as e:
                    print(e)
                    
        elif event == "winparamsfit":
            if not winparamsactive:
                fitmethod = value['fitmodel']
                fitmethod = FITMODELS[fitmethod]
                winparams = sg.Window("Guess for the solver", getInitialGuessParamsInputWindowLayout(fitmethod))
                winparamsactive = True
            else:
                try:
                    winparams.BringToFront()
                except Exception as e:
                    print(e)
                
        elif event == "fitbutton":
            fitCurrentMatrixWithCurrentFitModel(videoathand)
            qminfit = value['qminslider']
            qmaxfit = value['qmaxslider']
            
        elif event == "Fit and save all the matrices":
            fitcsvsavepathall = sg.PopupGetFile("Please choose a directory and a prefix for your fit",save_as=True)
            if not fitcsvsavepathall:
                continue
            for i_, videotofit in enumerate(computeddata):
                pb.UpdateBar(i_, len(computeddata))
                fitCurrentMatrixWithCurrentFitModel(videotofit)
                dtsest, qsest, ddmest, A, B, cumulants, f, plotfitmethod = fitteddata[videotofit]
                currentname    = basename(videotofit).replace(musthaves[0], '')
                prefix         = basename(fitcsvsavepathall).replace('.txt','').replace('.csv','')
                savedir        = dirname(fitcsvsavepathall)
                fitcsvsavepath = join(savedir, prefix+'_'+currentname)
                temperature, viscosity = readTempAndVisco()
                saveFitTextFile(fitcsvsavepath, qsest, A, B, cumulants, FITPARAMNAMES[plotfitmethod],\
                                temperature=temperature, viscosity=viscosity)
                
            pb.UpdateBar(1, 1)
            window.Refresh()
            
            
        elif event == "fitmodel":
            # make sure the init param windows is closed:
            if winparamsactive:
                winparams.Close()
                winparamsactive = False
        
        
        ##### watch the sliders:
        elif event == 'dtsminslider':
            dtsmin, dtsmax = value['dtsminslider'], value['dtsmaxslider']
            if dtsmin >= dtsmax:
                dtsmax = dtsmin + 1
                dtsmaxslider.Update(value=dtsmax)
        elif event == 'dtsmaxslider':
            dtsmin, dtsmax = value['dtsminslider'], value['dtsmaxslider']
            if dtsmin >= dtsmax:
                dtsmin = dtsmax - 1
                dtsminslider.Update(value=dtsmin)
        elif event == 'qminslider':
            qmin, qmax = value['qminslider'], value['qmaxslider']
            if qmin >= qmax:
                qmax = qmin + 1
                qmaxslider.Update(value=qmax)
        elif event == 'qmaxslider':
            qmin, qmax = value['qminslider'], value['qmaxslider']
            if qmin >= qmax:
                qmin = qmax-1
                qminslider.Update(value=qmin)
        
        ##### plot the fit params
        elif event == 'showfitparams':
            dtsest, qsest, ddmest, A, B, cumulants, f, plotfitmethod = fitteddata[videoathand]
            plotFitParams(qsest, cumulants, FITPARAMNAMES[plotfitmethod], title=basename(videoathand))
        ##### try to make a nice plot of some autocorrelation functions and their fits:
        elif event == "Plot some autocorrelation functions":
            ddmraw, dtsraw, qsraw = computeddata[videoathand]
            try:
                dtsest, qsest, ddmest, A, B, _, f, plotfitmethod = fitteddata[videoathand]
                plotSomefs(qmin=qminfit, qmax=qmaxfit,\
                           datasets=[[dtsraw, qsraw, ddmraw, None], [dtsest, qsest, ddmest, f, A, B]],\
                           title=basename(videoathand))
            except:
                print(traceback.format_exc())
                plotSomefs(qmin=value['qminslider'], qmax=value['qmaxslider'],\
                           datasets=[[dtsraw, qsraw, ddmraw, None]], \
                           title=basename(videoathand))
            
        ##### time dependant part, we will load the video, split it and
        ##### push the results into computeddata (with the appropriate keys)
        ##### like any other dataset. 
        elif event == "Split a video in N matrices (time dependant DDM)":
            if len(params)>1:
                vidstosplit = popupCheckboxListOfItems(list(params.keys()), default=False, title="Select the videos to be split:")
            else:
                vidstosplit = {next(iter(params)) : True}
            _success = False
            while not _success:
                try:
                    Npartitions    = int(sg.PopupGetText("In how many parts should the video(s) be split?"))
                    _success = True
                except:
                    _success = False
            for vidtosplit in [vidtosplit for vidtosplit in vidstosplit if vidstosplit[vidtosplit] == True]:
                freq           = float(params[vid]['framerate'])
                pixelsize      = float(params[vid]['pixelsize'])
                
               
                timedep = timeDependantDDM(freq, pixelsize, maxCouples, ptPerDecade, Npartitions)
                state.Update("Loading the video in the memory")
                loadthread = threading.Thread(target=timedep.loadVideo, args=[vidtosplit])
                loadthread.start()
                while not timedep.getStatus()[0]:
                    event, values = window.Read(timeout=200)
                    if event == 'Quit':
                        loadthread.join()
                        break
                state.Update("Computing the FFTs for all the subvideos")
                fftthread = threading.Thread(target=timedep.fftAllStacks)
                fftthread.start()
                while not timedep.getStatus()[1]:
                    event, values = window.Read(timeout=200)
                    if event == 'Quit':
                        fftthread.join()
                        break
                state.Update("Averaging in parallel")
                workthread = threading.Thread(target=timedep.ddmAllStacks)
                workthread.start()
                while not timedep.getStatus()[2]:
                    event, values = window.Read(timeout=200)
                    if event == 'Quit':
                        workthread.join()
                        break
                computeddata, displaycomputeddata = loadAnalyzedVideos(path)
                computedlist.Update(values=[key for key in displaycomputeddata.keys()],
                                    set_to_index=0)
                if len(computeddata) > 0:
                    videoathand = displaycomputeddata[value['computedlist']]
                    qsathand    = computeddata[videoathand][2]
                    dtsathand   = computeddata[videoathand][1]
                    qminslider.Update(range=(0,len(qsathand)-1))
                    qmaxslider.Update(range=(0,len(qsathand)-1))
                    dtsminslider.Update(range=(0, len(dtsathand)-1))
                    dtsmaxslider.Update(range=(0, len(dtsathand)-1))
                state.Update('Idle')
        
        elif event == "Save the DDM matrix as a text file":
            ddmcsvsavepath = sg.PopupGetFile("Please choose a directory and a title.",save_as=True)
            if not ddmcsvsavepath:
                continue
            saveMatrixCSV(ddmcsvsavepath, computeddata[videoathand])
        elif event == "Save the current fit parameters":
            fitcsvsavepath = sg.PopupGetFile("Please choose a directory and a title for your fit",
                                             save_as=True)
            if not fitcsvsavepath:
                continue
            dtsest, qsest, ddmest, A, B, cumulants, f, plotfitmethod = fitteddata[videoathand]
            temperature, viscosity = readTempAndVisco()
            saveFitTextFile(fitcsvsavepath, qsest, A, B, cumulants, FITPARAMNAMES[plotfitmethod],\
                            temperature=temperature, viscosity=viscosity)
        
        elif event == "Save the correlation functions":
            autocorcsvsavepath = sg.PopupGetFile("Please choose a directory and a title for your refined autocorrelation function",\
                                                 save_as=True)
            if not autocorcsvsavepath:
                continue
            dtsest, qsest, ddmest, A, B, cumulants, f, plotfitmethod = fitteddata[videoathand]
            saveAutocorrelationCSV(autocorcsvsavepath, computeddata[videoathand], A, B, qsest)
            
        elif event == "Concatenate videos":
            toastytoast = concatenateVideosGUI()
        elif event == "Rename the timestamps in a directory":
            renameTimestampGUI()

        elif event == 'About...':
            textcredits = """
                DDMSoft
                     
 Concept: Jérôme Crassous and Frédéric Dux
 Programming: Frédéric Dux
 Bugs and feature requests: duxfrederic@gmail.com
 """
            sg.Popup(textcredits, title="Credits")

        
        ##### now deal with the secondary windows
        if winparamsactive:
            evparams, val_params = winparams.Read(timeout=50)
            if not evparams == '__TIMEOUT__':
                print('second window event:', evparams)
            if evparams == 'submitinitialguess':
                FITPARAMGUESS[fitmethod] = readAllGuessParamsFromInputWindow(winparams, fitmethod)
            elif evparams == '_Close_' or evparams == None:
                winparamsactive  = False
                winparams.Close()
        if winCONTINactive:
            evparams, val_params = winCONTIN.Read(timeout=50)
            if not evparams == '__TIMEOUT__':
                print('CONTIN window event:', evparams)
            if evparams == "goCONTIN":
                try:
                    solution = runCONTIN(winCONTIN, videoathand)        
                    q_index = solution.q_index
                    q = computeddata[videoathand][2][q_index]
                    if videoathand in CONTINsolutions:
                        CONTINsolutions[videoathand][q] = solution
                    else:
                        CONTINsolutions[videoathand] = {}
                        CONTINsolutions[videoathand][q] = solution
                    plotCONTIN(solution)
                except Exception as e:
                    print(e)
                    sg.Popup("No data found", title="error")
            elif evparams == "plotCONTIN":
                q_index = int(float(winCONTIN.Element('continqindex').Get()))
                q = computeddata[videoathand][2][q_index]
                if videoathand in CONTINsolutions and q in CONTINsolutions[videoathand]:
                    solution = CONTINsolutions[videoathand][q]
                    addOrRemoveSizesCONTINSolution(solution)
                    plotCONTIN(solution)
                else:
                    sg.Popup("No CONTIN data yet", title='error')
            elif evparams == "saveCONTIN":
                q_index = int(float(winCONTIN.Element('continqindex').Get()))
                q = computeddata[videoathand][2][q_index]
                for vid in CONTINsolutions:
                    for q in CONTINsolutions[vid]:
                        addOrRemoveSizesCONTINSolution(CONTINsolutions[vid][q])
                savecontinpath = sg.PopupGetFile('Choose a file', save_as=True)
                saveCONTINfit(savecontinpath, CONTINsolutions, winCONTIN, videoathand, q)
            
            elif evparams == '_Close_' or evparams == None:
                winCONTINactive  = False
                winCONTIN.Close()            
                
        
    except Exception as e:
        sg.Popup('Error:', str(e), traceback.format_exc())
        pass        
        

window.Close()
if winparamsactive:
    winparamsactive  = False
    winparams.Close()
if winCONTINactive:
    winCONTINactive  = False
    winCONTIN.Close()

