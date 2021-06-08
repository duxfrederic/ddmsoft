#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A small GUI program to interface a custom DDM setup.

@author: Frédéric Dux, biosoft intern@IPC with Jerome Crassous
"""

from   os            import  name
from   os.path       import  expanduser
import PySimpleGUI   as      sg

sg.ChangeLookAndFeel('Reddit')
sg.SetOptions(element_padding=((2, 2), (6, 5)))
themes = sg.ListOfLookAndFeelValues()

if name != 'posix':
    sg.SetOptions(element_padding=((2, 2), (3, 3)))
    import ctypes
    import platform
    
    def make_dpi_aware():
        if int(platform.release()) >= 8:
            ctypes.windll.shcore.SetProcessDpiAwareness(True)
    make_dpi_aware()



menu_def = [['File', ['Open directory', 'Exit'  ]],      
            ['Tools', ["Concatenate videos" , 
                       "Split a video in N matrices (time dependant DDM)"]],
            ["Plotting", ["Show the DDM matrix (image)", "Plot some autocorrelation functions"]],
            [" Batch, export ", [
                          'Save the DDM matrix as a text file',
                          'Save the current fit parameters',
                          'Fit and save all the matrices',
                          'Save the correlation functions'
                          ]],
            ["More Fitting", ["CONTIN"]],
            ['Help', ['About...']], ]    

if name == "posix":
    DEFAULTDATADIRECTORY = expanduser('~')
    DEFAULTFONT          = ('Helvetica', 11)
    TITLEFONT            = ('Cambria', 12)
    slidersize           = 36
    multilinesize        = 83
    fitmenusize          = 40
    paramsize            = 77
    radiosize            = 40
    ddm_data_io_frame_layout = [
         # input directory
         [sg.Text('Directory containing the videos', size=(32, 1), auto_size_text=False, \
                  justification='right', font=DEFAULTFONT), 
          sg.InputText(DEFAULTDATADIRECTORY, size=(75,2), key='inputpath', font=DEFAULTFONT), 
          sg.FolderBrowse(initial_folder=DEFAULTDATADIRECTORY, font=DEFAULTFONT, key="browser")],
         # buttons
         [sg.Button('Load', font=DEFAULTFONT), sg.Quit(font=DEFAULTFONT)]        
    ]

    # frame with the big data crunching parameters
    ddm_compute_frame_layout = [
         # recompute toggles
         [sg.Radio(f'{"Keep existing matrices"}', "recompute", default=True, \
                   enable_events=True, key="notrecompute", font=DEFAULTFONT, size=(radiosize, 1)), 
          sg.Radio(f'{"Re-compute, overwrite existing":<40}', "recompute", \
                      enable_events=True, key='recompute', font=DEFAULTFONT, size=(radiosize, 1))], 
         # get the ncouples and ptperdecade params:
         [sg.Text(f"{'Max number of couples of frames during averaging'}",
                     font=DEFAULTFONT, size=(paramsize, 1)), 
          sg.InputText('300', size=(6,1), justification='right', key='maxcouples', enable_events=True,\
                       tooltip=" 50: fast, noisy. 300: okay slow, okay noise. 0: slow, sets the max possible value ",\
                       font=DEFAULTFONT)],
         [sg.Text(f'{"Number of considered lag times per decade of measured time"}', size=(paramsize,1),
                     font=DEFAULTFONT),
          sg.InputText('20',  size=(6,1), justification='right', key='ptperdecade',enable_events=True,\
                       tooltip=" Proportionally increases the calculation time. You will have more temporal rows in your DDM matrix though. ",
                       font=DEFAULTFONT)],
         [sg.Multiline(default_text='Please press load to display the videos and their parameters.', size=(multilinesize, 6), key='videodescriptions',\
                          enable_events=True, tooltip=" You can change the values in case they are not correct. ",
                          font=DEFAULTFONT)],
         [sg.Button('Process', tooltip=" Launch the data crunching of the videos displayed below ",
                    font=DEFAULTFONT),
          sg.Text("   (direction dependant dynamics: split in ", font=DEFAULTFONT), 
          sg.InputText("1", size=(3,1), key='Nangle', enable_events=True, font=DEFAULTFONT),
          sg.Text("parts.)", font=DEFAULTFONT)]
    ]

    ddm_fitting_frame_layout = [
            [sg.Button("Merge matrices", key="mergebutton", enable_events=True, 
                       tooltip=" Merge the matrices in the list below and creates a new matrix ",
                       font=DEFAULTFONT),
            sg.Button("Average matrices", key="averagebutton", enable_events=True, 
                       tooltip=" Average the matrices in the list below and creates one or more new matrices. ",
                       font=DEFAULTFONT)],
            [sg.Text("Matrices to analyse:",size=(fitmenusize+1,1),
                     font=DEFAULTFONT), 
             sg.Text("Model for the autocorrelation function:", 
                     size=(fitmenusize,1), font=DEFAULTFONT)],
            [sg.Combo([f'{"no video computed yet":<50}'], size=(fitmenusize,1), enable_events=True, 
                      key='computedlist', font=DEFAULTFONT),
             sg.Combo(["Single exponential decay", "Cumulants up to second order", \
                                 "Cumulants up to third order", "Stretched exponential",  \
                                 "Double exponential, second stretched", "Exponential and flow", 
                                 "Stretched exponential and flow", \
                                 "Dble exp (2nd stretched) with flow"], \
                                size=(fitmenusize,1), enable_events=True, key='fitmodel',
                                font=DEFAULTFONT, default_value="Single exponential decay")],
            [sg.Text('qmin:'), sg.Slider(range=(1,127), resolution=1, size=(slidersize+1,15), orientation='horizontal', 
                                         enable_events=True, key='qminslider', font=DEFAULTFONT),
             sg.Text('qmax:'), sg.Slider(range=(1,127), resolution=1, size=(slidersize,15), orientation='horizontal',
                                         enable_events=True, key='qmaxslider', font=DEFAULTFONT, default_value=127)],
            [sg.Text('tmin:'), sg.Slider(range=(0,100), resolution=1, size=(slidersize+1,15), orientation='horizontal', 
                                         enable_events=True, key='dtsminslider', font=DEFAULTFONT),
             sg.Text('tmax:'), sg.Slider(range=(0,100), resolution=1, size=(slidersize,15), orientation='horizontal',
                                         enable_events=True, key='dtsmaxslider', font=DEFAULTFONT, default_value=100)],
            [sg.Button("Initial guess for the fit", key="winparamsfit", enable_events=True, font=DEFAULTFONT),
             sg.Button("Fit the selected matrix", key="fitbutton", enable_events=True,\
                       font=DEFAULTFONT, tooltip=" Fit using the selected model. "),
             sg.Button("Show the fitted parameters", key="showfitparams", enable_events=True, font=DEFAULTFONT)]
    ]

    ddm_plotting_frame_layout = [
            [sg.Text('(Optional) temperature (°C)', font=DEFAULTFONT),
             sg.Input('', key='temperature', size=(15,1),font=DEFAULTFONT),
             sg.Text('(Optional) viscosity', font=DEFAULTFONT),   
             sg.Input('water', key='viscosity', size=(17,1), font=DEFAULTFONT)],
             
            [sg.Button("Plot the matrix and the fit", key='mainplotbutton', enable_events=True,
                       font=DEFAULTFONT), 
             sg.Button("Plot the amplitude, the noise, the diffusion", key="noiseandamp", enable_events=True,
                       font=DEFAULTFONT)]
    ]
            
    ddm_timedependant_frame_layout = [
            [sg.Text("Please select and load a directory containing a single video before using this.",
                     font=DEFAULTFONT)],
            [sg.Text("Number of portions to split the video in:",
                     font=DEFAULTFONT), 
             sg.Input(10, key="npartitions", enable_events=True, size=(20,1),
                      font=DEFAULTFONT),
             sg.Button("Compute the matrices", key='dependantddmccomputebutton', enable_events=True,
                       font=DEFAULTFONT)]
    ]

    # Layout the design of the GUI
    framewidth  = 200
    BORDERWIDTH = 0
    layout = [
         [sg.Menu(menu_def, background_color="#f8f8ff")], 
         
         # frame with all the computation parameters:
         [sg.Frame('DDM matrix computation', ddm_compute_frame_layout, font=TITLEFONT, title_color='blue',\
                   size=(framewidth,200), border_width=BORDERWIDTH)],
         [sg.Text(' '*80)],          
         [sg.Frame('DDM matrix fitting', ddm_fitting_frame_layout, font=TITLEFONT, title_color='blue',\
                   size=(framewidth,200), border_width=BORDERWIDTH)],
         [sg.Text(' '*80)],  
         [sg.Frame('Plotting', ddm_plotting_frame_layout, font=TITLEFONT, title_color='blue',\
                   border_width=BORDERWIDTH)],
     
         [sg.Text(' '*80)],  
         [sg.Text("Idle", key="state", size=(41,3), tooltip=" What the backend is currently doing. ",
                  font=DEFAULTFONT), 
          sg.ProgressBar(100, orientation='h', size=(39, 20), key='progressbar')],
    ]
else:
    DEFAULTDATADIRECTORY =  expanduser('~')
    DEFAULTFONT          = ('Segoe UI', 9)
    TITLEFONT            = ('Segoe UI', 10)

    ddm_data_io_frame_layout = [
         # input directory
         [sg.Text('Directory containing the videos', size=(27, 1), auto_size_text=False, \
                  justification='right', font=DEFAULTFONT), 
          sg.InputText(DEFAULTDATADIRECTORY, size=(90,1),key='inputpath', font=DEFAULTFONT), 
          sg.FolderBrowse(initial_folder=DEFAULTDATADIRECTORY, font=DEFAULTFONT, key="browser")],
         # buttons
         [sg.Button('Load', font=DEFAULTFONT), sg.Quit(font=DEFAULTFONT)]        
    ]

    # frame with the big data crunching parameters
    ddm_compute_frame_layout = [
         # recompute toggles
         [sg.Radio(f'{"Keep existing matrices":<50}', "recompute", default=True, \
                   enable_events=True, key="notrecompute", font=DEFAULTFONT), 
          sg.Radio(f'{"Re-compute DDM matrices (overwrite existing)":<50}', "recompute", \
                      enable_events=True, key='recompute', font=DEFAULTFONT)], 
         # get the ncouples and ptperdecade params:
         [sg.Text(f"{'Max number of couples of frames during averaging (0 for max possible value)':<90}",
                     font=DEFAULTFONT, size=(84,1)), 
          sg.InputText('300', size=(7,1), justification='right', key='maxcouples', enable_events=True,\
                       tooltip=" 50: very fast, very noisy. 600: okay slow, okay noise. 0: slow, extract the most out of your data. ",\
                       font=DEFAULTFONT)],
         [sg.Text(f'{"Number of considered lag times per decade of measured time (the default works nicely)":<90}',
                     font=DEFAULTFONT, size=(84,1)),
          sg.InputText('20',  size=(7,1), justification='right', key='ptperdecade',enable_events=True,\
                       tooltip=" Proportionally increases the calculation time. You will have more temporal rows in your DDM matrix though. ",
                       font=DEFAULTFONT)],
                        # information about the videos
         [sg.Text('Make sure the information below is correct. Correct it if no. (Might want to check your config files)',
                  font=DEFAULTFONT)],
         [sg.Multiline(default_text='Please press load to display the videos and their parameters.', 
                          size=(83, 6), 
                          key='videodescriptions',
                          enable_events=True, tooltip=" You can change the values in case they are not correct. ",
                          font=("Lucida Console", 9))],
         [sg.Button('Process', tooltip=" Launch the data crunching of the videos displayed below ",
                    font=DEFAULTFONT),
          sg.Text("   (direction dependant dynamics: split in ", font=DEFAULTFONT), 
          sg.InputText("1", size=(3,1), key='Nangle', enable_events=True, font=DEFAULTFONT,\
                       tooltip=" In how many directions should the dynamics be probed? \n The bigger the number, the lesser the quality of the statistics for each direction. "),
          sg.Text("parts.)", font=DEFAULTFONT)]
    ]
    slider_size = 31.6
    qslidertooltip = " Use \"Plot the matrix and the fit\" below to assess \n in which range of wavenumber q the correlation function looks good. "
    tslidertooltip = " Which times should be considered for the fit ? "
    ddm_fitting_frame_layout = [
            [sg.Button("Merge matrices", key="mergebutton", enable_events=True, 
                       tooltip=" Merge the matrices in the list below and creates a new matrix ",
                       font=DEFAULTFONT),
             sg.Button("Average matrices", key="averagebutton", enable_events=True, 
                       tooltip=" Average the matrices in the list below and creates one or more new matrices. ",
                       font=DEFAULTFONT)],
            [sg.Text("Matrices to analyse:",size=(45,1), font=DEFAULTFONT),
             sg.Text("Fitting model of the autocorrelation function:", 
             size=(43,1), font=DEFAULTFONT)],
            [sg.Combo([f'{"no video computed yet":<40}'], size=(44,6), enable_events=True, 
                      key='computedlist', font=DEFAULTFONT),
             sg.Combo(["Single exponential decay", "Cumulants up to second order", \
                                 "Cumulants up to third order", "Stretched exponential",  \
                                 "Double exponential, second stretched", "Exponential and anisotropic flow", 
                                 "Stretched exponential and anisotropic flow",\
                                 "Double exponential (2nd stretched) with flow"], \
                                size=(44,6), enable_events=True, key='fitmodel',
                                font=DEFAULTFONT, default_value='Single exponential decay')],
            [sg.Text('qmin:', font=('Lucida Console',10)), sg.Slider(range=(1,127), resolution=1, size=(slider_size,15), orientation='horizontal', 
                                         enable_events=True, key='qminslider', font=DEFAULTFONT, default_value=1, tooltip=qslidertooltip),
             sg.Text('qmax:', font=('Lucida Console',10)), sg.Slider(range=(1,127), resolution=1, size=(slider_size-1,15), orientation='horizontal',
                                         enable_events=True, key='qmaxslider', font=DEFAULTFONT, default_value=127, tooltip=qslidertooltip)],
            [sg.Text('tmin:', font=('Lucida Console',10)), sg.Slider(range=(0,50), resolution=1, size=(slider_size,15), orientation='horizontal', 
                                         enable_events=True, key='dtsminslider', font=DEFAULTFONT, default_value=0, tooltip=tslidertooltip),
             sg.Text('tmax:', font=('Lucida Console',10)), sg.Slider(range=(0,50), resolution=1, size=(slider_size-1,15), orientation='horizontal',
                                         enable_events=True, key='dtsmaxslider', font=DEFAULTFONT, default_value=50, tooltip=tslidertooltip)],
            [sg.Button("Initial guess for the fit", key="winparamsfit", enable_events=True, font=DEFAULTFONT),
             sg.Button("Fit the selected matrix", key="fitbutton", enable_events=True,  font=DEFAULTFONT, tooltip=" Fit using the selected model. "),
             sg.Button("Show the fitted parameters", key="showfitparams", enable_events=True, font=DEFAULTFONT)]
    ]

    ddm_plotting_frame_layout = [
            [sg.Text('(Optional) temperature (°C)', font=DEFAULTFONT),
             sg.Input('', key='temperature', size=(25,1),font=DEFAULTFONT,
                      tooltip=" used to convert the diffusion coefficients to equivalent hydrodynamic radii. Leave empty to keep the diffusion coefficient. "),
             sg.Text('(Optional) viscosity', font=DEFAULTFONT),   
             sg.Input('water', key='viscosity', size=(25,1), font=DEFAULTFONT,
                      tooltip=" If the solvent is not water, enter here the dynamic viscosity in SI units. ")],
             
            [sg.Button("Plot the matrix and the fit", key='mainplotbutton', enable_events=True,
                       font=DEFAULTFONT), 
             sg.Button("Plot the amplitude, the noise, the diffusion", key="noiseandamp", enable_events=True,
                       font=DEFAULTFONT)]
    ]
            
    ddm_timedependant_frame_layout = [
            [sg.Text("Please select and load a directory containing a single video before using this.",
                     font=DEFAULTFONT)],
            [sg.Text("Number of portions to split the video in:",
                     font=DEFAULTFONT), 
             sg.Input(10, key="npartitions", enable_events=True, size=(20,1),
                      font=DEFAULTFONT),
             sg.Button("Compute the matrices", key='dependantddmccomputebutton', enable_events=True,
                       font=DEFAULTFONT)]
    ]

    # Layout the design of the GUI
    framewidth = 150
    BORDERWIDTH = 0
    layout = [
         [sg.Menu(menu_def, )], 
         
         # frame with all the computation parameters:
         [sg.Frame('DDM matrix computation', ddm_compute_frame_layout, font=TITLEFONT, title_color='blue',\
                   size=(framewidth,200), border_width=BORDERWIDTH)],
          
         [sg.Frame('DDM matrix fitting', ddm_fitting_frame_layout, font=TITLEFONT, title_color='blue',\
                   size=(framewidth,200), border_width=BORDERWIDTH)],
         [sg.Frame('Plotting', ddm_plotting_frame_layout, font=TITLEFONT, title_color='blue',
                   border_width=BORDERWIDTH)],
  
         [sg.Text("Idle", key="state", size=(40,3), tooltip=" What the backend is currently doing. ",
                  font=DEFAULTFONT), 
          sg.ProgressBar(100, orientation='h', size=(32, 20), key='progressbar')],
    ]


if __name__ == "__main__":
    window = sg.Window('DDMSoft', layout, return_keyboard_events=True, icon='logo.png')
    window.set_icon("icon.png")
    while 1:
        event, value = window.Read(timeout=100)
        if event in ('Quit', 'Exit', None):
            break
    window.Close()
