##### a list of the available fitting models and their abbreviation counterpart:
FITMODELS =     {"Single exponential decay": "cumulant_1",
                 "Cumulants up to second order":"cumulant_2", 
                 "Cumulants up to third order":"cumulant_3", 
                 "Stretched exponential":"stretch",
                 "Double exponential, second stretched":"dblexp_2ndstretched",
                 "Exponential and flow":"expcos",
                 "Stretched exponential and flow":"expcosstretch",
                 "Dble exp (2nd stretched) with flow":"dblexpcosstretch"}
##### for each fitting model, a list of the names of the parameters:
FITPARAMNAMES = {"cumulant_1" : ["Diffusion coefficient", "A", "B"],
                 "cumulant_2" : ["Diffusion coefficient", "2nd cumulant", "A", "B"],
                 "cumulant_3" : ["Diffusion coefficient", "2nd cumulant", "3rd cumulant", "A", "B"],
                 "stretch"    : ["Diffusion coefficient", "stretch parameter", "A", "B"],
                 "dblexp_2ndstretched": ["Diffusion coefficient 1", "Diffusion coefficient 2",\
                                         "stretch coefficient (2)", "weighting parameter", "A", "B"],
                 "expcos"     : ["Diffusion coefficient", "effective flow speed", "A", "B"],
                 "expcosstretch" : ["Diffusion coefficient", "effective flow speed", "stretch parameter", "A", "B"],
                 "dblexpcosstretch": ["Diffusion coefficient 1", "Diffusion coefficient 2",\
                                       "stretch coefficient (2)", "weighting parameter", \
                                        "effective flow speed", "A", "B"]}
         
                 
##### for each fitting model and for each relevant param, a default value:
FITPARAMDEFAULTS={"cumulant_1" : [2e-12, '', ''],
                  "cumulant_2" : [2e-12, 0, '', ''],
                  "cumulant_3" : [2e-12, 0, 0, '', ''],
                  "stretch"    : [2e-12, 1, '', ''],
                  "dblexp_2ndstretched": [2e-12, 2e-13, 1, 1, '', ''],
                  "expcos"     : [2e-12, 0, '', ''],
                  "expcosstretch" : [2e-12, 0, 1, '', ''],
                  "dblexpcosstretch": [2e-12, 2e-13, 1, 1, 0, '', '']}

