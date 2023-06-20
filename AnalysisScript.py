import FemtoAnalysis as FA
import FemtoDreamReader as FDR

#
#   All options for UFFA:
#
#   settings = {
#           "function":     None,       # Function Name: 'cf', 'tf'
#           "pair":         None,       # No effect yet
#           "path":         "",         # system's directory of 'file'
#           "file":         None,       # name of 'file', can also be the full path to the file
#           "fileTDir":     "",         # TDirectory/TList inside of 'file'
#           "newfile":      None,       # 'new', 'recreate', 'update'
#           "mc":           None,       # bool to include mc data / file of mc data
#           "mcTDir":       "",         # if 'mc'/'templates' is a file, TDirectory/TList of the file
#           "outDir":       "",         # system's directory for output file
#           "rename":       None,       # rename output file
#           "bins":         None,       # bins for diff analysis (list) / template fits (int/list)
#           "rebin":        None,       # rebin factor (int) / factors (list)
#           "atype":        None,       # analysis type: "diff", "int"
#           "htype":        None,       # histogram type: "kstar", "kmult", "kmt"
#           "tftype":       None,       # template fit type: "dca", "cpa" (cpa not fully tested)
#           "data":         None,       # No effect yet
#           "templates":    None,       # MonteCarlo templates, 'mcTDir' is used to set the directory
#           "namelist":     None,       # list of names for mc templates
#           "fitrange":     None,       # fit range for the template fits
#           "normalize":    None,       # normalization range for the correlation function
#           "debug":        False       # debug option
#           }
#


ipath = ""
opath = ""

filename = "AnalysisResults_LHC22q_pass3.root"

templates = ""

# f1 = FDR.FemtoDreamReader(templates, "_base-tempFit")
# mc_plots = []
# mc_plots.append(f1.getHisto("Tracks_one_MC/hDCAxy_Primary"))
# mc_plots.append(f1.getHisto("Tracks_one_MC/hDCAxy_DaughterLambda"))
# mc_plots.append(f1.getHisto("Tracks_one_MC/hDCAxy_DaughterSigmaplus"))
# mc_plots.append(f1.getHisto("Tracks_one_MC/hDCAxy_Fake"))

namelist = ['prim', 'lam', 'sig', 'fake']

settings_cf = {
        "function":     'cf',
        "file":         ipath + filename,
        "fileTDir":     "",
        "newfile":      "new",
        "outDir":       opath,
        "mc":           None,
        "mcTDir":       "",
        "rename":       None,
        "atype":        'int',
        "htype":        'kmult',
        "bins":         [0, 20, 40, 60],
        "rebin":        None,
        "normalize":    [0.2, 0.3],
        "debug":        True
    }

settings_tf = {
        "function":     'tf',
        "file":         ipath + filename,
        "fileTDir":     "_base-tempFit",
        "newfile":      "new",
        "templates":    templates,
        #"templates":    mc_plots,
        "mcTDir":       "_base-tempFit",
        "outDir":       opath,
        "rename":       None,
        "bins":         [0.500, 0.678, 0.856, 1.033, 1.211, 1.388, 1.567, 1.743, 1.921, 2.099, 2.453, 2.986, 4.051],
        "rebin":        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 4],
        "tftype":       'dca',
        "namelist":     namelist,
        "fitrange":     0.3,
        "debug":        False
        }

FA.UFFA(settings_cf)

