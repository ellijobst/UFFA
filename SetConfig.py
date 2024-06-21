import ROOT
import FileUtils as FU


# generates list with rebin factors
def bin2list(rebin):
    rebin_list = []
    if type(rebin) == int or type(rebin) == str:
        rebin = [rebin]
    elif type(rebin) != list:
        return None
    rebin_list.extend(rebin)
    return rebin_list

# generates the proper settings dictionary
def config(dic_conf):
    """
    This function sets up all the configurable options
    so that it is consitent and expandable.

    Current options include:
            "function":     'cf', 'syst', 'tf' -> for correlation function, systematics, template fits
            "pair":         'pp', 'pl' -> for q&a plots relevant for individual analyses
            "path":         "string" -> full path to the root file, might include ~/ for home directory
            "file":         "string" -> name of the root file
            "fullpath":     "string" -> full path and file name equal to "path" + "file"
            "outDir":       "string" -> output directory
            "rename":       "string" -> rename output file
            "fileTDir":     "string" -> root file directory: path to directory inside the root file
            "SE_path":       "string" -> path + name of the se plot inside the provided "fileTDir" if given
            "ME_path":       "string" -> same as SE_path but for the ME distribution
            "newfile":      'new', 'recreate', 'update' -> same option as in ROOT, 'new' will rename if file already exists
            "mc":           'true', 'false' -> save monte carlo data from provided root file
            "mcTDir":       "string" -> root file directory for the monte carlo data
            "bins":         [list of floats] -> binning for differential analysis
            "diff3d":       'mt', 'mult' -> project 3D plots first in mt/mult 2D and after in mult/mt 1D or vice versa
            "bins3d":       [list of floats] -> binning for 3D plots to 2D plots
            "yield":        [GeV, Deviation] -> integrated analysis: include systematics inside deviation for the GeV range
            "rebin":        int or [list of ints] -> rebin output plots
                                for tf: int -> all dca/cpa rebinned with int
                                or      [list of ints] -> each int will correspond to one range of the binning if provided
            "atype":        'int', 'dif' -> integrated analysis or differential analysis
            "htype":        'k', 'mt', 'mult', 'mt3d', 'mult3d', 'mtmult', 'rew3d'
                                'k'     -> k* - relative pair momentum distribution
                                'mt'  -> mt vs k* distribution
                                'mult'  -> multiplicity vs k* distribution
                                'mt3d' -> mt vs k* from 3D distribution but integrated in mult
                                'mult3d' -> mult vs k* from 3D distribution but integrated in mt
                                'mtmult' -> mult vs mt vs k* 3D distribution
                                'rew3d' -> mult vs mt vs k* 3D differentially in mt and reweighted in mult
            "tftype":       'dca', 'cpa' -> option for the template fit plots
            "templates":    [list of th1 plots] -> list of dca/cpa plots for fitting
            "temp_init":   list of values to initialize fitting parameters
            "temp_limits": list of limits of the fitting parameters
            "namelist":     [list of strings] -> names of dca/cpa plots for fitting
            "fitrange":     float -> fitrange for the template fitter
            "normalize":    [float, float] -> normalization range for the correlation function
            "percentile":   [low, upper] or int -> lower and upper edges or 0 to int for percentile cut
            "include":      "string" or [list of strings] -> include these variations in the systematics
            "exclude":      "string" or [list of strings] -> exclude these variations in the systematics
            "interactive":  'True', 'False' -> include/exclude interactively variations in terminal
            "debug":        'True', 'False' -> debug information in console
            "print":        'True', 'False' -> print canvas as png
    """

    # settings dictionary skeleton
    dic = {
            "function":         None,
            "pair":             None,
            "path":             "",
            "file":             None,
            "fullpath":         None,
            "fileTDir":         "",
            "SE_path":          "",
            "SE_path":          "",
            "newfile":          None,
            "mc":               None,
            "mcTDir":           "",
            "outDir":           "",
            "rename":           None,
            "bins":             None,
            "bins3d":           None,
            "diff3d":           "",
            "diff3d2":          "",
            "yield":            None,
            "rebin":            None,
            "atype":            None,
            "htype":            None,
            "tftype":           None,
            "data":             None,
            "templates":        None,
            "temp_init":        None,
            "temp_limits":      None,
            "temp_fraction":    None,
            "namelist":         None,
            "fitrange":         None,
            "percentile":       None,
            "rewrange":         None,
            "normalize":        None,
            "include":          None,
            "exclude":          None,
            "debug":            False,
            "print":            False,
            "interactive":      False,
            "dimension":        None,    
            "kstar_axis":       None, 
            "mult_axis":        None,
            "higher_dimension_axes": None,
            "reweight":         False,
        }

    # keys to set values
    keys_k      = ['k', 'kstar']
    keys_mult   = ['mult', 'kmult']
    keys_mult3d = ['mult3d', 'kmult3d']
    keys_rew3d  = ['rew3d', 'rewmult']
    keys_mt     = ['mt', 'kmt']
    keys_mt3d   = ['mt3d', 'kmt3d']
    keys_mtmult = ['mtmult','kmtmult']
    keys_4d     = ['perc', '4d', '4dim', '4dims']
    keys_rew4d  = ['rew4d', 'rewperc']

    keys_int    = ['int', 'integrated']
    keys_dif    = ['diff', 'dif', 'differential']

    # initialize values
    entries = ['function',      # function to be used
               'path',          # path to input file
               'fileTDir',      # root folder where getter functions are used
               'mc',            # path to mc root file
               'mcTDir',        # mc folder inside root file
               'rename',        # rename output file
               'templates',     # list of template histos
               'namelist',      # list of template histos names
               'temp_init',     # initialize template fitting values
               'temp_limits',   # limits for template fitting parameters
               'temp_fraction', # set fraction with dictionary {'name', 'temp_init', 'temp_limits'}
               'fitrange',      # fit range for templates
               'signalrange',   # signal range to evaluate template fractions
               'normalize',     # range to normalize cf
               'data',          # not used
               'bins3d',        # bins to split 3d histo
               'bins',          # bins to split 2d histo
               'SE_path',       # root folder with SE histo
               'ME_path',       # root folder with ME histo
               "filepath",      # file path of input root file
               "atype",
               "dimension",     # dimension of the input histogram
               "kstar_axis",    # index of the axis corresponding to k-star
               "mult_axis",     # index of the axis corresponding to multiplicity
               "higher_dimension_axes", # list of [[index, [range]]] for higher dimensions
               "reweight",      # use reweighting or not
               ]
    for entry in entries:
        if entry in dic_conf:
            dic[entry] = dic_conf[entry]

    # type of particle pair
    if 'pair' in dic_conf:
        if dic_conf['pair']:
            dic['pair'] = dic_conf['pair'].lower()

    # file name, file directory
    if 'file' in dic_conf:
        if dic_conf['file']:
            path_name = dic_conf['file'].rsplit('/', 1)
            if len(path_name) == 1:
                dic['file']  = path_name[0]
            else:
                dic['path'] = FU.path_expand(path_name[0]) + '/'
                dic['file']  = path_name[1]
            dic['fullpath'] = dic['path'] + dic['file']

    # output directory
    if 'outDir' in dic_conf:
        if dic_conf['outDir'] != "" and dic_conf['outDir']:
            dic['outDir'] = FU.path_expand(dic_conf['outDir'])
            if ROOT.gSystem.AccessPathName(dic['outDir']):
                print("output directory \"" + dic['outDir'] + "\" does not exist!")
                exit()
    else:
        dic['outDir'] = dic['path']

    # create file
    if 'newfile' in dic_conf:
        if dic_conf['newfile'] in [1, "new"]:
            dic['newfile'] = "new"
        if dic_conf['newfile'] in [2, "recreate"]:
            dic['newfile'] = "recreate"
        if dic_conf['newfile'] in [3, "update"]:
            dic['newfile'] = "update"

    # histogram and analysis type #TODO: was ist das? WArum resettet es unsere eingaben?
    if 'type' in dic_conf:
        dic_conf['atype'] = dic_conf['type'][0]
        dic_conf['htype'] = dic_conf['type'][1]
        if len(dic_conf['type']) > 2:
            dic_conf['diff3d'] = dic_conf['type'][2]

    # analysis type
    if 'atype' in dic_conf:
        atype = dic_conf['atype']
        if type(atype) == str:
            atype = atype.lower()
        if atype in keys_int:
            dic['atype'] = 'int'
        elif atype in keys_dif:
            dic['atype'] = 'dif'

    # histogram type
    if 'htype' in dic_conf:
        htype = dic_conf['htype']
        if type(htype) == str:
            htype = htype.lower()
        if htype in keys_4d:
            dic['htype'] = '4d'
        elif htype in keys_rew4d:
            dic['htype'] = 'rew4d'
            dic['diff3d'] = 'mt'
        elif htype in keys_k:
            dic['htype'] = 'k'
        elif htype in keys_mult:
            dic['htype'] = 'mult'
        elif htype in keys_mult3d:
            dic['htype'] = 'mult3d'
        elif htype in keys_rew3d:
            dic['htype'] = 'rew3d'
            dic['diff3d'] = 'mt'
        elif htype in keys_mt:
            dic['htype'] = 'mt'
        elif htype in keys_mt3d:
            dic['htype'] = 'mt3d'
        elif htype in keys_mtmult:
            dic['htype'] = 'mtmult'

    # template fit type
    if 'tftype' in dic_conf:
        if dic_conf['tftype']:
            tftype = dic_conf['tftype'].lower()
            if tftype == 'dca':
                dic['tftype'] = 'dca'
            elif tftype == 'cpa':
                dic['tftype'] = 'cpa'

    # which axis to be used for the first split in a 3D analysis
    if 'diff3d' in dic_conf:
        diff3d = dic_conf['diff3d']
        if type(diff3d) == str:
            diff3d = diff3d.lower()
        if diff3d in keys_mult:
            dic['diff3d'] = 'mult'
            dic['diff3d2'] = 'mt'
        elif diff3d in keys_mt:
            dic['diff3d'] = 'mt'
            dic['diff3d2'] = 'mult'

    # yield setting to exclude systematic variations below a value of GeV that vary by a given percentage
    # input: [GeV, %]
    if 'yield' in dic_conf:
        dic['yield'] = dic_conf['yield']
        if dic['yield'] and len(dic['yield']) != 2:
            print("'yield' accepts [GeV, deviation], where GeV defines [0, GeV) and deviation the percentage!")

    # rebin factor/s
    if 'rebin' in dic_conf:
        dic['rebin'] = bin2list(dic_conf['rebin'])

    # rewrange
    if 'rewrange' in dic_conf:
        if dic_conf['rewrange'] and type(dic_conf['rewrange']) != list:
            print("'rewrange' accepts [min, max]!")
        dic['rewrange'] = dic_conf['rewrange']

    # percentile range
    if 'percentile' in dic_conf:
        if type(dic_conf['percentile']) == int:
            dic['percentile'] = [0, dic_conf['percentile']]
        elif type(dic_conf['percentile']) == list:
            if len(dic_conf['percentile']) == 1:
                dic['percentile'] = [0, dic_conf['percentile'][0]]
            else:
                dic['percentile'] = dic_conf['percentile']

    # include variations
    if 'include' in dic_conf:
        dic['include'] = bin2list(dic_conf['include'])

    # exclude variations
    if 'exclude' in dic_conf:
        dic['exclude'] = bin2list(dic_conf['exclude'])

    if 'debug' in dic_conf:
        dic['debug'] = bool(dic_conf['debug'])

    if 'print' in dic_conf:
        dic['print'] = bool(dic_conf['print'])

    if 'interactive' in dic_conf:
        dic['interactive'] = bool(dic_conf['interactive'])

    return dic

