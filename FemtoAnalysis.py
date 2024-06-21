import ROOT
import FileUtils as FU
import FemtoDreamSaver as FDS
import FemtoDreamReader as FDR
import CorrelationHandler as CH
import CombinedTemplateFit as TF


# class that handles the retrieving of histos and computing of correlation functions
class cf_handler():
    def __init__(self, FileReader, conf):
        self._file  = FileReader
        self._pair  = conf['pair']
        self._atype = None                         # analysis type
        # self._htype = conf['htype']             # histo type
        self._mc    = conf['mc']                # bool monte carlo data
        self._bins  = conf['bins']              # bin range for differential
        # self._diff3d = conf['diff3d']           # which axis to split first in a 3D analysis
        # self._bins3d = conf['bins3d']   # bin range for the first differential split in case of a 3D analysis
        self._rebin_factors = conf['rebin_factors']             # rebin factors for all se, me, cf plots
        self._norm  = conf['normalize']         # normalization range
        self._perc  = conf['percentile']        # percentile range
        self._rew_range = conf['rewrange']      # reweighting range
        self._name_se = conf['SE_path']
        self._name_me = conf['ME_path']

        self._dimension = conf["dimension"]
        self._kstar_axis = conf["kstar_axis"]
        self._reweight_axis = conf["rew_axis"]
        self._reweighting_bins = conf["rew_bins"]
        self._rew_range = conf['rew_range']  
        self._projection_axis = conf["projection_axes"]
        
        self._se = None
        self._me = None
        self._se_mc = None
        self._me_mc = None
        self._event = None
        self._tracks = None
        self._tracks_mc = None
        self._v0 = None
        self._get_histos()

    # retrieves histos from the provided file reader
    def _get_histos(self):
        if self._name_se and self._name_me:
            self._se = self._file.get_histo(self._name_se)
            self._me = self._file.get_histo(self._name_me)
        else:
            print("No paths to SE and ME histograms were provided in config. Exiting...")
            exit()

        if self._mc:
            self._tracks_mc = self._file.get_tracks_mc()

        if self._pair == 'pp':
            self._event = self._file.get_event()
            self._tracks = self._file.get_tracks()
        elif self._pair == 'pl':
            self._v0 = self._file.get_v0()

    #TODO
    def get_histos(self):
        """
        computes the cf for integrated or differential analysis and for mc data
        and returns the histos for all the different options:
        [histos, histos_unw, histos_mc, histos_unw_mc, self._event, self._tracks, self._tracks_mc]
        """
        histos = []
        histos_mc = []
        histos_unw = []
        histos_unw_mc = []

        if self._dimension == 1:  # integrated analysis (1D)
            self._atype = "int"
            histos, histos_unw = AnalysisUtils.getIntegrated(self._se, self._me, self._htype, self._rebin_factors, self._norm, self._rew_range)
            if self._mc:
                histos_mc, histos_unw_mc = AnalysisUtils.getIntegrated(self._se_mc, self._me_mc, self._htype, self._rebin_factors, self._norm, self._rew_range)
        elif self._dimension >= 2:  # differential analysis
            self._atype = "dif"  
            if self._dimension == 3: 
                se2d, me2d = self.Project3Dto2D()
            elif self._dimension >= 4: 
                se2d, me2d = self.ProjectNDto2D()
            elif self._dimension == 2:
                se2d, me2d = self._se, self._me
            else:
                print("Invalid input histogram dimension selected, expecting integer. Exiting...")
                exit()
            
            if self._reweight:
                rew_histos = self.Reweight2D(se2d, me2d,)
                self.GetDifferential(rew_histos)
            else:
                self.GetDifferential(se2d, me2d)
            # get differntial
            


            #     histos, histos_unw = getDiffReweight3D(self._se, self._me, self._bins3d, self._bins, self._rebin_factors, self._norm, self._rew_range)
            # else:
            #     histos = getDifferential3D(self._se, self._me, self._diff3d, self._bins3d, self._bins, self._rebin_factors, self._norm)
            # #4D analysis
            
            #     # se3d, me3d = getProj4d(self._se, self._me, self._perc, self._projection_axis)
            #     se2d, me2d = self.ProjectNDto2D()
            #     if self._reweight:
            #         histos, histos_unw = getDiffReweight3D(se3d, me3d, self._bins3d, self._bins, self._rebin_factors, self._norm, self._rew_range)
            #     else:
            #         histos = getDifferential3D(se3d, me3d, self._diff3d, self._bins3d, self._bins, self._rebin_factors, self._norm)
            # # 2D
            # else:
            #     histos = getDifferential(self._se, self._me, self._htype, self._bins, self._rebin_factors, self._norm)
            #     if self._mc:
            #         histos_mc = getDifferential(self._se_mc, self._me_mc, self._htype, self._bins, self._rebin_factors, self._norm)

        return [histos, histos_unw, histos_mc, histos_unw_mc, self._event, self._tracks, self._tracks_mc, self._v0]

    # returns a list of cf and their rebinned version
    # [[cf, [rebin 1, rebin 2, ...]], [bin2...], ...] same for unweighted if integrated analysis
    #TODO
    def get_cf(self):
        histos = []
        histos_unw = []
        cf_list = []
        cf_list_unw = []

        # integrated analysis
        if self._atype == 'int':
            histos, histos_unw = AnalysisUtils.getIntegrated(self._se, self._me, self._htype, self._rebin_factors, self._norm)
            if self._htype == 'mult':
                cf_list_unw.append(histos_unw[1])
                cf_list_unw.append([])
        # differential analysis
        elif self._atype == 'dif':
            histos = AnalysisUtils.GetDifferential(self._se, self._me, self._htype, self._bins, self._rebin_factors, self._norm)
        cf_list.append([histos[1][2], []])                          # cf, for differential 1st bin

        # rebinned entries appended to the empty list for the first bin
        if self._rebin_factors:
            for n in range(len(self._rebin_factors)):
                cf_list[0][1].append(histos[1][3][n][2])            # rebinned cf
                if self._atype == 'int' and self._htype == 'mult':
                    cf_list_unw[1].append(histos_unw[2][n][1])      # rebinned unw cf for integrated

        # repeat for the rest of the bins in case of differential analysis
        if self._atype == 'dif':
            for n in range(2, len(self._bins)):
                cf_list.append([histos[n][2], []])
                if self._rebin_factors:
                    for nn in range(len(self._rebin_factors)):
                        cf_list[n][1].append(histos[n][3][nn][2])    # rebinned cf appended to rebin list
        return [cf_list, cf_list_unw]

    # returns a list of se and their rebinned version
    #TODO
    def get_se(self):
        histos = []
        se_list = []

        # integrated analysis
        if self._atype == 'int':
            histos, histos_unw = AnalysisUtils.getIntegrated(self._se, self._me, self._htype, self._rebin_factors, self._norm)
        # differential analysis
        elif self._atype == 'dif':
            histos = AnalysisUtils.GetDifferential(self._se, self._me, self._htype, self._bins, self._rebin_factors, self._norm)
        se_list.append([histos[1][0], []])                          # se for differential 1st bin

        # rebinned entries appended to the empty list for the first bin
        if self._rebin_factors:
            for n in range(len(self._rebin_factors)):
                se_list[0][1].append(histos[1][3][n][0])            # rebinned se

        # repeat for the rest of the bins in case of differential analysis
        if self._atype == 'dif':
            for n in range(1, len(self._bins) - 1):
                se_list.append([histos[n][0], []])
                if self._rebin_factors:
                    for nn in range(len(self._rebin_factors)):
                        se_list[n][1].append(histos[n][3][nn][0])    # rebinned se appended to rebin list
        return se_list

    # returns all the cf's for a 3D mt/mult histo
    # [[[bin1-1 cf, [rebin cf]], [bin1-2 cf, [rebin cf]], ...], [[bin2-1 cf, [rebin cf]], [bin2-2 cf, [rebin cf]], ...], ...]
    # not used for cf
    def get_cf_3d(self):
        cf_list = []

        se = self._se
        me = self._me
        if self._htype in ['4d', 'rew4d']:
            se, me = AnalysisUtils.getProj4d(self._se, self._me, self._perc)

        if self._htype in ['rew3d', 'rew4d']:
            histos, histos_unw = AnalysisUtils.getDiffReweight3D(se, me, self._bins3d, self._bins, self._rebin_factors, self._norm, self._rew_range)
        else:
            histos = AnalysisUtils.getDifferential3D(se, me, self._diff3d, self._bins3d, self._bins, self._rebin_factors, self._norm)

        histos = histos[1:]     # remove TH3 histos
        for n, bin1 in enumerate(histos):
            cf_list.append([])
            bin1 = bin1[1:]     # remove TH2 histos
            for nn, th1 in enumerate(bin1):
                cf_list[n].append([th1[2], []])
                if self._rebin_factors:
                    for nnn in range(len(self._rebin_factors)):
                        cf_list[n][nn][1].append(th1[3][nnn][2])

        return cf_list

    # returns all the cf's for a 3D mt/mult histo
    #is not used for cf
    def get_se_3d(self):
        se_list = []

        se = self._se
        me = self._me
        if self._htype in ['4d', 'rew4d']:
            se, me = AnalysisUtils.getProj4d(self._se, self._me, self._perc)

        if self._htype in ['rew3d', 'rew4d']:
            histos, histos_unw = AnalysisUtils.getDiffReweight3D(se, me, self._bins3d, self.bins, self._rebin_factors, self._nrm, self._rew_range)
        else:
            histos = AnalysisUtils.getDifferential3D(se, me, self._diff3d, self._bins3d, self._bins, self._rebin_factors, self._norm)

        histos = histos[1:]     # remove TH3 histos
        for n, bin1 in enumerate(histos):
            se_list.append([])
            bin1 = bin1[1:]     # remove TH2 histos
            for nn, th1 in enumerate(bin1):
                se_list[n].append([th1[0], []])
                if self._rebin_factors:
                    for nnn in range(len(self._rebin_factors)):
                        se_list[n][nn][1].append(th1[3][nnn][0])

        return se_list
    
    def GetDifferential(self, rew_histos, title = None):
        """
        for 2d histogram
        returns [[iSE, iME], [se, me, cf]] for a list of mt or mult ranges
        [[iSE, iME], [se, me, cf, [rebin: [...], [...], ...], [bin 2 [rebin]], ...]
        """
        histos = []
        conf = "" if not title else title + " "         # append to given name
        norm = self._normalization_range
        #htype should be mult or mt
        # assert self._reweight_axis != None, "getDifferential: no differential axis input!"

        # bins=self._reweight_bins

        # # axis for reweighting is set to y axis automatically when projecting down from higher dimensions
        # if self._dimension >= 3:
        #     axis = se2d.GetYaxis()
        # elif self._dimension == 2:
        #     if self._reweight_axis == 0:
        #         axis = se2d.GetXaxis()
        #     elif self._reweight_axis == 1:
        #         axis = se2d.GetYaxis()
        #     else:
        #         print("Invalid reweighting axid index selected. Cannot be higher than 1 for a 2D histogram. Exiting...")
        #         exit()

        # htype = axis.GetTitle()
        # conf += "htype"+": "
        # histos.append([se2d.Clone(f"SE k{htype}"), me2d.Clone(f"ME k{htype}")]) 

        # print(htype)
        # print(bins)

        # divide in bins
        # mt_histos = self.Get2DHistosFromBins(se2d, me2d, htype)

        for n, [name, se, me] in enumerate(rew_histos, 1):
            histos.append(AnalysisUtils.getCorrelation(se, me, name, conf + name, norm))
            histos[n].append([])
            if self._rebin_factors:       # append a list of rebinned [se, me, cf] in the original [se, me, cf, []]
                for factor in self._rebin_factors:
                    se_rebin = AnalysisUtils.rebin_hist(se, factor)
                    me_rebin = AnalysisUtils.rebin_hist(me, factor)
                    rebin_conf = " rebin: " + str(factor)
                    histos[n][3].append(AnalysisUtils.getCorrelation(se_rebin, me_rebin, name, conf + name + rebin_conf, norm))
        return histos
    
    #move to analysis class
    def Get2DHistosFromBins(self, se, me, label):
        """
        this function divides a 2D histogram into the reweighting bins
        returns list of [name, se, me] tuples for each bin
        """
        bins = self._reweight_bins

        if self._dimension >= 3:
            axis = se.GetYaxis()
        elif self._dimension == 2:
            if self._reweight_axis == 0:
                axis = se.GetXaxis()
            elif self._reweight_axis == 1:
                axis = se.GetYaxis()
            else:
                print("Invalid reweighting axid index selected. Cannot be higher than 1 for a 2D histogram. Exiting...")
                exit()

        if type(bins) == list:
            limits = []
            for value_diff in bins:
                value_bin = axis.FindBin(float(value_diff))
                limits.append(value_bin)
        else:
            print("Get2DHistosFromBins: bin input \"" + str(bins) + "\" not a list of ranges!")
            exit()

        histos = []
        #loop through the reweight bins
        for n in range(1, len(limits)):

            diff_low = bins[n - 1]
            diff_up  = bins[n]
            print(diff_low, diff_up)

            bin_low = axis.FindBin(diff_low)
            bin_up  = axis.FindBin(diff_up)
            print(bin_low, bin_up)
            print("bin_up low edge:", axis.GetBinLowEdge(bin_up))
            print("bin_up - 1 =", bin_up - 1)
            print("bin_up upper edge:", axis.GetBinLowEdge(bin_up))

            if diff_up == axis.GetBinLowEdge(bin_up):
                bin_up -= 1

            name = label + ": [%.2f-%.2f)" % (bins[n-1], bins[n])
            axis.SetRange(bin_low, bin_up)
            axis.SetRange(bin_low, bin_up)
            histos.append([name, se.Clone(name), me.Clone(name)])

        return histos
    
    #move to analysis class
    def ProjectNDto2D(self):
        # clone SE,ME distribution
        seND = self._se.Clone("Ndim_se")
        meND = self._se.Clone("Ndim_me")

        for proj_axis in self._projection_axes:
            axis_index = proj_axis[0]
            proj_range = proj_axis[1]

            axis = seND.GetAxis(axis_index)

            # set the projection range
            perc_low = proj_range[0]
            perc_up  = proj_range[1]

            # get bins
            bin_low = axis.FindBin(perc_low)
            bin_up  = axis.FindBin(perc_up)

            # lower the upper bin number if it is right on the next bin's lower edge
            if perc_up == axis.GetBinLowEdge(bin_up):
                bin_up -= 1

            # adjust the axis ranges
            seND.GetAxis(axis_index).SetRange(bin_low, bin_up)
            meND.GetAxis(axis_index).SetRange(bin_low, bin_up)

        se = seND.Projection(self._kstar_axis, self._reweight_axis)
        me = meND.Projection(self._kstar_axis, self._reweight_axis)

        return [se, me]
    
    #move to analysis class
    def Project3Dto2D(self):
        """
        this function projects a 3D histogram to a 2D histogram.
        The output 2D histogram will have k-star as its XAxis!
        Use this when reweighting is True 
        """
        # clone SE,ME distribution
        se3D = self._se.Clone("3dim_se")
        me3D = self._se.Clone("3dim_me")

        assert len(self._projection_axis) == 1, "Project2Dto3D: Too many projection axes specified for a 3D histogram!"
        
        axis_index = self._projection_axis[0]
        proj_range = self._projection_axis[1]

        assert axis_index != self._kstar_axis, "Project3Dto2D: Cannot project axis with same index as k* axis!"

        if axis_index == 0:
            se_axis = se3D.GetXaxis()
            me_axis = me3D.GetXaxis()
            if self._kstar_axis == 1:
                proj_option = "zy"
            elif self._kstar_axis == 2:
                proj_option = "yz"
            else:
                print("Project2Dto3D: Invalid k* axis selected. Exiting...")
                exit()
        elif axis_index == 1:
            se_axis = se3D.GetYaxis()
            me_axis = me3D.GetYaxis()
            if self._kstar_axis == 0:
                proj_option = "zx"
            elif self._kstar_axis == 2:
                proj_option = "xz"
            else:
                print("Project2Dto3D: Invalid k* axis selected. Exiting...")
                exit()
        elif axis_index == 2:
            se_axis = se3D.GetZaxis()
            me_axis = me3D.GetZaxis()
            if self._kstar_axis == 0:
                proj_option = "yx"
            elif self._kstar_axis == 1:
                proj_option = "xy"
            else:
                print("Project2Dto3D: Invalid k* axis selected. Exiting...")
                exit()
        else:
            print("Project2Dto3D: Axis index out of range for a 3D histogram. Exiting...")
            exit()

        # set the projection range
        perc_low = proj_range[0]
        perc_up  = proj_range[1]

        # get bins
        bin_low = se_axis.FindBin(perc_low)
        bin_up  = se_axis.FindBin(perc_up)

        # lower the upper bin number if it is right on the next bin's lower edge
        if perc_up == se_axis.GetBinLowEdge(bin_up):
            bin_up -= 1

        # adjust the axis ranges
        se_axis.SetRange(bin_low, bin_up)
        me_axis.SetRange(bin_low, bin_up)

        se = se3D.Project3D(proj_option)
        me = me3D.Project3D(proj_option)

        return [se, me]

    #move to analysis class
    def Reweight2D(self, se2d, me2d):
        """
        this function reweights the 2D histogram in each reweighting bin. 
        The output is a list [name, se, me_reweighted] for each bin.
        """
        label= se2d.GetYaxis().GetTitle() #TODO: check if this works!
        rew_range = self._rew_range
        histos = self.Get2DHistosFromBins(se2d, me2d, label)
        out = []

        for hist in histos:
            name = hist[0]
            se = hist[1]
            me = hist[2]
            out.append([name, se, AnalysisUtils.reweight(se, me, rew_range)[3].Clone(name)])        # append the reweighted th2 ME distribution

        return out


class AnalysisUtils():
    # # splits th2 in section based on provided bins -> Get2DHistosFromBins
    # def getBinRangeHistos(iSE, iME, bins, projection_axis):
    #     """
    #     This function splits 2D histograms in the ranges
    #     defined in the option 'bins'.

    #     The output is a list of ["range", SE, ME] for each bin:
    #         [[name, SE, ME], [bin2], ...]
    #     where the name is a string containing the limits.
    #     """
    #     yAxis = iSE.GetAxis(projection_axis)
    #     remaining_axis = 0 if projection_axis == 1 else 1
    #     xAxis = iSE.GetAxis(remaining_axis)

    #     if type(bins) == list:
    #         limits = []
    #         for value_diff in bins:
    #             value_bin = yAxis.FindBin(float(value_diff))
    #             limits.append(value_bin)
    #     else:
    #         print("Error in getBinRangeHistos: bin input \"" + str(bins) + "\" not a list of ranges!")
    #         exit()


    #     histos = []
    #     for n in range(1, len(bins)):
    #         diff_low = bins[n - 1]
    #         diff_up  = bins[n]

    #         bin_low = xAxis.FindBin(diff_low)
    #         bin_up  = xAxis.FindBin(diff_up)

    #         if diff_up == xAxis.GetBinLowEdge(bin_up):
    #             bin_up -= 1

    #         name = "[%.2f-%.2f)" % (bins[n - 1], bins[n])
    #         se = iSE.ProjectionX("se_k", bin_low, bin_up)
    #         me = iME.ProjectionX("me_k", bin_low, bin_up)
    #         histos.append([name, se.Clone(), me.Clone()])

    #     return histos

    # # splits th3 in section based on provided bins -> Project2Dto3D and Get2DHistosFromBins
    # def getBinRangeHistos3D(iSE, iME, diff3d, bins3d, mult_axis):
    #     """
    #     This function takes as input 3D SE and ME plots
    #     and splits them into 2D plots in mt/mult according to 'diff3d'
    #     in the ranges defined in 'bins3d'.

    #     The output is a list of ["range", SE mt/mult vs k*, ME mt/mult vs k*] for each bin:
    #         [[name, SE, ME], [bin2], ...]
    #     where the name is a string containing the limits and SE, ME are 2D plots.
    #     """
    #     if diff3d == 'mt':#TODO
    #         diffAxisSE = iSE.GetYaxis()
    #         diffAxisME = iME.GetYaxis()
    #         projOpt = "zx"
    #     elif diff3d == 'mult':
    #         if mult_axis == 1:
    #             diffAxisSE = iSE.GetYaxis()
    #             diffAxisME = iME.GetYaxis()
    #             projOpt = "zx"
    #         elif mult_axis == 2:
    #             diffAxisSE = iSE.GetZaxis()
    #             diffAxisME = iME.GetZaxis()
    #             projOpt = "yx"
    #         else:
    #             print("no multiplicity axis specified.")
    #     else:
    #         print("Error in getBinRangeHistos: diff3d axis not known. Please choose either 'mt' or 'mult'")
    #         exit()

    #     if type(bins3d) == list:
    #         limits = []
    #         for value_diff in bins3d:
    #             value_bin = diffAxisSE.FindBin(float(value_diff))
    #             limits.append(value_bin)
    #     else:
    #         print("Error in getBinRangeHistos: bin input \"" + str(bins3d) + "\" not a list of ranges!")
    #         exit()

    #     histos = []
    #     for n in range(1, len(limits)):

    #         diff_low = bins3d[n - 1]
    #         diff_up  = bins3d[n]
    #         print(diff_low, diff_up)

    #         bin_low = diffAxisSE.FindBin(diff_low)
    #         bin_up  = diffAxisSE.FindBin(diff_up)
    #         print(bin_low, bin_up)
    #         print("bin_up low edge:", diffAxisSE.GetBinLowEdge(bin_up))
    #         print("bin_up - 1 =", bin_up - 1)
    #         print("bin_up upper edge:", diffAxisSE.GetBinLowEdge(bin_up))

    #         if diff_up == diffAxisSE.GetBinLowEdge(bin_up):
    #             bin_up -= 1

    #         name = diff3d + ": [%.2f-%.2f)" % (bins3d[n - 1], bins3d[n])
    #         diffAxisSE.SetRange(bin_low, bin_up)
    #         diffAxisME.SetRange(bin_low, bin_up)
    #         se = iSE.Project3D("SE_"+projOpt+"_"+name)
    #         me = iME.Project3D("ME_"+projOpt+"_"+name)
    #         histos.append([name, se.Clone(), me.Clone()])

    #     return histos

    #keep
    def getCorrelation(se, me, name, conf, norm = None):
        """
        helper function for the correlation function
        outputs tuple [SE, ME, CF]
        """
        ch = CH.CorrelationHandler(name, se, me)

        ch.make_cf()
        minmax = norm if norm else [0.24, 0.34]
        
        ch.normalize_cf(minmax[0], minmax[1])
        se = ch.get_se().Clone("SE")
        me = ch.get_me().Clone("ME")
        cf = ch.get_cf().Clone("CF")
        se.SetTitle(conf)
        me.SetTitle(conf)
        cf.SetTitle(conf)
        del ch
        return [se, me, cf]

    #moved
    def getDifferential(iSE, iME, htype, bins, rebin, norm, mult_axis, mt_axis, title = None):
        """
        for 2d histogram
        returns [[iSE, iME], [se, me, cf]] for a list of mt or mult ranges
        [[iSE, iME], [se, me, cf, [rebin: [...], [...], ...], [bin 2 [rebin]], ...]
        """
        histos = []
        conf = "" if not title else title + " "         # append to given name
    
        #htype should be mult or mt
        assert htype in ["mult", "mt"], "getDifferential: no kmT or kmult input!"

        conf += "htype"+": "
        histos.append([iSE.Clone(f"SE k{htype}"), iME.Clone(f"ME k{htype}")]) 

        print(htype)
        print(bins)

        # divide in bins
        projection_axis = mult_axis if htype == "mult" else mt_axis
        mt_histos = AnalysisUtils.getBinRangeHistos(iSE, iME, bins, projection_axis)

        for n, [name, se, me] in enumerate(mt_histos, 1):
            histos.append(AnalysisUtils.getCorrelation(se, me, name, conf + name, norm))
            histos[n].append([])
            if rebin:       # append a list of rebinned [se, me, cf] in the original [se, me, cf, []]
                for factor in rebin:
                    se_rebin = AnalysisUtils.rebin_hist(se, factor)
                    me_rebin = AnalysisUtils.rebin_hist(me, factor)
                    rebin_conf = " rebin: " + str(factor)
                    histos[n][3].append(AnalysisUtils.getCorrelation(se_rebin, me_rebin, name, conf + name + rebin_conf, norm))
        return histos

    # # [[iSE, iME], [[1st proj SE, 1st proj ME], [se, me, cf, [rebin], [bin 2 [rebin]]]], ...]
    # def getDifferential3D(iSE, iME, diff3d, bins3d, bins, rebin, norm, mult_axis):
    #     """
    #     This function takes as input 3D mult-mt-k* plots
    #     and splits them first in mt/mult according to 'diff3d' in the limits defined in 'bins3d'.
    #     The 2D mt/mult-k* plots are then projected in k* in the limits defined in 'bins'.

    #     The output is a list of the mt bins, each bin is then a list of the mult bins
    #     which include the SE, ME, CF and the rebinned plots:
    #         [ [ [SE, ME, CF, [rebinned SE, ME, CF]], [mult/mt bin2], ... ], [mt/mult bin2], ... ]
    #     """
    #     histos = []
    #     histos.append([iSE.Clone("SE kmTmult"), iME.Clone("ME kmTmult")])

    #     diff3d_histos = AnalysisUtils.getBinRangeHistos3D(iSE, iME, diff3d, bins3d, mult_axis)

    #     htypeSplit2 = ""
    #     if diff3d == 'mult':
    #         htypeSplit2 = "mt"
    #     elif diff3d == 'mt':
    #         htypeSplit2 = "mult"

    #     for title, se, me in diff3d_histos:
    #         histos.append(AnalysisUtils.getDifferential(se, me, htypeSplit2, bins, rebin, norm, title))

    #     return histos

    # # [[iSE, iME], [[1st proj SE, 1st proj ME], [se, me, cf, [rebin], [bin 2 [rebin]]]], ...]
    # def getDiffReweight3D(iSE, iME, bins3d, bins, rebin, norm, rew_range):
    #     """
    #     This function takes as input 3D mult-mt-k* plots
    #     and splits them first in mt according to the limits defined in 'bins3d'.

    #     The 2D mult-k* plots are then reweighted in mult and splits according
    #     to the limits defined in 'bins'.

    #     The output is a list of the mt bins, each bin is then a list of the mult bins
    #     which include the SE, ME, CF and the rebinned plots:
    #         [ [ [SE, ME, CF, [rebinned SE, ME, CF]], [mult bin2], ... ], [mt bin2], ... ]
    #     """
    #     histos = []
    #     histos.append([iSE.Clone("SE kmTmult"), iME.Clone("ME kmTmult")])

    #     histos_unw = AnalysisUtils.getDifferential3D(iSE, iME, "mt", bins3d, bins, rebin, norm)

    #     histos_diff3d = AnalysisUtils.reweight3D(iSE, iME, bins3d, rew_range)

    #     for title, se, me in histos_diff3d:
    #         histos.append(AnalysisUtils.getDifferential(se, me, "mult", bins, rebin, norm, title))

    #     return histos, histos_unw

    # returns a list of [[iSE, iME], [se, me, cf]] for rel pair k* input
    # or reweights and returns ([[iSE, iME], [se, me, cf, [rebin]]], [me_unw, cf_unw, [rebin]]) for kmult
    # and [[iSE, iME], [se, me, cf, [rebin]]] for kmT
    def getIntegrated(iSE, iME, htype, rebin, norm, rew_range): #TODO
        histos = []
        histos_unw = []
        if htype == 'k':      # k* input
            histos.append([iSE.Clone("SE kstar"), iME.Clone("ME kstar")])
            se = iSE
            me = iME
        elif htype == 'mult':    # kmult input
            histos.append([iSE.Clone("SE kmult"), iME.Clone("ME kmult")])
            hReweight = AnalysisUtils.reweight(iSE, iME, rew_range)
            se          = hReweight[0]
            me          = hReweight[1]
            me_unw      = hReweight[2]
            se_mult     = hReweight[4]
            me_mult     = hReweight[5]
            me_mult_unw = hReweight[6]
            histos[0].append(se_mult.Clone("SE mult"))
            histos[0].append(me_mult.Clone("ME mult"))
            histos[0].append(me_mult_unw.Clone("ME mult unw"))
        elif htype == 'mt':    # kmT input
            histos.append([iSE.Clone("SE kmT"), iME.Clone("ME kmT")])
            hReweight = AnalysisUtils.reweight(iSE, iME, rew_range)
            se = hReweight[0]
            me = hReweight[2]   # unweighted me, i.e. normal me projection of the kmT histo

        histos.append(AnalysisUtils.getCorrelation(se, me, "cf", "", norm))
        if rebin:               # append rebinned histos to list of histos
            histos_rebin = []
            for factor in rebin:
                se_rebin = AnalysisUtils.rebin_hist(se, factor)
                me_rebin = AnalysisUtils.rebin_hist(me, factor)
                rebin_conf = " rebin: " + str(factor)
                histos_rebin.append(AnalysisUtils.getCorrelation(se_rebin, me_rebin, "rebin: " + str(factor), rebin_conf, norm))
            histos[1].append(histos_rebin)
        if htype == 'mult':      # 2nd list with unweighted histos
            se, me, cf = AnalysisUtils.getCorrelation(se, me_unw, "cf_unw", "unweighted", norm)
            histos_unw.append(me.Clone("ME unw"))
            histos_unw.append(cf.Clone("CF unw"))
            histos_unw.append([])
            if rebin:           # append rebinned histos to list of histos
                histos_rebin = []
                for factor in rebin:
                    se_rebin = AnalysisUtils.rebin_hist(se, factor)
                    me_rebin = AnalysisUtils.rebin_hist(me, factor)
                    rebin_conf = " rebin: " + str(factor)
                    se_rebin, me_rebin, cf_rebin = AnalysisUtils.getCorrelation(se_rebin, me_rebin, "rebin: " + str(factor), rebin_conf, norm)
                    histos_unw[2].append([me_rebin.Clone("ME unw"), cf_rebin.Clone("CF unw")])

        return histos, histos_unw


    # needed
    def reweight(iSE, iME, rew_range):
        """
        This function takes as input a 2D mt/mult vs k* SE and ME distribution
        and reweights the ME distribution in each bin projection of mt/mult.

        The output is a list that includes all plots that can be generated:
            [0] SE 1D k*
            [1] ME 1D k* reweighted
            [2] ME 1D k* unweighted
            [3] ME 2D mt/mult vs k* reweighted
            [4] SE 1D mt/mult
            [5] ME 1D mt/mult reweighted
            [6] ME 1D mt/mult unweighted
        """

        me = iME.Clone("ME kmult reweighted")
        me.Reset("ICESM")
        me_axis = me.GetYaxis()

        se_k = iSE.ProjectionX("se_k")
        me_k = iME.ProjectionX("me_k")

        if rew_range:
            int_min = se_k.FindBin(rew_range[0])
            int_max = se_k.FindBin(rew_range[1])
        else:
            int_min = 0
            int_max = se_k.GetNbinsX()

        se_int = se_k.Integral(int_min, int_max)
        me_int = me_k.Integral(int_min, int_max)

        se_mult = iSE.ProjectionY("se_mult")
        me_mult = iME.ProjectionY("me_mult")

        me_k_unw = iME.ProjectionX("me_k_unw")
        me_mult_unw = iME.ProjectionY("me_mult_unw")

        me_k.Reset("ICESM")
        me_mult.Reset("ICESM")

        # loop for the projection of each multiplicity slice
        for ybin in range(1, iSE.GetNbinsY() + 1):
            se_n = iSE.ProjectionX("se_bin", ybin, ybin)
            me_n = iME.ProjectionX("me_bin", ybin, ybin)

            se_int = 1
            me_int = 1

            if se_int:
                se_ratio = se_n.Integral(int_min, int_max) / se_int
            if me_int:
                me_ratio = me_n.Integral(int_min, int_max) / me_int

            if me_ratio > 0. and se_ratio > 0.:
                me_n.Scale(se_ratio / me_ratio)
                me_mult.SetBinContent(ybin, me_n.Integral(int_min, int_max))
                me_k.Add(me_n)
                for xbin in range(1, me_n.GetNbinsX() + 1):        # fill th2 reweighted ME
                    #me.Fill(me_n.GetBinContent(xbin), me_axis.GetBinCenter(ybin))
                    #me.Fill(me_axis.GetBinCenter(ybin), me_n.GetBinContent(xbin))
                    me.SetBinContent(xbin, ybin, me_n.GetBinContent(xbin))

        return [se_k, me_k, me_k_unw, me, se_mult, me_mult, me_mult_unw]

    
    # def reweight3D(iSE, iME, bins3d, rew_range):
    #     """
    #     This function takes as input a 3D SE and ME distribution
    #     and splits them in mt by the provided binning.
    #     The resulting 2D mult/k* plots are then reweighted.

    #     The output is a list for the individual mt bins with the name (mt limits), SE, ME:
    #         [[name, SE mult/k*, ME mult/k* reweighted], [bin 2], ...]
    #     """
    #     histos = AnalysisUtils.getBinRangeHistos3D(iSE, iME, "mt", bins3d)        # split th3 in mt and get a list of mult-k histos
    #     out = []

    #     for hist in histos:
    #         out.append([hist[0], hist[1], AnalysisUtils.reweight(hist[1], hist[2], rew_range)[3].Clone(hist[0])])        # append the reweighted th2 ME distribution

    #     return out

    # # 4d percentile histos
    # def getProj4d(iSE, iME, perc_range, axis_index):
    #     """
    #     Project a 4D histogram to 3 dimension
    #     --
    #     iSE: Same Event Distribution
    #     iME: Mixed Event Distribution
    #     perc_range: percentile range for projection
    #     axis_index: axis that should be projected
    #     """
    #     se4d = iSE.Clone("4d_perc_se")
    #     me4d = iME.Clone("4d_perc_me")

    #     axis = se4d.GetAxis(axis_index)

    #     # set the projection range
    #     perc_low = perc_range[0]
    #     perc_up  = perc_range[1]

    #     # get bins
    #     bin_low = axis.FindBin(perc_low)
    #     bin_up  = axis.FindBin(perc_up)

    #     # lower the upper bin number if it is right on the next bin's lower edge
    #     if perc_up == axis.GetBinLowEdge(bin_up):
    #         bin_up -= 1

    #     se4d.GetAxis(axis_index).SetRange(bin_low, bin_up)
    #     me4d.GetAxis(axis_index).SetRange(bin_low, bin_up)

    #     # find remaining axes
    #     dim = [0,1,2,3].remove(axis_index)
    #     a,b,c = dim[0], dim[1], dim[2]

    #     se = se4d.Projection(a, b, c)
    #     me = me4d.Projection(a, b, c)

    #     return [se, me]
     
    
    def rebin_hist(input_histo, binning):
        """
        returns rebinned copy of histo
        """
        histo = input_histo.Clone()
        histo = histo.Rebin(binning)
        return histo

