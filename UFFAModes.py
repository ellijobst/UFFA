import ROOT
import FemtoAnalysis as FA 
import FemtoDreamReader as FDR
import FemtoDreamSaver as FDS
import CombinedTemplateFit as TF
from FemtoAnalysis import Systematics

class UFFA():
    def UFFA(settings):
        conf = FA.config(settings)
        if conf['function'] == 'cf':
            UFFA.UFFA_cf(conf)
        elif conf['function'] == 'tf':
            UFFA.UFFA_tf(conf)
        elif conf['function'] == 'tf2d':
            UFFA.UFFA_tf2d(conf)
        elif conf['function'] == 'ctf':
            UFFA.UFFA_ctf(conf)
        elif conf['function'] == 'syst':
            if conf['htype'] in ['mtmult', 'rew3d', '4d', 'rew4d']:
                UFFA.UFFA_syst_3d(conf)
            else:
                UFFA.UFFA_syst(conf)

    # correlation function
    def UFFA_cf(settings):
        conf = FA.config(settings)
        fdr = FDR.FemtoDreamReader(conf['fullpath'], conf['fileTDir'])
        ch = FA.cf_handler(fdr, conf)
        fds = FDS.FemtoDreamSaver(conf, ch.get_histos())

    # template fits
    def UFFA_tf(settings):
        conf = FA.config(settings)
        if conf['file']:
            fdr1 = FDR.FemtoDreamReader(conf['fullpath'], conf['fileTDir'])
            dca_data = fdr1.get_dca()
        elif conf['data']:
            dca_data = conf['data']
        else:
            print('UFFA_tf: Missing input data!')
        if conf['templates']:
            if type(conf['templates']) == str:
                fdr2 = FDR.FemtoDreamReader(conf['templates'], conf['mcTDir'])
                dca_mcplots = fdr2.get_dca_mc()
            else:
                dca_mcplots = conf['templates']
        else:
            dca_mcplots = fdr1.get_dca_mc()

        fds = FDS.FemtoDreamSaver(settings)
        ofile = fds.getFile()

        TF.TemplateFit(ofile, dca_data, dca_mcplots, conf['tftype'], conf['namelist'], conf['fitrange'], conf['signalrange'], conf['bins'], conf['rebin'], conf['outDir'], conf['temp_init'], conf['temp_limits'], conf['temp_fraction'], conf['print'])

    # template fits 2d
    def UFFA_tf2d(settings):
        conf = FA.config(settings)
        dca_data = conf['data']
        dca_mcplots = conf['templates']

        fds = FDS.FemtoDreamSaver(settings)
        ofile = fds.getFile()

        TF.TemplateFit2D(ofile, dca_data, dca_mcplots, conf['namelist'], conf['fitrange'], conf['signalrange'], conf['bins'], conf['rebin'], conf['outDir'], conf['temp_init'], conf['temp_limits'], conf['temp_fraction'], conf['print'], conf['debug'])

    # combined template fits
    def UFFA_ctf(settings):
        conf = FA.config(settings)
        if conf['file']:
            fdr1 = FDR.FemtoDreamReader(conf['fullpath'], conf['fileTDir'])
            dca_data = fdr1.get_dca()
        elif conf['data']:
            dca_data = conf['data']
        else:
            print('UFFA_tf: Missing input data!')
        if conf['templates']:
            if type(conf['templates']) == str:
                fdr2 = FDR.FemtoDreamReader(conf['templates'], conf['mcTDir'])
                dca_mcplots = fdr2.get_dca_mc()
            else:
                dca_mcplots = conf['templates']
        else:
            dca_mcplots = fdr1.get_dca_mc()

        fds = FDS.FemtoDreamSaver(settings)
        ofile = fds.getFile()

        TF.CombinedFit(ofile, conf['outDir'], dca_data, dca_mcplots, conf['namelist'], conf['fitrange'], conf['signalrange'], conf['bins'], conf['rebin'], conf['temp_init'], conf['temp_limits'], conf['temp_fraction'], conf['print'])

    # systematics
    def UFFA_syst(settings):
        conf = FA.config(settings)
        fdr = FDR.FemtoDreamReader(conf['fullpath'], conf['fileTDir'])

        # default cf
        ch = FA.cf_handler(fdr, conf)
        cf, cf_unw = ch.get_cf()                                # [[cf, [rebins]], [bin2...], ...], [[cf unw, [rebins]], [bin2...], ...]

        # input same event for yield filtering
        if conf['yield']:
            se = fdr.get_se()
            pair_num_se = se.Integral(se.FindBin(0), se.FindBin(conf['yield'][0]))
        if conf['debug']:
            se_all = ch.get_se()

        cf_list = []
        if conf['rebin']:
            len_rebin = len(conf['rebin'])

        if conf['atype'] == 'int':                              # integrated
            ck, ck_rebin = cf[0]
            cf_list.append([ck, ck_rebin])
            syst = [[Systematics(ck), []]]                          # [[syst cf, [rebins]]]
            if conf['rebin']:
                for i in range(len_rebin):
                    syst[0][1].append(Systematics(ck_rebin[i]))
        elif conf['atype'] == 'dif':                            # differential
            syst = []
            for n, [ck, ck_rebin] in enumerate(cf):
                cf_list.append([ck, ck_rebin])
                syst.append([Systematics(ck), []])                  # [[syst cf, [rebins]], [bin2...], ...]
                if conf['rebin']:
                    for i in range(len_rebin):
                        syst[n][1].append(Systematics(ck_rebin[i]))

        # loop over data variations in file and calculate the cf for each
        # which is then saved in a th2 from which the systematic error is computed and saved in a th1
        file_dir = fdr.get_dir()
        fdr.cd(0)           # class method of FileSaver to return to root of file
        folders = fdr.get_folder_names()
        for folder in folders:
            fdr.cd(folder)

            # allows to include/exclude specific variations
            if conf['exclude'] and folder in conf['exclude']:
                continue
            elif conf['include']:
                if folder in conf['include']:
                    pass
                else:
                    continue
            elif folder.rsplit('_')[-1][:3] != "Var":
                continue

            ch_var = FA.cf_handler(fdr, conf)
            cf_var, cf_var_unw = ch_var.get_cf()

            if conf['debug']:
                print("Variation: \"" + folder + "\"")
            if conf['yield']:
                se_var = fdr.get_se()
                pair_num_var = se_var.Integral(se_var.FindBin(0), se_var.FindBin(conf['yield'][0]))
                deviation = abs(pair_num_se - pair_num_var) / pair_num_se
                if deviation > conf['yield'][1]:
                    if conf['debug']:
                        dev = deviation * 100
                        print("Integrated yield k*: [0, " + str(conf['yield'][0]) + ") differs by " + f"{dev:.1f} %")
                        if deviation > conf['yield'][1]:
                            print("Variation: Excluded!\n")
                            continue
            if conf['debug'] and conf['htype'] != 'k':
                se_var_all = ch_var.get_se()
                tab = '\t'
                print("Differential yield:")
                for n, bin1 in enumerate(se_var_all):
                    yield_all = se_all[n][0].Integral()
                    yield_all_var = se_var_all[n][0].Integral()
                    deviation = (abs(yield_all - yield_all_var) / yield_all) * 100
                    print(f"{tab}{conf['htype']:s}:  [{conf['bins'][n]:.2f}, {conf['bins'][n + 1]:.2f}) {tab} {deviation:5.2f} %")
                print()

            for n, [ck_var, ck_var_rebin] in enumerate(cf_var):
                syst[n][0].AddVar(ck_var)
                if conf['rebin']:
                    for i in range(len_rebin):
                        syst[n][1][i].AddVar(ck_var_rebin[i])
            del ch_var

        # generate th2 plots for systematics
        for n in range(len(syst)):
            syst[n][0].GenSyst()
            if conf['rebin']:
                for i in range(len_rebin):
                    syst[n][1][i].GenSyst()

        syst_plots = []                                         # [[[cf, diff, syst, dev], [rebins]], [bin2...], ...]
        for n in range(len(syst)):
            syst_plots.append([syst[n][0].GetAll(), []])
            if conf['rebin']:
                for i in range(len_rebin):
                    syst_plots[n][1].append(syst[n][1][i].GetAll())

        # generates the graphs with the systematic errors for the cf and the rebinned entries
        tgraphs = []
        for n, (hist, hist_rebin) in enumerate(cf_list):
            tgraphs.append([ROOT.TGraphErrors(), []])
            for i in range(1, hist.GetNbinsX() + 1):
                tgraphs[n][0].SetName("CF_syst_graph")
                tgraphs[n][0].SetPoint(i - 1, hist.GetBinCenter(i), hist.GetBinContent(i))
                tgraphs[n][0].SetPointError(i - 1, 0, syst_plots[n][0][2].GetBinContent(i))
            if conf['rebin']:
                for i in range(len_rebin):
                    tgraphs[n][1].append(ROOT.TGraphErrors())
                    for j in range(1, hist.GetNbinsX() + 1):
                        tgraphs[n][1][i].SetName("CF_syst_graph")
                        tgraphs[n][1][i].SetPoint(j - 1, hist_rebin[i].GetBinCenter(j), hist_rebin[i].GetBinContent(j))
                        tgraphs[n][1][i].SetPointError(j - 1, 0, syst_plots[n][1][i][2].GetBinContent(j))

        histos = (cf_list, syst_plots, tgraphs)
        fds = FDS.FemtoDreamSaver(conf, histos)

    # systematics
    def UFFA_syst_3d(settings):
        conf = FA.config(settings)
        fdr = FDR.FemtoDreamReader(conf['fullpath'], conf['fileTDir'])

        # default cf
        ch = FA.cf_handler(fdr, conf)
        histos = ch.get_cf_3d()                                # [[cf, [rebins]], [bin2...], ...], [[cf unw, [rebins]], [bin2...], ...]

        # input same event for yield filtering
        if conf['yield']:
            se = fdr.get_se()
            pair_num_se = se.Integral(se.FindBin(0), se.FindBin(conf['yield'][0]))
        if conf['debug']:
            se_all = ch.get_se_3d()

        syst = []
        syst_plots = []
        cf_raw = []

        if conf['rebin']:
            len_rebin = len(conf['rebin'])

        # create systematic object for all entries
        for n, bin1 in enumerate(histos):
            syst.append([])
            for nn, [cf, cf_rebin] in enumerate(bin1):
                syst[n].append([Systematics(cf), []])
                if conf['rebin']:
                    for nnn in range(len_rebin):
                        syst[n][nn][1].append(Systematics(cf_rebin[nnn]))

        # loop over data variations in file and calculate the cf for each
        # which is then saved in a th2 from which the systematic error is computed and saved in a th1
        file_dir = fdr.get_dir()
        fdr.cd(0)                               # class method of FileSaver to return to root of file
        folders = fdr.get_folder_names()
        folder_counter = -1
        for folder in folders:
            fdr.cd(folder)

            # include/exclude specific variations
            if conf['exclude'] and folder in conf['exclude']:
                continue
            elif conf['include']:
                if folder in conf['include']:
                    pass
                else:
                    continue
            elif folder.rsplit('_')[-1][:3] != "Var":
                continue

            ch_var = FA.cf_handler(fdr, conf)
            histos_var = ch_var.get_cf_3d()

            if conf['debug']:
                print("Variation: \"" + folder + "\"")
            # compare integrated yields in given range
            if conf['yield']:
                se_var = fdr.get_se()
                pair_num_var = se_var.Integral(se_var.FindBin(0), se_var.FindBin(conf['yield'][0]))
                deviation = abs(pair_num_se - pair_num_var) / pair_num_se
                if deviation > conf['yield'][1]:
                    if conf['debug']:
                        dev = deviation * 100
                        print("Integrated yield k*: [0, " + str(conf['yield'][0]) + ") differs by " + f"{dev:.1f} %")
                        if deviation > conf['yield'][1]:
                            print("Variation: Excluded!\n")
                            continue
            if conf['debug']:
                se_var_all = ch_var.get_se_3d()
                tab = '\t'
                for n, bin1 in enumerate(se_var_all):
                    print(f"Differential yield {conf['diff3d']:s}: [{conf['bins3d'][n]:.2f}, {conf['bins3d'][n + 1]:.2f})")
                    for nn, bin2 in enumerate(bin1):
                        yield_all = se_all[n][nn][0].Integral()
                        yield_all_var = se_var_all[n][nn][0].Integral()
                        deviation = (abs(yield_all - yield_all_var) / yield_all) * 100
                        print(f"{tab}{conf['diff3d2']:s}:  [{conf['bins'][nn]:.2f}, {conf['bins'][nn + 1]:.2f}) {tab} {deviation:5.2f} %")
                    print()
                if conf['interactive']:
                    option = input("Include [Y/n] ")
                    if option and option.lower()[0] == 'n':
                        print("\"" + folder + "\" excluded!\n")
                        continue
            folder_counter += 1

            cf_raw.append([])   # add entry for folder
            # add rebinned variations
            for n, bin1 in enumerate(histos_var):
                cf_raw[folder_counter].append([])
                for nn, [cf, cf_rebin] in enumerate(bin1):
                    cf_raw[folder_counter][n].append([cf.Clone("CF_" + folder.rsplit('_')[-1]), []])
                    syst[n][nn][0].AddVar(cf)
                    if conf['rebin']:
                        for nnn in range(len_rebin):
                            cf_raw[folder_counter][n][nn][1].append(cf_rebin[nnn].Clone("CF_" + folder.rsplit('_')[-1]))
                            syst[n][nn][1][nnn].AddVar(cf_rebin[nnn])
            del ch_var

        # generate th2 plots for systematics
        for n, bin1 in enumerate(syst):
            syst_plots.append([])
            for nn, bin2 in enumerate(bin1):
                syst[n][nn][0].GenSyst()
                syst_plots[n].append([syst[n][nn][0].GetAll(), []])
                if conf['rebin']:
                    for nnn in range(len_rebin):
                        syst[n][nn][1][nnn].GenSyst()
                        syst_plots[n][nn][1].append(syst[n][nn][1][nnn].GetAll())

        # generates the graphs with the systematic errors for the cf and the rebinned entries
        tgraphs = []
        for n, bin1 in enumerate(histos):
            tgraphs.append([])
            for nn, [hist, hist_rebin] in enumerate(bin1):
                tgraphs[n].append([ROOT.TGraphErrors(), []])
                for nnn in range(1, hist.GetNbinsX() + 1):
                    tgraphs[n][nn][0].SetName("CF syst graph")
                    tgraphs[n][nn][0].SetPoint(nnn - 1, hist.GetBinCenter(nnn), hist.GetBinContent(nnn))
                    tgraphs[n][nn][0].SetPointError(nnn - 1, 0, syst_plots[n][nn][0][2].GetBinContent(nnn))
                if conf['rebin']:
                    for nnn in range(len_rebin):
                        tgraphs[n][nn][1].append(ROOT.TGraphErrors())
                        for nnnn in range(1, hist.GetNbinsX() + 1):
                            tgraphs[n][nn][1][nnn].SetName("CF syst graph")
                            tgraphs[n][nn][1][nnn].SetPoint(nnnn - 1, hist_rebin[nnn].GetBinCenter(nnnn), hist_rebin[nnn].GetBinContent(nnnn))
                            tgraphs[n][nn][1][nnn].SetPointError(nnnn - 1, 0, syst_plots[n][nn][1][nnn][2].GetBinContent(nnnn))

        all_histos = (histos, syst_plots, tgraphs, cf_raw)
        fds = FDS.FemtoDreamSaver(conf, all_histos)
