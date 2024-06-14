import ROOT
import unittest
import FemtoAnalysis as FA
import FemtoDreamReader as FDR
import os
import sys
import numpy as np

# NOTE: Update the config in run_UFFA method for more specific changes
TestingConfiguration = {
            "nbins" : 1000,
            "kstar_range" : (0, 0.5),
            "normalization_range" : (0.28, 0.4),
            "rebin_factors" : [2,5,10],
            "TDir" : "femto-dream-pair-task-track-track",
            "input_file" : "TestUFFA2.root",
            "precision" : 0.001,
            "precision_rebin" : 0.01
        }

class TestCF(unittest.TestCase):
    def _setUp(self, nbins, kstar_range, normalization_range):
        self.nbins=nbins
        self.kstar_min = kstar_range[0]
        self.kstar_max = kstar_range[1]
        self.normalize_min = normalization_range[0]
        self.normalize_max = normalization_range[1]


        # create histograms for UFFA input
        SEHisto = ROOT.TH1F("relPairDist", "Same Event Histogram", self.nbins, self.kstar_min, self.kstar_max)
        MEHisto = ROOT.TH1F("relPairDist", "Mixed Event Histogram", self.nbins, self.kstar_min, self.kstar_max)

        # functions that describe the SE/ME distributions
        SE_function = ROOT.TF1("SEgauss", "1000*(gaus(0)+1)", self.kstar_min, self.kstar_max)
        SE_function.SetParameter(0, 1)
        SE_function.SetParameter(1, 0.1)
        SE_function.SetParameter(2, 0.05)


        ME_function = ROOT.TF1("MEgauss", "1000000*(gaus(0)+1)", self.kstar_min, self.kstar_max)
        ME_function.SetParameter(0, 1)
        ME_function.SetParameter(1, 0.125)
        ME_function.SetParameter(2, 0.05)

        # fix bin edges
        self.bin_edges = np.linspace(self.kstar_min, self.kstar_max, self.nbins+1)
        self.bin_centers = (self.bin_edges[:-1] + self.bin_edges[1:]) / 2

        # get list of SE/ME distribution values
        self.SE_list = [int(SE_function.Eval(x)) for x in self.bin_centers]
        self.ME_list = [int(ME_function.Eval(x)) for x in self.bin_centers]

        # fill histograms for UFFA input
        for i in range(0, self.nbins, 1):
            SEHisto.SetBinContent(i+1, self.SE_list[i])
            MEHisto.SetBinContent(i+1, self.ME_list[i])

        # save UFFA input to file
        tempFile = ROOT.TFile(TestingConfiguration["input_file"], "RECREATE")
        folder = tempFile.mkdir(TestingConfiguration["TDir"])
        folder.cd()

        SE = folder.mkdir("SameEvent")
        SE.cd()
        SEHisto.Write()

        ME = folder.mkdir("MixedEvent")
        ME.cd()
        MEHisto.Write()
        tempFile.Close()

    def rebinDistribution(self, distribution_list, rebin_factor):
        # check if rebinning is possible
        assert len(distribution_list)%rebin_factor == 0, "Rebin not possible. Binnumber is not a multiple of the rebin factor!"

        nbins_rebinned = int(len(distribution_list)/rebin_factor)
        distribution_list_rebinned = []
        for i in range(0, nbins_rebinned, 1):
            new_bin_indices = [rebin_factor*i+j for j in range(0, rebin_factor, 1)]
            new_bin_content = sum(distribution_list[k] for k in new_bin_indices)
            distribution_list_rebinned.append(new_bin_content)

        return distribution_list_rebinned

    def CalculateCF(self, SE_distribution, ME_distribution):
        assert len(SE_distribution)==len(ME_distribution), "SE and ME distributions do not have same binnumber."

        CF_vals = [SE_distribution[i]/ME_distribution[i] for i in range(len(SE_distribution))]
        return CF_vals

    def update_bins(self, rebin_factor):
        self.nbins = int(TestingConfiguration["nbins"]/rebin_factor)
        self.bin_edges = np.linspace(self.kstar_min, self.kstar_max, self.nbins+1)
        self.bin_centers = (self.bin_edges[:-1] + self.bin_edges[1:]) / 2

    def normalizeCF(self, CF):
        integral_indizes = [i for i in range(len(CF)) if self.normalize_min <= self.bin_centers[i] < self.normalize_max]
        integral = sum([CF[i] for i in integral_indizes])*(self.bin_edges[1]-self.bin_edges[0])
        N = (self.normalize_max-self.normalize_min)/integral

        CF_normalized = [N*x for x in CF]
        return CF_normalized
            
    def tearDown(self):
        # # delete files
        # if os.path.exists("TestUFFA.root"):
        #     os.remove("TestUFFA.root")
        # if os.path.exists("UFFA_pd.root"):
        #     os.remove("UFFA_pd.root")
        pass

    def run_CFtest(self, CFhisto_UFFA, CFHisto_predicted, index, delta):
        self.assertAlmostEqual(
                CFhisto_UFFA.GetBinContent(index), 
                CFHisto_predicted.GetBinContent(index), 
                msg=f"CFs differ too much in bin {index}.\n CF-predicted: {CFHisto_predicted.GetBinContent(index)}, CF-UFFA: {CFhisto_UFFA.GetBinContent(index)}",
                delta=delta
            )
        
    def run_UFFA(self):
        DataDir = "./"
        InputFile = TestingConfiguration["input_file"]

        function = "cf"
        TDir=TestingConfiguration["TDir"]
        Normalize = TestingConfiguration["normalization_range"]
        Rebin = TestingConfiguration["rebin_factors"]
        mtBins = [0,4.5]
        multBins = [0,200] # integrate in charge track multiplicity
        multPercentileBins = [0,100] # multiplicity percentile
        #reweight=[0,0.4]

        OutputPath = "./" #"./UFFA_test/"
        # if not os.path.exists(OutputPath):
        #     os.mkdir(OutputPath)

        BaseName="pd"
        
        config = {
            "function":     'cf',
            "file":         DataDir + InputFile,
            "fileTDir":     TDir,
            "newfile":      "recreate",
            "outDir":       OutputPath,
            "rename":       BaseName,
            "atype":        "int",
            "htype":        "kstar",
            "bins3d":       mtBins,
            "bins":         multBins,
            "rebin":        Rebin,
            #"rewrange":     reweight,
            "percentile":   multPercentileBins,
            "normalize":    Normalize,
            "debug":        False
        }

        FA.UFFA(config)

    def test_CF(self):
        # setup
        self._setUp(nbins=TestingConfiguration["nbins"], kstar_range=TestingConfiguration["kstar_range"], normalization_range=TestingConfiguration["normalization_range"])
        
        self.run_UFFA()

        tempFile = ROOT.TFile.Open(TestingConfiguration["input_file"], "UPDATE")
        CF_Histo = ROOT.TH1F("CF", "CF normalized Histogram", self.nbins, self.kstar_min, self.kstar_max)

        # CF
        CF = self.CalculateCF(self.SE_list, self.ME_list)
        CF_normalized = self.normalizeCF(CF)
        for i in range(self.nbins):
            CF_Histo.SetBinContent(i+1, CF_normalized[i])
        
        # save for comparison
        self.CF_Histo = CF_Histo.Clone()
        self.CF_Histo.SetDirectory(0)
        CF_Histo.Write()
        tempFile.Close()

        UFFA_output_file = ROOT.TFile.Open("UFFA_pd.root", "READ")
        CF_UFFAHisto = UFFA_output_file.Get("femto-dream-pair-task-track-track_std").Get("CF")

        for i in range(1, self.nbins+1, 1):
            with self.subTest(i=i):
                self.run_CFtest(CF_UFFAHisto, self.CF_Histo, index=i, delta=TestingConfiguration["precision"])
    
    def test_rebinnedCF(self):
        # setup
        self._setUp(nbins=TestingConfiguration["nbins"], kstar_range=TestingConfiguration["kstar_range"], normalization_range=TestingConfiguration["normalization_range"])
        
        self.run_UFFA()

        for rebin_factor in TestingConfiguration["rebin_factors"]:
            SE_rebinned = self.rebinDistribution(self.SE_list, rebin_factor)
            ME_rebinned = self.rebinDistribution(self.ME_list, rebin_factor)

            # important: for rebinned we need to update our bins
            self.update_bins(rebin_factor)
            
            # ame as above
            tempFile = ROOT.TFile.Open(TestingConfiguration["input_file"], "UPDATE")
            CF_Histo = ROOT.TH1F(f"CF_rebin{rebin_factor}", f"CF normalized Histogram rebin_{rebin_factor}", self.nbins, self.kstar_min, self.kstar_max)

            # CF
            CF = self.CalculateCF(SE_rebinned, ME_rebinned)
            CF_normalized = self.normalizeCF(CF)


            for i in range(self.nbins):
                CF_Histo.SetBinContent(i+1, CF_normalized[i])

            self.CF_Histo = CF_Histo.Clone()
            self.CF_Histo.SetDirectory(0)
            CF_Histo.Write()
            tempFile.Close()

            UFFA_output_file = ROOT.TFile.Open("UFFA_pd.root", "READ")
            CF_UFFAHisto = UFFA_output_file.Get(f"femto-dream-pair-task-track-track_std/rebin_{rebin_factor}").Get("CF")

            for i in range(1, self.nbins+1, 1):
                with self.subTest(i=i):
                    self.run_CFtest(CF_UFFAHisto, self.CF_Histo, index=i, delta=TestingConfiguration["precision_rebin"])

        

