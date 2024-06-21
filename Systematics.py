import ROOT

class Systematics():
    """
    class that returns the systematics of a CF
    ---
    - add variations with AddVar(var) before calling GenSyst()
    - GetAll() : returns [th2 cf, th2 difference, th1 systematics, th1 std dev]
    """
    counter = 0 #TODO: put limits in config and not here!
    ybins = 1200
    def __init__(self, cf):
        self._cf = cf
        self._xaxis = cf.GetXaxis()
        self._xbins = cf.GetNbinsX()

        self._var = ROOT.TH2D("CF_th2", "CF_th2", self._xbins, self._xaxis.GetXmin(), self._xaxis.GetXmax(), Systematics.ybins, 0, 10)
        self._dif = ROOT.TH2D("diff_th2", "diff_th2", self._xbins, self._xaxis.GetXmin(), self._xaxis.GetXmax(), Systematics.ybins, -5, 5)
        self._sys = ROOT.TH1D("syst_th1", "syst_th1", self._xbins, self._xaxis.GetXmin(), self._xaxis.GetXmax())
        self._dev = ROOT.TH1D("dev_th1", "dev_th1", self._xbins, self._xaxis.GetXmin(), self._xaxis.GetXmax())
        Systematics.counter = Systematics.counter + 1

    def AddVar(self, cf_var):
        for i in range(1, self._xbins + 1):
            if self._cf.GetBinCenter(i) > 3:    # break over 3GeV
                break
            self._var.Fill(cf_var.GetBinCenter(i), cf_var.GetBinContent(i))     # fill th2 cf histo with variation

    def GenSyst(self):
        for i in range(1, self._xbins + 1):
            if self._cf.GetBinCenter(i) > 3:    # break over 3GeV
                break
            cf_proj = self._var.ProjectionY("cf_xbin" + str(i), i, i)
            dev = cf_proj.GetStdDev()
            self._dev.SetBinContent(i, dev)
            cont_def = self._cf.GetBinContent(i)
            var_min = cf_proj.GetBinCenter(cf_proj.FindFirstBinAbove(0))
            var_max = cf_proj.GetBinCenter(cf_proj.FindLastBinAbove(0))
            self._dif.SetBinContent(i, self._dif.GetYaxis().FindBin(var_min), 1)
            self._dif.SetBinContent(i, self._dif.GetYaxis().FindBin(var_max), 1)
            dif_proj = self._dif.ProjectionY("diff_xbin" + str(i), i, i)
            proj_min = dif_proj.GetBinCenter(dif_proj.FindFirstBinAbove(0))     # minimum value of difference
            proj_max = dif_proj.GetBinCenter(dif_proj.FindLastBinAbove(0))      # maximum value of difference
            self._sys.SetBinContent(i, (proj_max - proj_min) / (12**0.5))       # assume a square distribution
        self._var.SetDirectory(0)
        self._dif.SetDirectory(0)
        self._sys.SetDirectory(0)
        self._dev.SetDirectory(0)

    def SetBinning(self, n):
        Systematics.ybins = n

    def GetVar(self):
        return self._var

    def GetDiff(self):
        return self._dif

    def GetSyst(self):
        return self._sys

    def GetDev(self):
        return self._dev

    def GetAll(self):
        return [self._var, self._dif, self._sys, self._dev]
