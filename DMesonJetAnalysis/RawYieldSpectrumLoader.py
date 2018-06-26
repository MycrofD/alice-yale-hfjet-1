#!/usr/bin/env python
# load raw spectrum

import ROOT
import numpy
import DMesonJetFDCorrection
import DMesonJetUtils
import HistogramNormalizator

class RawYieldSpectrumLoader:

    def __init__(self, input_path=None, train=None, ana_name=None):
        self.fInputPath = input_path
        self.fTrain = train
        self.fAnalysisName = ana_name
        self.fJetPtBins = [5, 6, 8, 10, 14, 20, 30]
        self.fUseReflections = True
        self.fReflectionFit = "DoubleGaus"
        self.fReflectionRoS = 0
        self.fEvents = None
        self.fDMeson = None
        self.fJetType = None
        self.fJetRadius = None
        self.fVariableName = None
        self.fKinematicCuts = None
        self.fRawYieldMethod = None
        self.fDataSpectrumList = None
        self.fDataFile = None
        self.fMultiTrialSubDir = None
        self.fDataJetList = None
        self.fDataMesonList = None
        self.fFDConfig = None
        self.fTrigger = "AnyINT"
        self.fSuffix = ""

    def LoadSpectrumConfig(self, spectrum):
        iSpectrum = spectrum["name"].find("Spectrum") + 1
        if iSpectrum < 1:
            print("RawYieldSpectrumLoader.LoadSpectrumConfig: Failed to extract the variable name from spectrum name '{}'".format(spectrum["name"]))
            return
        self.fVariableName = spectrum["name"][0:iSpectrum - 1]

        if spectrum["type"] == "inv_mass_fit":
            self.fRawYieldMethod = "InvMassFit"
        elif spectrum["type"] == "side_band":
            self.fRawYieldMethod = "SideBand"
        else:
            print("RawYieldSpectrumLoader.LoadSpectrumConfig: Method '{}' for spectrum '{}' unknown".format(spectrum["type"], spectrum["name"]))
            return

        if len(spectrum["name"]) > iSpectrum + 8:  # check whether the length of the name is longer than 'var' + 'Spectrum'
            self.fKinematicCuts = spectrum["name"][iSpectrum + 8:]
        else:
            self.fKinematicCuts = ""

        self.fSuffix = spectrum["suffix"].replace(self.fRawYieldMethod, "")
        if self.fSuffix and self.fSuffix[0] == "_": self.fSuffix = self.fSuffix[1:]

    def LoadNumberOfEvents(self):
        if not self.fDataMesonList: self.LoadDataListFromDMesonJetAnalysis()
        histEvents = self.fDataMesonList.FindObject("histEvents")
        self.fEvents = histEvents.GetBinContent(histEvents.GetXaxis().FindBin("Normalized"))
        return self.fEvents

    def LoadDataFileFromDMesonJetAnalysis(self):
        self.fDataFile = ROOT.TFile("{0}/{1}/{2}.root".format(self.fInputPath, self.fTrain, self.fAnalysisName))
        return self.fDataFile

    def LoadDataListFromDMesonJetAnalysis(self):
        if not self.fDataFile: self.LoadDataFileFromDMesonJetAnalysis()
        if self.fVariableName:
            var = self.fVariableName.replace("z", "Z")
        else:
            print("No variable name provided!")
            exit(1)
        if self.fUseReflections:
            if not self.fReflectionRoS == 0:
                print("****Attention Attention Attention****")
                print("You asked for reflections with RoS != 0. This is not available in DMesonJetAnalysis!")
                exit(1)
            refl = self.fReflectionFit
        else:
            refl = None
        spectrumName = "{}Spectrum".format(var)
        if self.fDMeson:
            if self.fTrigger:
                dataMesonListName = "{}_{}".format(self.fTrigger, self.fDMeson)
            else:
                dataMesonListName = self.fDMeson
            self.fDataMesonList = self.fDataFile.Get(dataMesonListName)
            if not self.fDataMesonList:
                print("Could not find list {0} in file {1}". format(dataMesonListName, self.fDataFile.GetName()))
                exit(1)
            if self.fJetType and self.fJetRadius:
                dataJetListName = "_".join([self.fJetType, self.fJetRadius])
                self.fDataJetList = self.fDataMesonList.FindObject(dataJetListName)
                if not self.fDataJetList:
                    print("Could not find list {0}/{1} in file {2}". format(self.fDMeson, dataJetListName, self.fDataFile.GetName()))
                    self.fDataMesonList.Print()
                    exit(1)
                if var:
                    dataListName = "_".join([s for s in [self.fDMeson, dataJetListName, spectrumName, self.fKinematicCuts, self.fRawYieldMethod, refl] if s])
                    self.fDataSpectrumList = self.fDataJetList.FindObject(dataListName)
                    if not self.fDataSpectrumList:
                        print("Could not find list {0}/{1}/{2} in file {3}". format(self.fDMeson, dataJetListName, dataListName, self.fDataFile.GetName()))
                        self.fDataJetList.Print()
                        exit(1)
            if var and not (self.fJetType and self.fJetRadius):
                dataListName = "_".join([s for s in [self.fDMeson, spectrumName, self.fKinematicCuts, self.fRawYieldMethod, refl] if s])
                self.fDataSpectrumList = self.fDataMesonList.FindObject(dataListName)
                if not self.fDataSpectrumList:
                    print("Could not find list {0}/{1} in file {2}". format(self.fDMeson, dataListName, self.fDataFile.GetName()))
                    self.fDataMesonList.Print()
                    exit(1)
        return self.fDataSpectrumList

    def GetDefaultSpectrumFromDMesonJetAnalysis(self, fd, fd_error_band, ry_error_band):
        if self.fVariableName:
            var = self.fVariableName.replace("z", "Z")
        else:
            print("No variable name provided!")
            exit(1)
        if self.fUseReflections:
            if not self.fReflectionRoS == 0:
                print("****Attention Attention Attention****")
                print("You asked for reflections with RoS != 0. This is not available in DMesonJetAnalysis!")
                exit(1)
            refl = self.fReflectionFit
        else:
            refl = None
        if not self.fDataSpectrumList: self.LoadDataListFromDMesonJetAnalysis()
        spectrumName = "{}Spectrum".format(var)
        inputSpectrumName = "_".join([s for s in [self.fDMeson, self.fJetType, self.fJetRadius, spectrumName, self.fKinematicCuts, self.fRawYieldMethod, refl] if s])
        h_orig = self.fDataSpectrumList.FindObject(inputSpectrumName)
        if not h_orig:
            print("Could not find histogram {0} in list {1}". format(inputSpectrumName, self.fDataSpectrumList.GetName()))
            self.fDataSpectrumList.Print()
            exit(1)
        h = h_orig.Clone(h_orig.GetName())
        if var == "JetPt":
            h.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})")
        elif var == "JetZ":
            h.GetXaxis().SetTitle("#it{z}_{||}^{ch jet}")
        h.GetYaxis().SetTitle("raw yield")

        if fd: self.ApplyFDCorrection(h, fd_error_band)
        if ry_error_band <> 0: self.ApplyRawYieldSystFromMultiTrial(h, ry_error_band)

        return h

    def GetDefaultSpectrumInvMassFitFromMultiTrial(self):
        if not self.fVariableName == "JetPt":
            print("Fatal error: inv. mass fit multi-trial available only for jet pt. Requested for: {}".format(self.fVariableName))
            exit(1)
        else:
            var = self.fVariableName
        spectrumName = "{}Spectrum".format(var)
        inputSpectrumName = "_".join([s for s in [self.fDMeson[3:], spectrumName, self.fKinematicCuts, self.fRawYieldMethod, self.fSuffix] if s])

        if self.fUseReflections:
            self.fMultiTrialSubDir = "RawYieldUnc_refl_{0}".format(self.fReflectionFit)
            if not self.fReflectionRoS == 0: self.fMultiTrialSubDir += "_{0}".format(self.fReflectionRoS)
        else:
            self.fMultiTrialSubDir = "RawYieldUnc"

        h = ROOT.TH1D("TrialExpoFreeSigmaInvMassFit", "Trial Expo Free Sigma Inv. Mass Fit", len(self.fJetPtBins) - 1, numpy.array(self.fJetPtBins, dtype=numpy.float64))
        for ibin, (ptmin, ptmax) in enumerate(zip(self.fJetPtBins[:-1], self.fJetPtBins[1:])):
            fname = "{}/{}/{}/{}/{}_RawYieldVariations_Dzero_InvMassFit_{:.1f}to{:.1f}.root".format(self.fInputPath, self.fTrain, self.fAnalysisName, self.fMultiTrialSubDir, inputSpectrumName, ptmin, ptmax)
            file = ROOT.TFile(fname)
            if not file or file.IsZombie():
                print("Could not open file {0}".format(fname))
                file.ls()
                exit(1)
            else:
                print("File {0} open successfully".format(fname))
            hry = file.Get("hRawYieldTrialExpoFreeS")
            if not hry:
                print("Could not find histogram {0} in file {1}".format("hRawYieldTrialExpoFreeS", fname))
                file.ls()
                exit(1)
            ry = hry.GetBinContent(1)  # this should be with the normal inv mass binning
            ery = hry.GetBinError(1)
            h.SetBinContent(ibin + 1, ry)
            h.SetBinError(ibin + 1, ery)
        h.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})")
        h.GetYaxis().SetTitle("raw yield")
        return h

    def GetDefaultSpectrumSideBandFromMultiTrial(self):
        if self.fVariableName:
            var = self.fVariableName.replace("Z", "z")
        else:
            print("No variable name provided!")
            exit(1)
        spectrumName = "{}Spectrum".format(var.replace("z", "Z"))
        inputSpectrumName = "_".join([s for s in [self.fDMeson[3:], spectrumName, self.fKinematicCuts, self.fRawYieldMethod, self.fSuffix] if s])

        if self.fUseReflections:
            self.fMultiTrialSubDir = "RawYieldUnc_refl_{0}".format(self.fReflectionFit)
            if not self.fReflectionRoS == 0: self.fMultiTrialSubDir += "_{0}".format(self.fReflectionRoS)
        else:
            self.fMultiTrialSubDir = "RawYieldUnc"

        fname = "{input_path}/{train}/{name}/{subdir}/{spectrum_name}_TrialExpoFreeS_{var}_Dzero_SideBand.root".format(input_path=self.fInputPath, train=self.fTrain, name=self.fAnalysisName, subdir=self.fMultiTrialSubDir, var=var, spectrum_name=inputSpectrumName)
        file = ROOT.TFile(fname)
        if not file or file.IsZombie():
            print("Could not open file {0}".format(fname))
            file.ls()
            exit(1)
        else:
            print("File {0} open successfully".format(fname))
        hname = "f{var}SpectrSBDef".format(var=var)
        h = file.Get(hname)
        if not h:
            print("Could not find histogram {0} in file {1}".format(hname, fname))
            file.ls()
            exit(1)
        h_copy = h.Clone("TrialExpoFreeSigmaSideBand")
        h_copy.SetTitle("Trial Expo Free Sigma Inv. Mass Fit")
        if var == "JetPt":
            h_copy.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})")
        elif var == "Jetz":
            h_copy.GetXaxis().SetTitle("#it{z}_{||}^{ch jet}")
        h_copy.GetYaxis().SetTitle("raw yield")
        return h_copy

    def GetDefaultSpectrumFromMultiTrial(self, fd, fd_error_band, ry_error_band):
        if self.fRawYieldMethod == "SideBand": h = self.GetDefaultSpectrumSideBandFromMultiTrial()
        elif self.fRawYieldMethod == "InvMassFit": h = self.GetDefaultSpectrumInvMassFitFromMultiTrial()
        else:
            print("Method {0} unknown!".format(self.fRawYieldMethod))
            exit(1)
        if fd: self.ApplyFDCorrection(h, fd_error_band)
        if ry_error_band <> 0: self.ApplyRawYieldSystFromMultiTrial(h, ry_error_band)
        return h

    def GetAverageSpectrumFromMultiTrial(self, fd, fd_error_band):
        if self.fVariableName:
            var = self.fVariableName.replace("Z", "z")
        else:
            print("No variable name provided!")
            exit(1)
        spectrumName = "{}Spectrum".format(var)
        inputSpectrumName = "_".join([s for s in [self.fDMeson[3:], spectrumName, self.fKinematicCuts, self.fRawYieldMethod, self.fSuffix] if s])

        if self.fRawYieldMethod == "InvMassFit" and not var == "JetPt":
            print("Fatal error: multi-trial available only for jet pt. Requested for: {}".format(self.fVariableName))
            exit(1)
        if self.fUseReflections:
            self.fMultiTrialSubDir = "RawYieldUnc_refl_{0}".format(self.fReflectionFit)
            if not self.fReflectionRoS == 0: self.fMultiTrialSubDir += "_{0}".format(self.fReflectionRoS)
        else:
            self.fMultiTrialSubDir = "RawYieldUnc"

        fname = "{input_path}/{train}/{name}/{subdir}/{spectrum_name}_FinalRawYieldCentralPlusSystUncertainty_{var}_Dzero_{meth}.root".format(input_path=self.fInputPath, train=self.fTrain, name=self.fAnalysisName, subdir=self.fMultiTrialSubDir, var=var, meth=self.fRawYieldMethod, spectrum_name=inputSpectrumName)
        file = ROOT.TFile(fname)
        if not file or file.IsZombie():
            print("Could not open file {0}".format(fname))
            file.ls()
            exit(1)
        hname = "{var}RawYieldCentral".format(var=var)
        h = file.Get(hname)
        if not h:
            print("Could not find histogram {0} in file {1}".format(hname, fname))
            file.ls()
            exit(1)
        h_copy = h.Clone("{0}_copy".format(h.GetName()))
        if fd: self.ApplyFDCorrection(h_copy, fd_error_band)
        return h_copy

    def ApplyFDCorrection(self, h, error_band):
        print("Applying FD correction with error band point: {0} (0 = central point, +/- = upper/lower error band)".format(error_band))
        if not isinstance(error_band, basestring): error_band = -error_band
        fdHist = self.GetFDCorrection(error_band)
        h.Add(fdHist, -1)

    def GetFDCorrection(self, error_band=0):
        if not self.fFDConfig:
            print("ERROR!!!!!!!!!")
            print("The FD correction configuration was not provided.")
            exit(1)

        if self.fVariableName:
            var = self.fVariableName.replace("z", "Z")
        else:
            print("No variable name provided!")
            exit(1)
        spectrumName = "{}Spectrum".format(var)
        if self.fKinematicCuts: spectrumName += "_{0}".format(self.fKinematicCuts)
        fdCorrection = DMesonJetFDCorrection.DMesonJetFDCorrection(self.fFDConfig, spectrumName, self.fInputPath, "D0", "Charged", "R040")
        fdHist = fdCorrection.fFDHistogram.Clone("FD")
        fdSyst = dict()
        fdSyst[1] = fdCorrection.fFDUpSystUncHistogram.Clone("FDUpSystUnc")
        fdSyst[-1] = fdCorrection.fFDLowSystUncHistogram.Clone("FDLowSystUnc")
        fdSystGraph = fdCorrection.fFDSystUncGraph.Clone("FDSystUnc")
        fdSyst["tot_up"] = fdCorrection.fFDTotUpSystUncHistogram.Clone("FDTotUpSystUnc")
        fdSyst["tot_low"] = fdCorrection.fFDTotLowSystUncHistogram.Clone("FDTotLowSystUnc")
        fdTotSystGraph = fdCorrection.fFDTotSystUncGraph.Clone("FDTotSystUnc")

        if self.fEvents is None: self.LoadNumberOfEvents()
        norm = self.fEvents / HistogramNormalizator.CROSS_SECTION * HistogramNormalizator.BRANCHING_RATIO

        if isinstance(error_band, basestring) and "graph" in error_band:
        # a graph was requested
            if "tot" in error_band:
                graph = fdTotSystGraph
            else:
                graph = fdSystGraph
            for ibin in xrange(0, graph.GetN()):
                graph.SetPoint(ibin, graph.GetX()[ibin], graph.GetY()[ibin] * norm * fdHist.GetXaxis().GetBinWidth(ibin + 1))
                graph.SetPointEYlow(ibin, graph.GetErrorYlow(ibin) * norm * fdHist.GetXaxis().GetBinWidth(ibin + 1))
                graph.SetPointEYhigh(ibin, graph.GetErrorYhigh(ibin) * norm * fdHist.GetXaxis().GetBinWidth(ibin + 1))
            return graph
        elif error_band in fdSyst:
        # error band requested: can be the partial systematic, including only the POWHEG scale and quark mass variations,
        # or the total systematic which includes also other uncertainties added in quadrature, such as PDF variation
            if error_band == -1 or error_band == "tot_low": s = -1
            else: s = 1
            fdHist.Add(fdSyst[error_band], s)
        elif error_band == 0:
        # central points were requested
            pass
        else:
        # a specific variation was requested
            fdHist = fdCorrection.GetFDHistogram(error_band)

        for ibin in xrange(1, fdHist.GetNbinsX() + 1):
            fdHist.SetBinContent(ibin, fdHist.GetBinContent(ibin) * norm * fdHist.GetXaxis().GetBinWidth(ibin))
            fdHist.SetBinError(ibin, fdHist.GetBinError(ibin) * norm * fdHist.GetXaxis().GetBinWidth(ibin))

        return fdHist

    def GetRawYieldSystFromMultiTrial(self):
        if self.fVariableName:
            var = self.fVariableName.replace("Z", "z")
        else:
            print("No variable name provided!")
            exit(1)
        spectrumName = "{}Spectrum".format(var)
        inputSpectrumName = "_".join([s for s in [self.fDMeson[3:], spectrumName, self.fKinematicCuts, self.fRawYieldMethod, self.fSuffix] if s])

        if self.fUseReflections:
            self.fMultiTrialSubDir = "RawYieldUnc_refl_{0}".format(self.fReflectionFit)
            if not self.fReflectionRoS == 0: self.fMultiTrialSubDir += "_{0}".format(self.fReflectionRoS)
        else:
            self.fMultiTrialSubDir = "RawYieldUnc"
        fname = "{input_path}/{train}/{name}/{subdir}/{spectrum_name}_FinalRawYieldUncertainty_{var}_Dzero_{meth}.root".format(input_path=self.fInputPath, train=self.fTrain, name=self.fAnalysisName, subdir=self.fMultiTrialSubDir, var=var, meth=self.fRawYieldMethod, spectrum_name=inputSpectrumName)
        file = ROOT.TFile(fname)
        if not file or file.IsZombie():
            print("Could not open file {0}".format(fname))
            file.ls()
            exit(1)
        rySyst = file.Get("{}RawYieldUncert".format(var))
        if not rySyst:
            print("Could not find histogram {0} in file {1}".format("JetRawYieldUncert", fname))
            file.ls()
            exit(1)
        return rySyst

    def ApplyRawYieldSystFromMultiTrial(self, h, error_band):
        print("Applying raw yield systematic with error band point: {0} (0 = central point, +/- = upper/lower error band)".format(error_band))
        rySyst = self.GetRawYieldSystFromMultiTrial()
        h.Add(rySyst, error_band)
