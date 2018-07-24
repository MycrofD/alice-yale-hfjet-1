#!/usr/bin/env python
# python program to perform a D meson jet unfolding

import math
import os
import collections
import array
import numpy
import ROOT

import DMesonJetUtils
import DMesonJetCompare
import UnfoldingResponseMatrix
import RawYieldSpectrumLoader

globalList = []

class DMesonJetUnfoldingEngine:

    def __init__(self, input_path, dataTrain, analysisName, config, fd_error_band=0, ry_error_band=0):
        # input path
        self.fInputPath = input_path

        # response matrix
        self.fDMesonResponse = config["d_meson_response"]
        self.fSpectrumResponseName = config["spectrum_response"]
        if "d_meson_truth" in config: self.fDMesonTruth = config["d_meson_truth"]
        else: self.fDMesonTruth = None
        if "spectrum_truth" in config: self.fSpectrumTruthName = config["spectrum_truth"]
        else: self.fSpectrumTruthName = None

        # data
        self.fDataTrain = dataTrain
        self.fDataList = None
        self.fAnalysisName = analysisName
        self.fDMeson = config["d_meson"]
        self.fJetType = config["jet_type"]
        self.fJetRadius = config["jet_radius"]
        if "jet_reco_scheme" in config:
            self.fJetRecoScheme = config["jet_reco_scheme"]
        else:
            self.fJetRecoScheme = "pt_scheme"
        self.fVariableName = config["variable"]
        self.fSpectrumName = "{}Spectrum".format(self.fVariableName)
        self.fKinematicCuts = config["kinematic_cuts"]
        self.fRawYieldMethod = config["raw_yield_method"]
        self.fFullSpectrumName = "_".join([s for s in [self.fSpectrumName, self.fKinematicCuts, self.fRawYieldMethod] if s])
        self.fUseMultiTrialInput = True
        if "FD" in config:
            self.fFDConfig = config["FD"]
        else:
            self.fFDConfig = None
            if fd_error_band != 0:
                print("Error: no FD configuration provided, but FD error band requested!")
                exit(1)

        if "rebin_input" in config:
            self.fRebinInput = config["rebin_input"]
        else:
            self.fRebinInput = False

        # feed-down systematic uncertainty
        self.fFDErrorBand = fd_error_band

        # raw yield extraction systematic uncertainty
        self.fRYErrorBand = ry_error_band

        # use raw yield with reflections
        self.fUseReflections = True

        # reflection systematic uncertainties
        self.fReflectionFit = "DoubleGaus"
        self.fReflectionRoS = 0

        # unfolding settings
        self.fPriors = config["priors"]
        self.fUnfoldingConfig = config["methods"]
        self.fDefaultPrior = config["default_prior"]
        self.fDefaultMethod = config["default_method"]
        self.fDoErrorCrossChecks = True

        # name
        self.fName = config["name"]
        if isinstance(self.fFDErrorBand, basestring):
            self.fName += "_FD{}".format(self.fFDErrorBand)
        elif self.fFDErrorBand > 0:
            self.fName += "_FDUpperBand"
        elif self.fFDErrorBand < 0:
            self.fName += "_FDLowerBand"
        if self.fRYErrorBand > 0:
            self.fName += "_RYUpperBand"
        elif self.fRYErrorBand < 0:
            self.fName += "_RYLowerBand"

        self.fNumberOfEvents = None

        # self unfolding option (for debugging)
        if self.fName == "__self_unfolding__":
            self.fNumberOfEvents = config["self_unfolding_events"]
            self.fEvents = ROOT.TH1D("Events", "Events", 1, 0, 1)
            self.fEvents.SetBinContent(1, self.fNumberOfEvents)
            self.fBins = config["self_unfolding_bins"]
            self.fSelfUnfoldingScaling = config["self_unfolding_scaling"]

        # variables computed during execution
        self.fResponseMatrices = collections.OrderedDict()
        self.fUnfoldedSpectra = collections.OrderedDict()
        self.fUnfoldedSpectraErrors = collections.OrderedDict()
        self.fRefoldedSpectra = collections.OrderedDict()
        self.fCovarianceMatrices = collections.OrderedDict()
        self.fSvdDvectors = collections.OrderedDict()
        self.fPearsonMatrices = collections.OrderedDict()
        self.fBinByBinCorrectionFactors = collections.OrderedDict()
        self.fCustomPriors = dict()
        self.fCanvases = []

    def SavePlots(self, path, format):
        for c in self.fCanvases:
            c.SaveAs("{0}/{1}.{2}".format(path, c.GetName(), format))
            if "Ratio" in c.GetName() and ("UnfoldingRegularization" in c.GetName() or "UnfoldingMethod" in c.GetName() or "UnfoldingPrior" in c.GetName()):
                if c.GetUymax() > 2:
                    for obj in c.GetListOfPrimitives():
                        if isinstance(obj, ROOT.TH1):
                            obj.SetMaximum(2)
                            c.SaveAs("{0}/{1}_zoom.{2}".format(path, c.GetName(), format))
                            break

    def GenerateRootList(self):
        rlist = ROOT.TList()
        rlist.SetName(self.fName)
        rlist.Add(self.fEvents)
        rlist.Add(self.fInputSpectrum)
        if self.fTruthSpectrum:
            rlist.Add(self.fTruthSpectrum)
        rlist.Add(self.fDetectorResponse)
        rlist.Add(self.fDetectorTrainTruth)
        for r in self.fResponseMatrices.itervalues():
            rlist.Add(r.GenerateRootList())
        unfoldingHistos = dict()
        for method in self.fUnfoldingConfig.iterkeys():
            rMethodList = ROOT.TList()
            rMethodList.SetName(method)
            rlist.Add(rMethodList)
            unfoldingHistos[method] = rMethodList
        for (method, reg, prior), h in self.fUnfoldedSpectra.iteritems():
            unfoldingHistos[method].Add(h)
        for (method, reg, prior), h in self.fRefoldedSpectra.iteritems():
            unfoldingHistos[method].Add(h)
        for (method, reg, prior), h in self.fCovarianceMatrices.iteritems():
            unfoldingHistos[method].Add(h)
        for h in self.fSvdDvectors.itervalues():
            unfoldingHistos["Svd"].Add(h)
        for (method, reg, prior), h in self.fPearsonMatrices.iteritems():
            unfoldingHistos[method].Add(h)
        for (method, reg, prior), (h1, h2, h3) in self.fUnfoldedSpectraErrors.iteritems():
            unfoldingHistos[method].Add(h1)
            if h2: unfoldingHistos[method].Add(h2)
            if h3: unfoldingHistos[method].Add(h3)
        for h in self.fBinByBinCorrectionFactors.itervalues():
            unfoldingHistos["BinByBin"].Add(h)
        return rlist

    def LoadData(self, dataFile, responseFile, eff, use_overflow):
        self.fUseOverflow = use_overflow
        self.fDataFile = dataFile
        self.fResponseFile = responseFile
        print("Now loading data")
        res = self.LoadResponse(eff)
        if not res: return False
        res = self.LoadInputSpectrum()
        if not res: return False
        res = self.LoadTruthSpectrum()
        return res

    def LoadResponse(self, eff):
        responseListName = "{}_Jet_AKT{}{}_{}_{}".format(self.fDMesonResponse, self.fJetType, self.fJetRadius, self.fJetRecoScheme, self.fSpectrumResponseName)
        responseList = self.fResponseFile.Get(responseListName)
        if not responseList:
            print("Could not find list {0} in file {1}". format(responseListName, self.fResponseFile.GetName()))
            exit(1)
        detectorResponseName = "{}_Jet_AKT{}{}_{}_{}_DetectorResponse".format(self.fDMesonResponse, self.fJetType, self.fJetRadius, self.fJetRecoScheme, self.fSpectrumResponseName)
        detectorResponse = responseList.FindObject(detectorResponseName)
        if not detectorResponse:
            print("Could not find histogram {0} in list {1} in file {2}". format(detectorResponseName, responseListName, self.fResponseFile.GetName()))
            exit(1)
        self.fDetectorResponse = detectorResponse.Clone("{0}_DetectorResponse".format(self.fName))
        self.fDetectorResponse.SetTitle("{0} Detector Response".format(self.fName))

        if eff:
            detTrainTruthName = "{}_Jet_AKT{}{}_{}_{}_Truth".format(self.fDMesonResponse, self.fJetType, self.fJetRadius, self.fJetRecoScheme, self.fSpectrumResponseName)
            detTrainTruth = responseList.FindObject(detTrainTruthName)
            if not detTrainTruth:
                print("Could not find histogram {0} in list {1} in file {2}". format(detTrainTruthName, responseListName, self.fResponseFile.GetName()))
                exit(1)
            self.fDetectorTrainTruth = detTrainTruth.Clone("{0}_ResponseTruth".format(self.fName))
        else:
            self.fDetectorTrainTruth = self.fDetectorResponse.ProjectionY("{0}_ResponseTruth".format(self.fName), 0, -1)
        self.fDetectorTrainTruth.SetTitle("{0} Response Truth".format(self.fName))
        return True

    def LoadTruthSpectrum(self):
        if self.fName == "__self_unfolding__":
            self.fTruthSpectrum = self.fDetectorTrainTruth.Rebin(len(self.fBins) - 1, "{0}_TruthSpectrum".format(self.fName), array.array('d', self.fBins))
            self.fTruthSpectrum.SetTitle("{0} Truth Spectrum".format(self.fName))
            self.fTruthSpectrum.Scale(self.fSelfUnfoldingScaling * self.fNumberOfEvents)
            return True

        if not self.fDMesonTruth:
            self.fTruthSpectrum = None
            return True

        dataTruthMesonList = self.fDataFile.Get(self.fDMesonTruth)
        if not dataTruthMesonList:
            print("Could not find list {0} in file {1}". format(self.fDMesonTruth, self.fDataFile.GetName()))
            exit(1)
        dataTruthJetListName = "_".join([self.fJetType, self.fJetRadius])
        dataTruthJetList = dataTruthMesonList.FindObject(dataTruthJetListName)
        if not dataTruthJetList:
            print("Could not find list {0}/{1} in file {2}". format(self.fDMesonTruth, dataTruthJetListName, self.fDataFile.GetName()))
            exit(1)
        dataTruthListName = "{0}_{1}_{2}".format(self.fDMesonTruth, dataTruthJetListName, self.fSpectrumTruthName)
        dataTruthList = dataTruthJetList.FindObject(dataTruthListName)
        if not dataTruthList:
            print("Could not find list {0}/{1}/{2} in file {3}". format(self.fDMesonTruth, dataTruthJetListName, dataTruthListName, self.fDataFile.GetName()))
            exit(1)
        truthSpectrumName = "{0}_{1}_{2}".format(self.fDMesonTruth, dataTruthJetListName, self.fSpectrumTruthName)
        truthSpectrum = dataTruthList.FindObject(truthSpectrumName)
        if not truthSpectrum:
            print("Could not find histogram {0} in list {1} in file {2}". format(truthSpectrumName, dataTruthListName, self.fDataFile.GetName()))
            exit(1)
        self.fTruthSpectrum = truthSpectrum.Clone("{0}_TruthSpectrum".format(self.fName))
        self.fTruthSpectrum.SetTitle("{0} Truth Spectrum".format(self.fName))

        return True

    def LoadInputSpectrum(self):
        if self.fName == "__self_unfolding__":
            self.fInputSpectrum = self.fDetectorResponse.ProjectionX("temp", 1, self.fDetectorResponse.GetNbinsY()).Rebin(len(self.fBins) - 1, "{0}_InputSpectrum".format(self.fName), array.array('d', self.fBins))
            self.fInputSpectrum.SetTitle("{0} Input Spectrum".format(self.fName))
            self.fInputSpectrum.Scale(self.fSelfUnfoldingScaling * self.fNumberOfEvents)
            return True

        wrap = RawYieldSpectrumLoader.RawYieldSpectrumLoader(self.fInputPath, self.fDataTrain, self.fAnalysisName)

        if self.fFDConfig:
            useFD = True
        else:
            useFD = False

        wrap.fFDConfig = self.fFDConfig
        wrap.fDataFile = self.fDataFile
        wrap.fDMeson = self.fDMeson
        wrap.fJetType = self.fJetType
        wrap.fJetRadius = self.fJetRadius
        wrap.fJetRecoScheme = self.fJetRecoScheme
        wrap.fVariableName = self.fVariableName
        wrap.fKinematicCuts = self.fKinematicCuts
        wrap.fDataSpectrumList = self.fDataList
        wrap.fUseReflections = self.fUseReflections
        wrap.fReflectionFit = self.fReflectionFit
        wrap.fRawYieldMethod = self.fRawYieldMethod

        self.fNumberOfEvents = wrap.LoadNumberOfEvents()
        self.fEvents = ROOT.TH1D("Events", "Events", 1, 0, 1)
        self.fEvents.SetBinContent(1, self.fNumberOfEvents)

        wrap.fReflectionRoS = self.fReflectionRoS

        self.fDataList = wrap.fDataSpectrumList

        if self.fUseMultiTrialInput:
            print("Getting input spectrum from the multi-trial analysis, with method {0}".format(self.fRawYieldMethod))
            inputSpectrum = wrap.GetDefaultSpectrumFromMultiTrial(useFD, self.fFDErrorBand, self.fRYErrorBand)
        else:
            print("Getting input spectrum from DMesonJetAnalysis, with method {0}".format(self.fRawYieldMethod))
            inputSpectrum = wrap.GetDefaultSpectrumFromDMesonJetAnalysis(useFD, self.fFDErrorBand, self.fRYErrorBand)

        if not inputSpectrum:
            print("Could not find input spectrum!")
            exit(1)

        if self.fRebinInput:
            self.fInputSpectrum = DMesonJetUtils.Rebin1D_fromBins(inputSpectrum, "{0}_InputSpectrum".format(self.fName), len(self.fRebinInput) - 1, numpy.array(self.fRebinInput, dtype=numpy.float32))
        else:
            self.fInputSpectrum = inputSpectrum.Clone("{0}_InputSpectrum".format(self.fName))
        self.fInputSpectrum.SetTitle("{0} Input Spectrum".format(self.fName))

        return True

    def Start(self, doPlotting=True):
        self.GenerateResponse()
        self.GenerateMeasuredUncertainty()
        for method, config in self.fUnfoldingConfig.iteritems():
            if method == "Svd":
                self.UnfoldSvd()
            elif method == "Bayes":
                self.UnfoldBayes(config["iter_min"], config["iter_max"], config["iter_step"])
            elif method == "BinByBin":
                self.UnfoldBinByBin()
            else:
                print("Unfolding method {0} not known!".format(method))
        self.GeneratePearsonMatrices()
        if doPlotting: self.Plot()

    def GenerateMeasuredUncertainty(self):
        self.fInputSpectrumErrors = ROOT.TH1D("{0}_Errors".format(self.fInputSpectrum.GetName()), "{0} Errors".format(self.fInputSpectrum.GetTitle()), self.fInputSpectrum.GetNbinsX(), self.fInputSpectrum.GetXaxis().GetXbins().GetArray())
        self.fInputSpectrumErrors.GetXaxis().SetTitle(self.fInputSpectrum.GetXaxis().GetTitle())
        self.fInputSpectrumErrors.GetYaxis().SetTitle("relative statistical uncertainty")
        for ibin in range(self.fInputSpectrum.GetNbinsX() + 2):
            if self.fInputSpectrum.GetBinContent(ibin) == 0:
                continue
            self.fInputSpectrumErrors.SetBinContent(ibin, self.fInputSpectrum.GetBinError(ibin) / self.fInputSpectrum.GetBinContent(ibin))

    def GeneratePearsonMatrices(self):
        for id, cov in self.fCovarianceMatrices.iteritems():
            pearson = cov.Clone(cov.GetName().replace("CovMat", "Pearson"))
            pearson.Reset()
            for x in range(1, pearson.GetNbinsX() + 1):
                if not cov.GetBinContent(x, x) > 0:
                    continue
                sigmax = math.sqrt(cov.GetBinContent(x, x))
                for y in range(1, pearson.GetNbinsY() + 1):
                    if not cov.GetBinContent(y, y) > 0:
                        continue
                    sigmay = math.sqrt(cov.GetBinContent(y, y))
                    pearson.SetBinContent(x, y, cov.GetBinContent(x, y) / sigmax / sigmay)
            self.fPearsonMatrices[id] = pearson

    def Plot(self):
        # ROOT.gStyle.SetPalette(55,ROOT.nullptr)
        if len(self.fUnfoldingConfig) > 1:
            self.CompareMethods()
        if len(self.fPriors) > 1:
            self.ComparePriors()
        self.UnfoldingSummary()
        self.PlotPearsonMatrices()
        self.PlotSvdDvectors()
        self.PlotResponses()
        self.PlotPriors()
        self.PlotBinByBinCorrectionFactors()
        self.PlotUnfoldingErrors()

    def PlotBinByBinCorrectionFactors(self):
        if not "BinByBin" in self.fUnfoldingConfig.keys():
            return
        baselineId = self.fDefaultPrior
        baseline = self.fBinByBinCorrectionFactors[baselineId].Clone("{0}_copy".format(self.fBinByBinCorrectionFactors[baselineId].GetName()))
        baseline.SetTitle(self.fDefaultPrior)
        globalList.append(baseline)
        spectra = []
        for priorIt, corrFact in self.fBinByBinCorrectionFactors.iteritems():
            if baselineId == priorIt:
                continue
            cf = corrFact.Clone("{0}_copy".format(corrFact.GetName()))
            cf.SetTitle(priorIt)
            spectra.append(cf)
            globalList.append(cf)
        if len(spectra) < 1:
            return
        comp = DMesonJetCompare.DMesonJetCompare("{0}_BinByBinCorrectionFactors".format(self.fName))
        comp.fOptSpectrum = "hist"
        comp.fOptRatio = "hist"
        comp.fYaxisRatio = "Prior x / {0}".format(self.fDefaultPrior)
        comp.fDoSpectraPlot = "lineary"
        comp.fDoRatioPlot = "lineary"
        r = comp.CompareSpectra(baseline, spectra)
        for obj in r:
            globalList.append(obj)
            if isinstance(obj, ROOT.TCanvas):
                self.fCanvases.append(obj)

    def PlotPriors(self):
        default_reg = self.GetDefaultRegularization(self.fDefaultMethod, self.fDefaultPrior)
        default_unfolded = self.fUnfoldedSpectra[self.fDefaultMethod, default_reg, self.fDefaultPrior]
        if self.fTruthSpectrum:
            baselineId = None
            baseline = self.fTruthSpectrum.Clone("{0}_copy".format(self.fTruthSpectrum.GetName()))
            baseline.SetTitle("MC Truth")
            yaxisRatio = "Prior / Truth"
        elif default_unfolded:
            baselineId = None
            baseline = default_unfolded.Clone("{0}_copy".format(default_unfolded.GetName()))
            yaxisRatio = "Prior / Unfolded"
        else:
            baselineId = self.fDefaultPrior
            baseline = self.fResponseMatrices[baselineId].fTruth.Clone("{0}_copy".format(self.fResponseMatrices[baselineId].fTruth.GetName()))
            baseline.SetTitle(self.fDefaultPrior)
            yaxisRatio = "Prior x / {0}".format(self.fDefaultPrior)
        baseline.Scale(1. / baseline.Integral())
        globalList.append(baseline)
        spectra = []
        for priorIt, resp in self.fResponseMatrices.iteritems():
            if baselineId == priorIt:
                continue
            spectrum = resp.fTruth.Clone("{0}_copy".format(resp.fTruth.GetName()))
            spectrum.Scale(1. / spectrum.Integral())
            spectrum.SetTitle(priorIt)
            spectra.append(spectrum)
            globalList.append(spectrum)
        if len(spectra) < 1:
            return
        comp = DMesonJetCompare.DMesonJetCompare("{0}_Priors".format(self.fName))
        if self.fName.startswith("JetZ"):
            comp.fX1LegRatio = 0.15
            comp.fX2LegRatio = 0.50
        comp.fOptSpectrum = "hist"
        comp.fOptRatio = "hist"
        comp.fYaxisRatio = yaxisRatio
        if self.fVariableName == "JetZ": comp.fDoSpectraPlot = "lineary"
        else: comp.fDoSpectraPlot = "logy"
        comp.fDoRatioPlot = "lineary"
        r = comp.CompareSpectra(baseline, spectra)
        for obj in r:
            globalList.append(obj)
            if isinstance(obj, ROOT.TCanvas):
                self.fCanvases.append(obj)

    def PlotResponses(self):
        for name, resp in self.fResponseMatrices.iteritems():
            cname = "{0}_Response_Prior{1}".format(self.fName, name)
            c = ROOT.TCanvas(cname, cname, 850, 750)
            self.fCanvases.append(c)
            globalList.append(c)
            c.Divide(2, 2)

            pad = c.cd(1)
            pad.SetLogz()
            pad.SetLeftMargin(0.16)
            pad.SetRightMargin(0.24)
            pad.SetBottomMargin(0.15)
            rm = resp.fNormResponse.DrawCopy("colz")
            rm.GetXaxis().SetTitleFont(43)
            rm.GetXaxis().SetTitleSize(19)
            rm.GetXaxis().SetTitleOffset(2.6)
            rm.GetXaxis().SetLabelFont(43)
            rm.GetXaxis().SetLabelSize(16)
            rm.GetYaxis().SetTitleFont(43)
            rm.GetYaxis().SetTitleSize(19)
            rm.GetYaxis().SetTitleOffset(2.6)
            rm.GetYaxis().SetLabelFont(43)
            rm.GetYaxis().SetLabelSize(16)
            rm.GetZaxis().SetTitleFont(43)
            rm.GetZaxis().SetTitleSize(16)
            rm.GetZaxis().SetTitleOffset(2.8)
            rm.GetZaxis().SetLabelFont(43)
            rm.GetZaxis().SetLabelSize(15)
            rm.SetMaximum(1)
            globalList.append(rm)

            pad = c.cd(2)
            pad.SetLeftMargin(0.14)
            pad.SetBottomMargin(0.15)
            self.PlotEfficiency(resp.fTotalEfficiency)

            pad = c.cd(3)
            pad.SetLeftMargin(0.14)
            pad.SetBottomMargin(0.15)
            self.PlotEfficiency(resp.fRecoEfficiency)

            pad = c.cd(4)
            pad.SetLeftMargin(0.14)
            pad.SetBottomMargin(0.15)
            self.PlotEfficiency(resp.fKineEfficiency)

    def PlotEfficiency(self, eff):
        eff_copy = eff.DrawCopy("hist")
        if eff_copy.GetMinimum() > 0.99:  # fix ugly plotting when the efficiency is 1
            eff_copy.SetMaximum(1.2)
            eff_copy.SetMinimum(0.8)
        eff_copy.GetXaxis().SetTitleFont(43)
        eff_copy.GetXaxis().SetTitleSize(19)
        eff_copy.GetXaxis().SetTitleOffset(2.6)
        eff_copy.GetXaxis().SetLabelFont(43)
        eff_copy.GetXaxis().SetLabelSize(16)
        eff_copy.GetYaxis().SetTitleFont(43)
        eff_copy.GetYaxis().SetTitleSize(19)
        eff_copy.GetYaxis().SetTitleOffset(2.6)
        eff_copy.GetYaxis().SetLabelFont(43)
        eff_copy.GetYaxis().SetLabelSize(16)
        globalList.append(eff_copy)

    def GetDefaultRegularization(self, method, prior):
        if isinstance(self.fUnfoldingConfig[method]["default_reg"], dict):
            reg = self.fUnfoldingConfig[method]["default_reg"][prior]
        else:
            reg = self.fUnfoldingConfig[method]["default_reg"]
        return reg

    def FilterObjectsFromDict(self, objects, method, reg, prior):
        for (methodIt, regIt, priorIt), cov in objects.iteritems():
            if method and not methodIt == method:
                continue
            if prior and not priorIt == prior:
                continue
            if reg and not regIt == reg:
                continue
            yield (methodIt, regIt, priorIt), cov

    def PlotUnfoldingErrors(self):
        for method, config in self.fUnfoldingConfig.iteritems():
            for prior in self.fPriors:
                self.PlotUnfoldingErrorsForMethod(method, prior)
        self.PlotUnfoldingErrorsMethodComparison()

    def PlotUnfoldingErrorsMethodComparison(self):
        baseline = self.fInputSpectrumErrors.Clone("{0}_copy".format(self.fInputSpectrumErrors.GetName()))
        baseline.SetTitle("Input spectrum")
        globalList.append(baseline)

        spectra = []
        for method in self.fUnfoldingConfig.iterkeys():
            reg = self.GetDefaultRegularization(method, self.fDefaultPrior)
            id = (method, reg, self.fDefaultPrior)
            (hErr, hCov, hToyMC) = self.fUnfoldedSpectraErrors[id]
            h = hErr.Clone("{0}_copy".format(hErr.GetName()))
            h.SetTitle(method)
            spectra.append(h)
            globalList.append(h)
        comp = DMesonJetCompare.DMesonJetCompare("{0}_UnfoldingStatisticalUncertaintyCompareMethods".format(self.fName))
        comp.fOptSpectrum = "hist"
        comp.fOptSpectrumBaseline = "hist"
        comp.fDoSpectraPlot = "lineary"
        comp.fDoRatioPlot = False
        r = comp.CompareSpectra(baseline, spectra)
        for obj in r:
            globalList.append(obj)
            if isinstance(obj, ROOT.TCanvas):
                self.fCanvases.append(obj)

    def PlotUnfoldingErrorsForMethod(self, method, prior):
        # first compare different regularization strengths
        baseline = self.fInputSpectrumErrors.Clone("{0}_copy".format(self.fInputSpectrumErrors.GetName()))
        baseline.SetTitle("Input spectrum")
        globalList.append(baseline)
        spectra = []
        for (methodIt, regIt, priorIt), unfolded in self.fUnfoldedSpectra.iteritems():
            if not (methodIt == method and  priorIt == prior):
                continue
            id = (method, regIt, prior)
            (hErr, hCov, hToyMC) = self.fUnfoldedSpectraErrors[id]
            h = hErr.Clone("{0}_copy".format(hErr.GetName()))
            h.SetTitle("Reg={0}".format(regIt))
            spectra.append(h)
            globalList.append(h)

        if len(spectra) > 0:
            comp = DMesonJetCompare.DMesonJetCompare("{0}_UnfoldingStatisticalUncertainty_{1}_Prior{2}".format(self.fName, method, prior))
            comp.fOptSpectrum = "hist"
            comp.fOptSpectrumBaseline = "hist"
            comp.fDoSpectraPlot = "lineary"
            comp.fDoRatioPlot = False
            r = comp.CompareSpectra(baseline, spectra)
            for obj in r:
                globalList.append(obj)
                if isinstance(obj, ROOT.TCanvas):
                    self.fCanvases.append(obj)

        # now compare different error estimation methods
        if self.fDoErrorCrossChecks:
            reg = self.GetDefaultRegularization(method, prior)
            baselineId = (method, reg, prior)
            (baselineErr, baselineCov, baselineToyMC) = self.fUnfoldedSpectraErrors[baselineId]
            baseline = baselineErr.Clone("{0}_copy".format(baselineErr.GetName()))
            baseline.SetTitle("kError")
            globalList.append(baseline)

            spectra = []
            hcov = baselineCov.Clone("{0}_copy".format(baselineCov.GetName()))
            hcov.SetTitle("kCovariance")
            globalList.append(hcov)
            spectra.append(hcov)
            htoy = baselineToyMC.Clone("{0}_copy".format(baselineToyMC.GetName()))
            htoy.SetTitle("kCovToy")
            globalList.append(htoy)
            spectra.append(htoy)

            comp = DMesonJetCompare.DMesonJetCompare("{0}_UnfoldingStatisticalUncertaintyStrategy_{1}_Reg{2}_Prior{3}".format(self.fName, method, reg, prior))
            comp.fOptSpectrum = "hist"
            comp.fOptSpectrumBaseline = "hist"
            comp.fDoSpectraPlot = "lineary"
            comp.fDoRatioPlot = False
            r = comp.CompareSpectra(baseline, spectra)
            for obj in r:
                globalList.append(obj)
                if isinstance(obj, ROOT.TCanvas):
                    self.fCanvases.append(obj)

    def PlotPearsonMatrices(self):
        for method, config in self.fUnfoldingConfig.iteritems():
            for prior in self.fPriors:
                self.PlotPearsonMatricesForMethod(method, prior)

    def PlotPearsonMatricesForMethod(self, method, prior):
        pearsons = []
        filteredPearsons = self.FilterObjectsFromDict(self.fPearsonMatrices, method, None, prior)
        for (methodIt, reg, priorIt), cov in filteredPearsons:
            cov_copy = cov.Clone()
            if reg:
                cov_copy.SetTitle("reg={0}".format(reg))
            else:
                cov_copy.SetTitle("No regularization".format(reg))
            pearsons.append(cov_copy)
            globalList.append(cov_copy)
        if len(pearsons) < 1:
            return

        cname = "{0}_Pearson_{1}_Prior{2}".format(self.fName, method, prior)
        rows = int(math.floor(math.sqrt(len(pearsons))))
        cols = int(math.ceil(float(len(pearsons)) / rows))
        c = ROOT.TCanvas(cname, cname, 360 * cols, 350 * rows)
        self.fCanvases.append(c)
        globalList.append(c)
        c.Divide(cols, rows)

        for i, cov in enumerate(pearsons):
            c.cd(i + 1)
            cov.GetXaxis().SetTitleFont(43)
            cov.GetXaxis().SetTitleSize(16)
            cov.GetXaxis().SetTitleOffset(2)
            cov.GetXaxis().SetLabelFont(43)
            cov.GetXaxis().SetLabelSize(15)
            cov.GetYaxis().SetTitleFont(43)
            cov.GetYaxis().SetTitleSize(16)
            cov.GetYaxis().SetTitleOffset(2)
            cov.GetYaxis().SetLabelFont(43)
            cov.GetYaxis().SetLabelSize(15)
            if i % cols == cols - 1:
                cov.Draw("colz")
            else:
                cov.Draw("col")
            cov.GetZaxis().SetRangeUser(-1, 1)
            pave = ROOT.TPaveText(0, 1, 1, 0.9, "NB NDC")
            pave.SetBorderSize(0)
            pave.SetFillStyle(0)
            pave.SetTextFont(43)
            pave.SetTextSize(17)
            pave.SetTextAlign(22)
            pave.AddText(cov.GetTitle())
            pave.Draw()
            globalList.append(pave)

    def PlotSvdDvectors(self):
        for prior, d in self.fSvdDvectors.iteritems():
            cname = d.GetName()
            c = ROOT.TCanvas(cname, cname)
            self.fCanvases.append(c)
            c.cd()
            c.SetLogy()
            globalList.append(c)
            d_copy = d.DrawCopy()
            globalList.append(d_copy)

    def UnfoldingSummary(self):
        for method, config in self.fUnfoldingConfig.iteritems():
            reg = self.GetDefaultRegularization(method, self.fDefaultPrior)
            self.UnfoldingSummaryForMethod(method, reg, self.fDefaultPrior)
            for prior in self.fPriors:
                self.RegularizationComparisonForMethod(method, prior)

        if self.fTruthSpectrum:
            ratioTruthOverMeas = self.fTruthSpectrum.Clone("{0}_TruthOverMeasured".format(self.fName))
            cnameTOM = "{0}_UnfoldingSummary_TruthOverMeasured".format(self.fName)
            ratioTruthOverMeas.GetYaxis().SetTitle("Truth / Measured")
            globalList.append(ratioTruthOverMeas)
            cTOM = ROOT.TCanvas(cnameTOM, cnameTOM)
            self.fCanvases.append(cTOM)
            cTOM.cd()
            globalList.append(cTOM)
            ratioTruthOverMeas.Divide(self.fInputSpectrum)
            ratioTruthOverMeas.SetFillColor(ROOT.kRed - 10)
            ratioTruthOverMeas.SetFillStyle(1001)
            ratioTruthOverMeas.SetLineColor(ROOT.kRed + 2)
            ratioTruthOverMeas.SetLineWidth(2)
            ratioTruthOverMeas.SetMarkerColor(ROOT.kRed + 2)
            ratioTruthOverMeas.SetMarkerStyle(ROOT.kFullCircle)
            ratioTruthOverMeas.SetMarkerSize(0.9)
            ratioTruthOverMeas.Draw()

    def RegularizationComparisonForMethod(self, method, prior):
        if self.fTruthSpectrum:
            baselineId = None
            baseline = self.fTruthSpectrum.Clone("{0}_copy".format(self.fTruthSpectrum.GetName()))
            baseline.SetTitle("MC Truth")
            yaxisRatio = "Unfolded / Truth"
        else:
            reg = self.GetDefaultRegularization(method, prior)
            baselineId = (method, reg, prior)
            baseline = self.fUnfoldedSpectra[baselineId].Clone("{0}_copy".format(self.fUnfoldedSpectra[baselineId].GetName()))
            baseline.SetTitle("Reg = {0}".format(reg))
            yaxisRatio = "Reg = x / Reg = {0}".format(reg)
        globalList.append(baseline)
        spectra = []
        for (methodIt, regIt, priorIt), unfolded in self.fUnfoldedSpectra.iteritems():
            if not (methodIt == method and  priorIt == prior):
                continue
            if baselineId == (methodIt, regIt, priorIt):
                continue
            spectrum = unfolded.Clone("{0}_copy".format(unfolded.GetName()))
            spectrum.SetTitle("Reg = {0}".format(regIt))
            spectra.append(spectrum)
            globalList.append(spectrum)
        if (self.fTruthSpectrum and len(spectra) < 2) or len(spectra) < 1:
            return
        comp = DMesonJetCompare.DMesonJetCompare("{0}_UnfoldingRegularization_{1}_Prior{2}".format(self.fName, method, prior))
        comp.fYaxisRatio = yaxisRatio
        comp.fOptRatio = "hist"
        if self.fVariableName == "JetZ": comp.fDoSpectraPlot = "lineary"
        else: comp.fDoSpectraPlot = "logy"
        comp.fDoRatioPlot = "lineary"
        r = comp.CompareSpectra(baseline, spectra)
        for obj in r:
            globalList.append(obj)
            if isinstance(obj, ROOT.TCanvas):
                self.fCanvases.append(obj)

    def UnfoldingSummaryForMethod(self, method, reg, prior):
        cname = "{0}_UnfoldingSummary_{1}".format(self.fName, method)
        c = ROOT.TCanvas(cname, cname)
        self.fCanvases.append(c)
        c.cd()
        globalList.append(c)

        leg = ROOT.TLegend(0.55, 0.70, 0.85, 0.87)
        globalList.append(leg)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(43)
        leg.SetTextSize(16)

        if self.fTruthSpectrum:
            truth = self.fTruthSpectrum.DrawCopy("hist")
            globalList.append(truth)
            truth.SetName("{0}_{1}_{2}_{3}".format(self.fUnfoldedSpectra[(method, reg, prior)].GetName(), method, reg, prior))
            truth.SetMarkerStyle(ROOT.kOpenSquare)
            truth.SetMarkerColor(ROOT.kOrange + 2)
            truth.SetMarkerSize(1.1)
            truth.SetFillColor(ROOT.kOrange - 9)
            truth.SetFillStyle(1001)
            truth.SetLineColor(ROOT.kOrange + 2)
            leg.AddEntry(truth, "Truth", "f")
            if self.fTruthSpectrum.GetMinimum() > 0 and self.fVariableName != "JetZ": c.SetLogy()
        else:
            truth = None

        if truth:
            meas = self.fInputSpectrum.DrawCopy("same")
        else:
            meas = self.fInputSpectrum.DrawCopy()
            if self.fInputSpectrum.GetMinimum() > 0 and self.fVariableName != "JetZ": c.SetLogy()
        globalList.append(meas)
        meas.SetName("{0}_{1}_{2}_{3}".format(meas, method, reg, prior))
        meas.SetMarkerStyle(ROOT.kFullCircle)
        meas.SetMarkerColor(ROOT.kRed + 2)
        meas.SetMarkerSize(0.9)
        meas.SetLineColor(ROOT.kRed + 2)
        leg.AddEntry(meas, "Measured", "pe")

        unfolded = self.fUnfoldedSpectra[(method, reg, prior)].DrawCopy("same")
        globalList.append(unfolded)
        unfolded.SetName("{0}_{1}_{2}_{3}".format(unfolded.GetName(), method, reg, prior))
        unfolded.SetMarkerStyle(ROOT.kFullSquare)
        unfolded.SetMarkerColor(ROOT.kBlue + 2)
        unfolded.SetMarkerSize(0.9)
        unfolded.SetLineColor(ROOT.kBlue + 2)
        leg.AddEntry(unfolded, "{0}, reg={1}, prior={2}".format(method, reg, prior), "pe")

        refolded = self.fRefoldedSpectra[(method, reg, prior)].DrawCopy("same")
        globalList.append(refolded)
        refolded.SetName("{0}_{1}_{2}_{3}".format(refolded.GetName(), method, reg, prior))
        refolded.SetMarkerStyle(ROOT.kOpenCircle)
        refolded.SetMarkerColor(ROOT.kGreen + 2)
        refolded.SetMarkerSize(1.1)
        refolded.SetLineColor(ROOT.kGreen + 2)
        leg.AddEntry(refolded, "Refolded", "l")

        leg.Draw()

        # ratios
        # unfolded/measured
        ratioUnfoldedOverMeas = unfolded.Clone("{0}_UnfoldedOverMeasured".format(meas.GetName()))
        cnameUOM = "{0}_UnfoldingSummary_{1}_UnfoldedOverMeasured".format(self.fName, method)
        ratioUnfoldedOverMeas.GetYaxis().SetTitle("Unfolded / Measured")
        globalList.append(ratioUnfoldedOverMeas)
        cUOM = ROOT.TCanvas(cnameUOM, cnameUOM)
        self.fCanvases.append(cUOM)
        cUOM.cd()
        globalList.append(cUOM)
        ratioUnfoldedOverMeas.Divide(meas)
        ratioUnfoldedOverMeas.SetLineWidth(2)
        ratioUnfoldedOverMeas.Draw("hist")

        # refolded/measured and unfolded/truth
        cnameROM = "{0}_UnfoldingSummary_{1}_RefoldedOverMeasured".format(self.fName, method)
        cROM = ROOT.TCanvas(cnameROM, cnameROM)
        self.fCanvases.append(cROM)
        cROM.cd()
        globalList.append(cROM)
        ratioRefoldedOverMeas = refolded.Clone("{0}_RefoldedOverMeasured".format(meas.GetName()))
        globalList.append(ratioRefoldedOverMeas)
        ratioRefoldedOverMeas.Divide(meas)
        ratioRefoldedOverMeas.SetLineWidth(2)
        ratioRefoldedOverMeas.Draw("hist")
        if truth:
            ratioRefoldedOverMeas.GetYaxis().SetTitle("ratio")
            ratioUnfoldedOverTruth = unfolded.Clone("{0}_UnfoldedOverTruth".format(meas.GetName()))
            globalList.append(ratioUnfoldedOverTruth)
            ratioUnfoldedOverTruth.Divide(truth)
            minY = min(ratioUnfoldedOverTruth.GetMinimum(), ratioRefoldedOverMeas.GetMinimum()) * 0.8
            maxY = max(ratioUnfoldedOverTruth.GetMaximum(), ratioRefoldedOverMeas.GetMaximum()) * 1.2
            ratioRefoldedOverMeas.GetYaxis().SetRangeUser(minY, maxY)
            ratioUnfoldedOverTruth.Draw("same")
            legROM = ROOT.TLegend(0.55, 0.78, 0.85, 0.87)
            globalList.append(legROM)
            legROM.SetFillStyle(0)
            legROM.SetBorderSize(0)
            legROM.SetTextFont(43)
            legROM.SetTextSize(16)
            legROM.AddEntry(ratioRefoldedOverMeas, "Refolded / Measured", "l")
            legROM.AddEntry(ratioUnfoldedOverTruth, "Unfolded / Truth", "pe")
            legROM.Draw()
        else:
            ratioRefoldedOverMeas.GetYaxis().SetTitle("Refolded / Measured")

    def CompareMethods(self):
        if self.fTruthSpectrum:
            baselineId = None
            baseline = self.fTruthSpectrum.Clone("{0}_copy".format(self.fTruthSpectrum.GetName()))
            baseline.SetTitle("MC Truth")
            yaxisRatio = "Unfolded / Truth"
        else:
            reg = self.GetDefaultRegularization(self.fDefaultMethod, self.fDefaultPrior)
            baselineId = (self.fDefaultMethod, reg, self.fDefaultPrior)
            baseline = self.fUnfoldedSpectra[baselineId].Clone("{0}_copy".format(self.fUnfoldedSpectra[baselineId].GetName()))
            yaxisRatio = "Other Method / {0}".format(self.fDefaultMethod)
        globalList.append(baseline)
        spectra = []
        for method, config in self.fUnfoldingConfig.iteritems():
            reg = self.GetDefaultRegularization(method, self.fDefaultPrior)
            id = (method, reg, self.fDefaultPrior)
            if id == baselineId:
                continue
            h = self.fUnfoldedSpectra[id].Clone("{0}_copy".format(self.fUnfoldedSpectra[id].GetName()))
            spectra.append(h)
            globalList.append(h)

        comp = DMesonJetCompare.DMesonJetCompare("{0}_UnfoldingMethod".format(self.fName))
        comp.fYaxisRatio = yaxisRatio
        comp.fOptRatio = "hist"
        if self.fVariableName == "JetZ": comp.fDoSpectraPlot = "lineary"
        else: comp.fDoSpectraPlot = "logy"
        comp.fDoRatioPlot = "lineary"
        r = comp.CompareSpectra(baseline, spectra)
        for obj in r:
            globalList.append(obj)
            if isinstance(obj, ROOT.TCanvas):
                self.fCanvases.append(obj)

    def ComparePriors(self):
        for method in self.fUnfoldingConfig.iterkeys():
            self.ComparePriorsForMethod(method)

    def ComparePriorsForMethod(self, method):
        config = self.fUnfoldingConfig[method]
        if self.fTruthSpectrum:
            baselineId = None
            baseline = self.fTruthSpectrum.Clone("{0}_copy".format(self.fTruthSpectrum.GetName()))
            baseline.SetTitle("MC Truth")
            yaxisRatio = "Unfolded / Truth"
        else:
            reg = self.GetDefaultRegularization(method, self.fDefaultPrior)
            baselineId = (method, reg, self.fDefaultPrior)
            baseline = self.fUnfoldedSpectra[baselineId].Clone("{0}_copy".format(self.fUnfoldedSpectra[baselineId].GetName()))
            yaxisRatio = "Prior = x / Prior = {0}".format(self.fDefaultPrior)
        globalList.append(baseline)
        spectra = []
        for prior in self.fPriors:
            reg = self.GetDefaultRegularization(method, prior)
            id = (method, reg, prior)
            if id == baselineId:
                continue
            h = self.fUnfoldedSpectra[id].Clone("{0}_copy".format(self.fUnfoldedSpectra[id].GetName()))
            spectra.append(h)
            globalList.append(h)

        comp = DMesonJetCompare.DMesonJetCompare("{0}_UnfoldingPrior_{1}".format(self.fName, method))
        if self.fName.startswith("JetZ"):
            comp.fX1LegRatio = 0.15
            comp.fX2LegRatio = 0.50
        comp.fYaxisRatio = yaxisRatio
        comp.fOptRatio = "hist"
        comp.fColors = comp.fColors[1:]
        comp.fLines = comp.fLines[1:]
        comp.fMarkers = comp.fMarkers[1:]
        if self.fVariableName == "JetZ": comp.fDoSpectraPlot = "lineary"
        else: comp.fDoSpectraPlot = "logy"
        comp.fDoRatioPlot = "lineary"
        r = comp.CompareSpectra(baseline, spectra)
        for obj in r:
            globalList.append(obj)
            if isinstance(obj, ROOT.TCanvas):
                self.fCanvases.append(obj)

    def UnfoldSvd(self):
        max_reg = self.fInputSpectrum.GetNbinsX()
        # if self.fUseOverflow: max_reg += 2
        for prior, resp in self.fResponseMatrices.iteritems():
            for reg in range(1, max_reg + 1):
                print("Unfolding {0}, reg={1}, prior={2}".format("SVD", reg, prior))
                unfold = ROOT.RooUnfoldSvd(resp.fRooUnfoldResponse, self.fInputSpectrum, reg)
                unfolded = unfold.Hreco()
                if not prior in self.fSvdDvectors.keys():
                    dvector = unfold.Impl().GetD().Clone("{0}_SvdDvector_Prior{1}".format(self.fName, prior))
                    dvector.GetXaxis().SetTitle("k")
                    dvector.GetYaxis().SetTitle("d")
                    self.fSvdDvectors[prior] = dvector

                unfolded.SetName("{0}_UnfoldedSpectrum_{1}_Reg{2}_Prior{3}".format(self.fName, "Svd", reg, prior))
                unfolded.SetTitle("{0}, k={1}, prior={2}".format("Svd", reg, prior))
                self.fUnfoldedSpectra["Svd", reg, prior] = unfolded

                # refolded
                refolded = resp.fRooUnfoldResponse.ApplyToTruth(unfolded)
                refolded.SetName("{0}_RefoldedSpectrum_{1}_Reg{2}_Prior{3}".format(self.fName, "Svd", reg, prior))
                refolded.SetTitle("{0}, k={1}, prior={2}".format("Svd", reg, prior))
                self.fRefoldedSpectra["Svd", reg, prior] = refolded

                # covariance
                cov = ROOT.TH2D(unfold.Ereco())
                cov.SetName("{0}_CovMat_{1}_Reg{2}_Prior{3}".format(self.fName, "Svd", reg, prior))
                cov.GetXaxis().SetTitle("bin number")
                cov.GetYaxis().SetTitle("bin number")
                self.fCovarianceMatrices["Svd", reg, prior] = cov

                # errors
                unfoldingErr = DMesonJetUtils.V2TH1(unfold.ErecoV(ROOT.RooUnfold.kErrors))
                unfoldingErr.SetName("{0}_UnfoldErr_{1}_Reg{2}_Prior{3}".format(self.fName, "Svd", reg, prior))
                unfoldingErr.SetTitle("{0} Unfolding Errors {1} Reg={2} Prior={3}".format(self.fName, "Svd", reg, prior))
                unfoldingErr.GetXaxis().SetTitle(unfolded.GetXaxis().GetTitle())
                unfoldingErr.GetXaxis().Set(unfolded.GetNbinsX(), unfolded.GetXaxis().GetXbins().GetArray())
                unfoldingErr.GetYaxis().SetTitle("relative statistical uncertainty")
                DMesonJetUtils.DivideNoErrors(unfoldingErr, unfolded)
                if self.fDoErrorCrossChecks:
                    diagCov = DMesonJetUtils.V2TH1(unfold.ErecoV(ROOT.RooUnfold.kCovariance))
                    diagCov.SetName("{0}_DiagCov_{1}_Reg{2}_Prior{3}".format(self.fName, "Svd", reg, prior))
                    diagCov.SetTitle("{0} Unfolding Errors (covariance) {1} Reg={2} Prior={3}".format(self.fName, "Svd", reg, prior))
                    diagCov.GetXaxis().SetTitle(unfolded.GetXaxis().GetTitle())
                    diagCov.GetXaxis().Set(unfolded.GetNbinsX(), unfolded.GetXaxis().GetXbins().GetArray())
                    diagCov.GetYaxis().SetTitle("relative statistical uncertainty")
                    DMesonJetUtils.DivideNoErrors(diagCov, unfolded)
                    errToy = DMesonJetUtils.V2TH1(unfold.ErecoV(ROOT.RooUnfold.kCovToy))
                    errToy.SetName("{0}_MCToyErr_{1}_Reg{2}_Prior{3}".format(self.fName, "Svd", reg, prior))
                    errToy.SetTitle("{0} Unfolding Errors (MC-toy) {1} Reg={2} Prior={3}".format(self.fName, "Svd", reg, prior))
                    errToy.GetXaxis().SetTitle(unfolded.GetXaxis().GetTitle())
                    errToy.GetXaxis().Set(unfolded.GetNbinsX(), unfolded.GetXaxis().GetXbins().GetArray())
                    errToy.GetYaxis().SetTitle("relative statistical uncertainty")
                    DMesonJetUtils.DivideNoErrors(errToy, unfolded)
                else:
                    diagCov = None
                    errToy = None
                self.fUnfoldedSpectraErrors["Svd", reg, prior] = unfoldingErr, diagCov, errToy

    def UnfoldBayes(self, iter_min, iter_max, iter_step):
        for prior, resp in self.fResponseMatrices.iteritems():
            for reg in range(iter_min, iter_max + 1, iter_step):
                print("Unfolding {0}, reg={1}, prior={2}".format("Bayes", reg, prior))
                unfold = ROOT.RooUnfoldBayes(resp.fRooUnfoldResponse, self.fInputSpectrum, reg)
                unfolded = unfold.Hreco()

                unfolded.SetName("{0}_UnfoldedSpectrum_{1}_Reg{2}_Prior{3}".format(self.fName, "Bayes", reg, prior))
                unfolded.SetTitle("{0}, iter={1}, prior={2}".format("Bayes", reg, prior))
                self.fUnfoldedSpectra["Bayes", reg, prior] = unfolded

                # refolded
                refolded = resp.fRooUnfoldResponse.ApplyToTruth(unfolded)
                refolded.SetName("{0}_RefoldedSpectrum_{1}_Reg{2}_Prior{3}".format(self.fName, "Bayes", reg, prior))
                refolded.SetTitle("{0}, iter={1}, prior={2}".format("Bayes", reg, prior))
                self.fRefoldedSpectra["Bayes", reg, prior] = refolded

                # covariance
                cov = ROOT.TH2D(unfold.Ereco())
                cov.SetName("{0}_CovMat_{1}_{2}_{3}".format(self.fName, "Bayes", reg, prior))
                cov.GetXaxis().SetTitle("bin number")
                cov.GetYaxis().SetTitle("bin number")
                self.fCovarianceMatrices["Bayes", reg, prior] = cov

                # errors
                unfoldingErr = DMesonJetUtils.V2TH1(unfold.ErecoV(ROOT.RooUnfold.kErrors))
                unfoldingErr.SetName("{0}_UnfoldErr_{1}_Reg{2}_Prior{3}".format(self.fName, "Bayes", reg, prior))
                unfoldingErr.SetTitle("{0} Unfolding Errors {1} Reg={2} Prior={3}".format(self.fName, "Bayes", reg, prior))
                unfoldingErr.GetXaxis().SetTitle(unfolded.GetXaxis().GetTitle())
                unfoldingErr.GetXaxis().Set(unfolded.GetNbinsX(), unfolded.GetXaxis().GetXbins().GetArray())
                unfoldingErr.GetYaxis().SetTitle("relative statistical uncertainty")
                DMesonJetUtils.DivideNoErrors(unfoldingErr, unfolded)
                if self.fDoErrorCrossChecks:
                    diagCov = DMesonJetUtils.V2TH1(unfold.ErecoV(ROOT.RooUnfold.kCovariance))
                    diagCov.SetName("{0}_DiagCov_{1}_Reg{2}_Prior{3}".format(self.fName, "Bayes", reg, prior))
                    diagCov.SetTitle("{0} Unfolding Errors (covariance) {1} Reg={2} Prior={3}".format(self.fName, "Bayes", reg, prior))
                    diagCov.GetXaxis().SetTitle(unfolded.GetXaxis().GetTitle())
                    diagCov.GetXaxis().Set(unfolded.GetNbinsX(), unfolded.GetXaxis().GetXbins().GetArray())
                    diagCov.GetYaxis().SetTitle("relative statistical uncertainty")
                    DMesonJetUtils.DivideNoErrors(diagCov, unfolded)
                    errToy = DMesonJetUtils.V2TH1(unfold.ErecoV(ROOT.RooUnfold.kCovToy))
                    errToy.SetName("{0}_MCToyErr_{1}_Reg{2}_Prior{3}".format(self.fName, "Bayes", reg, prior))
                    errToy.SetTitle("{0} Unfolding Errors (MC-toy) {1} Reg={2} Prior={3}".format(self.fName, "Bayes", reg, prior))
                    errToy.GetXaxis().SetTitle(unfolded.GetXaxis().GetTitle())
                    errToy.GetXaxis().Set(unfolded.GetNbinsX(), unfolded.GetXaxis().GetXbins().GetArray())
                    errToy.GetYaxis().SetTitle("relative statistical uncertainty")
                    DMesonJetUtils.DivideNoErrors(errToy, unfolded)
                else:
                    diagCov = None
                    errToy = None
                self.fUnfoldedSpectraErrors["Bayes", reg, prior] = unfoldingErr, diagCov, errToy

    def UnfoldBinByBin(self):
        for prior, resp in self.fResponseMatrices.iteritems():
            print("Unfolding {0}, prior={1}".format("BinByBin", prior))
            binBybinCorrFactors = resp.fResponse.ProjectionY("{0}_BinByBinCorrectionFactors_Prior{1}".format(self.fName, prior))
            binBybinCorrFactors.SetTitle("{0} Correction Factors, prior={1}".format("BinByBin", prior))
            binBybinCorrFactors.Divide(resp.fResponse.ProjectionX())
            self.fBinByBinCorrectionFactors[prior] = binBybinCorrFactors

            unfold = ROOT.RooUnfoldBinByBin(resp.fRooUnfoldResponse, self.fInputSpectrum)
            unfolded = unfold.Hreco()

            unfolded.SetName("{0}_UnfoldedSpectrum_{1}_Prior{2}".format(self.fName, "BinByBin", prior))
            unfolded.SetTitle("{0}, prior={1}".format("BinByBin", prior))
            self.fUnfoldedSpectra["BinByBin", None, prior] = unfolded

            # refolded
            refolded = resp.fRooUnfoldResponse.ApplyToTruth(unfolded)
            # refolded = unfolded.Clone("refolded")
            # refolded.Divide(binBybinCorrFactors)
            refolded.SetName("{0}_RefoldedSpectrum_{1}_Prior{2}".format(self.fName, "BinByBin", prior))
            refolded.SetTitle("{0}, prior={1}".format("BinByBin", prior))
            self.fRefoldedSpectra["BinByBin", None, prior] = refolded

            # covariance
            cov = ROOT.TH2D(unfold.Ereco())
            cov.SetName("{0}_CovMat_{1}_{2}".format(self.fName, "BinByBin", prior))
            cov.GetXaxis().SetTitle("bin number")
            cov.GetYaxis().SetTitle("bin number")
            self.fCovarianceMatrices["BinByBin", None, prior] = cov

            # errors
            unfoldingErr = DMesonJetUtils.V2TH1(unfold.ErecoV(ROOT.RooUnfold.kErrors))
            unfoldingErr.SetName("{0}_UnfoldErr_{1}_Prior{2}".format(self.fName, "BinByBin", prior))
            unfoldingErr.SetTitle("{0} Unfolding Errors {1} Prior={2}".format(self.fName, "BinByBin", prior))
            unfoldingErr.GetXaxis().SetTitle(unfolded.GetXaxis().GetTitle())
            unfoldingErr.GetXaxis().Set(unfolded.GetNbinsX(), unfolded.GetXaxis().GetXbins().GetArray())
            unfoldingErr.GetYaxis().SetTitle("relative statistical uncertainty")
            DMesonJetUtils.DivideNoErrors(unfoldingErr, unfolded)
            if self.fDoErrorCrossChecks:
                diagCov = DMesonJetUtils.V2TH1(unfold.ErecoV(ROOT.RooUnfold.kCovariance))
                diagCov.SetName("{0}_DiagCov_{1}_Prior{2}".format(self.fName, "BinByBin", prior))
                diagCov.SetTitle("{0} Unfolding Errors (covariance) {1} Prior={2}".format(self.fName, "BinByBin", prior))
                diagCov.GetXaxis().SetTitle(unfolded.GetXaxis().GetTitle())
                diagCov.GetXaxis().Set(unfolded.GetNbinsX(), unfolded.GetXaxis().GetXbins().GetArray())
                diagCov.GetYaxis().SetTitle("relative statistical uncertainty")
                DMesonJetUtils.DivideNoErrors(diagCov, unfolded)
                errToy = DMesonJetUtils.V2TH1(unfold.ErecoV(ROOT.RooUnfold.kCovToy))
                errToy.SetName("{0}_MCToyErr_{1}_Prior{2}".format(self.fName, "BinByBin", prior))
                errToy.SetTitle("{0} Unfolding Errors (MC-toy) {1} Prior={2}".format(self.fName, "BinByBin", prior))
                errToy.GetXaxis().SetTitle(unfolded.GetXaxis().GetTitle())
                errToy.GetXaxis().Set(unfolded.GetNbinsX(), unfolded.GetXaxis().GetXbins().GetArray())
                errToy.GetYaxis().SetTitle("relative statistical uncertainty")
                DMesonJetUtils.DivideNoErrors(errToy, unfolded)
            else:
                diagCov = None
                errToy = None
            self.fUnfoldedSpectraErrors["BinByBin", None, prior] = unfoldingErr, diagCov, errToy

    def GenerateResponse(self):
        for prior in self.fPriors:
            print("Generating response for prior {0}".format(prior))
            (detResp, trainTruth), (smallBinDetResp, smallBinTrainTruth) = self.NormalizeAndRebinResponseMatrix(prior)
            detResp.GetZaxis().SetTitle("counts")
            detRespForRooUnfold = detResp.Clone()
            if self.fUseOverflow:
                # if using overflow bins, need to reset the content of the measured axis overflow bins (not measured)
                # not to be done in detResp, since there the overflow bins are needed to calculate the kinematic efficiency
                for ybin in range(0, detRespForRooUnfold.GetNbinsY() + 2):
                    detRespForRooUnfold.SetBinContent(0, ybin, 0)
                    detRespForRooUnfold.SetBinError(0, ybin, 0)
                    detRespForRooUnfold.SetBinContent(detRespForRooUnfold.GetNbinsX() + 1, ybin, 0)
                    detRespForRooUnfold.SetBinError(detRespForRooUnfold.GetNbinsX() + 1, ybin, 0)
            rooUnfoldResp = ROOT.RooUnfoldResponse(ROOT.nullptr, trainTruth, detRespForRooUnfold)
            rooUnfoldResp.UseOverflow(self.fUseOverflow)
            detRespNorm = self.Normalize2D(detResp, self.GenerateFlatPrior(False))
            detRespNorm.GetZaxis().SetTitle("Probability Density")
            respMatrix = UnfoldingResponseMatrix.ResponseMatrix("Prior{0}".format(prior), rooUnfoldResp, detResp, trainTruth, smallBinDetResp, smallBinTrainTruth, detRespNorm)
            self.fResponseMatrices[prior] = respMatrix
            if self.fName == "__self_unfolding__" and not self.fInputSpectrum:
                respMatrix.SelfUnfoldingData(self)

    def NormalizeAndRebinResponseMatrix(self, prior):
        (normResp, priorHist) = self.NormalizeResponseMatrix(prior)
        coarseResp = self.Rebin(normResp)
        coarsePrior = self.Rebin(priorHist)
        return (coarseResp, coarsePrior), (normResp, priorHist)

    def Rebin(self, hist):
        if isinstance(hist, ROOT.TH2):
            return DMesonJetUtils.Rebin2D(hist, self.fInputSpectrum.GetXaxis(), self.fInputSpectrum.GetXaxis(), True)
        elif isinstance(hist, ROOT.TH1):
            return DMesonJetUtils.Rebin1D(hist, self.fInputSpectrum.GetXaxis())
        else:
            print("Object {0} of unrecognized type!".format(hist))
            return None  # will fail

    def ComputePrior(self, prior, scale_bin_width=True):
        if "_multiply_" in prior:
            priorHistos = [self.ComputePrior(prior_mult, False) for prior_mult in prior.split("_multiply_")]
            priorHist = priorHistos[0][0].Clone("myprior")
            axisCompare = priorHistos[0][1]
            for h, a in priorHistos[1:]:
                priorHist.Multiply(h)
                if a != DMesonJetUtils.AxisCompare.Identical:
                    if axisCompare == DMesonJetUtils.AxisCompare.Identical:
                        axisCompare = a
                    elif axisCompare != DMesonJetUtils.AxisCompare.Identical:
                        print("ComputePrior, Error in binning while computing prior '{}'".format(prior))
                        exit(1)
            if scale_bin_width and not "ResponseTruth" in prior:
                for ibin in range(0, priorHist.GetNbinsX() + 2):
                    priorHist.SetBinContent(ibin, priorHist.GetBinContent(ibin) * priorHist.GetXaxis().GetBinWidth(ibin))
        elif prior == "ResponseTruth":
            priorHist = self.fDetectorTrainTruth
            axisCompare = DMesonJetUtils.AxisCompare.Identical
        elif prior == "Flat":
            priorHist = self.GenerateFlatPrior(scale_bin_width)
            axisCompare = DMesonJetUtils.AxisCompare.Identical
        elif "PowerLaw" in prior:
            priorHist = self.GeneratePowerLawPrior(-int(prior.replace("PowerLaw_", "")), 3, scale_bin_width)
            axisCompare = DMesonJetUtils.AxisCompare.Identical
        elif "pol" in prior:
            pars_string = prior.replace("pol(", "")
            pars_string = pars_string.replace(")", "")
            pars = [float(v) for v in pars_string.split(",")]
            priorHist = self.GeneratePolinomialPrior(pars, scale_bin_width)
            axisCompare = DMesonJetUtils.AxisCompare.Identical
        else:
            priorHist = self.GetCustomPrior(prior)
            axisCompare = DMesonJetUtils.AxisCompare.CheckConsistency(priorHist.GetXaxis(), self.fDetectorTrainTruth.GetXaxis())

        return priorHist, axisCompare

    def NormalizeResponseMatrix(self, prior):
        if prior == "ResponseTruth": return self.fDetectorResponse, self.fDetectorTrainTruth

        priorHist, axisCompare = self.ComputePrior(prior)

        if axisCompare == DMesonJetUtils.AxisCompare.Identical:
            detectorResponse = self.fDetectorResponse
            detectorTrainTruth = self.fDetectorTrainTruth

        elif axisCompare == DMesonJetUtils.AxisCompare.IsContainedSameBinning or \
        axisCompare == DMesonJetUtils.AxisCompare.IsContained or \
        (axisCompare == DMesonJetUtils.AxisCompare.SameLimits and priorHist.GetXaxis().GetNbins() <= self.fDetectorTrainTruth.GetNbins()):
            print("The custom prior {0} has a different binning compared to the response matrix. Trying to rebin the response matrix, hopefully things will go well (and you know what you are doing!)".format(prior))
            print("Prior axis: {}".format(list(priorHist.GetXaxis().GetXbins())))
            print("Response matrix axis: {}".format(list(self.fDetectorResponse.GetYaxis().GetXbins())))
            detectorResponse = DMesonJetUtils.Rebin2D(self.fDetectorResponse, self.fDetectorResponse.GetXaxis(), priorHist.GetXaxis(), False)
            detectorTrainTruth = DMesonJetUtils.Rebin1D(self.fDetectorTrainTruth, priorHist.GetXaxis())

        elif axisCompare == DMesonJetUtils.AxisCompare.ContainsSameBinning or \
        axisCompare == DMesonJetUtils.AxisCompare.Contains or \
        (axisCompare == DMesonJetUtils.AxisCompare.SameLimits and priorHist.GetXaxis().GetNbins() > self.fDetectorTrainTruth.GetNbins()):
            print("The custom prior {0} has a different binning compared to the response matrix. Trying to rebin the custom prior, hopefully things will go well (and you know what you are doing!)".format(prior))
            print("Prior axis: {}".format(list(priorHist.GetXaxis().GetXbins())))
            print("Response matrix axis: {}".format(list(self.fDetectorResponse.GetYaxis().GetXbins())))
            priorHist = DMesonJetUtils.Rebin1D(priorHist, self.fDetectorTrainTruth.GetXaxis())
            detectorResponse = self.fDetectorResponse
            detectorTrainTruth = self.fDetectorTrainTruth

        else:
            print("DMesonJetUnfolding, NormalizeResponseMatrix: this case not implemented. Aborting.")
            exit(1)

        # will zero entries in the prior corresponding to y slices of the response matrix that have less than 1% of the total integral
        # commenting this out .... not sure why I did this!
#         totint = detectorResponse.Integral(0, detectorResponse.GetXaxis().GetNbins() + 1, 0, detectorResponse.GetYaxis().GetNbins() + 1)
#         for ybin in range(0, detectorResponse.GetYaxis().GetNbins() + 2):
#             inty = detectorResponse.Integral(0, detectorResponse.GetXaxis().GetNbins() + 1, ybin, ybin)
#             if inty < totint * 1e-2:
#                 priorHist.SetBinContent(ybin, 0)
#                 priorHist.SetBinError(ybin, 0)

        priorHist_Integral = priorHist.Integral()
        if priorHist_Integral != 0:
            priorHist.Scale(detectorTrainTruth.Integral() / priorHist_Integral)
        else:
            print("Error (DMesonJetUnfolding, NormalizeResponseMatrix): prior histogram '{}' is empty!".format(priorHist.GetName()))
            exit(1)
        priorEffHist = priorHist.Clone("priorEffHist")
        priorEffHist.Multiply(detectorResponse.ProjectionY())
        priorEffHist.Divide(detectorTrainTruth)

        resp = self.Normalize2D(detectorResponse, priorEffHist)

        return resp, priorHist

    def SetCustomPrior(self, name, hist):
        self.fCustomPriors[name] = hist

    def GetCustomPrior(self, name):
        if not name in self.fCustomPriors:
            print("Need to set the custom prior {0} by using SetCustomPrior".format(name))
            exit(1)
        return self.fCustomPriors[name]

    def GeneratePowerLawPrior(self, a, peak, scale_bin_width):
        print("Generating power law prior with index {0}".format(a))
        priorHist = self.fDetectorTrainTruth.Clone("myprior")
        priorHist.Reset()
        # this is a modified power law, where at low pT an exponential dumps the function to zero, to avoid infinities
        # other protections are added via the [2] parameter to avoid intermediate infinities in the function evaluation
        # the resulting modified power law has a local maximum, which is set to 3 GeV/c by default
        f = ROOT.TF1("f", "(x>[2])*TMath::Max([2],x)^[0]*TMath::Exp([1]/TMath::Max([2],x))", priorHist.GetXaxis().GetBinLowEdge(0), priorHist.GetXaxis().GetBinUpEdge(priorHist.GetXaxis().GetNbins() + 1))
        f.SetParameter(0, a)
        f.SetParameter(1, a * peak)
        f.SetParameter(2, 1e-6)
        for ibin in range(0, priorHist.GetNbinsX() + 2):
            if scale_bin_width:
                priorHist.SetBinContent(ibin, f.Eval(priorHist.GetXaxis().GetBinCenter(ibin)) * priorHist.GetXaxis().GetBinWidth(ibin))
            else:
                priorHist.SetBinContent(ibin, f.Eval(priorHist.GetXaxis().GetBinCenter(ibin)))
        return priorHist

    def GeneratePolinomialPrior(self, pars, scale_bin_width):
        priorHist = self.fDetectorTrainTruth.Clone("myprior")
        priorHist.Reset()
        f = ROOT.TF1("f", "pol{}".format(len(pars)), priorHist.GetXaxis().GetBinLowEdge(0), priorHist.GetXaxis().GetBinUpEdge(priorHist.GetXaxis().GetNbins() + 1))
        for i, p in enumerate(pars): f.SetParameter(i, p)
        for ibin in range(0, priorHist.GetNbinsX() + 2):
            if scale_bin_width:
                priorHist.SetBinContent(ibin, f.Eval(priorHist.GetXaxis().GetBinCenter(ibin)) * priorHist.GetXaxis().GetBinWidth(ibin))
            else:
                priorHist.SetBinContent(ibin, f.Eval(priorHist.GetXaxis().GetBinCenter(ibin)))
        return priorHist

    def GenerateFlatPrior(self, scale_bin_width):
        priorHist = self.fDetectorTrainTruth.Clone("myprior")
        for ibin in range(0, priorHist.GetNbinsX() + 2):
            priorHist.SetBinError(ibin, 0)
            if scale_bin_width:
                priorHist.SetBinContent(ibin, priorHist.GetXaxis().GetBinWidth(ibin))
            else:
                priorHist.SetBinContent(ibin, 1)
        return priorHist

    def Normalize2D(self, hist, norm):
        resp = hist.Clone("myresp")
        for ybin in range(0, resp.GetYaxis().GetNbins() + 2):
            inty = resp.Integral(0, resp.GetXaxis().GetNbins() + 1, ybin, ybin)
            if inty == 0:
                continue
            scaling = norm.GetBinContent(ybin) / inty
            for xbin in range(0, resp.GetXaxis().GetNbins() + 2):
                binValue = resp.GetBinContent(xbin, ybin) * scaling
                binErr = resp.GetBinError(xbin, ybin) * scaling
                resp.SetBinContent(xbin, ybin, binValue)
                resp.SetBinError(xbin, ybin, binErr)
        return resp


class DMesonJetUnfolding:

    def __init__(self, name, input_path, dataTrain, data, responseTrain, response, mt, refl, refl_fit, refl_ros):
        self.fName = name
        if mt: self.fName += "_mt"
        if refl:
            self.fName += "_refl_{0}".format(refl_fit)
            if not refl_ros == 0: self.fName += "_{0}".format(refl_ros)
        self.fInputPath = input_path
        self.fDataTrain = dataTrain
        self.fAnalysisName = data
        self.fResponseTrain = responseTrain
        self.fResponseAnalysisName = response
        self.fUseMultiTrial = mt
        self.fUseReflections = refl
        self.fReflectionFit = refl_fit
        self.fReflectionRoS = refl_ros
        self.fUnfoldingEngine = []
        self.OpenFiles()

    def OpenFiles(self):
        self.fDataFile = ROOT.TFile("{0}/{1}/{2}.root".format(self.fInputPath, self.fDataTrain, self.fAnalysisName))
        self.fResponseFile = ROOT.TFile("{0}/{1}/{2}.root".format(self.fInputPath, self.fResponseTrain, self.fResponseAnalysisName))

    def StartUnfolding(self, config, eff, use_overflow, fd_error_band, ry_error_band, plot):
        eng = DMesonJetUnfoldingEngine(self.fInputPath, self.fDataTrain, self.fAnalysisName, config, fd_error_band, ry_error_band)
        self.fUnfoldingEngine.append(eng)
        eng.fUseMultiTrialInput = self.fUseMultiTrial
        eng.fUseReflections = self.fUseReflections
        eng.fReflectionFit = self.fReflectionFit
        eng.fReflectionRoS = self.fReflectionRoS
        r = eng.LoadData(self.fDataFile, self.fResponseFile, eff, use_overflow)
        if r: eng.Start(plot)

    def SaveRootFile(self):
        path = "{0}/{1}".format(self.fInputPath, self.fName)
        if not os.path.isdir(path):
            os.makedirs(path)
        fname = "{0}/{1}.root".format(path, self.fName)
        print("Storing results in file {0}".format(fname))
        file = ROOT.TFile(fname, "recreate")
        file.cd()
        for eng in self.fUnfoldingEngine:
            rlist = eng.GenerateRootList()
            rlist.Write(rlist.GetName(), ROOT.TObject.kSingleKey)
        file.Close()

    def SavePlots(self, format):
        path = "{0}/{1}/{2}".format(self.fInputPath, self.fName, format)
        if not os.path.isdir(path):
            os.makedirs(path)
        for eng in self.fUnfoldingEngine:
            eng.SavePlots(path, format)
