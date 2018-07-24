#!/usr/bin/env python
# python base classes and utilities for D Meson jet analysis

import ROOT
import DMesonJetUtils
import Axis
from DMesonJetBase import AnalysisType


class Spectrum:

    def __init__(self, config, dmeson, jtype, jradius, jscheme, jtitle, binSet, effWeight):
        self.fConfig = config
        self.fDMeson = dmeson
        self.fJetType = jtype
        self.fJetRadius = jradius
        self.fJetRecoScheme = jscheme
        self.fJetTitle = jtitle
        self.fReflFitterName = ""
        if "suffix" in config:
            self.fSuffix = config["suffix"]
        else:
            self.fSuffix = None
        self.fSimpleName = config["name"]
        self.fName = '_'.join(obj for obj in [self.fDMeson, self.fJetType, self.fJetRadius, self.fJetRecoScheme, self.fSimpleName, self.fSuffix] if obj)
        self.fBinSet = binSet
        self.fHistogram = None
        self.fNormHistogram = None
        self.fUncertainty = None
        self.fMass = None
        self.fMassWidth = None
        self.fBackground = None
        self.fAxis = []
        self.fSideBandWindowInvMassHistos = dict()
        self.fSignalWindowInvMassHistos = dict()
        self.fSkipBins = None
        self.fEfficiencyWeight = effWeight

        # S-B analysis
        self.fSideBandHistograms = None
        self.fSideBandLeftHistograms = None
        self.fSideBandRightHistograms = None
        self.fSignalHistograms = None
        self.fSideBandWindowTotalHistogram = None
        self.fSignalWindowTotalHistogram = None

        # L-S analysis
        self.fLikeSignSubtractedBinSet = None
        self.fLikeSignHistograms = None
        self.fUnlikeSignHistograms = None
        self.fLikeSignTotalHistogram = None
        self.fUnlikeSignTotalHistogram = None

        self.fTitle = config["title"]
        if binSet.fTitle:
            if self.fTitle:
                self.fTitle += ", "
            self.fTitle += binSet.fTitle

        self.fComparisonDone = False
        if "compare" in config:
            self.fCompare = list(config["compare"])
        else:
            self.fCompare = None

        if "comp_titles" in config:
            self.fComparisonTitles = list(config["comp_titles"])
        else:
            self.fComparisonTitles = None

        if config["type"] == "side_band":
            self.fAnalysisType = AnalysisType.SideBand
            self.fSideBandMinSigmas = config["side_band"]["min_sigmas"]
            self.fSideBandMaxSigmas = config["side_band"]["max_sigmas"]
            self.fBinCountSignalSigmas = config["side_band"]["max_signal_sigmas"]
            self.fBackupSigma = config["side_band"]["backup_sigma"]
            self.fBackupMean = config["side_band"]["backup_mean"]
            if "skip_bins" in config["side_band"]:
                self.fSkipBins = config["side_band"]["skip_bins"]
        elif config["type"] == "like_sign":
            if config["like_sign"]["mode"] == "bin_count":
                self.fAnalysisType = AnalysisType.LikeSign
            elif config["like_sign"]["mode"] == "fit":
                self.fAnalysisType = AnalysisType.LikeSignFit
            else:
                print("Like Sign mode {0} not recognized!!".format(config["like_sign"]["mode"]))
                exit(1)
            self.fSideBandMinSigmas = config["like_sign"]["min_sigmas"]
            self.fSideBandMaxSigmas = config["like_sign"]["max_sigmas"]
            self.fBinCountSignalSigmas = config["like_sign"]["max_signal_sigmas"]
            self.fBinCountNormSignalSigmas = config["like_sign"]["max_signal_norm_sigmas"]
            self.fBackupSigma = config["like_sign"]["backup_sigma"]
            self.fBackupMean = config["like_sign"]["backup_mean"]
            self.fLikeSignTree = config["like_sign"]["name"]
            if "skip_bins" in config["like_sign"]:
                self.fSkipBins = config["like_sign"]["skip_bins"]
        elif config["type"] == "inv_mass_fit":
            self.fAnalysisType = AnalysisType.InvMassFit
        elif config["type"] == "truth":
            self.fAnalysisType = AnalysisType.Truth
        else:
            print("Analysis type {0} not recognized!".format(config["type"]))

        if "axis" in config:
            if len(config["axis"]) > 2:
                print("Error: cannot do bin counting spectra (e.g. side band) with more than 2 axis. Spectrum {}".format(config["name"]))
                exit(1)
            for axis_name, axis_bins in config["axis"].iteritems():
                self.fAxis.append(Axis.Axis(axis_name, axis_bins, "", (jtype != "Full")))
        else:
            for axis in binSet.fAxis:
                self.fAxis.append(axis)

        print("Spectrum {0} with analysis type {1} added".format(self.fName, self.fAnalysisType.name))

    def Clone(self):
        s_new = Spectrum(self.fConfig, self.fDMeson, self.fJetType, self.fJetRadius, self.fJetRecoScheme, self.fJetTitle, self.fBinSet, self.fEfficiencyWeight)
        return s_new

    def GenerateRootList(self):
        rlist = ROOT.TList()
        rlist.SetName(self.fName)
        if self.fHistogram:
            rlist.Add(self.fHistogram)
        if self.fNormHistogram:
            rlist.Add(self.fNormHistogram)
        if self.fUncertainty:
            rlist.Add(self.fUncertainty)
        if self.fMass:
            rlist.Add(self.fMass)
        if self.fMassWidth:
            rlist.Add(self.fMassWidth)
        if self.fBackground:
            rlist.Add(self.fBackground)
        if self.fAnalysisType == AnalysisType.SideBand and self.fSideBandHistograms and self.fSignalHistograms:
            SBlist = ROOT.TList()
            SBlist.SetName("SideBandAnalysis")
            for h in self.fSideBandWindowInvMassHistos.itervalues():
                SBlist.Add(h)
            for h in self.fSignalWindowInvMassHistos.itervalues():
                SBlist.Add(h)
            for h in self.fSideBandHistograms:
                SBlist.Add(h)
            for h in self.fSideBandLeftHistograms:
                SBlist.Add(h)
            for h in self.fSideBandRightHistograms:
                SBlist.Add(h)
            for h in self.fSignalHistograms:
                SBlist.Add(h)
            SBlist.Add(self.fSideBandWindowTotalHistogram)
            SBlist.Add(self.fSignalWindowTotalHistogram)
            rlist.Add(SBlist)
        if self.fAnalysisType == AnalysisType.LikeSign and self.fLikeSignSubtractedBinSet:
            LSlist = self.fLikeSignSubtractedBinSet.GenerateInvMassRootList()
            LSlist.SetName("LikeSignAnalysis")
            for h in self.fLikeSignHistograms:
                LSlist.Add(h)
            for h in self.fUnlikeSignHistograms:
                LSlist.Add(h)
            LSlist.Add(self.fLikeSignTotalHistogram)
            LSlist.Add(self.fUnlikeSignTotalHistogram)
            rlist.Add(LSlist)
        return rlist

    def GenerateNormalizedSpectrum(self, events, weighted=False):
        if self.fHistogram:
            self.fNormHistogram = self.GenerateNormalizedSpectrumForHistogram(self.fHistogram, events, weighted)

    def GenerateNormalizedSpectrumForHistogram(self, hist, events, weighted):
        hname = "{0}_Normalized".format(hist.GetName())
        result = hist.Clone(hname)
        result.SetTitle(hname)
        if len(self.fAxis) == 1:
            if weighted:
                axisTitle = "#frac{{d#sigma}}{{d{var}}}".format(var=self.fAxis[0].GetVariableName())
                if self.fAxis[0].GetVariableUnits():
                    axisTitle += " [mb ({0})^{{-1}}]".format(self.fAxis[0].GetVariableUnits())
                else:
                    axisTitle += " (mb)"
            else:
                axisTitle = "#frac{{1}}{{#it{{N}}_{{evt}}}} #frac{{d#it{{N}}}}{{d{var}}}".format(var=self.fAxis[0].GetVariableName())
                if self.fAxis[0].GetVariableUnits():
                    axisTitle += " ({0})^{{-1}}".format(self.fAxis[0].GetVariableUnits())
            result.GetYaxis().SetTitle(axisTitle)
        elif len(self.fAxis) == 2:
            units1 = self.fAxis[0].GetVariableUnits()
            units2 = self.fAxis[1].GetVariableUnits()
            if weighted:
                axisTitle = "#frac{{d^{{2}}#sigma}}{{d{var1} d{var2}}}".format(var1=self.fAxis[0].GetVariableName(), var2=self.fAxis[1].GetVariableName())
                if units1 != units2:
                    axisTitle += " [mb"
                    if units1:
                        axisTitle += " ({0})^{{-1}}".format(units1)
                    if units2:
                        axisTitle += " ({0})^{{-1}}".format(units2)
                    axisTitle += "]"
                else:
                    if units1:
                        axisTitle += " [mb ({0})^{{-2}}]".format(units1)
            else:
                axisTitle = "#frac{{1}}{{#it{{N}}_{{evt}}}} #frac{{d^{{2}}#it{{N}}}}{{d{var1} d{var2}}}".format(var1=self.fAxis[0].GetVariableName(), var2=self.fAxis[1].GetVariableName())
                if units1 != units2:
                    if units1:
                        axisTitle += " ({0})^{{-1}}".format(units1)
                    if units2:
                        axisTitle += " ({0})^{{-1}}".format(units2)
                else:
                    if units1:
                        axisTitle += " ({0})^{{-2}}".format(units1)
            result.GetZaxis().SetTitle(axisTitle)

        if not weighted and events > 0:
            result.Scale(1. / events, "width")
        else:
            result.Scale(1, "width")
        return result

    def BuildHistograms(self):
        if len(self.fAxis) == 1:
            self.fHistogram = DMesonJetUtils.BuildHistogram(self.fAxis, self.fName, "counts")
            self.fUncertainty = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_Unc".format(self.fName), "relative statistical uncertainty")
            if self.fAnalysisType == AnalysisType.SideBand:
                self.fSideBandHistograms = []
                self.fSideBandLeftHistograms = []
                self.fSideBandRightHistograms = []
                self.fSignalHistograms = []
                self.fMass = None
                self.fMassWidth = None
                self.fBackground = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_Bkg".format(self.fName), "background |#it{{m}} - <#it{{m}}>| < {0}#sigma".format(int(self.fBinCountSignalSigmas)))
                self.fSideBandWindowTotalHistogram = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_SideBandWindowTotal".format(self.fName), "counts")
                self.fSignalWindowTotalHistogram = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_SignalWindowTotal".format(self.fName), "counts")
            elif self.fAnalysisType == AnalysisType.LikeSign:
                self.fLikeSignHistograms = []
                self.fUnlikeSignHistograms = []
                self.fMass = None
                self.fMassWidth = None
                self.fBackground = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_Bkg".format(self.fName), "background |#it{{m}} - <#it{{m}}>| < {0}#sigma".format(int(self.fBinCountSignalSigmas)))
                self.fLikeSignTotalHistogram = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_LikeSignTotal".format(self.fName), "counts")
                self.fUnlikeSignTotalHistogram = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_UnlikeSignTotal".format(self.fName), "counts")
            elif self.fAnalysisType == AnalysisType.InvMassFit or self.fAnalysisType == AnalysisType.LikeSignFit:
                self.fMass = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_Mass".format(self.fName), "D^{0} mass (GeV/#it{c}^{2})")
                self.fMassWidth = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_MassWidth".format(self.fName), "D^{0} mass width (GeV/#it{c}^{2})")
                self.fBackground = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_Bkg".format(self.fName), "background |#it{m} - <#it{m}>| < 2#sigma")
            elif self.fAnalysisType == AnalysisType.Truth and self.fDMeson in self.fBinSet.fNeedInvMass:
                self.fMass = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_Mass".format(self.fName), "D^{0} mass (GeV/#it{c}^{2})")
                self.fMassWidth = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_MassWidth".format(self.fName), "D^{0} mass width (GeV/#it{c}^{2})")
        elif len(self.fAxis) == 2:
            self.fHistogram = DMesonJetUtils.BuildHistogram(self.fAxis, self.fName, "counts")
            self.fUncertainty = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_Unc".format(self.fName), "relative statistical uncertainty")
            if self.fAnalysisType == AnalysisType.SideBand:
                self.fSideBandHistograms = []
                self.fSignalHistograms = []
                self.fMass = None
                self.fMassWidth = None
                self.fBackground = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_Bkg".format(self.fName), "background |#it{{m}} - <#it{{m}}>| < {0}#sigma".format(int(self.fBinCountSignalSigmas)))
                self.fSideBandWindowTotalHistogram = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_SideBandWindowTotal".format(self.fName), "counts")
                self.fSignalWindowTotalHistogram = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_SignalWindowTotal".format(self.fName), "counts")
            elif self.fAnalysisType == AnalysisType.InvMassFit:
                self.fMass = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_Mass".format(self.fName), "mass")
                self.fBackground = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_Bkg".format(self.fName), "background |#it{m} - <#it{m}>| < 3#sigma")
        else:
            print("Cannot handle spectra with more than two axis!")
            return

    def Fill(self, dmeson, jet, w):
        values = []
        for axis in self.fAxis:
            if axis.fName == "jet_pt":
                values.append(jet.fPt)
            elif axis.fName == "jet_corrpt":
                values.append(jet.fCorrPt)
            elif axis.fName == "jet_bkgpt":
                values.append(jet.fPt - jet.fCorrPt)
            elif axis.fName == "jet_eta":
                values.append(jet.fEta)
            elif axis.fName == "jet_phi":
                values.append(jet.fPhi)
            elif axis.fName == "jet_n":
                values.append(jet.fN)
            elif axis.fName == "d_z":
                values.append(jet.fZ)
            elif axis.fName == "d_corrz":
                values.append(jet.fCorrZ)
            elif axis.fName == "d_pt":
                values.append(dmeson.fPt)
            elif axis.fName == "d_eta":
                values.append(dmeson.fEta)
            elif axis.fName == "d_phi":
                values.append(dmeson.fPhi)
        if len(values) == 1:
            self.fHistogram.Fill(values[0], w)
        elif len(values) == 2:
            self.fHistogram.Fill(values[0], values[1], w)
        elif len(values) == 3:
            self.fHistogram.Fill(values[0], values[1], values[2], w)
        else:
            print("Cannot handle histograms with more than two axis!")
