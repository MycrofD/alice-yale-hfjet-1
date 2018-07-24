#!/usr/bin/env python
# python base classes and utilities for D Meson jet analysis

import array
import copy
import collections
import numpy
import math

import ROOT

import DetectorResponseLoader
import DMesonJetCuts
import Axis
import Spectrum
from DMesonJetBase import AnalysisType

class BinMultiSet:

    def __init__(self):
        self.fBinSets = collections.OrderedDict()

    def GenerateInvMassRootLists(self):
        for binSet in self.fBinSets.itervalues():
            yield binSet.GenerateInvMassRootList()

    def AddBinSet(self, binSet):
        self.fBinSets[binSet.fBinSetName] = binSet

    def Initialize(self, dmeson, jtype, jradius, scheme, jtitle, reflections, inputPath):
        self.fDMeson = dmeson
        self.fJetType = jtype
        self.fJetRadius = jradius
        self.fJetRecoScheme = scheme
        self.fJetTitle = jtitle
        for k, binSet in self.fBinSets.items():
            if not dmeson in binSet.fActiveMesons:
                del self.fBinSets[k]
                continue
            r = binSet.Initialize(dmeson, jtype, jradius, scheme, jtitle, reflections, inputPath)
            if not r:
                del self.fBinSets[k]

    def SetWeightEfficiency(self, weight):
        for binSet in self.fBinSets.itervalues():
            binSet.fWeightEfficiency = weight

    def FindBin(self, dmeson, jet, dmesonDef):
        for binSet in self.fBinSets.itervalues():
            if not dmesonDef in binSet.fNeedInvMass:
                continue
            bins = binSet.FindBin(dmeson, jet)
            for bin, w in bins:
                yield bin, w

    def FindSpectra(self, dmeson, jet):
        for binSet in self.fBinSets.itervalues():
            if not binSet.fCuts.ApplyCuts(dmeson, jet):
                continue
            for s in binSet.fSpectra.itervalues():
                if s.fAnalysisType != AnalysisType.Truth:
                    continue
                yield s, binSet.fWeightEfficiency.GetEfficiencyWeight(dmeson, jet)

    def FindAllSpectra(self):
        for binSet in self.fBinSets.itervalues():
            for s in binSet.fSpectra.itervalues():
                yield s

    def ResetAllComaprisonDone(self):
        for binSet in self.fBinSets.itervalues():
            for s in binSet.fSpectra.itervalues():
                s.fComparisonDone = False

    def GetAllComparisonGroups(self):
        groups = set()
        for binSet in self.fBinSets.itervalues():
            for s in binSet.fSpectra.itervalues():
                if s.fCompare: 
                    groups.update(s.fCompare)
        return groups

class BinSet:

    def __init__(self, name, title, active_mesons, need_inv_mass, limitSetList, spectra, axis, cutList, efficiency, refl_templ_names, fitOptions):
        self.fBinSetName = name
        self.fTitle = title
        self.fBins = []
        self.fSpectraConfigs = spectra
        self.fSpectra = collections.OrderedDict()
        self.fAxis = axis
        self.fCuts = DMesonJetCuts.DMesonJetCuts(cutList)
        self.fFitOptions = fitOptions
        self.fEfficiency = efficiency
        self.fNeedInvMass = need_inv_mass
        self.fActiveMesons = active_mesons
        self.fBinCountSpectraAxis = dict()
        self.fLimitSetList = limitSetList
        self.fReflectionTemplates = None
        self.fReflOverSign = 0
        self.fReflectionTemplateNames = refl_templ_names

    def AmIJetDependent(self):
        for axis in self.fAxis:
            if "jet" in axis.fName or axis.fName == "d_z" or axis.fName == "d_corrz": return True
        for axis in self.fBinCountSpectraAxis:
            for a in axis:
                if "jet" in a.fName or a.fName == "d_z" or a.fName == "d_corrz": return True
        for cut in self.fCuts.fCuts:
            if cut["object"] == "jet": return True
        return False

    def LoadEfficiency(self, inputPath, dmeson, jetName, eff, axis):
        dmeson = dmeson.replace("_kSignalOnly", "")  # for backward compatibility
        dmeson = dmeson.replace("_SignalOnly", "")  # for >= vAN-20180222
        dmeson = dmeson.replace("_WrongPID", "")  # for backward compatibility
        dmeson = dmeson.replace("_D0Reflection", "")  # for >= vAN-20180222
        if eff and "file_name" in eff:
            eff["file_name"] = "{0}/{1}".format(inputPath, eff["file_name"])
        if axis:
            weightEfficiency = DetectorResponseLoader.DMesonJetEfficiency.fromConfigAdvanced(eff, "Prompt", dmeson, jetName, axis.GetNbins(), axis.GetBinsArray())
        else:
            weightEfficiency = DetectorResponseLoader.DMesonJetEfficiency.fromConfig(eff, "Prompt", dmeson, jetName)
        return weightEfficiency

    def LoadReflectionsForTemplateName(self, reflFileName, refl_fit):
        reflFileName = "{}_fitted_{}.root".format(reflFileName, refl_fit)
        print("Loading reflections from file '{}' for bin set '{}'".format(reflFileName, self.fBinSetName))
        file = ROOT.TFile(reflFileName, "read")
        if not file or file.IsZombie():
            print("BinSet.LoadReflectionsForTemplateName: Could not open file {}".format(reflFileName))
            exit(1)
        result = []
        for ibin in range(0, len(self.fBins)):
            reflHistoname = "histRflFitted{fit}_ptBin{bin}".format(fit=refl_fit, bin=ibin)  # name of template histo
            reflHisto = file.Get(reflHistoname)
            if not reflHisto:
                print("Histogram '{}' could not be loaded from file '{}'".format(reflHistoname, reflFileName))
                exit(1)
            MCSigHistoname = "histSgn_{0}".format(ibin)  # name of template histo
            MCSigHisto = file.Get(MCSigHistoname)
            if not MCSigHisto:
                print("Histogram '{}' could not be loaded from file '{}'".format(MCSigHistoname, reflFileName))
                exit(1)

            result.append((reflHisto, MCSigHisto))
        return result

    def LoadReflections(self, dmeson, jtype, jradius, reflections):
        # hard coded conversion from D meson name to D meson cut name
        dmeson_cuts = dmeson[3:]

        if "DPtBins" in self.fBinSetName:
            varname = "DPt"
        elif "JetPtBins" in self.fBinSetName:
            varname = "JetPt"
        elif "JetZBins" in self.fBinSetName:
            varname = "JetZ"

        if jtype and jradius:
            jet_def = "_{}_{}".format(jtype, jradius)
        else:
            jet_def = ""

        if self.fEfficiency:
            is_eff_corrected = "_efficiency"
        else:
            is_eff_corrected = ""

        reflFileName = "../rawYieldUnc/reflTemp/{refl_name}{is_eff_corrected}_{cuts}_{var}{jet_def}_{binlist_name}".format(refl_name=reflections,
                                                                                                             cuts=dmeson_cuts, jet_def=jet_def, binlist_name=self.fBinSetName,
                                                                                                             var=varname, is_eff_corrected=is_eff_corrected)

        self.fReflectionTemplates = dict()
        for refl_fit in self.fReflectionTemplateNames:
            self.fReflectionTemplates[refl_fit] = self.LoadReflectionsForTemplateName(reflFileName, refl_fit)

    def Initialize(self, dmeson, jtype, jradius, scheme, jtitle, reflections, inputPath):
        self.fName = "_".join(obj for obj in [dmeson, jtype, jradius, scheme, self.fBinSetName] if not obj == None)
        jetDep = self.AmIJetDependent()
        if jtype or jradius or scheme:
            if not jetDep: return False
        else:
            if jetDep: return False

        if jtype == "Full":
            for a in self.fAxis:
                a.fChargedJet = False

        if jtype and jradius and scheme:
            jetName = "Jet_AKT{0}{1}_{2}".format(jtype, jradius, scheme)
        else:
            # This is a temporary hack. The detector response analysis does not have an option for "no jet"
            jetName = "Jet_AKTChargedR040_pt_scheme"

        if "MCTruth" in dmeson:
            self.fWeightEfficiency = DetectorResponseLoader.DMesonJetEfficiency(None)
        else:
            self.fWeightEfficiency = self.LoadEfficiency(inputPath, dmeson, jetName, self.fEfficiency, None)

        for s in self.fSpectraConfigs:
            # skip spectra that do not have this D meson active
            if not dmeson in s["active_mesons"]: continue

            # add bin counting axis
            if "axis" in s:
                if len(s["axis"]) > 2:
                    print("Error: cannot do bin counting spectra (e.g. side band) with more than 2 axis. Spectrum {}".format(s["name"]))
                    exit(1)
                bin_count_spectra = []
                axis_name = s["axis"].keys()[0]
                axis_bins = s["axis"].values()[0]
                bcaxis = Axis.Axis(axis_name, axis_bins, "", (jtype != "Full"))
                bin_count_spectra.append(bcaxis)
                self.fBinCountSpectraAxis[s["name"]] = bin_count_spectra

            # configure efficiency
            if "efficiency" in s and not "MCTruth" in dmeson:
                effWeight = self.LoadEfficiency(inputPath, dmeson, jetName, s["efficiency"], self.fAxis[0])
            else:
                effWeight = DetectorResponseLoader.DMesonJetEfficiency(None)

            spectrum = Spectrum.Spectrum(s, dmeson, jtype, jradius, scheme, jtitle, self, effWeight)

            self.fSpectra[spectrum.fName] = spectrum

        limits = dict()
        self.AddBinsRecursive(self.fLimitSetList, limits)

        if reflections:
            if scheme != "pt_scheme":
                print("BinSet.Initialize: Error! Reflections have only been implemented for the pt scheme!")
                print("************* CONTINUE AT YOUR OWN RISK! *******************")
                #exit(1)
            self.LoadReflections(dmeson, jtype, jradius, reflections)

        return True

    def GenerateInvMassRootList(self):
        rlist = ROOT.TList()
        rlist.SetName(self.fName)
        for bin in self.fBins:
            if bin.fInvMassHisto: rlist.Add(bin.fInvMassHisto)
            if bin.fSelectionHisto: rlist.Add(bin.fSelectionHisto)
            for fitter in bin.fMassFitters.itervalues(): rlist.Add(fitter)
            for _, h in bin.fBinCountSpectra.itervalues(): rlist.Add(h)
        return rlist

    def AddBinsRecursive(self, limitSetList, limits):
        if len(limitSetList) > 0:
            (limitSetName, limitSet) = limitSetList[0]
            for min, max in zip(limitSet[:-1], limitSet[1:]):
                limits[limitSetName] = min, max
                self.AddBinsRecursive(limitSetList[1:], limits)
        else:
            bin = BinLimits(limits)
            bin.fBinCountSpectra = dict()
            for n, a in self.fBinCountSpectraAxis.iteritems():
                bin.fBinCountSpectra[n] = (a, None)
            self.fBins.append(bin)

    def FindBin(self, dmeson, jet):
        if not self.fCuts.ApplyCuts(dmeson, jet):
            return
        else:
            w = self.fWeightEfficiency.GetEfficiencyWeight(dmeson, jet)
        for bin in self.fBins:
            if bin.IsInBinLimits(dmeson, jet):
                yield bin, w


class BinLimits:

    def __init__(self, limits=dict()):
        self.fLimits = copy.deepcopy(limits)
        self.fInvMassHisto = None
        self.fSelectionHisto = None
        self.fMassFitter = None
        self.fMassFitters = dict()
        self.fBinCountSpectra = None
        self.fCounts = 0
        self.fSumw2 = 0
        self.fEntries = 0

    def IsWeighted(self):
        if math.fabs(self.fEntries - self.fCounts) < 1e-6:
            return False
        else:
            return True

    def Fill(self, dmeson, jet, w):
        self.fEntries += 1
        self.fCounts += w
        self.fSumw2 += w * w
        if self.fInvMassHisto: self.fInvMassHisto.Fill(dmeson.fInvMass, w)
        if self.fSelectionHisto: self.fSelectionHisto.Fill(dmeson.fSelectionType, w)

        for axis, hist in self.fBinCountSpectra.itervalues():
            obsVal = []
            for a in axis:
                if a.fName == "jet_pt":
                    obsVal.append(jet.fPt)
                elif a.fName == "d_z":
                    obsVal.append(jet.fZ)
                elif a.fName == "d_corrz":
                    obsVal.append(jet.fCorrZ)
                elif a.fName == "jet_n":
                    obsVal.append(jet.fN)
                elif a.fName == "jet_corrpt":
                    obsVal.append(jet.fCorrPt)
                elif a.fName == "jet_bkgpt":
                    obsVal.append(jet.fPt - jet.fCorrPt)
                elif a.fName == "d_pt":
                    obsVal.append(dmeson.fPt)

            if hist.GetDimension() == 2:
                hist.Fill(dmeson.fInvMass, obsVal[0], w)
            elif hist.GetDimension() == 3:
                hist.Fill(dmeson.fInvMass, obsVal[0], obsVal[1], w)

    def AddFromAxis(self, axis, binIndex):
        if binIndex == 0:
            min = axis.fBins[0]
            max = axis.fBins[-1]
        else:
            min = axis.fBins[binIndex - 1]
            max = axis.fBins[binIndex]
        self.fLimits[axis.fName] = min, max

    def SetMassFitter(self, fitter):
        self.fMassFitter = fitter
        if not fitter in self.fMassFitters:
            self.AddMassFitter(fitter)

    def AddMassFitter(self, fitter, name="default"):
        self.fMassFitters[name] = fitter

    def IsSameOf(self, bin):
        for name, (min, max) in self.fLimits.iteritems():
            if bin.fLimits.has_key(name):
                if not bin.fLimits[name] == (min, max):
                    return False
            else:
                return False
        return True

    def SetJetCorrPtLimits(self, min, max):
        self.fLimits["jet_corrpt"] = min, max

    def SetJetPtLimits(self, min, max):
        self.fLimits["jet_pt"] = min, max

    def SetDPtLimits(self, min, max):
        self.fLimits["d_pt"] = min, max

    def SetDZLimits(self, min, max):
        self.fLimits["d_z"] = min, max

    def SetDCorrZLimits(self, min, max):
        self.fLimits["d_corrz"] = min, max

    def SetJetEtaLimits(self, min, max):
        self.fLimits["jet_eta"] = min, max

    def SetDEtaLimits(self, min, max):
        self.fLimits["d_eta"] = min, max

    def IsInBinLimits(self, dmeson, jet):
        for name, (min, max) in self.fLimits.iteritems():
            if not min < max:
                continue
            if name == "d_pt" and (dmeson.fPt < min or dmeson.fPt >= max):
                return False
            elif name == "jet_pt" and (jet.fPt < min or jet.fPt >= max):
                return False
            elif name == "jet_corrpt" and (jet.fCorrPt < min or jet.fCorrPt >= max):
                return False
            elif name == "d_eta" and (dmeson.fEta < min or dmeson.fEta >= max):
                return False
            elif name == "jet_eta" and (jet.fEta < min or jet.fEta >= max):
                return False
            elif name == "d_z" and (jet.fZ < min or jet.fZ >= max):
                return False
            elif name == "d_corrz" and (jet.fCorrZ < min or jet.fCorrZ >= max):
                return False

        return True

    def GetBinCenter(self, axis):
        if axis in self.fLimits:
            (min, max) = self.fLimits[axis]
            return (min + max) / 2
        else:
            return -1

    def GetName(self):
        name = ""
        for varName, (min, max) in self.fLimits.iteritems():
            if varName == "d_pt":
                name += "DPt_{0}_{1}_".format(int(min * 100), int(max * 100))
            elif varName == "jet_pt":
                name += "JetPt_{0}_{1}_".format(int(min * 100), int(max * 100))
            elif varName == "jet_corrpt":
                name += "JetCorrPt_{0}_{1}_".format(int(min * 100), int(max * 100)).replace("-", "N")
            elif varName == "d_eta":
                name += "DEta_{0}_{1}_".format(int(min * 10), int(max * 10)).replace("-", "N")
            elif varName == "jet_eta":
                name += "JetEta_{0}_{1}_".format(int(min * 10), int(max * 10)).replace("-", "N")
            elif varName == "d_z":
                name += "JetZ_{0}_{1}_".format(int(min * 100), int(max * 100))
            elif varName == "d_corrz":
                name += "JetCorrZ_{0}_{1}_".format(int(min * 100), int(max * 100))

        # remove last "_"
        if name: name = name[:-1]

        return name

    def GetTitle(self):
        title = ""

        for varName, (min, max) in self.fLimits.iteritems():
            if varName == "d_pt":
                title += "{0:.1f} < #it{{p}}_{{T,D}} < {1:.1f} GeV/#it{{c}}, ".format(min, max)
            elif varName == "jet_pt":
                title += "{0:.0f} < #it{{p}}_{{T,jet}} < {1:.0f} GeV/#it{{c}}, ".format(min, max)
            elif varName == "jet_corrpt":
                title += "{0:.0f} < #it{{p}}_{{T,jet}}^{{sub}} < {1:.0f} GeV/#it{{c}}, ".format(min, max)
            elif varName == "d_eta":
                title += "{0:.1f} < #it{{#eta}}_{{D}} < {1:.1f}, ".format(min, max)
            elif varName == "jet_eta":
                title += "{0:.0f} < #it{{#eta}}_{{jet}} < {1:.0f}, ".format(min, max)
            elif varName == "d_z":
                title += "{0:.1f} < #it{{z}}_{{||, D}} < {1:.1f}, ".format(min, max)
            elif varName == "d_corrz":
                title += "{0:.1f} < #it{{z}}_{{||, D}}^{{sub}} < {1:.1f}, ".format(min, max)

        # remove last ", "
        if title: title = title[:-2]
        return title

    def Print(self):
        print(self.GetTitle())

    def CreateHistograms(self, trigger, DMesonDef, xAxis, yAxis, nMassBins, minMass, maxMass):
        self.CreateInvMassHisto(trigger, DMesonDef, xAxis, yAxis, nMassBins, minMass, maxMass)
        self.CreateQAHistos(trigger, DMesonDef, yAxis)

    def CreateQAHistos(self, trigger, DMesonDef, yAxis):
        xAxis = "Selection Type"
        if trigger:
            hname = "Selection_{0}_{1}_{2}".format(trigger, DMesonDef, self.GetName())
            htitle = "{0} - {1} Selection Type: {2};{3};{4}".format(trigger, DMesonDef, self.GetTitle(), xAxis, yAxis)
        else:
            hname = "Selection_{0}_{1}".format(DMesonDef, self.GetName())
            htitle = "{0} Selection Type: {1};{2};{3}".format(DMesonDef, self.GetTitle(), xAxis, yAxis)
        self.fSelectionHisto = ROOT.TH1D(hname, htitle, 3, 0.5, 3.5)

    def CreateInvMassHisto(self, trigger, DMesonDef, xAxis, yAxis, nMassBins, minMass, maxMass):
        if trigger:
            hname = "InvMass_{0}_{1}_{2}".format(trigger, DMesonDef, self.GetName())
            htitle = "{0} - {1} Invariant Mass: {2};{3};{4}".format(trigger, DMesonDef, self.GetTitle(), xAxis, yAxis)
        else:
            hname = "InvMass_{0}_{1}".format(DMesonDef, self.GetName())
            htitle = "{0} Invariant Mass: {1};{2};{3}".format(DMesonDef, self.GetTitle(), xAxis, yAxis)

        if not "MCTruth" in DMesonDef:
            self.fInvMassHisto = ROOT.TH1D(hname, htitle, nMassBins, minMass, maxMass)
            self.fInvMassHisto.Sumw2()
            self.fInvMassHisto.SetMarkerSize(0.9)
            self.fInvMassHisto.SetMarkerStyle(ROOT.kFullCircle)
            self.fInvMassHisto.SetMarkerColor(ROOT.kBlue + 2)
            self.fInvMassHisto.SetLineColor(ROOT.kBlue + 2)

        massAxis = array.array('d', numpy.linspace(minMass, maxMass, nMassBins * 100 + 1, True))
        for sname, (axis, _) in self.fBinCountSpectra.iteritems():
            if trigger: hnameSB = "InvMassBinCounting_{}_{}_{}_{}".format(trigger, DMesonDef, self.GetName(), sname)
            else: hnameSB = "InvMassBinCounting_{}_{}_{}".format(DMesonDef, self.GetName(), sname)
            if len(axis) == 2:
                if trigger: htitleSB = "{} - {} Invariant Mass: {};{};{};{}".format(trigger, DMesonDef, self.GetTitle(), xAxis, axis[0].GetTitle(), axis[1].GetTitle())
                else: htitleSB = "{} Invariant Mass: {};{};{};{}".format(DMesonDef, self.GetTitle(), xAxis, axis[0].GetTitle(), axis[1].GetTitle())
                h = ROOT.TH3D(hnameSB, htitleSB, nMassBins * 100, massAxis, len(axis[0].fBins) - 1, array.array('d', axis[0].fBins), len(axis[1].fBins) - 1, array.array('d', axis[1].fBins))
                h.Sumw2()
            elif len(axis) == 1:
                if trigger: htitleSB = "{} - {} Invariant Mass: {};{};{}".format(trigger, DMesonDef, self.GetTitle(), xAxis, axis[0].GetTitle())
                else: htitleSB = "{} Invariant Mass: {};{};{}".format(DMesonDef, self.GetTitle(), xAxis, axis[0].GetTitle())
                h = ROOT.TH2D(hnameSB, htitleSB, nMassBins * 100, massAxis, len(axis[0].fBins) - 1, array.array('d', axis[0].fBins))
                h.Sumw2()
            else:
                print("Error!! Bin counting histograms only allowed for a maximum of 2 axis. Spectrum {}".format(sname))
                exit(1)
            self.fBinCountSpectra[sname] = (axis, h)
