#!/usr/bin/env python
# python program to perform a D meson jet analysis

import copy
import collections
import os
import numpy

import ROOT

import DMesonJetCompare
import DMesonJetUtils
import DMesonJetCuts
from DMesonJetBase import AnalysisType
import BinSet
import Axis

globalList = []

class MatchingSpectrum:
    def __init__(self, trigger, dmeson, s):
        self.fName = "_".join([obj for obj in [trigger, dmeson, s["name"]] if obj])
        self.fTitle = s["title"]
        self.fVariable = s["variable"]
        if "bins" in s:
            binsx = s["bins"]
            binsy = s["bins"]
        else:
            binsx = s["binsx"]
            binsy = s["binsy"]
        self.fHistogram = ROOT.TH2D(self.fName, self.fTitle, len(binsx) - 1, numpy.array(binsx, dtype=numpy.float), len(binsy) - 1, numpy.array(binsy, dtype=numpy.float))
        self.fHistogram.GetXaxis().SetTitle(s["x_axis_title"])
        self.fHistogram.GetYaxis().SetTitle(s["y_axis_title"])
        self.fHistogram.GetZaxis().SetTitle(s["z_axis_title"])
        self.fHistogram.GetZaxis().SetRangeUser(s["zmin"], s["zmax"])
        self.fCuts = DMesonJetCuts.DMesonJetCuts(s["cuts"])
        self.fLogScale = s["log_scale"]

    def Fill(self, dmeson, jet1, jet2, w):
        if self.fCuts.ApplyCuts(dmeson, jet1):
            v1 = getattr(jet1, self.fVariable)
            v2 = getattr(jet2, self.fVariable)
            self.fHistogram.Fill(v1, v2, w)

    def MakeProjections(self):
        self.fProjectionX = self.fHistogram.ProjectionX("{}_xproj".format(self.fHistogram.GetName()))
        self.fProjectionY = self.fHistogram.ProjectionY("{}_yproj".format(self.fHistogram.GetName()))

    def CalculateInterval(self, th):
        tot = self.fProjectionX.Integral(1, self.fProjectionX.GetNbinsX())
        partial = 0.0
        for ibin in range(1, self.fProjectionY.GetNbinsX() + 1):
            partial += self.fProjectionY.GetBinContent(ibin)
            if partial / tot > th:
                break
        return self.fProjectionY.GetXaxis().GetBinLowEdge(1), self.fProjectionY.GetXaxis().GetBinUpEdge(ibin)

class DMesonJetMatchingContainer:

    def __init__(self, trigger, DMesonDef, spectra_definitions, jet_definitions):
        self.fDMesonDef = DMesonDef
        self.fTrigger = trigger
        self.fSpectraDefinitions = spectra_definitions
        self.fJetDefinitions = jet_definitions
        self.GenerateJetBranchNames()
        self.GenerateHistograms()

    def GenerateHistograms(self):
        self.fSpectra = []
        for s in self.fSpectraDefinitions:
            if not self.fDMesonDef in s["active_mesons"]:
                continue
            self.fSpectra.append(MatchingSpectrum(self.fTrigger, self.fDMesonDef, s))

    def GetSpectra(self):
        for s in self.fSpectra:
            yield s.fHistogram

    def GetJetBranches(self):
        return self.fJetBranches

    def GenerateJetBranchNames(self):
        if len(self.fJetDefinitions) != 2:
            print("Must provide 2 and only 2 jet definitions.")
            exit(1)
        self.fJetBranches = []
        for jet_def in self.fJetDefinitions:
            jetName = "Jet_AKT{}{}_pt_scheme".format(jet_def["type"], jet_def["radius"])
            self.fJetBranches.append(jetName)

    def Fill(self, event, eventWeight):
        dmeson = event.DmesonJet
        jet1 = getattr(event, self.fJetBranches[0])
        jet2 = getattr(event, self.fJetBranches[1])
        for spectrum in self.fSpectra:
            spectrum.Fill(dmeson, jet1, jet2, eventWeight)

class DMesonJetAnalysisEngine:

    def __init__(self, collision, trigger, dmeson, spectra_definitions, jets, projector):
        self.fCollision = collision
        self.fTrigger = trigger
        self.fDMeson = dmeson
        self.fJetDefinitions = jets
        self.fProjector = projector
        self.fCanvases = []
        self.fEvents = 0
        self.fIsWeighted = not (self.fProjector.fWeight == 1)
        self.fSpectraDefinitions = spectra_definitions

    def SaveRootFile(self, myfile):
        rlist = ROOT.TList()
        if self.fTrigger:
            rlist.SetName("{}_{}".format(self.fTrigger, self.fDMeson))
        else:
            rlist.SetName(self.fDMeson)
        if self.fHistEvents:
            rlist.Add(self.fHistEvents)

        for s in self.fDataContainer.fSpectra:
            print("Adding '{}'".format(s.fHistogram.GetName()))
            rlist.Add(s.fHistogram)
            print("Adding '{}'".format(s.fProjectionX.GetName()))
            rlist.Add(s.fProjectionX)
            print("Adding '{}'".format(s.fProjectionY.GetName()))
            rlist.Add(s.fProjectionY)

        if rlist.GetEntries() > 0:
            myfile.cd()
            rlist.Write(rlist.GetName(), ROOT.TObject.kSingleKey)

    def SavePlots(self, path, myformat):
        for c in self.fCanvases:
            if c: c.SaveAs("{0}/{1}.{2}".format(path, c.GetName(), myformat))

    def DoProjections(self):
        self.fDataContainer = DMesonJetMatchingContainer(self.fTrigger, self.fDMeson, self.fSpectraDefinitions, self.fJetDefinitions)
        self.fProjector.StartProjection(self.fTrigger, self.fDMeson, None, self.fDataContainer)

        self.fHistEvents = self.fProjector.fHistEvents
        if self.fHistEvents:
            self.fEvents = self.fHistEvents.GetBinContent(self.fHistEvents.GetXaxis().FindBin("Normalized"))

        th = 0.95
        for s in self.fDataContainer.fSpectra:
            s.MakeProjections()
            xmin, xmax = s.CalculateInterval(th)
            print("{}: {}%% in the interval: [{}, {}]".format(s.fName, th*100, xmin, xmax))
            print("The mean of the x axis is {}".format(s.fProjectionX.GetMean()))
            print("The mean of the y axis is {}".format(s.fProjectionY.GetMean()))
            print("The ratio is {}".format(s.fProjectionY.GetMean() / s.fProjectionX.GetMean()))

    def PlotSpectra(self):
        for s in self.fDataContainer.fSpectra:
            print("Plotting {}".format(s.fName))
            self.PlotSpectrum2D(s.fHistogram, s.fLogScale)
            self.PlotSpectrum1D(s.fProjectionX, s.fLogScale)
            self.PlotSpectrum1D(s.fProjectionY, s.fLogScale)

    def PlotSpectrum1D(self, s, logscale):
        cname = "{0}_canvas".format(s.GetName())

        c = ROOT.TCanvas(cname, cname)
        c.SetLeftMargin(0.12)
        if logscale:
            c.SetLogy()
        self.fCanvases.append(c)
        c.cd()
        h = s.DrawCopy()
        h.GetXaxis().SetTitleOffset(1)
        h.GetYaxis().SetTitleOffset(1)

        h.Scale(1.0 / h.Integral(1, h.GetNbinsX()), "width")
        h.GetYaxis().SetTitle("probability density")

        globalList.append(c)
        globalList.append(h)

    def PlotSpectrum2D(self, s, logscale):
        cname = "{0}_canvas".format(s.GetName())

        c = ROOT.TCanvas(cname, cname)
        c.SetRightMargin(0.18)
        c.SetLogz()
        self.fCanvases.append(c)
        c.cd()
        h = s.DrawCopy("colz")
        h.GetZaxis().SetTitleOffset(1.6)

        globalList.append(c)
        globalList.append(h)

class DMesonJetMatchingAnalysis:

    def __init__(self, name):
        self.fName = name
        self.fAnalysisEngine = []
        self.fCanvases = []

    def SetProjector(self, projector):
        self.fProjector = projector

    def StartAnalysis(self, config):
        self.fCollision = config["collision_system"]
        self.fJets = config["jets"]

        for trigger in config["trigger"]:
            for d_meson in config["d_meson"]:
                if trigger:
                    print("Projecting trigger {0}, D meson {1}".format(trigger, d_meson))
                else:
                    print("Projecting D meson {0}".format(d_meson))
                eng = DMesonJetAnalysisEngine(self.fCollision, trigger, d_meson,
                                                config["spectra"],
                                                self.fJets, self.fProjector)
                self.fAnalysisEngine.append(eng)
                eng.DoProjections()
                eng.PlotSpectra()

    def SaveRootFile(self, path):
        file = ROOT.TFile.Open("{0}/{1}.root".format(path, self.fName), "recreate")
        for eng in self.fAnalysisEngine:
            eng.SaveRootFile(file)
        file.Close()

    def SavePlots(self, path, format):
        fullPath = "{0}/{1}/{2}".format(path, self.fName, format)
        if not os.path.isdir(fullPath):
            os.makedirs(fullPath)
        for eng in self.fAnalysisEngine:
            eng.SavePlots(fullPath, format)

        for c in self.fCanvases:
            c.SaveAs("{0}/{1}.{2}".format(fullPath, c.GetName(), format))
