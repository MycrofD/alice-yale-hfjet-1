#!/usr/bin/env python
# python program to perform a D meson jet analysis

import math
import copy
import collections
import os

import ROOT

import DMesonJetProjectors
import DMesonJetCompare
import DMesonJetUtils
from DMesonJetBase import AnalysisType
import BinSet
import Axis

globalList = []

class DMesonJetVariable:
    def __init__(self, name, varname, bins, low, high):
        self.fName = name
        self.fVariableName = varname
        self.fNBins = bins
        self.fLow = low
        self.fHigh = high
        self.fHistogram = ROOT.TH1F(self.fName, self.fName, self.fNBins, self.fLow, self.fHigh)
        self.fHistogram.Sumw2()

    def FillStd(self, dmeson, w):
        v = getattr(dmeson, self.fVariableName)
        self.fHistogram.Fill(v, w)

    def FillAbs(self, dmeson, w):
        v = math.fabs(getattr(dmeson, self.fVariableName))
        self.fHistogram.Fill(v, w)

    def Analyze(self):
        integral = self.fHistogram.Integral(0, -1)
        if integral == 0: return
        self.fHistogram.Scale(1. / integral)
        self.fCutEfficiency = self.fHistogram.Clone("{}_CutEfficiency".format(self.fName))
        self.fCutEfficiency.Reset()
        self.CalculateCutEfficiency()

    def CalculateLeftCutEfficiency(self):
        partialIntErr = 0.
        partialInt = 0.
        for ibin in xrange(0, self.fHistogram.GetNbinsX() + 2):
            partialIntErr = math.sqrt(partialIntErr ** 2 + self.fHistogram.GetBinError(ibin) ** 2)
            partialInt += self.fHistogram.GetBinContent(ibin)
            self.fCutEfficiency.SetBinContent(ibin, partialInt)
            self.fCutEfficiency.SetBinError(ibin, partialIntErr)

    def CalculateRightCutEfficiency(self):
        partialIntErr = 0.
        partialInt = 0.
        for ibin in reversed(xrange(0, self.fHistogram.GetNbinsX() + 2)):
            partialIntErr = math.sqrt(partialIntErr ** 2 + self.fHistogram.GetBinError(ibin) ** 2)
            partialInt += self.fHistogram.GetBinContent(ibin)
            self.fCutEfficiency.SetBinContent(ibin, partialInt)
            self.fCutEfficiency.SetBinError(ibin, partialIntErr)

    @classmethod
    def DCA(cls):
        obj = cls("DCA", "fDCA", 30, 0, 0.3)
        obj.Fill = obj.FillStd
        obj.CalculateCutEfficiency = obj.CalculateLeftCutEfficiency
        return obj

    @classmethod
    def CosThetaStar(cls):
        obj = cls("CosThetaStar", "fCosThetaStar", 50, 0, 1)
        obj.Fill = obj.FillAbs
        obj.CalculateCutEfficiency = obj.CalculateLeftCutEfficiency
        return obj

    @classmethod
    def d0d0(cls):
        obj = cls("d0d0", "fd0d0", 100, -5e-5, 1e-4)
        obj.Fill = obj.FillStd
        obj.CalculateCutEfficiency = obj.CalculateLeftCutEfficiency
        return obj

    @classmethod
    def MaxNormd0(cls):
        obj = cls("MaxNormd0", "fMaxNormd0", 1000, 0, 100)
        obj.Fill = obj.FillAbs
        obj.CalculateCutEfficiency = obj.CalculateLeftCutEfficiency
        return obj

    @classmethod
    def CosPointing(cls):
        obj = cls("CosPointing", "fCosPointing", 50, 0, 1)
        obj.Fill = obj.FillStd
        obj.CalculateCutEfficiency = obj.CalculateRightCutEfficiency
        return obj

DMesonJetVariable.fgVariableList = [DMesonJetVariable.DCA, DMesonJetVariable.CosThetaStar, DMesonJetVariable.d0d0, DMesonJetVariable.MaxNormd0, DMesonJetVariable.CosPointing]

class DMesonJetTopoContainer:
    def __init__(self, jet_pt_bins, jet_type, jet_radius):
        self.fJetPtBins = jet_pt_bins
        self.fVariables = dict()
        self.fJetType = jet_type
        self.fJetRadius = jet_radius
        for jetPtBin in self.fJetPtBins:
            self.fVariables[jetPtBin] = self.GenerateVariables()

    def AddVariable(self, variables, v):
        variables[v.fName] = v

    def GenerateVariables(self):
        variables = dict()

        for var in DMesonJetVariable.fgVariableList:
            self.AddVariable(variables, var())

        return variables

    def Fill(self, event, eventWeight):
        dmeson = event.DmesonJet

        if self.fJetType or self.fJetRadius:
            jetName = "Jet_AKT{0}{1}_pt_scheme".format(self.fJetType, self.fJetRadius)
            jet = getattr(event, jetName)
        else:
            jet = None

        for (jetPtMin, jetPtMax) in self.fJetPtBins:
            if jet.fPt < jetPtMin or jet.fPt >= jetPtMax: continue
            for var in self.fVariables[(jetPtMin, jetPtMax)].itervalues():
                var.Fill(dmeson, eventWeight)

class DMesonJetTopoAnalysis:
    def __init__(self, name, trigger, dmeson, jet_pt_bins, projector):
        self.fName = name
        self.fTrigger = trigger
        self.fDMeson = dmeson
        self.fProjector = projector
        self.fJetPtBins = jet_pt_bins

    def SaveRootFile(self, file):
        rlist = ROOT.TList()
        if self.fTrigger:
            rlist.SetName("{}_{}".format(self.fTrigger, self.fDMeson))
        else:
            rlist.SetName(self.fDMeson)

        for (jetPtMin, jetPtMax), histograms in self.fHistograms.iteritems():
            hListName = "Histograms_{}_{}".format(jetPtMin, jetPtMax)
            hList = ROOT.TList()
            hList.SetName(jetName)
            for h in histograms:
                hList.Add(h)
            rlist.Add(hList)

        for (jetPtMin, jetPtMax), histograms in self.fCutEfficiency.iteritems():
            hListName = "CutEfficiency_{}_{}".format(jetPtMin, jetPtMax)
            hList = ROOT.TList()
            hList.SetName(jetName)
            for h in histograms:
                hList.Add(h)
            rlist.Add(hList)

        if rlist.GetEntries() > 0:
            file.cd()
            rlist.Write(rlist.GetName(), ROOT.TObject.kSingleKey)

    def DoProjections(self):
        data_container = DMesonJetTopoContainer(self.fJetPtBins, "Charged", "R040")
        self.fProjector.StartProjection(self.fTrigger, self.fDMeson, None, data_container)
        self.fVariables = data_container.fVariables

    def Start(self):
        self.DoProjections()
        self.Analyze()

    def Analyze(self):
        for (jetPtMin, jetPtMax), variables in self.fVariables.iteritems():
            for v in variables.itervalues():
                v.Analyze()

class DMesonJetTopoAnalysisManager:
    def __init__(self, dmeson):
        self.fDMeson = dmeson
        self.fAnalysisList = dict()
        self.fJetPtBins = [(5, 15), (15, 30)]

    def AddAnalysis(self, name, trigger, projector, dmesonSuffix):
        self.fAnalysisList[name] = DMesonJetTopoAnalysis(name, trigger, "{}_{}".format(self.fDMeson, dmesonSuffix), self.fJetPtBins, projector)

    def StartAnalysis(self):
        for ana in self.fAnalysisList.itervalues():
            ana.Start()
        self.Plot()

    def SaveRootFile(self, path):
        pass

    def SavePlots(self, path, format):
        for c in self.fCanvases:
            if c: c.SaveAs("{0}/{1}.{2}".format(path, c.GetName(), format))

    def Plot(self):
        self.fCanvases = []
        self.fCompareObjects = []

        for (jetPtMin, jetPtMax) in self.fJetPtBins:
            for variableName in [varFunct().fName for varFunct in DMesonJetVariable.fgVariableList]:
                hVariable = []
                hEfficiency = []
                for ana in self.fAnalysisList.itervalues():
                    hVar = ana.fVariables[(jetPtMin, jetPtMax)][variableName].fHistogram
                    hVar.SetTitle(ana.fName)
                    hVariable.append(hVar)
                    hEff = ana.fVariables[(jetPtMin, jetPtMax)][variableName].fCutEfficiency
                    hEff.SetTitle(ana.fName)
                    hEfficiency.append(hEff)

                cname = "{}_{}_{}".format(hEfficiency[0].GetName(), jetPtMin, jetPtMax)
                comp = DMesonJetCompare.DMesonJetCompare(cname)
                comp.fDoSpectraPlot = "lineary"
                comp.fDoRatioPlot = False
                comp.CompareSpectra(hEfficiency[0], hEfficiency[1:])
                self.fCompareObjects.append(comp)
                self.fCanvases.append(comp.fCanvasSpectra)

                cname = "{}_{}_{}".format(hVariable[0].GetName(), jetPtMin, jetPtMax)
                comp = DMesonJetCompare.DMesonJetCompare(cname)
                comp.fDoSpectraPlot = "logy"
                comp.fDoRatioPlot = False
                comp.CompareSpectra(hVariable[0], hVariable[1:])
                self.fCompareObjects.append(comp)
                self.fCanvases.append(comp.fCanvasSpectra)
