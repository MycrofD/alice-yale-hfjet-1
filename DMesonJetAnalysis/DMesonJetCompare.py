#!/usr/bin/env python
# python utilities for D Meson jet analysis

import ROOT
import math
import DMesonJetUtils
import random

class DMesonJetCompare:
    def __init__(self, name):
        self.fName = name
        self.fBaselineHistogram = None
        self.fHistograms = None
        self.fRatios = []
        self.fFitResults = []
        self.fOptSpectrumBaseline = ""
        self.fOptSpectrum = ""
        self.fOptRatio = ""
        self.fYaxisRatio = "Ratio"
        self.fDoSpectraPlot = "logy"
        self.fDoRatioPlot = "lineary"
        self.fCanvasSpectra = None
        self.fCanvasRatio = None
        self.fLegendSpectra = None
        self.fLegendRatio = None
        self.fBaselineRatio = None
        self.fColors = [ROOT.kBlack, ROOT.kBlue + 2, ROOT.kRed + 2, ROOT.kGreen + 2, ROOT.kOrange + 2, ROOT.kAzure + 2, ROOT.kMagenta + 2, ROOT.kCyan + 2, ROOT.kPink + 1, ROOT.kTeal + 2, ROOT.kYellow + 2]
        self.fMarkers = [ROOT.kOpenCircle, ROOT.kFullCircle, ROOT.kFullSquare, ROOT.kFullTriangleUp, ROOT.kFullTriangleDown, ROOT.kFullDiamond, ROOT.kFullStar, ROOT.kStar, ROOT.kFullCross, ROOT.kMultiply, ROOT.kPlus]
        self.fLines = [1, 2, 9, 5, 7, 10, 4, 3, 6, 8, 9]
        self.fFills = [3001, 3002, 3003, 3004, 3005, 3006, 3007]
        self.fMainHistogram = None
        self.fMainRatioHistogram = None
        self.fMaxSpectrum = None
        self.fMinSpectrum = None
        self.fMaxRatio = None
        self.fMinRatio = None
        self.fRatioRelativeUncertainty = None
        self.fResults = None
        self.fNColsLegRatio = 1
        self.fNColsLegSpectrum = 1
        self.fX1LegRatio = 0.55
        self.fY1LegRatio = 0.87
        self.fX1LegSpectrum = 0.55
        self.fX2LegSpectrum = 0.90
        self.fY1LegSpectrum = 0.87
        self.fLogUpperSpace = 10  # this factor will be used to adjust the y axis in log scale
        self.fLogLowerSpace = 2  # this factor will be used to adjust the y axis in log scale
        self.fLinUpperSpace = 0.9  # this factor will be used to adjust the y axis in linear scale
        self.fLinLowerSpace = 0.2  # this factor will be used to adjust the y axis in linear scale
        self.fLegTextSize = 20
        self.fLegLineHeight = 0.04
        self.fBaselineForRatio = None
        self.fSeparateBaselineUncertainty = False
        self.fNoErrorInBaseline = False
        self.fRatioRelativeUncertaintyTitle = "Rel. Unc."
        self.fGridyRatio = False
        self.fFitFunction = "expo(0)+expo(2)"
        self.fDoSpectrumLegend = True
        self.fDoRatioLegend = True

    def SetRatioRelativeUncertaintyFromHistogram(self, hist):
        self.fRatioRelativeUncertainty = hist.Clone("{0}_unc".format(hist.GetName()))
        self.fRatioRelativeUncertainty.SetTitle(self.fRatioRelativeUncertaintyTitle)
        for ibin in range(0, self.fRatioRelativeUncertainty.GetNbinsX() + 2):
            if self.fRatioRelativeUncertainty.GetBinContent(ibin) != 0:
                self.fRatioRelativeUncertainty.SetBinError(ibin, self.fRatioRelativeUncertainty.GetBinError(ibin) / self.fRatioRelativeUncertainty.GetBinContent(ibin))
                self.fRatioRelativeUncertainty.SetBinContent(ibin, 1)
            else:
                self.fRatioRelativeUncertainty.SetBinError(ibin, 0)
                self.fRatioRelativeUncertainty.SetBinContent(ibin, 0)

    def PrepareSpectraCanvas(self):
        if not self.fCanvasSpectra:
            print("Creating new canvas {0}".format(self.fName))
            self.fCanvasSpectra = ROOT.TCanvas(self.fName, self.fName)

        if self.fDoSpectraPlot == "logy":
            self.fCanvasSpectra.SetLogy()

        if self.fDoSpectrumLegend:
            if self.fLegendSpectra:
                y1 = self.fLegendSpectra.GetY1() - self.fLegLineHeight * (len(self.fHistograms) + 1) / self.fNColsLegSpectrum
                if y1 < 0.2: y1 = 0.2
                self.fLegendSpectra.SetY1(y1)
            else:
                y1 = self.fY1LegSpectrum - self.fLegLineHeight * (len(self.fHistograms) + 1) / self.fNColsLegSpectrum
                if y1 < 0.2: y1 = 0.2
                self.fLegendSpectra = ROOT.TLegend(self.fX1LegSpectrum, y1, self.fX2LegSpectrum, self.fY1LegSpectrum)
                self.fLegendSpectra.SetName("{0}_legend".format(self.fCanvasSpectra.GetName()))
                self.fLegendSpectra.SetNColumns(self.fNColsLegSpectrum)
                self.fLegendSpectra.SetFillStyle(0)
                self.fLegendSpectra.SetBorderSize(0)
                self.fLegendSpectra.SetTextFont(43)
                self.fLegendSpectra.SetTextSize(self.fLegTextSize)

        if "hist" in self.fOptSpectrumBaseline:
            self.fBaselineHistogram.SetLineColor(self.fColors[0])
            self.fBaselineHistogram.SetLineWidth(2)
            self.fBaselineHistogram.SetLineStyle(self.fLines[0])
            if self.fDoSpectrumLegend: self.fLegendSpectra.AddEntry(self.fBaselineHistogram, self.fBaselineHistogram.GetTitle(), "l")
        elif "e2" in self.fOptSpectrumBaseline:
            self.fBaselineHistogram.SetLineColor(self.fColors[0])
            self.fBaselineHistogram.SetFillColor(self.fColors[0])
            self.fBaselineHistogram.SetLineWidth(1)
            self.fBaselineHistogram.SetLineStyle(self.fLines[0])
            self.fBaselineHistogram.SetFillStyle(self.fFills[0])
            if self.fDoSpectrumLegend: self.fLegendSpectra.AddEntry(self.fBaselineHistogram, self.fBaselineHistogram.GetTitle(), "f")
        else:
            self.fBaselineHistogram.SetMarkerColor(self.fColors[0])
            self.fBaselineHistogram.SetLineColor(self.fColors[0])
            self.fBaselineHistogram.SetMarkerStyle(self.fMarkers[0])
            self.fBaselineHistogram.SetMarkerSize(1.2)
            if self.fDoSpectrumLegend: self.fLegendSpectra.AddEntry(self.fBaselineHistogram, self.fBaselineHistogram.GetTitle(), "pe")

        print("Plotting histogram '{0}' with option '{1}'".format(self.fBaselineHistogram.GetName(), self.fOptSpectrumBaseline))
        self.fCanvasSpectra.cd()
        self.fBaselineHistogram.Draw(self.fOptSpectrumBaseline)

        if "frac" in self.fBaselineHistogram.GetYaxis().GetTitle():
            self.fCanvasSpectra.SetLeftMargin(0.12)
            self.fBaselineHistogram.GetYaxis().SetTitleOffset(1.4)
        if not "same" in self.fOptSpectrum:
            self.fOptSpectrum += "same"

        if not self.fMainHistogram:
            self.fMainHistogram = self.fBaselineHistogram
        self.fBaselineForRatio = self.fBaselineHistogram.Clone("{0}_copy".format(self.fBaselineHistogram.GetName()))

        m = DMesonJetUtils.FindMinimum(self.fBaselineHistogram, 0., not "hist" in self.fOptSpectrumBaseline)
        if not m is None:
            if self.fMinSpectrum is None:
                self.fMinSpectrum = m
            else:
                self.fMinSpectrum = min(self.fMinSpectrum, m)
        m = DMesonJetUtils.FindMaximum(self.fBaselineHistogram, 0., not "hist" in self.fOptSpectrumBaseline)
        if not m is None:
            if self.fMaxSpectrum is None:
                self.fMaxSpectrum = m
            else:
                self.fMaxSpectrum = max(self.fMaxSpectrum, m)

    def PrepareRatioCanvas(self):
        cname = "{0}_Ratio".format(self.fName)
        if not self.fCanvasRatio:
            self.fCanvasRatio = ROOT.TCanvas(cname, cname)
        self.fCanvasRatio.cd()
        if self.fDoRatioPlot == "logy":
            self.fCanvasRatio.SetLogy()

        n = len(self.fHistograms)
        if self.fRatioRelativeUncertainty or self.fSeparateBaselineUncertainty: n += 1

        if self.fDoRatioLegend:
            if self.fLegendRatio:
                y1 = self.fLegendRatio.GetY1() - self.fLegLineHeight * n / self.fNColsLegRatio
                if y1 < 0.2: y1 = 0.2
                self.fLegendRatio.SetY1(y1)
            else:
                y1 = self.fY1LegRatio - self.fLegLineHeight * n / self.fNColsLegRatio
                if y1 < 0.2: y1 = 0.2
                self.fLegendRatio = ROOT.TLegend(self.fX1LegRatio, y1, 0.9, self.fY1LegRatio)
                self.fLegendRatio.SetName("{0}_legend".format(self.fCanvasRatio.GetName()))
                self.fLegendRatio.SetNColumns(self.fNColsLegRatio)
                self.fLegendRatio.SetFillStyle(0)
                self.fLegendRatio.SetBorderSize(0)
                self.fLegendRatio.SetTextFont(43)
                self.fLegendRatio.SetTextSize(self.fLegTextSize)

        if self.fGridyRatio: self.fCanvasRatio.SetGridy()

        if self.fSeparateBaselineUncertainty or self.fNoErrorInBaseline:
            for ibin in range(0, self.fBaselineForRatio.GetNbinsX() + 2):
                self.fBaselineForRatio.SetBinError(ibin, 0)

        if self.fSeparateBaselineUncertainty:
            self.SetRatioRelativeUncertaintyFromHistogram(self.fBaselineHistogram)
            opt = "e2"
            if "same" in self.fOptRatio:
                opt += "same"
            h = self.fRatioRelativeUncertainty.DrawCopy(opt)
            h.SetFillColor(self.fColors[0])
            h.SetFillStyle(self.fFills[0])
            h.SetLineColor(self.fColors[0])
            h.GetYaxis().SetTitle(self.fYaxisRatio)
            if self.fDoRatioLegend: self.fLegendRatio.AddEntry(h, h.GetTitle(), "f")
            self.fResults.append(h)
            if not "same" in self.fOptRatio:
                self.fOptRatio += "same"
            if not self.fMainRatioHistogram:
                self.fMainRatioHistogram = h
            m = DMesonJetUtils.FindMinimum(h, 0., True)
            if not m is None:
                if self.fMinRatio is None:
                    self.fMinRatio = m
                else:
                    self.fMinRatio = min(self.fMinRatio, m)
            m = DMesonJetUtils.FindMaximum(h, 0., True)
            if not m is None:
                if self.fMaxRatio is None:
                    self.fMaxRatio = m
                else:
                    self.fMaxRatio = max(self.fMaxRatio, m)

    def PlotHistogram(self, color, marker, line, h):
        m = DMesonJetUtils.FindMinimum(h, 0., not "hist" in self.fOptSpectrum)
        if not m is None:
            if self.fMinSpectrum is None:
                self.fMinSpectrum = m
            else:
                self.fMinSpectrum = min(self.fMinSpectrum, m)
        m = DMesonJetUtils.FindMaximum(h, 0., not "hist" in self.fOptSpectrum)
        if not m is None:
            if self.fMaxSpectrum is None:
                self.fMaxSpectrum = m
            else:
                self.fMaxSpectrum = max(self.fMaxSpectrum, m)

        print("Plotting histogram '{0}' with option '{1}'".format(h.GetName(), self.fOptSpectrum))
        self.fCanvasSpectra.cd()
        h.Draw(self.fOptSpectrum)

        if "hist" in self.fOptSpectrum:
            h.SetLineColor(color)
            h.SetLineWidth(2)
            h.SetLineStyle(line)
            if self.fDoSpectrumLegend: self.fLegendSpectra.AddEntry(h, h.GetTitle(), "l")
        else:
            h.SetMarkerColor(color)
            h.SetLineColor(color)
            h.SetMarkerStyle(marker)
            h.SetMarkerSize(1.2)
            if self.fDoSpectrumLegend: self.fLegendSpectra.AddEntry(h, h.GetTitle(), "pe")

    def FitAndMakeConsistent(self, h, templateH):
        fit_func = ROOT.TF1("{0}_fit".format(h.GetName()), self.fFitFunction, h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
        fit_func.SetParameter(0, 1)
        fit_func.SetParameter(1, -0.5)
        fit_func.SetParameter(2, 0)
        fit_func.SetParameter(3, 0)
        fit_func.SetParameter(0, 1. / fit_func.Eval(h.GetXaxis().GetBinCenter(1)))
        fit_func.SetParameter(2, 1)
        fit_func.SetParameter(3, -0.2)
        fit_func.SetParameter(2, 1. / fit_func.Eval(h.GetXaxis().GetBinCenter(int(h.GetNbinsX() / 2))))
        fitR = h.Fit(fit_func, "NS")
        fitOk = int(fitR)
        if not fitOk == 0: return None
        h_fit = DMesonJetUtils.soft_clone(templateH, "{0}_fith".format(h.GetName()))
        for ibin in range(1, h_fit.GetNbinsX() + 1):
            valErr = fit_func.IntegralError(h_fit.GetXaxis().GetBinLowEdge(ibin), h_fit.GetXaxis().GetBinUpEdge(ibin)) / (h_fit.GetXaxis().GetBinUpEdge(ibin) - h_fit.GetXaxis().GetBinLowEdge(ibin))
            val = fit_func.Integral(h_fit.GetXaxis().GetBinLowEdge(ibin), h_fit.GetXaxis().GetBinUpEdge(ibin)) / (h_fit.GetXaxis().GetBinUpEdge(ibin) - h_fit.GetXaxis().GetBinLowEdge(ibin))
            print("integral = {0:.5f}, central = {1:.5f}".format(val, fit_func.Eval((h_fit.GetXaxis().GetBinCenter(ibin)))))
            h_fit.SetBinContent(ibin, val)
            h_fit.SetBinError(ibin, valErr)
        self.fFitResults.append(fit_func)
        self.fFitResults.append(h_fit)
        fit_func.SetLineColor(h.GetLineColor())
        self.fCanvasSpectra.cd()
        fit_func.Draw("same")
        return h_fit

    def PlotRatio(self, color, marker, line, h):
        if self.CheckConsistency(h, self.fBaselineForRatio):
            hRatio = h.Clone("{0}_Ratio".format(h.GetName()))
        else:
            print("Trying to fit histogram {0} with function {1}".format(h.GetName(), self.fFitFunction))
            hRatio = self.FitAndMakeConsistent(h, self.fBaselineForRatio)
            if not hRatio:
                print("Fin unsuccessfull!")
                return
        hRatio.GetYaxis().SetTitle(self.fYaxisRatio)
        if not self.fBaselineRatio:
            self.fBaselineRatio = hRatio
        if "hist" in self.fOptRatio:
            hRatio.SetLineColor(color)
            hRatio.SetLineWidth(3)
            hRatio.SetLineStyle(line)
            if self.fDoRatioLegend: self.fLegendRatio.AddEntry(hRatio, h.GetTitle(), "l")
        else:
            hRatio.SetMarkerColor(color)
            hRatio.SetLineColor(color)
            hRatio.SetMarkerStyle(marker)
            hRatio.SetMarkerSize(1.2)
            if self.fDoRatioLegend: self.fLegendRatio.AddEntry(hRatio, h.GetTitle(), "pe")
        self.fRatios.append(hRatio)
        hRatio.SetTitle("{0} Ratio".format(h.GetTitle()))
        hRatio.Divide(self.fBaselineForRatio)
        self.fCanvasRatio.cd()
        hRatio.Draw(self.fOptRatio)
        if not self.fMainRatioHistogram:
            self.fMainRatioHistogram = hRatio
        m = DMesonJetUtils.FindMinimum(hRatio, 0., not "hist" in self.fOptRatio)
        if not m is None:
            if self.fMinRatio is None:
                self.fMinRatio = m
            else:
                self.fMinRatio = min(self.fMinRatio, m)
        m = DMesonJetUtils.FindMaximum(hRatio, 0., not "hist" in self.fOptRatio)
        if not m is None:
            if self.fMaxRatio is None:
                self.fMaxRatio = m
            else:
                self.fMaxRatio = max(self.fMaxRatio, m)

        if not "same" in self.fOptRatio:
            self.fOptRatio += "same"

    def CompareSpectra(self, baseline, histos):
        while len(histos) + 1 > len(self.fColors): self.fColors += random.sample(self.fColors[1:], len(self.fColors) - 1)
        while len(histos) + 1 > len(self.fMarkers): self.fMarkers += random.sample(self.fMarkers[1:], len(self.fMarkers) - 1)
        while len(histos) + 1 > len(self.fLines): self.fLines += random.sample(self.fLines[1:], len(self.fLines) - 1)
        while len(histos) + 1 > len(self.fFills): self.fFills += random.sample(self.fFills[1:], len(self.fFills) - 1)
        self.fResults = []
        print("CompareSpectra: {0}".format(self.fName))
        self.fBaselineHistogram = baseline
        self.fHistograms = histos
        print("Baseline: {0}".format(self.fBaselineHistogram.GetName()))
        for s in self.fHistograms:
            print(s.GetName())

        if self.fDoSpectraPlot:
            self.PrepareSpectraCanvas()

        if self.fDoRatioPlot:
            self.PrepareRatioCanvas()

        for color, marker, line, h in zip(self.fColors[1:], self.fMarkers[1:], self.fLines[1:], self.fHistograms):
            if self.fDoSpectraPlot:
                self.PlotHistogram(color, marker, line, h)
            if self.fDoRatioPlot:
                self.PlotRatio(color, marker, line, h)
        self.AdjustYLimits()
        self.GenerateResults()
        return self.fResults

    def CompareUncertainties(self, baseline, histos):
        baseline_unc = DMesonJetUtils.GetRelativeUncertaintyHistogram(baseline)
        histos_unc = [DMesonJetUtils.GetRelativeUncertaintyHistogram(h) for h in histos]
        self.CompareSpectra(baseline_unc, histos_unc)
        self.fResults.append(baseline_unc)
        self.fResults.extend(histos_unc)
        return self.fResults

    def GenerateResults(self):
        self.fResults.extend(self.fRatios)
        self.fResults.extend(self.fFitResults)
        if self.fCanvasSpectra: self.fResults.append(self.fCanvasSpectra)
        if self.fCanvasRatio: self.fResults.append(self.fCanvasRatio)
        if self.fLegendSpectra: self.fResults.append(self.fLegendSpectra)
        if self.fLegendRatio: self.fResults.append(self.fLegendRatio)
        if self.fRatioRelativeUncertainty: self.fResults.append(self.fRatioRelativeUncertainty)

    def AdjustYLimits(self):
        if not self.fMaxRatio is None and not self.fMinRatio is None:
            if self.fDoRatioPlot == "logy":
                max = self.fMaxRatio * self.fLogUpperSpace
                min = self.fMinRatio / self.fLogLowerSpace
            else:
                max = self.fMaxRatio + (self.fMaxRatio - self.fMinRatio) * self.fLinUpperSpace
                min = self.fMinRatio - (self.fMaxRatio - self.fMinRatio) * self.fLinLowerSpace
            if min < 0 and self.fMinRatio >= 0: min = 0
            self.fMainRatioHistogram.SetMinimum(min)
            self.fMainRatioHistogram.SetMaximum(max)
            self.fCanvasRatio.cd()
            if self.fDoRatioLegend: self.fLegendRatio.Draw()

        if not self.fMaxSpectrum is None and not self.fMinSpectrum is None:
            if self.fDoSpectraPlot == "logy":
                max = self.fMaxSpectrum * self.fLogUpperSpace
                min = self.fMinSpectrum / self.fLogLowerSpace
            else:
                max = self.fMaxSpectrum + (self.fMaxSpectrum - self.fMinSpectrum) * self.fLinUpperSpace
                min = self.fMinSpectrum - (self.fMaxSpectrum - self.fMinSpectrum) * self.fLinLowerSpace
            if min < 0 and self.fMinSpectrum >= 0: min = 0
            self.fMainHistogram.SetMinimum(min)
            self.fMainHistogram.SetMaximum(max)
            self.fCanvasSpectra.cd()
            if self.fDoSpectrumLegend: self.fLegendSpectra.Draw()

    def CheckConsistency(self, h1, h2):
        if h1.GetNbinsX() != h2.GetNbinsX(): return False
        for ibin in range(0, h2.GetNbinsX() + 2):
            if h1.GetXaxis().GetBinLowEdge(ibin) != h2.GetXaxis().GetBinLowEdge(ibin): return False
        return True
