#!/usr/bin/env python

import subprocess
import IPython
import ROOT
import RawYieldSpectrumLoader

globalList = []

input_path = "/Volumes/DATA/ALICE/JetResults"

def PlotSBSpectra(pad, ptmin, ptmax, sbList, plotleg=False):
    pad.SetTicks(1, 1)
    # pad.SetLogy()
    pad.SetLeftMargin(0.22)
    pad.SetRightMargin(0.04)
    pad.SetTopMargin(0.04)
    pad.SetBottomMargin(0.19)

    sbHistL = sbList.FindObject("D0_Charged_R040_JetPtSpectrum_DPt_30_SideBand_SideBandWindowR_DPt_{0:.0f}_{1:.0f}".format(ptmin * 100, ptmax * 100))
    h = sbHistL.DrawCopy("axis")
    globalList.append(h)

    sbHistL_copy = sbHistL.DrawCopy("p0 same")
    sbHistL_copy.Scale(1. / sbHistL_copy.Integral())
    globalList.append(sbHistL_copy)
    sbHistL_copy.SetMarkerColor(ROOT.kBlue + 2)
    sbHistL_copy.SetLineColor(ROOT.kBlue + 2)
    sbHistL_copy.SetMarkerStyle(ROOT.kOpenCircle)
    sbHistL_copy.SetMarkerSize(0.9)

    sbHistR = sbList.FindObject("D0_Charged_R040_JetPtSpectrum_DPt_30_SideBand_SideBandWindowL_DPt_{0:.0f}_{1:.0f}".format(ptmin * 100, ptmax * 100))

    sbHistR_copy = sbHistR.DrawCopy("p0 same")
    sbHistR_copy.Scale(1. / sbHistR_copy.Integral())
    globalList.append(sbHistR_copy)
    sbHistR_copy.SetMarkerColor(ROOT.kRed + 2)
    sbHistR_copy.SetLineColor(ROOT.kRed + 2)
    sbHistR_copy.SetMarkerStyle(ROOT.kOpenSquare)
    sbHistR_copy.SetMarkerSize(0.9)

    h.GetYaxis().SetTitle("arb. units")
    h.GetXaxis().SetTitleFont(43)
    h.GetXaxis().SetTitleOffset(3.3)
    h.GetXaxis().SetTitleSize(19)
    h.GetXaxis().SetLabelFont(43)
    h.GetXaxis().SetLabelOffset(0.009)
    h.GetXaxis().SetLabelSize(18)
    h.GetYaxis().SetTitleFont(43)
    h.GetYaxis().SetTitleOffset(4.3)
    h.GetYaxis().SetTitleSize(19)
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelOffset(0.009)
    h.GetYaxis().SetLabelSize(23)

    # miny = min([DMesonJetUtils.FindMinimum(sbHist_copy), DMesonJetUtils.FindMinimum(sigHist_copy), DMesonJetUtils.FindMinimum(subHist_copy)])
    # maxy = max([DMesonJetUtils.FindMaximum(sbHist_copy), DMesonJetUtils.FindMaximum(sigHist_copy), DMesonJetUtils.FindMaximum(subHist_copy)])
    # miny /= 2
    # maxy *= 3
    diff = sbHistL_copy.GetMaximum() - sbHistL_copy.GetMinimum()
    miny = sbHistL_copy.GetMinimum() - 0.10 * diff
    maxy = sbHistL_copy.GetMaximum() + 0.3 * diff
    h.SetMaximum(maxy)
    h.SetMinimum(miny)

    if plotleg:
        leg = ROOT.TLegend(0.49, 0.70, 0.87, 0.90, "", "NB NDC")
        globalList.append(leg)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(43)
        leg.SetTextSize(19)
        leg.SetTextAlign(13)
        leg.AddEntry(sbHistL_copy, "Left SB", "p")
        leg.AddEntry(sbHistR_copy, "Right SB", "p")
        leg.Draw()

def SideBandPlot():
    loader = RawYieldSpectrumLoader.RawYieldSpectrumLoader(input_path, "Jets_EMC_pp_823_824_825_826", "LHC10_Train823_noRefl")
    loader.fDMeson = "D0"
    loader.fJetType = "Charged"
    loader.fJetRadius = "R040"
    loader.fVariableName = "JetPt"
    loader.fKinematicCuts = "DPt_30"
    loader.fRawYieldMethod = "SideBand"
    loader.LoadDataListFromDMesonJetAnalysis()
    dptbinList = loader.fDataJetList.FindObject("D0_Charged_R040_DPtBins_JetPt_5_30")
    spectrumList = loader.fDataSpectrumList
    sbList = spectrumList.FindObject("SideBandAnalysis")

    cname = "SideBandLeftVsRight"
    canvas = ROOT.TCanvas(cname, cname, 900, 900)
    globalList.append(canvas)
    canvas.Divide(3, 3)
    bins = [3, 4, 5, 6, 7, 8, 10, 12, 16, 30]
    leg = True
    for i, (minPt, maxPt) in enumerate(zip(bins[:-1], bins[1:])):
        PlotSBSpectra(canvas.cd(i + 1), minPt, maxPt, sbList, leg)
        leg = False

def main():
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    subprocess.call("make")
    ROOT.gSystem.Load("MassFitter.so")

    SideBandPlot()

    for obj in globalList:
        if isinstance(obj, ROOT.TCanvas):
            obj.SaveAs("{0}/{1}.pdf".format(input_path, obj.GetName()))

if __name__ == '__main__':
    main()

    IPython.embed()
