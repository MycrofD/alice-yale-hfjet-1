#!/usr/bin/env python

import IPython
import numpy
import ROOT
import DMesonJetUtils
import RawYieldSpectrumLoader
import HistogramNormalizator
import DMesonJetCompare

globalList = []

input_path = "/Volumes/DATA/ALICE/JetResults"

pass4DmesonAna = "HFPtSpectrum_Pass4_combinedFDForLow_mergeThr1.root"

antiPartNorm = 2.0  # D0 / D0bar
events = 0

def LoadEvents(file, hname):
    global events
    slash = hname.find("/")
    heventsName = hname[:slash + 1]
    heventsName += "Events"
    hevents = DMesonJetUtils.GetObject(file, heventsName)
    events = hevents.GetBinContent(1)
    print("The number of events is {0}".format(events))

def GetJetPtSpectrum():
    global events
    unfolding = "LHC10_Train823_LHC15i2_Train961_efficiency_mt_refl_DoubleGaus_15"
    fname = "{0}/{1}/{1}.root".format(input_path, unfolding)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    else:
        print("File {0} open successfully".format(fname))

    hname = "SideBand_DPt_30/Bayes/SideBand_DPt_30_UnfoldedSpectrum_Bayes_Reg4_PriorResponseTruth"

    LoadEvents(file, hname)
    h = DMesonJetUtils.GetObject(file, hname)
    print("The number of events is {0}".format(events))
    h.Scale(HistogramNormalizator.CROSS_SECTION / (events * HistogramNormalizator.BRANCHING_RATIO * antiPartNorm), "width")
    return h

def GetDPtSpectrum(kincuts=None, jet_radius=None, jet_type=None):
    global events
    loader = RawYieldSpectrumLoader.RawYieldSpectrumLoader(input_path, "Jets_EMC_pp_823_824_825_826", "LHC10_Train823_efficiency")
    loader.fUseReflections = False
    loader.fDMeson = "D0"
    loader.fVariableName = "DPt"
    loader.fJetType = jet_type
    loader.fJetRadius = jet_radius
    loader.fKinematicCuts = kincuts
    loader.fRawYieldMethod = "InvMassFit"
    h = loader.GetDefaultSpectrumFromDMesonJetAnalysis(True, 0, 0)
    h.Scale(HistogramNormalizator.CROSS_SECTION / (events * HistogramNormalizator.BRANCHING_RATIO * antiPartNorm), "width")
    return h

def GetTheoryJetPtCrossSection():
    fname = "{0}/PromptDJetsPrediction.root".format(input_path)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    hStat = DMesonJetUtils.GetObject(file, "default/JetPtSpectrum_DPt_30/GeneratorLevel_JetPtSpectrum")
    if not hStat:
        print("Cannot get theory cross section with statistical uncertainty!")
        exit(1)
    hSystUp = DMesonJetUtils.GetObject(file, "SystematicUncertainty/JetPtSpectrum_DPt_30//GeneratorLevel_JetPtSpectrum/GeneratorLevel_JetPtSpectrum_UpperSyst")
    if not hSystUp:
        print("Cannot get theory cross section upper systematic uncertainty!")
        exit(1)
    hSystLow = DMesonJetUtils.GetObject(file, "SystematicUncertainty/JetPtSpectrum_DPt_30//GeneratorLevel_JetPtSpectrum/GeneratorLevel_JetPtSpectrum_LowerSyst")
    if not hSystUp:
        print("Cannot get theory cross section lower systematic uncertainty!")
        exit(1)

    # scale for the bin width and the antiparticle factor
    hStat.Scale(0.5, "width")
    hSystUp.Scale(0.5, "width")
    hSystLow.Scale(0.5, "width")
    return hStat, hSystUp, hSystLow

def GetTheoryDPtCrossSection():
    fname = "{0}/PromptDJetsPrediction.root".format(input_path)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    hStat = DMesonJetUtils.GetObject(file, "default/DPtSpectrum/GeneratorLevel_DPtSpectrum")
    if not hStat:
        print("Cannot get theory cross section with statistical uncertainty!")
        exit(1)
    hSystUp = DMesonJetUtils.GetObject(file, "SystematicUncertainty/DPtSpectrum//GeneratorLevel_DPtSpectrum/GeneratorLevel_DPtSpectrum_UpperSyst")
    if not hSystUp:
        print("Cannot get theory cross section upper systematic uncertainty!")
        exit(1)
    hSystLow = DMesonJetUtils.GetObject(file, "SystematicUncertainty/DPtSpectrum//GeneratorLevel_DPtSpectrum/GeneratorLevel_DPtSpectrum_LowerSyst")
    if not hSystUp:
        print("Cannot get theory cross section lower systematic uncertainty!")
        exit(1)

    # scale for the bin width and the antiparticle factor
    hStat.Scale(0.5, "width")
    hSystUp.Scale(0.5, "width")
    hSystLow.Scale(0.5, "width")
    return hStat, hSystUp, hSystLow

def GetPass4AnalysisSpectrum():
    file = ROOT.TFile(pass4DmesonAna)
    h = file.Get("histoSigmaCorr")
    if not h:
        print("Cannot get pass 4 analysis spectrum!")
        exit(1)
    h_copy = h.Clone("pass4")
    file.Close()
    return h_copy


def CompareJetPtvsDPt():
    jetPtSpectrumHist = GetJetPtSpectrum()
    dptSpectrumHist = GetDPtSpectrum()
    jetPtSpectrumHist.SetTitle("Data, D^{0}-jet")
    dptSpectrumHist.SetTitle("Data, D^{0}")

    jetPtSpectrumHist_theory, jetPtSpectrumHist_theory_systup, jetPtSpectrumHist_theory_systdown = GetTheoryJetPtCrossSection()
    dptSpectrumHist_theory, dptSpectrumHist_theory_systup, dptSpectrumHist_theory_systdown = GetTheoryDPtCrossSection()
    jetPtSpectrumHist_theory.SetTitle("POWHEG, D^{0}-jet")
    dptSpectrumHist_theory.SetTitle("POWHEG, D^{0}")

    hpass4_old = GetPass4AnalysisSpectrum()
    hpass4 = hpass4_old.Clone("hpass4")
    hpass4.Scale(1e-9 / HistogramNormalizator.BRANCHING_RATIO)
    hpass4.SetTitle("Data, D^{0}")
    globalList.append(hpass4)

    canvasSpectrum = ROOT.TCanvas("Comparison_DPt_JetPt_Spectra", "Comparison_DPt_JetPt_Spectra")
    globalList.append(canvasSpectrum)
    canvasSpectrum.SetTicks(1, 1)
    canvasSpectrum.cd()
    h_axis = dptSpectrumHist_theory.DrawCopy("axis")
    h_axis.GetYaxis().SetTitle("#frac{d#sigma}{d#it{p}_{T}} [mb (GeV/#it{c})^{-1}]")
    h_axis.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")

    comp = DMesonJetCompare.DMesonJetCompare("Comparison_DPt_JetPt_Spectra")
    comp.fCanvasSpectra = canvasSpectrum
    comp.fMainHistogram = h_axis
    comp.fOptSpectrumBaseline = "same"
    comp.fOptSpectrum = "same"
    comp.fX1LegSpectrum = 0.58
    comp.fY1LegSpectrum = 0.52

    comp.fMarkers = [ROOT.kFullCircle, ROOT.kFullSquare]
    comp.fColors = [ROOT.kRed + 2, ROOT.kBlue + 2]
    r = comp.CompareSpectra(jetPtSpectrumHist, [hpass4])  # , [dptSpectrumHist])
    for obj in r:
        if not obj in globalList:
            globalList.append(obj)

    comp.fMarkers = [ROOT.kOpenCircle, ROOT.kOpenSquare]
    comp.fColors = [ROOT.kOrange + 2, ROOT.kGreen + 2]
    r = comp.CompareSpectra(jetPtSpectrumHist_theory, [dptSpectrumHist_theory])
    for obj in r:
        if not obj in globalList:
            globalList.append(obj)

    canvasSpectrum.cd()

#     h = jetPtSpectrumHist_theory.DrawCopy("same hist l")
#     h.Add(jetPtSpectrumHist_theory_systup, 1)
#     h.SetLineColor(ROOT.kOrange + 2)
#     globalList.append(h)
#     h = jetPtSpectrumHist_theory.DrawCopy("same hist l")
#     h.Add(jetPtSpectrumHist_theory_systdown, -1)
#     h.SetLineColor(ROOT.kOrange + 2)
#     globalList.append(h)
#
#     h = dptSpectrumHist_theory.DrawCopy("same hist l")
#     h.Add(dptSpectrumHist_theory_systup, 1)
#     h.SetLineColor(ROOT.kGreen + 2)
#     globalList.append(h)
#     h = dptSpectrumHist_theory.DrawCopy("same hist l")
#     h.Add(dptSpectrumHist_theory_systdown, -1)
#     h.SetLineColor(ROOT.kGreen + 2)
#     globalList.append(h)

    pave = ROOT.TPaveText(0.19, 0.80, 0.60, 0.90, "NB NDC")
    globalList.append(pave)
    pave.SetBorderSize(0)
    pave.SetFillStyle(0)
    pave.SetTextFont(43)
    pave.SetTextSize(20)
    pave.SetTextAlign(13)
    pave.AddText("pp, #sqrt{#it{s}} = 7 TeV")
    pave.AddText("Fully corrected spectra, statistical uncertainty only")
    pave.Draw()

    pave = ROOT.TPaveText(0.35, 0.55, 0.81, 0.75, "NB NDC")
    globalList.append(pave)
    pave.SetBorderSize(0)
    pave.SetFillStyle(0)
    pave.SetTextFont(43)
    pave.SetTextSize(20)
    pave.SetTextAlign(13)
    pave.AddText("Charged D^{0}-jets: Anti-#it{k}_{T}, #it{R}=0.4")
    pave.AddText("|#eta_{jet}| < 0.5 with D^{0}, #it{p}_{T,D} > 3 GeV/#it{c}")
    pave.AddText("D^{0}: |#eta_{D}| < 0.5")
    pave.Draw()

    dptBins = [3, 4, 5, 6, 7, 8, 12, 16]
    hpass4_copy = DMesonJetUtils.Rebin1D_fromBins(hpass4, "hpass4_copy", len(dptBins) - 1, numpy.array(dptBins, dtype=float))
    hpass4_copy.SetBinContent(6, hpass4_copy.GetBinContent(6) / 2)
    globalList.append(hpass4_copy)

    dptSpectrumHist_copy = DMesonJetUtils.Rebin1D_fromBins(dptSpectrumHist, "dptSpectrumHist_copy", len(dptBins) - 1, numpy.array(dptBins, dtype=float))
    dptSpectrumHist_copy.SetBinContent(6, dptSpectrumHist_copy.GetBinContent(6) / 2)
    globalList.append(dptSpectrumHist_copy)

    comp = DMesonJetCompare.DMesonJetCompare("Comparison_DPt_Pass4_Spectra")
    comp.fLinUpperSpace = 0.2
    hpass4_copy.SetTitle("D^{0} from pass4 analysis")
    dptSpectrumHist_copy.SetTitle("D^{0} from D^{0}-jet analysis")
    hpub = GetPublishedDmeson()
    globalList.append(hpub)
    hpub.Scale(1e-3)

    r = comp.CompareSpectra(hpass4_copy, [dptSpectrumHist_copy])
    # ratio = comp.fRatios[0]
    # ratio.Fit("pol0")
    for obj in r:
        if not obj in globalList:
            globalList.append(obj)

def GetPublishedDmeson():
    dptBins = [3, 4, 5, 6, 7, 8, 12, 16]
    yval = [59.7, 29.1, 12.5, 6.37, 3.07, 1.23, 0.215]
    yerr = [13.313526955694348, 6.168468205316454, 2.5495097567963922, 1.2870120434556935, 0.7083784299369935, 0.2469817807045694, 0.06280127387243033]
    h = ROOT.TH1D("published", "Published D^{0}", len(dptBins) - 1, numpy.array(dptBins, dtype=float))
    for ibin in range(1, h.GetNbinsX() + 1):
        h.SetBinContent(ibin, yval[ibin - 1])
        h.SetBinError(ibin, yerr[ibin - 1])
    return h

def main():
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    CompareJetPtvsDPt()

    for obj in globalList:
        if isinstance(obj, ROOT.TCanvas):
            obj.SaveAs("{0}/{1}.pdf".format(input_path, obj.GetName()))

if __name__ == '__main__':
    main()

    IPython.embed()
