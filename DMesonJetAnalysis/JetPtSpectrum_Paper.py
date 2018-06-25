#!/usr/bin/env python
# python script to make the jet pt spectrum and compare with theory

import IPython
import ROOT
import DMesonJetUtils

globalList = []

input_path = "/Volumes/DATA/ALICE/JetResults"

def GetMeasuredCrossSection():
    fname = "{0}/JetPtSpectrum_DPt_30_Systematics.root".format(input_path)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    hStat = DMesonJetUtils.GetObject(file, "FinalSpectrum/CentralPointsStatisticalUncertainty")
    if not hStat:
        print("Cannot get measured cross section with statistical uncertainty!")
        exit(1)
    hSyst = DMesonJetUtils.GetObject(file, "FinalSpectrum/CentralPointsSystematicUncertainty")
    if not hSyst:
        print("Cannot get measured cross section with systematic uncertainty!")
        exit(1)
    return hStat, hSyst

def GetTheoryCrossSection():
    # fname = "{0}/PromptDJetsPrediction_1505317519.root".format(input_path)
    fname = "{0}/PromptDJetsPrediction_1483386026.root".format(input_path)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    hStat = DMesonJetUtils.GetObject(file, "default/JetPtSpectrum_DPt_30_CrossSection/GeneratorLevel_JetPtSpectrum")
    if not hStat:
        print("Cannot get theory cross section with statistical uncertainty!")
        exit(1)
    hSystUp = DMesonJetUtils.GetObject(file, "SystematicUncertainty/JetPtSpectrum_DPt_30_CrossSection//GeneratorLevel_JetPtSpectrum/GeneratorLevel_JetPtSpectrum_UpperSyst")
    if not hSystUp:
        print("Cannot get theory cross section upper systematic uncertainty!")
        exit(1)
    hSystLow = DMesonJetUtils.GetObject(file, "SystematicUncertainty/JetPtSpectrum_DPt_30_CrossSection//GeneratorLevel_JetPtSpectrum/GeneratorLevel_JetPtSpectrum_LowerSyst")
    if not hSystUp:
        print("Cannot get theory cross section lower systematic uncertainty!")
        exit(1)
    hSyst = DMesonJetUtils.GetObject(file, "SystematicUncertainty/JetPtSpectrum_DPt_30_CrossSection//GeneratorLevel_JetPtSpectrum/GeneratorLevel_JetPtSpectrum_CentralAsymmSyst")
    if not hSyst:
        print("Cannot get theory cross section lower systematic uncertainty!")
        exit(1)

    return hStat, hSystUp, hSystLow, hSyst

def PlotCrossSections(dataStat, dataSyst, theoryStat, theorySyst):
    cname = "D0JetCrossSection_Paper"
    canvas = ROOT.TCanvas(cname, cname, 700, 750)
    globalList.append(canvas)
    canvas.Divide(1, 2)
    padMain = canvas.cd(1)
    padMain.SetPad(0, 0.35, 1, 1)
    padMain.SetTopMargin(0.05)
    padMain.SetBottomMargin(0)
    padMain.SetLeftMargin(0.17)
    padMain.SetRightMargin(0.05)
    padMain.SetTicks(1, 1)
    padRatio = canvas.cd(2)
    padRatio.SetPad(0, 0., 1, 0.35)
    padRatio.SetTopMargin(0)
    padRatio.SetBottomMargin(0.27)
    padRatio.SetLeftMargin(0.17)
    padRatio.SetRightMargin(0.05)
    padRatio.SetGridy()
    padRatio.SetTicks(1, 1)

    padMain.cd()
    padMain.SetLogy()
    h = ROOT.TH1I("axis", "axis", 1000, 0, 50)
    h.Draw("axis")
    globalList.append(h)
    h.GetXaxis().SetRangeUser(5, 30)
    h.GetYaxis().SetRangeUser(2e-5, 2e-2)
    h.GetYaxis().SetTitleFont(43)
    h.GetYaxis().SetTitleSize(26)
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelSize(22)
    h.GetYaxis().SetTitleOffset(1.9)
    h.GetYaxis().SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T,jet}^{ch}d#it{#eta}_{jet}} [mb (GeV/#it{c})^{-1}]")

    dataSyst_copy = dataSyst.Clone("{0}_copy".format(dataSyst.GetName()))
    dataSyst_copy.Draw("2")
    globalList.append(dataSyst_copy)

    dataStat_copy = dataStat.DrawCopy("same p e0 x0")
    dataStat_copy.SetMarkerSize(1.4)
    globalList.append(dataStat_copy)

    theoryStat_copy = theoryStat.DrawCopy("same p e0 x0")
    globalList.append(theoryStat_copy)
    theoryStat_copy.SetLineColor(ROOT.kBlue + 2)
    theoryStat_copy.SetMarkerColor(ROOT.kBlue + 2)
    theoryStat_copy.SetMarkerStyle(ROOT.kOpenCircle)
    theoryStat_copy.SetMarkerSize(1.2)

    theorySyst_copy = theorySyst.Clone("theorySyst_copy")
    globalList.append(theorySyst_copy)
    theorySyst_copy.Draw("2")
    theorySyst_copy.SetFillStyle(0)
    theorySyst_copy.SetLineWidth(2)
    theorySyst_copy.SetLineColor(ROOT.kBlue + 2)

    # Calculate ratio

    ratioDataSyst = dataSyst_copy.Clone("ratioDataSyst")
    globalList.append(ratioDataSyst)
    
    ratioDataStat = dataStat_copy.Clone("ratioDataStat")
    globalList.append(ratioDataStat)

    ratioTheorySyst = theorySyst_copy.Clone("ratioTheorySyst")
    globalList.append(ratioTheorySyst)
    
    ratioTheoryStat = theoryStat_copy.Clone("ratioTheorySyst")
    globalList.append(ratioTheoryStat)

    for ibin in range(1, ratioDataStat.GetNbinsX() + 1):
        ratioTheorySyst.SetPointEYlow(ibin - 1, ratioTheorySyst.GetErrorYlow(ibin - 1) / ratioDataStat.GetBinContent(ibin))
        ratioTheorySyst.SetPointEYhigh(ibin - 1, ratioTheorySyst.GetErrorYhigh(ibin - 1) / ratioDataStat.GetBinContent(ibin))
        ratioTheorySyst.SetPoint(ibin - 1, ratioTheorySyst.GetX()[ibin - 1], ratioTheorySyst.GetY()[ibin - 1] / ratioDataStat.GetBinContent(ibin))

        ratioDataSyst.SetPointEYlow(ibin - 1, ratioDataSyst.GetErrorYlow(ibin - 1) / ratioDataStat.GetBinContent(ibin))
        ratioDataSyst.SetPointEYhigh(ibin - 1, ratioDataSyst.GetErrorYhigh(ibin - 1) / ratioDataStat.GetBinContent(ibin))
        ratioDataSyst.SetPoint(ibin - 1, ratioDataSyst.GetX()[ibin - 1], 1.0)

        ratioTheoryStat.SetBinError(ibin, ratioTheoryStat.GetBinError(ibin) / ratioDataStat.GetBinContent(ibin))
        ratioTheoryStat.SetBinContent(ibin, ratioTheoryStat.GetBinContent(ibin) / ratioDataStat.GetBinContent(ibin))

        ratioDataStat.SetBinError(ibin, ratioDataStat.GetBinError(ibin) / ratioDataStat.GetBinContent(ibin))
        ratioDataStat.SetBinContent(ibin, 1.0)

    # Plotting ratio

    padRatio.cd()

    hRatio = ROOT.TH1I("axis", "axis", 1000, 0, 50)
    hRatio.GetXaxis().SetRangeUser(5, 30)
    hRatio.Draw("axis")
    globalList.append(hRatio)
    hRatio.GetYaxis().SetTitle("MC / Data")
    hRatio.GetXaxis().SetTitleFont(43)
    hRatio.GetXaxis().SetTitleSize(26)
    hRatio.GetXaxis().SetLabelFont(43)
    hRatio.GetXaxis().SetLabelSize(22)
    hRatio.GetYaxis().SetTitleFont(43)
    hRatio.GetYaxis().SetTitleSize(26)
    hRatio.GetYaxis().SetLabelFont(43)
    hRatio.GetYaxis().SetLabelSize(22)
    hRatio.GetYaxis().SetTitleOffset(1.4)
    hRatio.GetXaxis().SetTitleOffset(2.9)
    hRatio.GetYaxis().SetRangeUser(0, 2.99)
    hRatio.GetYaxis().SetNdivisions(509)
    hRatio.GetXaxis().SetTitle("#it{p}_{T,jet}^{ch} (GeV/#it{c})")

    ratioDataSyst.Draw("2")
    ratioDataStat.Draw("same p e0 x0")

    ratioTheorySyst.Draw("2")
    ratioTheoryStat.Draw("same p e0 x0")

    # Plotting labels

    padMain.cd()
    leg1 = ROOT.TLegend(0.19, 0.26, 0.54, 0.06, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(20)
    leg1.SetTextAlign(12)
    leg1.SetMargin(0.2)
    entry = leg1.AddEntry(None, "Data", "pf")
    entry.SetFillColor(dataSyst_copy.GetFillColor())
    entry.SetFillStyle(dataSyst_copy.GetFillStyle())
    entry.SetLineWidth(0)
    entry.SetMarkerColor(dataStat_copy.GetMarkerColor())
    entry.SetMarkerStyle(dataStat_copy.GetMarkerStyle())
    entry = leg1.AddEntry(None, "POWHEG (ccbar) + PYTHIA 6", "pf")
    entry.SetFillColor(theorySyst_copy.GetFillColor())
    entry.SetLineColor(theorySyst_copy.GetLineColor())
    entry.SetLineWidth(theorySyst_copy.GetLineWidth())
    entry.SetMarkerColor(theoryStat_copy.GetMarkerColor())
    entry.SetMarkerStyle(theoryStat_copy.GetMarkerStyle())
    entry = leg1.AddEntry(None, "POWHEG (dijet) + PYTHIA 6", "pf")
    entry.SetFillColor(ROOT.kGreen + 2)
    entry.SetLineColor(ROOT.kGreen + 2)
    entry.SetLineWidth(theorySyst_copy.GetLineWidth())
    entry.SetMarkerColor(ROOT.kGreen + 2)
    entry.SetMarkerStyle(ROOT.kOpenSquare)
    leg1.Draw()

    padMain.cd()
    paveALICE = ROOT.TPaveText(0.36, 0.72, 0.71, 0.92, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(22)
    paveALICE.SetTextAlign(13)
    paveALICE.AddText("ALICE, pp, #sqrt{#it{s}} = 7 TeV")
    paveALICE.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.4, |#it{#eta}_{jet}| < 0.5")
    paveALICE.AddText("with D^{0}, #it{p}_{T,D} > 3 GeV/#it{c}")
    paveALICE.Draw()

    padRatio.RedrawAxis("g")
    padRatio.RedrawAxis()
    padMain.RedrawAxis()

    return canvas

def main():
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    dataStat, dataSyst = GetMeasuredCrossSection()
    theoryStat, _, _, theorySyst = GetTheoryCrossSection()
    canvas = PlotCrossSections(dataStat, dataSyst, theoryStat, theorySyst)
    canvas.SaveAs("{0}/D0JetCrossSection_Paper.pdf".format(input_path))
    canvas.SaveAs("{0}/D0JetCrossSection_Paper.C".format(input_path))

if __name__ == '__main__':
    main()

    IPython.embed()
