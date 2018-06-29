#!/usr/bin/env python

import math
from scipy import stats
import IPython
import ROOT
import DMesonJetUtils

globalList = []

input_path = "/Volumes/DATA/ALICE/JetResults"

JET_PT_ACC = 10.0

def GetMeasuredCrossSection():
    fname = "{0}/JetZCrossSection_JetPt_5_15_Systematics.root".format(input_path)
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
    hStat.Scale(1.0 / JET_PT_ACC)
    for ipoint in range(0, hSyst.GetN()):
        hSyst.SetPointEYlow(ipoint, hSyst.GetErrorYlow(ipoint) / JET_PT_ACC)
        hSyst.SetPointEYhigh(ipoint, hSyst.GetErrorYhigh(ipoint) / JET_PT_ACC)
        hSyst.SetPoint(ipoint, hSyst.GetX()[ipoint], hSyst.GetY()[ipoint] / JET_PT_ACC)

    return hStat, hSyst

def GetTheoryCrossSection():
    fname = "{0}/PromptDJetsPrediction_1483386026.root".format(input_path)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    hStat = DMesonJetUtils.GetObject(file, "default/JetZSpectrum_DPt_20_JetPt_5_15_CrossSection/GeneratorLevel_JetZSpectrum")
    if not hStat:
        print("Cannot get theory cross section with statistical uncertainty!")
        exit(1)
    hSyst = DMesonJetUtils.GetObject(file, "SystematicUncertainty/JetZSpectrum_DPt_20_JetPt_5_15_CrossSection//GeneratorLevel_JetZSpectrum/GeneratorLevel_JetZSpectrum_CentralAsymmSyst")
    if not hSyst:
        print("Cannot get theory cross section lower systematic uncertainty!")
        exit(1)

    hStat.Scale(1.0 / JET_PT_ACC)
    for ipoint in range(0, hSyst.GetN()):
        hSyst.SetPointEYlow(ipoint, hSyst.GetErrorYlow(ipoint) / JET_PT_ACC)
        hSyst.SetPointEYhigh(ipoint, hSyst.GetErrorYhigh(ipoint) / JET_PT_ACC)
        hSyst.SetPoint(ipoint, hSyst.GetX()[ipoint], hSyst.GetY()[ipoint] / JET_PT_ACC)

    return hStat, hSyst

def CalculateMean(stat, syst):
    weighted_sum = 0
    sum_of_weights = 0
    mean = 0
    sum_of_deviations_squared_error_weighted = 0
    mean_error = 0

    for ibin in range(1, stat.GetNbinsX() + 1):
        weighted_sum += stat.GetBinContent(ibin) * stat.GetXaxis().GetBinCenter(ibin)
        sum_of_weights += stat.GetBinContent(ibin)

    mean = weighted_sum / sum_of_weights

    for ibin in range(1, stat.GetNbinsX() + 1):
        err2 = stat.GetBinError(ibin) ** 2 + syst.GetErrorY(ibin - 1) ** 2
        dev2 = (stat.GetXaxis().GetBinCenter(ibin) - mean) ** 2
        sum_of_deviations_squared_error_weighted += err2 * dev2

    mean_error = math.sqrt(sum_of_deviations_squared_error_weighted) / sum_of_weights

    return mean, mean_error

def CalculateChi2(stat1, syst1, stat2, syst2):
    chi2 = 0
    for ibin in range(1, stat1.GetNbinsX() + 1):
        tot_err2 = stat1.GetBinError(ibin) ** 2 + syst1.GetErrorY(ibin - 1) ** 2 + stat2.GetBinError(ibin) ** 2 + syst2.GetErrorY(ibin - 1) ** 2
        diff2 = (stat1.GetBinContent(ibin) - stat2.GetBinContent(ibin)) ** 2
        chi2 += diff2 / tot_err2
    p = 1 - stats.chi2.cdf(chi2, stat1.GetNbinsX())

    return chi2, stat1.GetNbinsX(), p

def PlotCrossSections(dataStat, dataSyst, theoryStat, theorySyst):
    mean_z_data, mean_z_err_data = CalculateMean(dataStat, dataSyst)
    mean_z_theory, mean_z_err_theory = CalculateMean(theoryStat, theorySyst)
    print("The mean of the data distribution is {} +/- {}".format(mean_z_data, mean_z_err_data))
    print("The mean of the theory distribution is {} +/- {}".format(mean_z_theory, mean_z_err_theory))

    chi2, ndf, p = CalculateChi2(dataStat, dataSyst, theoryStat, theorySyst)
    print("The chi2 is {}, the ndf = {} and the p-value is {}".format(chi2, ndf, p))

    cname = "JetZCrossSection_JetPt_5_15_Paper"
    canvas = ROOT.TCanvas(cname, cname, 700, 700)
    globalList.append(canvas)
    canvas.Divide(1, 2)
    padMain = canvas.cd(1)
    padMain.SetPad(0, 0.35, 1, 1)
    padMain.SetBottomMargin(0)
    padMain.SetTopMargin(0.05)
    padMain.SetLeftMargin(0.18)
    padMain.SetRightMargin(0.05)
    padMain.SetTicks(1, 1)
    padRatio = canvas.cd(2)
    padRatio.SetPad(0, 0., 1, 0.35)
    padRatio.SetTopMargin(0)
    padRatio.SetBottomMargin(0.27)
    padRatio.SetLeftMargin(0.18)
    padRatio.SetRightMargin(0.05)
    padRatio.SetGridy()
    padRatio.SetTicks(1, 1)

    padMain.cd()
    h = ROOT.TH1I("hAxis", "hAxis", 100, 0, 1)
    globalList.append(h)
    h.Draw("axis")
    h.GetYaxis().SetTitle("#frac{d^{3}#sigma}{d#it{z}_{||}^{ch}d#it{p}_{T,jet}^{ch}d#it{#eta}} (GeV/#it{c})^{-1}")
    h.GetYaxis().SetRangeUser(-0.002, 0.032)
    h.GetYaxis().SetTitleFont(43)
    h.GetYaxis().SetTitleSize(26)
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelSize(22)
    h.GetYaxis().SetTitleOffset(1.9)

    dataSyst_copy = dataSyst.Clone("{0}_copy".format(dataSyst.GetName()))
    dataSyst_copy.Draw("2")
    globalList.append(dataSyst_copy)

    dataStat_copy = dataStat.DrawCopy("same p e0 x0")
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
    
    for ipoint in range(0, ratioTheorySyst.GetN()):
        ratioTheorySyst.SetPointEYlow(ipoint, ratioTheorySyst.GetErrorYlow(ipoint) / ratioDataSyst.GetY()[ipoint])
        ratioTheorySyst.SetPointEYhigh(ipoint, ratioTheorySyst.GetErrorYhigh(ipoint) / ratioDataSyst.GetY()[ipoint])
        ratioTheorySyst.SetPoint(ipoint, ratioTheorySyst.GetX()[ipoint], ratioTheorySyst.GetY()[ipoint] / ratioDataSyst.GetY()[ipoint])

        ratioDataSyst.SetPointEYlow(ipoint, ratioDataSyst.GetErrorYlow(ipoint) / ratioDataSyst.GetY()[ipoint])
        ratioDataSyst.SetPointEYhigh(ipoint, ratioDataSyst.GetErrorYhigh(ipoint) / ratioDataSyst.GetY()[ipoint])
        ratioDataSyst.SetPoint(ipoint, ratioDataSyst.GetX()[ipoint], 1.0)
    
    for ibin in range(1, ratioDataStat.GetNbinsX() + 1):
        ratioTheoryStat.SetBinError(ibin, ratioTheoryStat.GetBinError(ibin) / ratioDataStat.GetBinContent(ibin))
        ratioTheoryStat.SetBinContent(ibin, ratioTheoryStat.GetBinContent(ibin) / ratioDataStat.GetBinContent(ibin))

        ratioDataStat.SetBinError(ibin, ratioDataStat.GetBinError(ibin) / ratioDataStat.GetBinContent(ibin))
        ratioDataStat.SetBinContent(ibin, 1.0)
    
    # Plotting ratio

    padRatio.cd()

    hRatio = ROOT.TH1I("axis", "axis", 100, 0, 1)
    hRatio.GetXaxis().SetRangeUser(0, 1)
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
    hRatio.GetYaxis().SetRangeUser(0, 1.99)
    hRatio.GetYaxis().SetNdivisions(509)
    hRatio.GetXaxis().SetTitle("#it{z}_{||}^{ch}")

    ratioDataSyst.Draw("2")
    ratioDataStat.Draw("same p e0 x0")

    ratioTheorySyst.Draw("2")
    ratioTheoryStat.Draw("same p e0 x0")

    # Plotting labels

    padMain.cd()
    leg1 = ROOT.TLegend(0.21, 0.52, 0.50, 0.65, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(23)
    leg1.SetTextAlign(12)
    leg1.SetMargin(0.2)
    entry = leg1.AddEntry(None, "Data", "pf")
    entry.SetMarkerColor(dataStat_copy.GetMarkerColor())
    entry.SetMarkerStyle(dataStat_copy.GetMarkerStyle())
    entry.SetLineColor(dataSyst_copy.GetFillColor())
    entry.SetFillColor(dataSyst_copy.GetFillColor())
    entry.SetFillStyle(dataSyst_copy.GetFillStyle())
    entry = leg1.AddEntry(theoryStat_copy, "POWHEG+PYTHIA 6", "pf")
    entry.SetLineWidth(theorySyst_copy.GetLineWidth())
    leg1.Draw()

    padMain.cd()
    paveALICE = ROOT.TPaveText(0.20, 0.70, 0.55, 0.95, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(22)
    paveALICE.SetTextAlign(13)
    paveALICE.AddText("ALICE, pp, #sqrt{#it{s}} = 7 TeV")
    paveALICE.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.4, |#it{#eta}_{jet}| < 0.5")
    paveALICE.AddText("5 < #it{p}_{T,jet}^{ch} < 15 GeV/#it{c}, with D^{0}, #it{p}_{T,D} > 2 GeV/#it{c}")
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
    theoryStat, theorySyst = GetTheoryCrossSection()
    canvas = PlotCrossSections(dataStat, dataSyst, theoryStat, theorySyst)
    canvas.SaveAs("{}/{}.pdf".format(input_path, canvas.GetName()))
    canvas.SaveAs("{}/{}.C".format(input_path, canvas.GetName()))
    canvas.SaveAs("{}/{}.eps".format(input_path, canvas.GetName()))

if __name__ == '__main__':
    main()

    IPython.embed()