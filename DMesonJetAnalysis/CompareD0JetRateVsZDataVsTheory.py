#!/usr/bin/env python
# python script to compare D0 jet rate over inclusive jet vs z with theory
# use it with
# D0JetRateVsZ_JetPt_5_15_CompareFullVsCharged.yaml
# D0JetRateVsZ_JetPt_5_15_CompareLO.yaml
# D0JetRateVsZ_JetPt_5_15_CompareMonteCarlo_Paper.yaml
# D0JetRateVsZ_JetPt_5_15_ComparePtDCut.yaml

# D0JetRateVsZ_JetPt_15_30_CompareFullVsCharged.yaml
# D0JetRateVsZ_JetPt_15_30_CompareLO.yaml
# D0JetRateVsZ_JetPt_15_30_CompareMonteCarlo_Paper.yaml
# D0JetRateVsZ_JetPt_15_30_ComparePtDCut.yaml

import argparse
import yaml
import IPython
import ROOT
import DMesonJetUtils
import LoadInclusiveJetSpectrum
import LoadTheoryCrossSections
import HistogramNormalizator

globalList = []

def GetD0JetCrossSection(input_path, file_name):
    fname = "{}/{}.root".format(input_path, file_name)
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

def GetInclJetCrossSection(bins):
    return LoadInclusiveJetSpectrum.GetCrossSection(bins)

def GetD0JetTheoryCrossSectionAll(config, axis):
    return LoadTheoryCrossSections.GetD0JetTheoryCrossSectionAll(config, axis)

def GetInclusiveJetTheoryCrossSectionAll(config):
    return LoadTheoryCrossSections.GetInclusiveJetTheoryCrossSectionAll(config)

def PlotCrossSections(dataStat, dataSyst, config):
    cname = "{}_{}".format(config["name_prefix"], config["name"])
    if "canvas_h" in config:
        canvas_h = config["canvas_h"]
    else:
        canvas_h = 900
    if "canvas_w" in config:
        canvas_w = config["canvas_w"]
    else:
        canvas_w = 700
    canvas = ROOT.TCanvas(cname, cname, canvas_w, canvas_h)
    globalList.append(canvas)
    canvas.Divide(1, 2)
    padMain = canvas.cd(1)
    padMain.SetPad(0, 0.35, 1, 1)
    padMain.SetBottomMargin(0)
    padMain.SetLeftMargin(config["left_margin"])
    padMain.SetRightMargin(0.05)
    padMain.SetTicks(1, 1)
    padRatio = canvas.cd(2)
    padRatio.SetPad(0, 0., 1, 0.35)
    padRatio.SetTopMargin(0)
    padRatio.SetBottomMargin(config["bottom_margin"])
    padRatio.SetLeftMargin(config["left_margin"])
    padRatio.SetRightMargin(0.05)
    padRatio.SetGridy()
    padRatio.SetTicks(1, 1)

    padMain.cd()
    if "minx" in config:
        minx = config["minx"]
    else:
        minx = dataStat.GetXaxis().GetXmin()

    if "maxx" in config:
        maxx = config["maxx"]
    else:
        maxx = dataStat.GetXaxis().GetXmax()

    h = ROOT.TH1I("h", "h", 100, 0, 1)
    h.Draw("axis")
    globalList.append(h)
    h.GetYaxis().SetRangeUser(config["miny"], config["maxy"])
    h.GetYaxis().SetTitleFont(43)
    h.GetYaxis().SetTitleSize(config["axis"]["font_size"])
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelSize(config["axis"]["font_size"] - 3)
    h.GetYaxis().SetTitleOffset(config["axis"]["y_offset"])
    h.GetYaxis().SetTitle(config["axis"]["y_title"])

    dataSyst_copy = dataSyst.Clone("{0}_copy".format(dataSyst.GetName()))
    dataSyst_copy.Draw("2")
    globalList.append(dataSyst_copy)

    dataStat_copy = dataStat.DrawCopy("same p e0 x0")
    globalList.append(dataStat_copy)

    if "data_minx" in config and "data_maxx" in config:
        points_to_be_removed = set()
        for ipoint in range(0, dataSyst_copy.GetN()):
            if dataSyst_copy.GetX()[ipoint] < config["data_minx"] or dataSyst_copy.GetX()[ipoint] > config["data_maxx"]:
                points_to_be_removed.add(ipoint)
        for ipoint in points_to_be_removed:
            dataSyst_copy.RemovePoint(ipoint)
        dataStat_copy.GetXaxis().SetRangeUser(config["data_minx"], config["data_maxx"])

    for t in config["theory"]:
        if not t["active"]: continue
        h = t["histogram"].Clone(t["gen"])
        h.Draw("same e0")
        globalList.append(h)
        h.SetLineColor(t["color"])
        h.SetLineStyle(t["line"])
        h.SetLineWidth(2)
        h.SetMarkerStyle(1)
        h.SetMarkerColor(t["color"])
        if "data_minx" in config and "data_maxx" in config:
            h.GetXaxis().SetRangeUser(config["data_minx"], config["data_maxx"])
        t["histogram_plot"] = h

    padRatio.cd()

    hRatio = ROOT.TH1C("hratio", "hratio", 100, 0, 1)
    hRatio.Draw("axis")
    hRatio.GetYaxis().SetTitle("MC / Data")
    hRatio.GetXaxis().SetTitleFont(43)
    hRatio.GetXaxis().SetTitleSize(config["axis"]["font_size"])
    hRatio.GetXaxis().SetLabelFont(43)
    hRatio.GetXaxis().SetLabelSize(config["axis"]["font_size"] - 3)
    hRatio.GetYaxis().SetTitleFont(43)
    hRatio.GetYaxis().SetTitleSize(config["axis"]["font_size"])
    hRatio.GetYaxis().SetLabelFont(43)
    hRatio.GetYaxis().SetLabelSize(config["axis"]["font_size"] - 3)
    hRatio.GetYaxis().SetTitleOffset(config["axis"]["y_offset"])
    hRatio.GetXaxis().SetTitleOffset(config["axis"]["x_offset"])
    hRatio.GetYaxis().SetRangeUser(config["minr"], config["maxr"])
    hRatio.GetYaxis().SetNdivisions(509)
    hRatio.GetXaxis().SetTitle(config["axis"]["x_title"])
    globalList.append(hRatio)

    ratioSyst = dataSyst_copy.Clone("ratioSyst")
    globalList.append(ratioSyst)
    ratioSyst.Draw("2")

    ratioStat = dataStat_copy.Clone("ratioStat")
    globalList.append(ratioStat)
    ratioStat.Draw("same p e0 x0")
    ratioStat.SetFillStyle(0)
    if "data_minx" in config and "data_maxx" in config:
        ratioStat.GetXaxis().SetRangeUser(config["data_minx"], config["data_maxx"])

    for ibin in range(0, ratioSyst.GetN()):
        ratioSyst.SetPointEYlow(ibin, ratioSyst.GetErrorYlow(ibin) / ratioSyst.GetY()[ibin])
        ratioSyst.SetPointEYhigh(ibin , ratioSyst.GetErrorYhigh(ibin) / ratioSyst.GetY()[ibin])
        ratioSyst.SetPoint(ibin, ratioSyst.GetX()[ibin], 1.0)

    for ibin in range(1, ratioStat.GetNbinsX()+1):
        ratioStat.SetBinError(ibin, ratioStat.GetBinError(ibin) / ratioStat.GetBinContent(ibin))
        ratioStat.SetBinContent(ibin, 1.0)

    for t in config["theory"]:
        if not t["active"]:
            continue
        r = t["histogram_plot"].Clone()
        for ibin in range(1, r.GetNbinsX() + 1):
            if t["histogram_plot"].GetBinContent(ibin) == 0:
                continue
            r.SetBinError(ibin, t["histogram_plot"].GetBinError(ibin) / t["histogram_plot"].GetBinContent(ibin))
            r.SetBinContent(ibin, t["histogram_plot"].GetBinContent(ibin) / dataStat_copy.GetBinContent(ibin))
        r.SetLineColor(t["color"])
        r.SetLineStyle(t["line"])
        r.SetLineWidth(2)
        r.SetMarkerStyle(20)
        r.SetMarkerSize(0)
        r.Draw("same e0")
        if "data_minx" in config and "data_maxx" in config:
            r.GetXaxis().SetRangeUser(config["data_minx"], config["data_maxx"])
        globalList.append(r)
        t["ratio"] = r

    padMain.cd()

    y1 = config["title"]["y"]
    y2 = y1 - 0.08 * len(config["title"]["text"])
    if "x" in config["title"]:
        x1 = config["title"]["x"]
    else:
        x1 = 0.19
    x2 = x1 + 0.36
    if x2 > 0.99:
        x2 = 0.99
    padMain.cd()
    paveALICE = ROOT.TPaveText(x1, y1, x2, y2, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(config["title"]["font_size"])
    paveALICE.SetTextAlign(12)
    for line in config["title"]["text"]: 
        paveALICE.AddText(line)
    paveALICE.Draw()

    if "theory_legend" in config and "n_columns" in config["theory_legend"]:
        n_leg_columns = config["theory_legend"]["n_columns"]
    else:
        max_length = 0
        for t in config["theory"]: 
            if not t["active"]:
                continue
            if len(t["title"]) > max_length:
                max_length = len(t["title"])

        if max_length > 40:
            n_leg_columns = 1
        else:
            n_leg_columns = 2

    if "theory_legend" in config and "y" in config["theory_legend"]:
        y1 = config["theory_legend"]["y"]
    else:
        y1 = y2 - 0.03
    active_t = len([t for t in config["theory"] if "histogram_plot" in t])
    y2 = y1 - 0.06 * active_t / n_leg_columns
    if "theory_legend" in config and "x" in config["theory_legend"]:
        x1 = config["theory_legend"]["x"]
    else:
        x1 = 0.16
    x2 = x1 + 0.50
    if x2 > 0.95:
        x2 = 0.95

    leg1 = ROOT.TLegend(x1, y1, x2, y2, "", "NB NDC")
    leg1.SetNColumns(n_leg_columns)
    globalList.append(leg1)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(config["theory_legend"]["font_size"])
    leg1.SetTextAlign(12)
    leg1.SetMargin(0.1)
    for t in config["theory"]:
        if not t["active"]:
            continue
        leg1.AddEntry(t["histogram_plot"], t["title"], "l")
    leg1.Draw()

    if "data_legend" in config and "y" in config["data_legend"]:
        y1 = config["data_legend"]["y"]
    else:
        y1 = y2 - 0.02
    y2 = y1 - 0.12
    if "data_legend" in config and "x" in config["data_legend"]:
        x1 = config["data_legend"]["x"]
    else:
        x1 = 0.16
    x2 = x1 + 0.30
    leg1 = ROOT.TLegend(x1, y1, x2, y2, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(config["data_legend"]["font_size"])
    leg1.SetTextAlign(12)
    leg1.SetMargin(0.2)
    leg1.AddEntry(dataStat_copy, "Data", "p")
    leg1.AddEntry(dataSyst_copy, "Syst. Unc. (data)", "f")
    leg1.Draw()

    padRatio.RedrawAxis("g")
    padRatio.RedrawAxis()
    padMain.RedrawAxis()

    return canvas

def NormalizeData(config, d0jet_stat, d0jet_syst, incl_stat, incl_syst):
    normalizator = HistogramNormalizator.Normalizator(d0jet_stat, "rate")
    normalizator.fNormalizationGraph = incl_syst
    normalizator.fNormalizationHistogram = incl_stat
    normalizator.fNormalizationXmin = config["min_jet_pt"]
    normalizator.fNormalizationXmax = config["max_jet_pt"]
    normalizator.NormalizeHistogram(d0jet_syst)
    d0jet_stat = normalizator.fNormalizedHistogram
    d0jet_syst = normalizator.fNormalizedGraph
    return d0jet_stat, d0jet_syst

def NormalizeTheory(config):
    for t in config["theory"]:
        if not t["active"]:
            continue
        if not "inclusive" in t or not "histogram" in t["inclusive"] or not "histogram" in t:
            continue
        normalizator = HistogramNormalizator.Normalizator(t["histogram"], "rate")
        normalizator.fNormalizationHistogram = t["inclusive"]["histogram"]
        normalizator.fNormalizationXmin = config["min_jet_pt"]
        normalizator.fNormalizationXmax = config["max_jet_pt"]
        normalizator.NormalizeHistogram()
        t["histogram"] = normalizator.fNormalizedHistogram

def main(config):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    d0jet_stat, d0jet_syst = GetD0JetCrossSection(config["input_path"], config["data"])
    incl_stat, incl_syst = GetInclJetCrossSection([config["min_jet_pt"], config["max_jet_pt"]])
    GetD0JetTheoryCrossSectionAll(config, d0jet_stat.GetXaxis())
    GetInclusiveJetTheoryCrossSectionAll(config)

    d0jet_stat, d0jet_syst = NormalizeData(config, d0jet_stat, d0jet_syst, incl_stat, incl_syst)
    NormalizeTheory(config)

    canvas = PlotCrossSections(d0jet_stat, d0jet_syst, config)
    canvas.SaveAs("{}/{}.pdf".format(config["input_path"], canvas.GetName()))
    canvas.SaveAs("{}/{}.C".format(config["input_path"], canvas.GetName()))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Jet pt spectrum theory comparison.')
    parser.add_argument('yaml', metavar='conf.yaml')

    args = parser.parse_args()

    f = open(args.yaml, 'r')
    yconfig = yaml.load(f)
    f.close()

    main(yconfig)

    IPython.embed()
