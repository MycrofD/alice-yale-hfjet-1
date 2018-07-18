#!/usr/bin/env python
# python script to compare measured  D0 and inclusive jet cross sections with theory
# us it with
# D0JetRateVsJetPt_CompareMonteCarlo_OnlyRatio_Paper.yaml
# D0JetRateVsJetPt_CompareMonteCarlo_SpectrumRatio_Paper.yaml
# D0JetRateVsJetPt_ComparePowheg_Paper.yaml
# D0JetRateVsJetPt_CompareFullVsCharged.yaml
# D0JetRateVsJetPt_ComparePtDCut.yaml

import argparse
import math
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

def GetInclJetCrossSection():
    return LoadInclusiveJetSpectrum.GetCrossSection()

def GetD0JetTheoryCrossSectionAll(config, axis):
    return LoadTheoryCrossSections.GetD0JetTheoryCrossSectionAll(config, axis)

def GetInclusiveJetTheoryCrossSectionAll(config):
    return LoadTheoryCrossSections.GetInclusiveJetTheoryCrossSectionAll(config)

def PlotCrossSections(d0jet_stat, d0jet_syst, incl_stat, incl_syst, inclusive_jet_cross_sections, config):
    no_box_sys = ["[]", "||", "0[]", "0||"]
    if not "data_systematics_style" in config:
        config["data_systematics_style"] = "2"
    config["data_box_systematics"] = not (config["data_systematics_style"] in no_box_sys)

    d0jet_syst_copy = d0jet_syst.Clone("{0}_copy".format(d0jet_syst.GetName()))
    d0jet_syst_copy.SetFillStyle(1001)
    d0jet_syst_copy.SetLineWidth(0)
    d0jet_syst_copy.SetFillColor(ROOT.kGray)
    globalList.append(d0jet_syst_copy)

    d0jet_stat_copy = d0jet_stat.Clone("d0jet_stat_copy")
    d0jet_stat_copy.SetLineColor(ROOT.kRed + 2)
    d0jet_stat_copy.SetMarkerColor(ROOT.kRed + 2)
    d0jet_stat_copy.SetMarkerStyle(ROOT.kFullCircle)
    d0jet_stat_copy.SetMarkerSize(1.2)
    globalList.append(d0jet_stat_copy)

    incl_syst_copy = incl_syst.Clone("incl_syst_copy")
    globalList.append(incl_syst_copy)
    incl_syst_copy.SetFillStyle(1001)
    incl_syst_copy.SetLineWidth(0)
    incl_syst_copy.SetFillColor(ROOT.kGray)

    incl_stat_copy = incl_stat.Clone("incl_stat_copy")
    globalList.append(incl_stat_copy)
    incl_stat_copy.SetLineColor(ROOT.kBlack)
    incl_stat_copy.SetMarkerColor(ROOT.kBlack)
    incl_stat_copy.SetMarkerStyle(ROOT.kFullSquare)
    incl_stat_copy.SetMarkerSize(1.2)

    for t in inclusive_jet_cross_sections.itervalues():
        incl_h = t["histogram"].Clone("_".join([t["gen"], t["proc"]]))
        globalList.append(incl_h)
        incl_h.SetLineColor(t["color"])
        incl_h.SetLineStyle(t["line"])
        incl_h.SetLineWidth(1)
        incl_h.SetMarkerColor(t["color"])
        t["histogram_plot"] = incl_h

    for t in config["theory"]:
        if not t["active"]:
            continue
        if not "type" in t:
            t["type"] = "stat-only"
        if t["type"] == "stat-only":
            d0jet_h = t["histogram"].Clone("_".join([t["gen"], t["proc"]]))
            globalList.append(d0jet_h)
            d0jet_h.SetLineColor(t["color"])
            d0jet_h.SetLineStyle(t["line"])
            d0jet_h.SetLineWidth(1)
            d0jet_h.SetMarkerColor(t["color"])
            t["histogram_plot"] = d0jet_h
        elif t["type"] == "stat+syst":
            if not "systematics_style" in t:
                t["systematics_style"] = "2"
            t["box_systematics"] = not (t["systematics_style"] in no_box_sys)
            
            hSyst = t["systematics"].Clone("{0}_copy".format(t["systematics"].GetName()))
            hSyst.SetLineColor(t["color"])
            hSyst.SetFillColor(t["color"])
            hSyst.SetLineWidth(1)
            if "line" in t:
                hSyst.SetLineStyle(t["line"])
            if "fill" in t:
                hSyst.SetFillStyle(t["fill"])
            else:
                hSyst.SetFillStyle(0)
            globalList.append(hSyst)
            t["systematics_plot"] = hSyst

            hStat = t["histogram"].Clone("{0}_copy".format(t["histogram"].GetName()))
            hStat.SetMarkerStyle(getattr(ROOT, t["marker"]))
            hStat.SetLineColor(t["color"])
            hStat.SetMarkerColor(t["color"])
            globalList.append(hStat)
            t["histogram_plot"] = hStat

    normalizator = HistogramNormalizator.Normalizator(d0jet_stat_copy, "ratio")
    normalizator.fNormalizationHistogram = incl_stat
    normalizator.fNormalizationGraph = incl_syst
    normalizator.NormalizeHistogram(d0jet_syst)
    ratioStat = normalizator.fNormalizedHistogram
    ratioSyst = normalizator.fNormalizedGraph

    globalList.append(ratioStat)
    globalList.append(ratioSyst)

    for t in config["theory"]:
        if not t["active"]:
            continue
        if not t["inclusive"]:
            continue

        if t["type"] == "stat-only":
            incl_h = t["inclusive"]["histogram_plot"]
            d0jet_h = t["histogram_plot"]
            normalizator = HistogramNormalizator.Normalizator(d0jet_h, "ratio")
            normalizator.fNormalizationHistogram = incl_h
            normalizator.NormalizeHistogram()
            hratio = normalizator.fNormalizedHistogram
            hratio.SetName("_".join([t["gen"], t["proc"], "over", t["inclusive"]["gen"], t["inclusive"]["proc"]]))
            t["ratio_histogram"] = hratio
            hStat = t["ratio_histogram"].Clone("{0}_copy".format(t["ratio_histogram"].GetName()))
            t["ratio_histogram_plot"] = hStat
        elif t["type"] == "stat+syst":
            hSyst = t["ratio_systematics"].Clone("{0}_copy".format(t["ratio_systematics"].GetName()))
            hSyst.SetLineColor(t["color"])
            hSyst.SetFillColor(t["color"])
            hSyst.SetLineWidth(1)
            if "line" in t:
                hSyst.SetLineStyle(t["line"])
            if "fill" in t:
                hSyst.SetFillStyle(t["fill"])
            else:
                hSyst.SetFillStyle(0)
            globalList.append(hSyst)
            t["ratio_systematics_plot"] = hSyst

            hStat = t["ratio_histogram"].Clone("{0}_copy".format(t["ratio_histogram"].GetName()))
            hStat.SetMarkerStyle(getattr(ROOT, t["marker"]))
            hStat.SetLineColor(t["color"])
            hStat.SetMarkerColor(t["color"])
            globalList.append(hStat)
            t["ratio_histogram_plot"] = hStat

    # Done with all the preparations, now plotting

    if "plot_ratio_only" in config and config["plot_ratio_only"]:
        canvas = DrawRatioCanvas(config, ratioSyst, ratioStat)
    else:
        canvas = DrawTwoPanelCanvas(config, d0jet_stat_copy, d0jet_syst_copy, incl_stat_copy, incl_syst_copy, ratioSyst, ratioStat, inclusive_jet_cross_sections)

    return canvas

def GenerateHistogramAxis(config, xmin, xmax):
    hAxis = ROOT.TH1I("axis", "axis", 1000, xmin, xmax)
    hAxis.GetYaxis().SetRangeUser(config["miny"], config["maxy"])
    hAxis.GetYaxis().SetTitleFont(43)
    hAxis.GetYaxis().SetTitleSize(config["axis"]["font_size"])
    hAxis.GetYaxis().SetTitleOffset(config["axis"]["y_offset"])
    hAxis.GetYaxis().SetLabelFont(43)
    hAxis.GetYaxis().SetLabelSize(config["axis"]["font_size"] - 3)
    hAxis.GetYaxis().SetTitle(config["axis"]["y_title"])
    return hAxis

def GenerateHistogramRatioAxis(config, xmin, xmax):
    hAxisRatio = ROOT.TH1I("axis", "axis", 1000, xmin, xmax)
    hAxisRatio.GetYaxis().SetTitle("#it{R}(#it{p}_{T,jet}^{ch})")
    hAxisRatio.GetXaxis().SetTitleFont(43)
    hAxisRatio.GetXaxis().SetTitleSize(config["axis"]["font_size"])
    hAxisRatio.GetXaxis().SetLabelFont(43)
    hAxisRatio.GetXaxis().SetLabelSize(config["axis"]["font_size"] - 3)
    hAxisRatio.GetYaxis().SetTitleFont(43)
    hAxisRatio.GetYaxis().SetTitleSize(config["axis"]["font_size"])
    hAxisRatio.GetYaxis().SetLabelFont(43)
    hAxisRatio.GetYaxis().SetLabelSize(config["axis"]["font_size"] - 3)
    hAxisRatio.GetYaxis().SetTitleOffset(config["axis"]["y_offset"])
    hAxisRatio.GetXaxis().SetTitleOffset(config["axis"]["x_offset"])
    hAxisRatio.GetYaxis().SetRangeUser(config["minr"], config["maxr"])
    hAxisRatio.GetYaxis().SetNdivisions(509)
    hAxisRatio.GetXaxis().SetTitle(config["axis"]["x_title"])
    return hAxisRatio

def DrawTwoPanelCanvas(config, d0jet_stat_copy, d0jet_syst_copy, incl_stat_copy, incl_syst_copy, ratioSyst, ratioStat, inclusive_jet_cross_sections):
    hAxis = GenerateHistogramAxis(config, d0jet_stat_copy.GetXaxis().GetBinLowEdge(1), d0jet_stat_copy.GetXaxis().GetBinUpEdge(d0jet_stat_copy.GetXaxis().GetNbins()))
    hAxisRatio = GenerateHistogramRatioAxis(config, d0jet_stat_copy.GetXaxis().GetBinLowEdge(1), d0jet_stat_copy.GetXaxis().GetBinUpEdge(d0jet_stat_copy.GetXaxis().GetNbins()))

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
    if config["logy"]: padMain.SetLogy()
    padRatio = canvas.cd(2)
    padRatio.SetPad(0, 0., 1, 0.35)
    padRatio.SetTopMargin(0)
    padRatio.SetBottomMargin(config["bottom_margin"])
    padRatio.SetLeftMargin(config["left_margin"])
    padRatio.SetRightMargin(0.05)
    padRatio.SetGridy()
    padRatio.SetTicks(1, 1)

    padMain.cd()
    hAxis_copy = hAxis.DrawCopy("axis")
    hAxis_copy.GetYaxis().SetTitleOffset(2.2)
    globalList.append(hAxis_copy)

    d0jet_syst_copy.Draw("2")
    d0jet_stat_copy.Draw("same p x0 e0")
    incl_syst_copy.Draw("2")
    incl_stat_copy.Draw("same p x0 e0")

    for t in inclusive_jet_cross_sections.itervalues():
        if "histogram_plot" in t:
            incl_h = t["histogram_plot"]
            incl_h.Draw("same e0")

    for t in config["theory"]:
        if not t["active"]:
            continue
        if t["type"] == "stat-only":
            d0jet_h = t["histogram_plot"]
            d0jet_h.Draw("same e0")
        elif t["type"] == "stat+syst":
            hSyst = t["systematics_plot"]
            hSyst.Draw(t["systematics_style"])

            hStat = t["histogram_plot"]
            hStat.Draw("same p e0 x0")

    padRatio.cd()
    hAxisRatio_copy = hAxisRatio.DrawCopy("axis")
    globalList.append(hAxisRatio_copy)

    ratioSyst.Draw("2")
    ratioStat.Draw("same p x0 e0")

    for t in config["theory"]:
        if not t["active"]:
            continue
        if t["type"] == "stat-only":
            hratio = t["ratio_histogram"]
            hratio.Draw("same")
        elif t["type"] == "stat+syst":
            hSyst = t["ratio_systematics_plot"]
            hSyst.Draw(t["systematics_style"])

            hStat = t["ratio_histogram_plot"]
            hStat.Draw("same p e0 x0")

    # Now plotting labels

    padMain.cd()

    if "y" in config["title"]:
        y1 = config["title"]["y"]
    else:
        y1 = 0.90
    if "x" in config["title"]:
        x1 = config["title"]["x"]
    else:
        x1 = 0.19
    y2 = y1 - 0.07 * len(config["title"]["text"])
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
    paveALICE.SetTextAlign(13)
    for line in config["title"]["text"]: 
        paveALICE.AddText(line)
    paveALICE.Draw()

    if "theory_legend" in config and "n_columns" in config["theory_legend"]:
        n_leg_columns = config["theory_legend"]["n_columns"]
    else:
        max_length = 0
        for t in config["theory"]: 
            if not "histogram_plot" in t:
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
    active_t += len([t for t in inclusive_jet_cross_sections.itervalues() if "histogram_plot" in t])
    y2 = y1 - 0.06 * active_t / n_leg_columns
    if "theory_legend" in config and "x" in config["theory_legend"]:
        x1 = config["theory_legend"]["x"]
    else:
        x1 = 0.16
    x2 = x1 + 0.50
    if x2 > 0.95:
        x2 = 0.95

    leg1 = ROOT.TLegend(x1, y1, x2, y2, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetNColumns(n_leg_columns)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(config["theory_legend"]["font_size"])
    leg1.SetTextAlign(12)
    leg1.SetMargin(0.10)
    for t in config["theory"]: 
        if not "histogram_plot" in t:
            continue
        if t["type"] == "stat-only":
            leg1.AddEntry(t["histogram_plot"], "D^{{0}} Jets, {}".format(t["title"]), "l")
        elif t["type"] == "stat+syst":
            if t["box_systematics"]:
                entry = leg1.AddEntry(None, "D^{{0}} Jets, {}".format(t["title"]), "pf")
                entry.SetLineColor(t["systematics_plot"].GetLineColor())
                entry.SetLineWidth(t["systematics_plot"].GetLineWidth())
                entry.SetFillColor(t["systematics_plot"].GetFillColor())
                entry.SetFillStyle(t["systematics_plot"].GetFillStyle())
                entry.SetMarkerColor(t["histogram_plot"].GetMarkerColor())
                entry.SetMarkerStyle(t["histogram_plot"].GetMarkerStyle())
            else:
                leg1.AddEntry(t["histogram_plot"], "D^{{0}} Jets, {}".format(t["title"]), "p")
    for t in inclusive_jet_cross_sections.itervalues():
        if "histogram_plot" in t: 
            leg1.AddEntry(t["histogram_plot"], "Inclusive Jets, {}".format(t["title"]), "l")

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
    if config["data_box_systematics"]:
        entry = leg1.AddEntry(None, "D^{0} Jets, #it{p}_{T,D} > 3 GeV/#it{c}", "pf")
        entry.SetFillStyle(d0jet_syst_copy.GetFillStyle())
        entry.SetFillColor(d0jet_syst_copy.GetFillColor())
        entry.SetLineColor(d0jet_syst_copy.GetFillColor())
        entry.SetMarkerColor(d0jet_stat_copy.GetMarkerColor())
        entry.SetMarkerStyle(d0jet_stat_copy.GetMarkerStyle())
    else:
        entry = leg1.AddEntry(None, "D^{0} Jets, #it{p}_{T,D} > 3 GeV/#it{c}", "p")
        entry.SetMarkerColor(d0jet_stat_copy.GetMarkerColor())
        entry.SetMarkerStyle(d0jet_stat_copy.GetMarkerStyle())

    entry = leg1.AddEntry(None, "Inclusive Jets", "pf")
    entry.SetFillStyle(incl_syst_copy.GetFillStyle())
    entry.SetFillColor(incl_syst_copy.GetFillColor())
    entry.SetLineColor(incl_syst_copy.GetFillColor())
    entry.SetMarkerColor(incl_stat_copy.GetMarkerColor())
    entry.SetMarkerStyle(incl_stat_copy.GetMarkerStyle())
    leg1.Draw()

    padRatio.RedrawAxis("g")
    padRatio.RedrawAxis()
    padMain.RedrawAxis()

    return canvas

def DrawRatioCanvas(config, ratioSyst, ratioStat):
    hAxisRatio = GenerateHistogramRatioAxis(config, ratioStat.GetXaxis().GetBinLowEdge(1), ratioStat.GetXaxis().GetBinUpEdge(ratioStat.GetXaxis().GetNbins()))
    hAxisRatio.GetXaxis().SetTitleOffset(config["axis"]["x_offset"] * 0.35)

    ratioSyst_copy = ratioSyst.Clone()
    ratioStat_copy = ratioStat.Clone()
    globalList.append(ratioStat_copy)
    globalList.append(ratioSyst_copy)

    cname = "{}_{}_Ratio".format(config["name_prefix"], config["name"])
    if "canvas_h" in config:
        canvas_h = config["canvas_h"]
    else:
        canvas_h = 700
    if "canvas_w" in config:
        canvas_w = config["canvas_w"]
    else:
        canvas_w = 700
    canvas_ratio = ROOT.TCanvas(cname, cname, canvas_w, canvas_h)
    globalList.append(canvas_ratio)
    canvas_ratio.SetTicks(1, 1)
    canvas_ratio.SetLeftMargin(config["left_margin"])
    canvas_ratio.SetBottomMargin(config["bottom_margin"] * 0.35)
    canvas_ratio.SetRightMargin(0.05)
    canvas_ratio.SetTopMargin(0.05)

    canvas_ratio.cd()
    hAxisRatio.Draw("axis")
    globalList.append(hAxisRatio)

    # Plot data systematic uncertainties, if they are plotted as a filled box
    if config["data_box_systematics"]:
        ratioSyst_copy.Draw(config["data_systematics_style"])

    # Plot theory systematic uncertainties (if they are plotted as a filled box)
    for t in config["theory"]:
        if not t["active"]:
            continue
        if t["type"] == "stat+syst" and t["box_systematics"]:
            hSyst = t["ratio_systematics_plot"]
            hSyst.Draw(t["systematics_style"])

    # Plot data statistical uncertainties
    ratioStat_copy.Draw("same p e0 x0")    

    # Plot theory statistical uncertainties and systematic uncertainties (if they are not plotted as filled boxes)
    for t in config["theory"]:
        if not t["active"]:
            continue
        hStat = t["ratio_histogram_plot"]
        if t["type"] == "stat-only":
            hStat.Draw("same e0")
        elif t["type"] == "stat+syst":
            hStat.Draw("same p e0 x0")
            if not t["box_systematics"]:
                hSyst = t["ratio_systematics_plot"]
                hSyst.Draw(t["systematics_style"])

    # Plot data systematic uncertainties, if they are not plotted as filled boxes
    if not config["data_box_systematics"]:
        for ipoint in range(0, ratioSyst_copy.GetN()):
            ratioSyst_copy.SetPointEXhigh(ipoint, 0)
            ratioSyst_copy.SetPointEXlow(ipoint, 0)
        ratioSyst_copy.SetLineColor(ratioStat_copy.GetLineColor())
        ratioSyst_copy.Draw(config["data_systematics_style"])

    # Now plotting labels

    if "y" in config["title"]:
        y1 = config["title"]["y"]
    else:
        y1 = 0.90
    if "x" in config["title"]:
        x1 = config["title"]["x"]
    else:
        x1 = 0.19
    y2 = y1 - 0.045 * len(config["title"]["text"])
    x2 = x1 + 0.36
    if x2 > 0.99:
        x2 = 0.99

    paveALICE = ROOT.TPaveText(x1, y1, x2, y2, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(config["title"]["font_size"])
    paveALICE.SetTextAlign(13)
    for line in config["title"]["text"]: 
        paveALICE.AddText(line)
    paveALICE.Draw()

    if "theory_legend" in config and "n_columns" in config["theory_legend"]:
        n_leg_columns = config["theory_legend"]["n_columns"]
    else:
        max_length = 0
        for t in config["theory"]: 
            if not "histogram_plot" in t:
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
    y2 = y1 - 0.04 * active_t / n_leg_columns
    if "theory_legend" in config and "x" in config["theory_legend"]:
        x1 = config["theory_legend"]["x"]
    else:
        x1 = 0.16
    x2 = x1 + 0.50
    if x2 > 0.95:
        x2 = 0.95

    leg1 = ROOT.TLegend(x1, y1, x2, y2, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetNColumns(n_leg_columns)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(config["theory_legend"]["font_size"])
    leg1.SetTextAlign(12)
    leg1.SetMargin(0.1)
    for t in config["theory"]: 
        if not "histogram_plot" in t:
            continue
        if t["type"] == "stat-only":
            leg1.AddEntry(t["histogram_plot"], t["title"], "l")
        elif t["type"] == "stat+syst":
            if t["box_systematics"]:
                entry = leg1.AddEntry(None, t["title"], "pf")
                entry.SetLineColor(t["ratio_systematics_plot"].GetLineColor())
                entry.SetLineWidth(t["ratio_systematics_plot"].GetLineWidth())
                entry.SetFillColor(t["ratio_systematics_plot"].GetFillColor())
                entry.SetFillStyle(t["ratio_systematics_plot"].GetFillStyle())
                entry.SetMarkerColor(t["ratio_histogram_plot"].GetMarkerColor())
                entry.SetMarkerStyle(t["ratio_histogram_plot"].GetMarkerStyle())
            else:
                leg1.AddEntry(t["ratio_histogram_plot"], t["title"], "p")

    leg1.Draw()

    if "data_legend" in config and "y" in config["data_legend"]:
        y1 = config["data_legend"]["y"]
    else:
        y1 = y2 - 0.02
    y2 = y1 - 0.05
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

    if config["data_box_systematics"]:
        entry = leg1.AddEntry(None, "Data", "pf")
        entry.SetFillStyle(ratioSyst_copy.GetFillStyle())
        entry.SetFillColor(ratioSyst_copy.GetFillColor())
        entry.SetLineColor(ratioSyst_copy.GetFillColor())
        entry.SetMarkerColor(ratioStat_copy.GetMarkerColor())
        entry.SetMarkerStyle(ratioStat_copy.GetMarkerStyle())
    else:
        entry = leg1.AddEntry(None, "Data", "p")
        entry.SetMarkerColor(ratioStat_copy.GetMarkerColor())
        entry.SetMarkerStyle(ratioStat_copy.GetMarkerStyle())
    leg1.Draw()

    canvas_ratio.RedrawAxis("g")
    canvas_ratio.RedrawAxis()

    return canvas_ratio

def main(config):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    d0jet_stat, d0jet_syst = GetD0JetCrossSection(config["input_path"], config["data"])
    incl_stat, incl_syst = GetInclJetCrossSection()
    GetD0JetTheoryCrossSectionAll(config, d0jet_stat.GetXaxis())
    inclusive_jet_cross_sections = GetInclusiveJetTheoryCrossSectionAll(config)
    canvas = PlotCrossSections(d0jet_stat, d0jet_syst, incl_stat, incl_syst, inclusive_jet_cross_sections, config)
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
