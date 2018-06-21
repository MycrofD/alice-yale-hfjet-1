#!/usr/bin/env python
# python script to compare measured  D0 and inclusive jet cross sections with theory
# us it with
# JetPtSpectrum_ComparePYTHIA6vs8.yaml
# JetPtSpectrum_CompareHerwig.yaml
# JetPtSpectrum_ComparePowheg.yaml
# JetPtSpectrum_CompareTheory.yaml
# JetPtSpectrum_CompareFullVsCharged.yaml
# JetPtSpectrum_ComparePtDCut.yaml

import argparse
import math
import yaml
import IPython
import ROOT
import DMesonJetUtils
import LoadInclusiveJetSpectrum
import LoadTheoryCrossSections

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
        incl_h.SetLineWidth(2)
        incl_h.SetMarkerColor(t["color"])
        t["histogram_plot"] = incl_h

    for t in config["theory"]:
        if not t["active"]: continue
        if not t["inclusive"]: continue
        d0jet_h = t["histogram"].Clone("_".join([t["gen"], t["proc"]]))
        globalList.append(d0jet_h)
        d0jet_h.SetLineColor(t["color"])
        d0jet_h.SetLineStyle(t["line"])
        d0jet_h.SetLineWidth(2)
        d0jet_h.SetMarkerColor(t["color"])
        t["histogram_plot"] = d0jet_h

    ratioSyst = d0jet_syst_copy.Clone("ratioSyst")
    globalList.append(ratioSyst)

    ratioStat = d0jet_stat_copy.Clone("ratioStat")
    globalList.append(ratioStat)

    xsec_incl_tot = 0.0
    stat_xsec_incl_tot2 = 0.0
    syst_xsec_incl_tot = 0.0

    xsec_d0_tot = 0.0
    stat_xsec_d0_tot2 = 0.0
    syst_xsec_d0_tot = 0.0

    avg_tot = 0.0
    stat_avg_tot = 0.0
    syst_avg_tot = 0.0

    avg_flat = 0.0
    sum_of_w_avg_flat = 0.0
    stat_avg_flat2 = 0.0
    syst_avg_flat2 = 0.0

    for ibin in range(0, ratioSyst.GetN()):
        ratioSyst.SetPoint(ibin, ratioSyst.GetX()[ibin], d0jet_syst_copy.GetY()[ibin] / incl_syst_copy.GetY()[ibin])
        syst_erry2 = ((d0jet_syst_copy.GetErrorY(ibin) / d0jet_syst_copy.GetY()[ibin]) ** 2 + (incl_syst_copy.GetErrorY(ibin) / incl_syst_copy.GetY()[ibin]) ** 2 - 2 * 0.035 ** 2) * ratioSyst.GetY()[ibin] ** 2
        syst_erry = math.sqrt(syst_erry2)
        ratioSyst.SetPointError(ibin, ratioSyst.GetErrorX(ibin), ratioSyst.GetErrorX(ibin), syst_erry, syst_erry)

        xsec_incl_tot += incl_stat_copy.GetBinContent(ibin + 1)
        stat_xsec_incl_tot2 += incl_stat_copy.GetBinError(ibin + 1) ** 2
        syst_xsec_incl_tot += incl_syst_copy.GetErrorY(ibin)  # take the weighted average of the rel unc

        ratioStat.SetBinContent(ibin + 1, d0jet_stat_copy.GetBinContent(ibin + 1) / incl_stat_copy.GetBinContent(ibin + 1))
        stat_err_y2 = ((d0jet_stat_copy.GetBinError(ibin + 1) / d0jet_stat_copy.GetBinContent(ibin + 1)) ** 2 + (incl_stat_copy.GetBinError(ibin + 1) / incl_stat_copy.GetBinContent(ibin + 1)) ** 2) * ratioStat.GetBinContent(ibin + 1) ** 2
        stat_err_y = math.sqrt(stat_err_y2)
        ratioStat.SetBinError(ibin + 1, stat_err_y)

        tot_err_y2 = stat_err_y2 + syst_erry2
        tot_err_y = math.sqrt(tot_err_y2)

        xsec_d0_tot += d0jet_stat_copy.GetBinContent(ibin + 1)
        stat_xsec_d0_tot2 += d0jet_stat_copy.GetBinError(ibin + 1) ** 2
        syst_xsec_d0_tot += d0jet_syst_copy.GetErrorY(ibin)  # take the weithed average of the rel unc

        if ratioSyst.GetX()[ibin] > 8:
            avg_flat += ratioStat.GetBinContent(ibin + 1) / tot_err_y
            sum_of_w_avg_flat += 1.0 / tot_err_y
            stat_avg_flat2 += stat_err_y2 / tot_err_y2
            syst_avg_flat2 += syst_erry2 / tot_err_y2

        print("Bin {}, x = {}, ratio = {} +/- {} (stat) +/- {} (syst)".format(ibin, ratioSyst.GetX()[ibin], ratioStat.GetBinContent(ibin + 1), stat_err_y, syst_erry))

    stat_xsec_incl_tot = math.sqrt(stat_xsec_incl_tot2)
    stat_xsec_d0_tot = math.sqrt(stat_xsec_d0_tot2)

    avg_tot = xsec_d0_tot / xsec_incl_tot
    stat_avg_tot = math.sqrt(stat_xsec_incl_tot2 / xsec_incl_tot ** 2 + stat_xsec_d0_tot2 / xsec_d0_tot ** 2) * avg_tot

    syst_xsec_incl_tot2 = syst_xsec_incl_tot ** 2
    syst_xsec_d0_tot2 = syst_xsec_d0_tot ** 2

    syst_avg_tot = math.sqrt(syst_xsec_incl_tot2 / xsec_incl_tot ** 2 + syst_xsec_d0_tot2 / xsec_d0_tot ** 2) * avg_tot

    avg_flat /= sum_of_w_avg_flat
    stat_avg_flat = math.sqrt(stat_avg_flat2) / sum_of_w_avg_flat
    syst_avg_flat = math.sqrt(syst_avg_flat2) / sum_of_w_avg_flat

    print("Average ratio in full range: {} +/- {} (stat) +/- {} (syst)".format(avg_tot, stat_avg_tot, syst_avg_tot))
    print("Average ratio in flat region (pt > 8 GeV/c): {} +/- {} (stat) +/- {} (syst)".format(avg_flat, stat_avg_flat, syst_avg_flat))

    for t in config["theory"]:
        if not t["active"]: continue
        if not t["inclusive"]: continue
        incl_h = t["inclusive"]["histogram_plot"]
        d0jet_h = t["histogram_plot"]
        hratio = d0jet_h.Clone("_".join([t["gen"], t["proc"], "over", t["inclusive"]["gen"], t["inclusive"]["proc"]]))
        hratio.Divide(incl_h)
        t["ratio_histogram"] = hratio

    # Done with all the preparations, now plotting

    # First plot the cross section + ratio panel
    canvas = DrawTwoPanelCanvas(config, d0jet_stat_copy, d0jet_syst_copy, incl_stat_copy, incl_syst_copy, ratioSyst, ratioStat, inclusive_jet_cross_sections)

    # Now plot the ratio only
    canvas_ratio = DrawRatioCanvas(config, ratioSyst, ratioStat)

    return canvas, canvas_ratio

def GenerateHistogramAxis(config, xmin, xmax):
    hAxis = ROOT.TH1I("axis", "axis", 1000, xmin, xmax)
    hAxis.GetYaxis().SetTitleFont(43)
    hAxis.GetYaxis().SetTitleSize(23)
    hAxis.GetYaxis().SetLabelFont(43)
    hAxis.GetYaxis().SetLabelSize(22)
    hAxis.GetYaxis().SetTitleOffset(1.7)
    hAxis.GetYaxis().SetRangeUser(config["D0JetRate"]["miny"], config["D0JetRate"]["maxy"])
    if "y_axis_title" in config:
        hAxis.GetYaxis().SetTitle(config["y_axis_title"])
    return hAxis

def GenerateHistogramRatioAxis(config, xmin, xmax):
    hAxisRatio = ROOT.TH1I("axis", "axis", 1000, xmin, xmax)
    hAxisRatio.GetXaxis().SetTitleFont(43)
    hAxisRatio.GetXaxis().SetTitleSize(23)
    hAxisRatio.GetXaxis().SetTitleOffset(2.9)
    hAxisRatio.GetXaxis().SetLabelFont(43)
    hAxisRatio.GetYaxis().SetTitle("#it{R}(#it{p}_{T,jet}^{ch})")
    hAxisRatio.GetYaxis().SetTitleFont(43)
    hAxisRatio.GetYaxis().SetTitleSize(23)
    hAxisRatio.GetYaxis().SetTitleOffset(1.9)
    hAxisRatio.GetYaxis().SetRangeUser(config["D0JetRate"]["minr"], config["D0JetRate"]["maxr"])
    hAxisRatio.GetXaxis().SetLabelSize(22)
    hAxisRatio.GetYaxis().SetLabelFont(43)
    hAxisRatio.GetYaxis().SetLabelSize(22)
    hAxisRatio.GetYaxis().SetNdivisions(509)
    if "x_axis_title" in config:
        hAxisRatio.GetXaxis().SetTitle(config["x_axis_title"])
    return hAxisRatio

def DrawTwoPanelCanvas(config, d0jet_stat_copy, d0jet_syst_copy, incl_stat_copy, incl_syst_copy, ratioSyst, ratioStat, inclusive_jet_cross_sections):
    hAxis = GenerateHistogramAxis(config, d0jet_stat_copy.GetXaxis().GetBinLowEdge(1), d0jet_stat_copy.GetXaxis().GetBinUpEdge(d0jet_stat_copy.GetXaxis().GetNbins()))
    hAxisRatio = GenerateHistogramRatioAxis(config, d0jet_stat_copy.GetXaxis().GetBinLowEdge(1), d0jet_stat_copy.GetXaxis().GetBinUpEdge(d0jet_stat_copy.GetXaxis().GetNbins()))

    cname = "{}_{}".format(config["D0JetRate"]["name_prefix"], config["name"])
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
    padMain.SetLeftMargin(0.15)
    padMain.SetTicks(1, 1)
    if config["logy"]: padMain.SetLogy()
    padRatio = canvas.cd(2)
    padRatio.SetPad(0, 0., 1, 0.35)
    padRatio.SetTopMargin(0)
    padRatio.SetBottomMargin(0.27)
    padRatio.SetLeftMargin(0.15)
    padRatio.SetGridy()
    padRatio.SetTicks(1, 1)

    padMain.cd()
    hAxis.Draw("axis")
    globalList.append(hAxis)

    d0jet_syst_copy.Draw("2")
    d0jet_stat_copy.Draw("same p x0 e0")
    incl_syst_copy.Draw("2")
    incl_stat_copy.Draw("same p x0 e0")

    for t in inclusive_jet_cross_sections.itervalues():
        incl_h = t["histogram_plot"]
        incl_h.Draw("same e0")

    for t in config["theory"]:
        if not t["active"]: continue
        if not t["inclusive"]: continue
        d0jet_h = t["histogram_plot"]
        d0jet_h.Draw("same e0")

    padRatio.cd()
    hAxisRatio.Draw("axis")
    globalList.append(hAxisRatio)

    ratioSyst.Draw("2")
    ratioStat.Draw("same p x0 e0")

    for t in config["theory"]:
        if not t["active"]: continue
        if not t["inclusive"]: continue
        hratio = t["ratio_histogram"]
        hratio.Draw("same")

    # Now plotting labels

    padMain.cd()

    y1 = 0.87
    y2 = y1 - 0.06 * len(config["D0JetRate"]["title"])
    paveALICE = ROOT.TPaveText(0.16, y1, 0.55, y2, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(22)
    paveALICE.SetTextAlign(12)
    for line in config["D0JetRate"]["title"]: paveALICE.AddText(line)
    paveALICE.Draw()

    if "theory_legend" in config["D0JetRate"] and "n_columns" in config["D0JetRate"]["theory_legend"]:
        n_leg_columns = config["D0JetRate"]["theory_legend"]["n_columns"]
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

    if "theory_legend" in config["D0JetRate"] and "y" in config["D0JetRate"]["theory_legend"]:
        y1 = config["D0JetRate"]["theory_legend"]["y"]
    else:
        y1 = y2 - 0.03
    active_t = len([t for t in config["theory"] if "histogram_plot" in t])
    active_t += len([t for t in inclusive_jet_cross_sections.itervalues() if "histogram_plot" in t])
    y2 = y1 - 0.06 * active_t / n_leg_columns
    if "theory_legend" in config["D0JetRate"] and "x" in config["D0JetRate"]["theory_legend"]:
        x1 = config["D0JetRate"]["theory_legend"]["x"]
    else:
        x1 = 0.16
    x2 = 0.85

    leg1 = ROOT.TLegend(x1, y1, x2, y2, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetNColumns(n_leg_columns)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(19)
    leg1.SetTextAlign(12)
    leg1.SetMargin(0.08)
    for t in config["theory"]:
        if "histogram_plot" in t:
            leg1.AddEntry(t["histogram_plot"], "D^{{0}} Jets, {}".format(t["title"]), "l")
    for t in inclusive_jet_cross_sections.itervalues():
        if "histogram_plot" in t: 
            leg1.AddEntry(t["histogram_plot"], "Inclusive Jets, {}".format(t["title"]), "l")
    leg1.Draw()

    if "data_legend" in config["D0JetRate"] and "y" in config["D0JetRate"]["data_legend"]:
        y1 = config["D0JetRate"]["data_legend"]["y"]
    else:
        y1 = y2 - 0.02
    y2 = y1 - 0.12
    if "data_legend" in config["D0JetRate"] and "x" in config["D0JetRate"]["data_legend"]:
        x1 = config["D0JetRate"]["data_legend"]["x"]
    else:
        x1 = 0.16
    x2 = x1 + 0.30
    leg1 = ROOT.TLegend(x1, y1, x2, y2, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(19)
    leg1.SetTextAlign(12)
    leg1.SetMargin(0.2)
    entry = leg1.AddEntry(None, "D^{0} Jets, #it{p}_{T,D} > 3 GeV/#it{c}", "pf")
    entry.SetFillStyle(d0jet_syst_copy.GetFillStyle())
    entry.SetFillColor(d0jet_syst_copy.GetFillColor())
    entry.SetLineColor(d0jet_syst_copy.GetFillColor())
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
    hAxisRatio.GetYaxis().SetTitleOffset(1.3)

    cname = "{}_{}_Ratio".format(config["D0JetRate"]["name_prefix"], config["name"])
    if "canvas_h" in config:
        canvas_h = config["canvas_h"]
    else:
        canvas_h = 600
    if "canvas_w" in config:
        canvas_w = config["canvas_w"]
    else:
        canvas_w = 600
    canvas_ratio = ROOT.TCanvas(cname, cname, canvas_w, canvas_h)
    globalList.append(canvas_ratio)
    canvas_ratio.SetTicks(1, 1)
    canvas_ratio.SetLeftMargin(0.10)
    canvas_ratio.SetTopMargin(0.05)
    canvas_ratio.SetRightMargin(0.05)

    canvas_ratio.cd()
    hAxisRatio.Draw("axis")
    hAxisRatio.GetXaxis().SetTitleOffset(1.2)
    hAxisRatio.GetYaxis().SetRangeUser(0, 0.075)
    globalList.append(hAxisRatio)

    ratioSyst.Draw("2")
    ratioStat_copy = ratioStat.DrawCopy("same p e0 x0")

    line_styles = [2, 4, 7, 3, 5, 6, 8, 9]
    for t, line_style in zip(config["theory"], line_styles):
        if not t["active"]: continue
        if not t["inclusive"]: continue
        hratio = t["ratio_histogram"]
        hratio_copy = hratio.DrawCopy("same")
        hratio_copy.SetLineStyle(line_style)
        hratio_copy.SetLineWidth(3)
        t["ratio_copy_histogram"] = hratio_copy

    # Now plotting labels

    y1 = 0.90
    y2 = y1 - 0.06 * (len(config["D0JetRate"]["title"]) + 1)
    paveALICE = ROOT.TPaveText(0.14, y1, 0.55, y2, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(22)
    paveALICE.SetTextAlign(12)
    for line in config["D0JetRate"]["title"]: 
        paveALICE.AddText(line)
    paveALICE.AddText("with D^{0}, #it{p}_{T,D} > 3 GeV/#it{c}")
    paveALICE.Draw()

    n_leg_columns = 1
    y1 = 0.25
    active_t = len([t for t in config["theory"] if "histogram_plot" in t])
    y2 = y1 - 0.04 * active_t / n_leg_columns
    x1 = 0.35
    x2 = 0.80

    leg1 = ROOT.TLegend(x1, y1, x2, y2, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetNColumns(n_leg_columns)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(22)
    leg1.SetTextAlign(12)
    leg1.SetMargin(0.08)
    for t in config["theory"]:
        if "histogram_plot" in t:
            leg1.AddEntry(t["ratio_copy_histogram"], t["title"], "l")
    leg1.Draw()

    y1 = 0.90
    y2 = y1 - 0.09
    x1 = 0.62
    x2 = x1 + 0.30
    leg1 = ROOT.TLegend(x1, y1, x2, y2, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(22)
    leg1.SetTextAlign(12)
    leg1.SetMargin(0.2)
    entry = leg1.AddEntry(ratioStat, "Data", "pe")
    entry = leg1.AddEntry(ratioSyst, "Systematic Uncertainty", "f")
    leg1.Draw()

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
    canvas, canvas_ratio = PlotCrossSections(d0jet_stat, d0jet_syst, incl_stat, incl_syst, inclusive_jet_cross_sections, config)
    canvas.SaveAs("{}/{}.pdf".format(config["input_path"], canvas.GetName()))
    canvas.SaveAs("{}/{}.C".format(config["input_path"], canvas.GetName()))
    canvas_ratio.SaveAs("{}/{}.pdf".format(config["input_path"], canvas_ratio.GetName()))
    canvas_ratio.SaveAs("{}/{}.C".format(config["input_path"], canvas_ratio.GetName()))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Jet pt spectrum theory comparison.')
    parser.add_argument('yaml', metavar='conf.yaml')

    args = parser.parse_args()

    f = open(args.yaml, 'r')
    yconfig = yaml.load(f)
    f.close()

    main(yconfig)

    IPython.embed()