#!/usr/bin/env python
# python script to compare D0 jet spectrum with theory
# use it with
# JetPtCrossSection_CompareFullVsCharged.yaml
# JetPtCrossSection_CompareLO.yaml
# JetPtCrossSection_CompareMonteCarlo_Paper.yaml
# JetPtCrossSection_ComparePowheg_Paper.yaml
# JetPtCrossSection_ComparePtDCut.yaml
# JetPtCrossSection_CompareSystematicsPowheg_ccbar.yaml

# JetZDistr_JetPt_5_15_CompareMonteCarlo.yaml
# JetZDistr_JetPt_5_15_ComparePowheg.yaml

# JetZDistr_JetPt_15_30_CompareMonteCarlo.yaml
# JetZDistr_JetPt_15_30_ComparePowheg.yaml

# JetZCrossSection_JetPt_5_15_CompareMonteCarlo_Paper.yaml
# JetZCrossSection_JetPt_5_15_ComparePowheg_Paper.yaml

# JetZCrossSection_JetPt_15_30_CompareMonteCarlo_Paper.yaml
# JetZCrossSection_JetPt_15_30_ComparePowheg_Paper.yaml

import argparse
import yaml
import IPython
import ROOT
import DMesonJetUtils
import LoadTheoryCrossSections

globalList = []

def GetMeasuredCrossSection(input_path, file_name, scale):
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
    hStat.Scale(scale)
    for ipoint in range(0, hSyst.GetN()):
        hSyst.SetPointEYlow(ipoint, hSyst.GetErrorYlow(ipoint) * scale)
        hSyst.SetPointEYhigh(ipoint, hSyst.GetErrorYhigh(ipoint) * scale)
        hSyst.SetPoint(ipoint, hSyst.GetX()[ipoint], hSyst.GetY()[ipoint] * scale)
    return hStat, hSyst

def GetTheoryCrossSectionAll(config, axis):
    return LoadTheoryCrossSections.GetD0JetTheoryCrossSectionAll(config, axis)

def PlotCrossSections(dataStat, dataSyst, config):
    no_box_sys = ["[]", "||", "0[]", "0||"]
    if not "data_systematics_style" in config:
        config["data_systematics_style"] = "2"
    config["data_box_systematics"] = not (config["data_systematics_style"] in no_box_sys)

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
    if config["logy"]:
        padMain.SetLogy()
    
    if "minx" in config:
        minx = config["minx"]
    else:
        minx = dataStat.GetXaxis().GetXmin()

    if "maxx" in config:
        maxx = config["maxx"]
    else:
        maxx = dataStat.GetXaxis().GetXmax()

    h = ROOT.TH1I("h", "h", 1000, minx, maxx)
    h.Draw("axis")
    globalList.append(h)
    h.GetYaxis().SetRangeUser(config["miny"], config["maxy"])
    h.GetYaxis().SetTitleFont(43)
    h.GetYaxis().SetTitleSize(config["axis"]["font_size"])
    h.GetYaxis().SetTitleOffset(config["axis"]["y_offset"])
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelSize(config["axis"]["font_size"] - 3)
    h.GetYaxis().SetTitle(config["axis"]["y_title"])

    dataSyst_copy = dataSyst.Clone("{0}_copy".format(dataSyst.GetName()))
    globalList.append(dataSyst_copy)
    dataStat_copy = dataStat.Clone("{0}_copy".format(dataStat.GetName()))
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
            if not t["active"]:
                continue
            t["histogram"].GetXaxis().SetRangeUser(config["data_minx"], config["data_maxx"])
            if "type" in t  and t["type"] == "stat+syst":
                hSyst = t["systematics"]
                for ipoint in range(0, hSyst.GetN()):
                    if hSyst.GetX()[ipoint] < config["data_minx"] or hSyst.GetX()[ipoint] > config["data_maxx"]:
                        points_to_be_removed.add(ipoint)
                for ipoint in points_to_be_removed:
                    hSyst.RemovePoint(ipoint)

    # Plot data systematic uncertainties, if they are plotted as a filled box
    if config["data_box_systematics"]:
        dataSyst_copy.Draw(config["data_systematics_style"])

    # Plot theory systematic uncertainties (if they are plotted as a filled box)
    for t in config["theory"]:
        if not t["active"]:
            continue
        if not "type" in t:
            t["type"] = "stat-only"
        if t["type"] == "stat+syst":
            if not "systematics_style" in t:
                t["systematics_style"] = "2"
            t["box_systematics"] = not (t["systematics_style"] in no_box_sys)

            if t["box_systematics"]:
                hSyst = t["systematics"].Clone("{0}_copy".format(t["systematics"].GetName()))
                hSyst.SetLineColor(t["color"])
                hSyst.SetFillColor(t["color"])
                if "line" in t:
                    hSyst.SetLineStyle(t["line"])
                if "fill" in t:
                    hSyst.SetFillStyle(t["fill"])
                else:
                    hSyst.SetFillStyle(0)
                hSyst.Draw(t["systematics_style"])
                globalList.append(hSyst)
                t["systematics_plot"] = hSyst

    # Plot data statistical uncertainties
    dataStat_copy.Draw("same p e0 x0")

    # Plot theory statistical uncertainties and systematic uncertainties (if they are not plotted as filled boxes)
    for t in config["theory"]:
        if not t["active"]:
            continue
        if t["type"] == "stat-only":
            h = t["histogram"].Clone(t["gen"])
            h.Draw("same e0")
            if "data_minx" in config and "data_maxx" in config:
                h.GetXaxis().SetRangeUser(config["data_minx"], config["data_maxx"])
            globalList.append(h)
            h.SetLineColor(t["color"])
            h.SetLineStyle(t["line"])
            h.SetLineWidth(3)
            h.SetMarkerStyle(1)
            h.SetMarkerColor(t["color"])
            t["histogram_plot"] = h
        elif t["type"] == "stat+syst":
            if not t["box_systematics"]:
                hSyst = t["systematics"].Clone("{0}_copy".format(t["systematics"].GetName()))
                for ipoint in range(0, hSyst.GetN()):
                    hSyst.SetPointEXhigh(ipoint, 0)
                    hSyst.SetPointEXlow(ipoint, 0)
                hSyst.SetLineColor(t["color"])
                hSyst.Draw(t["systematics_style"])
                globalList.append(hSyst)
                t["systematics_plot"] = hSyst

            hStat = t["histogram"].Clone("{0}_copy".format(t["histogram"].GetName()))
            hStat.SetMarkerStyle(getattr(ROOT, t["marker"]))
            hStat.SetLineColor(t["color"])
            hStat.SetMarkerColor(t["color"])
            hStat.Draw("same p e0 x0")
            globalList.append(hStat)
            t["histogram_plot"] = hStat

    # Plot data systematic uncertainties, if they are not plotted as filled boxes
    if not config["data_box_systematics"]:
        for ipoint in range(0, dataSyst_copy.GetN()):
            dataSyst_copy.SetPointEXhigh(ipoint, 0)
            dataSyst_copy.SetPointEXlow(ipoint, 0)
        dataSyst_copy.SetLineColor(dataStat_copy.GetLineColor())
        dataSyst_copy.Draw(config["data_systematics_style"])

    # RATIO

    ratioSyst = dataSyst_copy.Clone("ratioSyst")
    globalList.append(ratioSyst)

    ratioStat = dataStat_copy.Clone("ratioStat")
    ratioStat.SetFillStyle(0)
    globalList.append(ratioStat)

    for ibin in range(0, ratioSyst.GetN()):
        ratioSyst.SetPointEYlow(ibin, ratioSyst.GetErrorYlow(ibin) / ratioSyst.GetY()[ibin])
        ratioSyst.SetPointEYhigh(ibin , ratioSyst.GetErrorYhigh(ibin) / ratioSyst.GetY()[ibin])
        ratioSyst.SetPoint(ibin, ratioSyst.GetX()[ibin], 1.0)
    
    for ibin in range(1, ratioStat.GetNbinsX()+1):
        if ratioStat.GetBinContent(ibin) == 0:
            continue
        ratioStat.SetBinError(ibin, ratioStat.GetBinError(ibin) / ratioStat.GetBinContent(ibin))
        ratioStat.SetBinContent(ibin, 1.0)

    for t in config["theory"]:
        if not "histogram_plot" in t:
            continue
        r = t["histogram_plot"].Clone()
        globalList.append(r)
        for ibin in range(1, r.GetNbinsX() + 1):
            if dataStat_copy.GetBinContent(ibin) <= 0: 
                r.SetBinError(ibin, 0)
                r.SetBinContent(ibin, 0)
                continue
            r.SetBinError(ibin, t["histogram_plot"].GetBinError(ibin) / dataStat_copy.GetBinContent(ibin))
            r.SetBinContent(ibin, t["histogram_plot"].GetBinContent(ibin) / dataStat_copy.GetBinContent(ibin))
        t["ratio_histogram_plot"] = r
        if t["type"] == "stat+syst" and "systematics_plot" in t:
            rSyst = t["systematics_plot"].Clone()
            globalList.append(rSyst)
            for ibin in range(0, rSyst.GetN()):
                rSyst.SetPointEYlow(ibin, rSyst.GetErrorYlow(ibin) / dataSyst_copy.GetY()[ibin])
                rSyst.SetPointEYhigh(ibin , rSyst.GetErrorYhigh(ibin) / dataSyst_copy.GetY()[ibin])
                rSyst.SetPoint(ibin, rSyst.GetX()[ibin], rSyst.GetY()[ibin] / dataSyst_copy.GetY()[ibin])
            t["ratio_systematics_plot"] = rSyst

    padRatio.cd()

    hRatio = ROOT.TH1I("h", "h", 1000, minx, maxx)
    hRatio.Draw("axis")
    globalList.append(hRatio)
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

    # Plot data systematic uncertainties, if they are plotted as a filled box
    if config["data_box_systematics"]:
        ratioSyst.Draw(config["data_systematics_style"])

    # Plot theory systematic uncertainties (if they are plotted as a filled box)
    for t in config["theory"]:
        if not t["active"]:
            continue
        if t["type"] == "stat+syst" and t["box_systematics"]:
            if not "ratio_systematics_plot" in t:
                continue
            rSyst = t["ratio_systematics_plot"]
            rSyst.Draw(t["systematics_style"])

    # Plot data statistical uncertainties
    ratioStat.Draw("same p e0 x0")

    # Plot theory statistical uncertainties and systematic uncertainties (if they are not plotted as filled boxes)
    for t in config["theory"]:
        if not t["active"]:
            continue
        r = t["ratio_histogram_plot"]
        if t["type"] == "stat-only":
            r.SetLineWidth(3)
            r.Draw("same e0")
        elif t["type"] == "stat+syst":
            r.Draw("same p e0 x0")
            if not t["box_systematics"] and "ratio_systematics_plot" in t:
                rSyst = t["ratio_systematics_plot"]
                for ipoint in range(0, rSyst.GetN()):
                    rSyst.SetPointEXhigh(ipoint, 0)
                    rSyst.SetPointEXlow(ipoint, 0)
                rSyst.Draw(t["systematics_style"])

    # Plot data systematic uncertainties, if they are not plotted as filled boxes
    if not config["data_box_systematics"]:
        for ipoint in range(0, ratioSyst.GetN()):
            ratioSyst.SetPointEXhigh(ipoint, 0)
            ratioSyst.SetPointEXlow(ipoint, 0)
        ratioSyst.SetLineColor(ratioStat.GetLineColor())
        ratioSyst.Draw(config["data_systematics_style"])

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
    y2 = y1 - 0.06
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
        entry.SetFillStyle(dataSyst_copy.GetFillStyle())
        entry.SetFillColor(dataSyst_copy.GetFillColor())
        entry.SetLineColor(dataSyst_copy.GetFillColor())
        entry.SetMarkerColor(dataStat_copy.GetMarkerColor())
        entry.SetMarkerStyle(dataStat_copy.GetMarkerStyle())
    else:
        entry = leg1.AddEntry(None, "Data", "p")
        entry.SetMarkerColor(dataStat_copy.GetMarkerColor())
        entry.SetMarkerStyle(dataStat_copy.GetMarkerStyle())

    leg1.Draw()

    padRatio.RedrawAxis("g")
    padRatio.RedrawAxis()
    padMain.RedrawAxis()

    return canvas

def main(config):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    if "scale" in config:
        scale = config["scale"]
    else:
        scale = 1.0
    dataStat, dataSyst = GetMeasuredCrossSection(config["input_path"], config["data"], scale)
    GetTheoryCrossSectionAll(config, dataStat.GetXaxis())
    canvas = PlotCrossSections(dataStat, dataSyst, config)
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
