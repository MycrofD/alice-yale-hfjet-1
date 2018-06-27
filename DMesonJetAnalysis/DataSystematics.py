#!/usr/bin/env python
# python script to calculate systematic uncertainties

import argparse
import IPython
import ROOT
import array
import numpy
import math
import yaml
import DMesonJetUnfolding
import DMesonJetUtils
import DMesonJetCompare
import copy
from collections import OrderedDict
import HistogramNormalizator

globalList = []

xaxis_title = ""
yaxis_title = ""

def main(config):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    global xaxis_title
    global yaxis_title

    if "JetZ" in config["name"]:
        xaxis_title = "#it{z}_{||}^{ch}"
        if config["normalization"] == "cross_section":
            yaxis_title = "#frac{d^{2}#sigma}{d#it{z}_{||}^{ch}d#it{#eta}} (mb)"
        elif config["normalization"] == "distribution":
            yaxis_title = "Probability Density"
        else:
            print("Normalization option '{}' invalid".format(config["normalization"]))
            exit(1)
    elif "JetPtSpectrum" in config["name"]:
        xaxis_title = "#it{p}_{T,jet}^{ch} (GeV/#it{c})"
        if config["normalization"] == "cross_section":
            yaxis_title = "#frac{d^{2}#sigma}{d#it{p}_{T,jet}^{ch}d#it{#eta}} [mb (GeV/#it{c})^{-1}]"
        elif config["normalization"] == "distribution":
            yaxis_title = "Probability Density (GeV/#it{c})^{-1}"
        else:
            print("Normalization option '{}' invalid".format(config["normalization"]))
            exit(1)

    results = Start(config)

    PlotSystematicUncertaintySummary(config["name"], results)

    output_file_name = "{0}/{1}.root".format(config["input_path"], config["name"])
    print("Storing results in {0}".format(output_file_name))
    outputFile = ROOT.TFile(output_file_name, "recreate")
    if not outputFile or outputFile.IsZombie():
       print("Could not open output file {0}".format(output_file_name))
       exit(1)
    for name, r in results.iteritems():
       rlist = GenerateRootList(r, name)
       rlist.Write(rlist.GetName(), ROOT.TObject.kSingleKey)
    outputFile.Close()
    globalList.append(rlist)

    SaveCanvases(config["input_path"])

def Start(config):
    histograms = LoadHistograms(config)
    results = dict()
    results["Variations"] = CompareVariations(config, histograms)
    results["Uncertainties"] = GenerateUncertainties(config, histograms)
    results["FinalSpectrum"] = PlotSpectrumStatAndSyst(config["name"], results)
    return results

def LoadHistograms(config):
    files = OpenFiles(config)
    events = LoadEvents(files, config)

    histograms = dict()

    h = DMesonJetUtils.GetObject(files[config["default"]["input_name"]], config["default"]["histogram_name"])
    if not h: 
        exit(1)
    h.GetXaxis().SetTitle(xaxis_title)
    h.GetYaxis().SetTitle(yaxis_title)
    normalizator = HistogramNormalizator.DataNormalizator(h, config["normalization"], events[config["default"]["input_name"]])
    normalizator.fXmin = config["xmin"]
    normalizator.fXmax = config["xmax"]
    h = normalizator.NormalizeHistogram()
    histograms["default"] = h
    v_types = ["variations", "up_variations", "low_variations"]
    for s in config["sources"]:
        if not s["active"]: continue
        histograms[s["name"]] = dict()
        for v_type in v_types:
            if v_type in s:
                histograms[s["name"]][v_type] = []
                for v in s[v_type]:
                    h = DMesonJetUtils.GetObject(files[v["input_name"]], v["histogram_name"])
                    if not h:
                        exit(1)
                    normalizator = HistogramNormalizator.DataNormalizator(h, config["normalization"], events[v["input_name"]])
                    normalizator.fXmin = config["xmin"]
                    normalizator.fXmax = config["xmax"]
                    h_copy = normalizator.NormalizeHistogram()
                    h_copy.SetTitle(v["histogram_title"])
                    h_copy.GetXaxis().SetTitle(xaxis_title)
                    h_copy.GetYaxis().SetTitle(yaxis_title)
                    histograms[s["name"]][v_type].append(h_copy)
    return histograms

def LoadEvents(files, config):
    events = dict()
    for name, f in files.iteritems():
        slash = config["default"]["histogram_name"].find("/")
        heventsName = config["default"]["histogram_name"][:slash + 1]
        heventsName += "Events"
        hevents = DMesonJetUtils.GetObject(files[name], heventsName)
        events[name] = hevents.GetBinContent(1)
    return events

def OpenFiles(config):
    files = dict()

    fname = "{0}/{1}/{1}.root".format(config["input_path"], config["default"]["input_name"])
    f = ROOT.TFile(fname)
    if not f or f.IsZombie():
        print("Could not open file '{0}'".format(fname))
        exit(1)
    files[config["default"]["input_name"]] = f

    for s in config["sources"]:
        if not s["active"]: continue
        input_names = LoadInputNames(s)
        for input_name in input_names:
            if not input_name in files:
                fname = "{0}/{1}/{1}.root".format(config["input_path"], input_name)
                f = ROOT.TFile(fname)
                if not f or f.IsZombie():
                    print("Could not open file '{0}'".format(fname))
                    exit(1)
                files[input_name] = f
    return files

def LoadInputNames(s):
    input_names = []
    v_types = ["variations", "up_variations", "low_variations"]
    for v_type in v_types:
        if v_type in s:
            for v in s[v_type]: input_names.append(v["input_name"])
    return input_names


def CompareVariations(config, histograms):
    result = dict()
    variations = []
    h = histograms["default"]
    baseline = h.Clone("default")
    baseline.SetTitle(config["default"]["title"])

    baseline.GetXaxis().SetTitle(xaxis_title)
    baseline.GetYaxis().SetTitle(yaxis_title)
    result[baseline.GetName()] = baseline
    v_types = ["variations", "up_variations", "low_variations"]
    for s in config["sources"]:
        if not s["active"]: continue
        for v_type in v_types:
            if v_type in histograms[s["name"]]:
                for h in histograms[s["name"]][v_type]:
                    result[h.GetName()] = h
                    variations.append(h)

    globalList.extend(variations)
    globalList.append(baseline)
    comp = DMesonJetCompare.DMesonJetCompare("CompareVariations_{0}".format(config["name"]))
    comp.fOptRatio = "hist"
    comp.fX1LegRatio = 0.15
    comp.fX1LegSpectrum = 0.45
    comp.fLegLineHeight = 0.05
    comp.fGridyRatio = True
    if "JetPtCrossSection" in config["name"]:
        comp.fDoSpectraPlot = "logy"
        comp.fLogUpperSpace = 5
    elif "JetZ" in config["name"]:
        comp.fDoSpectraPlot = "lineary"
        comp.fLinUpperSpace = 1.1
    r = comp.CompareSpectra(baseline, variations)
    if "JetPtCrossSection" in config["name"]:
        comp.fMainHistogram.SetMinimum(1e-6)
    elif "JetPt_15_30" in config["name"]:
        comp.fMainHistogram.SetMinimum(2e-4)
    for obj in r:
        globalList.append(obj)
    return result

def CalculateFixSystematicUncertainty(config):
    fixed_unc2 = 0
    print("Source & Uncertainty (\\%) \\\\ \\hline")
    for u in config["fixed_uncertainties"]:
        if u["plot"]: continue
        fixed_unc2 += u["uncertainty"] ** 2
        print("{0} & {1:.1f} \\\\".format(u["title"], u["uncertainty"] * 100))
    fixed_unc = math.sqrt(fixed_unc2)
    print("\\hline")
    print("{0} & {1:.1f}".format("Total", fixed_unc * 100))
    return fixed_unc

def GenerateUncertainties(config, histograms):
    baseline = histograms["default"]
    fixed_syst_unc = CalculateFixSystematicUncertainty(config)
    partialRelSystUncUp = [DMesonJetUtils.soft_clone(baseline, "{0}_up".format(s["name"]), s["title"], "relative uncertainty") for s in config["sources"] if s["active"]] + [DMesonJetUtils.soft_clone(baseline, "{0}_up".format(s["name"]), s["title"], "relative uncertainty") for s in config["fixed_uncertainties"] if s["plot"]]
    partialRelSystUncLow = [DMesonJetUtils.soft_clone(baseline, "{0}_low".format(s["name"]), s["title"], "relative uncertainty") if not s["symmetrize"] else None for s in config["sources"] if s["active"]] + [None for s in config["fixed_uncertainties"] if s["plot"]]
    partialRelSystUncUp_smooth = [DMesonJetUtils.soft_clone(baseline, "{0}_up".format(s["name"]), s["title"], "relative uncertainty") for s in config["sources"] if s["active"]] + [DMesonJetUtils.soft_clone(baseline, "{0}_up".format(s["name"]), s["title"], "relative uncertainty") for s in config["fixed_uncertainties"] if s["plot"]]
    partialRelSystUncLow_smooth = [DMesonJetUtils.soft_clone(baseline, "{0}_low".format(s["name"]), s["title"], "relative uncertainty") if not s["symmetrize"] else None for s in config["sources"] if s["active"]] + [None for s in config["fixed_uncertainties"] if s["plot"]]
    totRelSystUncUp = DMesonJetUtils.soft_clone(baseline, "tot_rel_syst_unc_up", "Total Systematic", "relative uncertainty")
    totRelSystUncLow = DMesonJetUtils.soft_clone(baseline, "tot_rel_syst_unc_low", "Total Systematic", "relative uncertainty")
    statUnc = DMesonJetUtils.soft_clone(baseline, "stat_unc", "Statistical", "relative uncertainty")
    totUncUp = DMesonJetUtils.soft_clone(baseline, "tot_unc_up", "Total", "relative uncertainty")
    totUncLow = DMesonJetUtils.soft_clone(baseline, "tot_unc_low", "Total", "relative uncertainty")
    yerrup = []
    yerrlow = []
    y = []
    xerrup = []
    xerrlow = []
    x = []

    for ibin in range(1, baseline.GetNbinsX() + 1):
        if baseline.GetBinContent(ibin) == 0:
            continue
        ivar = 0
        for s in config["sources"]:
            if not s["active"]: continue
            low_unc2 = 0
            up_unc2 = 0
            sym_unc2 = 0

            for h in histograms[s["name"]]["variations"]:
                diff = baseline.GetBinContent(ibin) - h.GetBinContent(ibin)
                if s["combine_strategy"] == "sum_in_quadrature":
                    if s["symmetrize"]:
                        sym_unc2 += diff ** 2
                    else:
                        if diff > 0:
                            up_unc2 += diff ** 2
                        else:
                            low_unc2 += diff ** 2
                elif s["combine_strategy"] == "envelope":
                    if s["symmetrize"]:
                        sym_unc2 = max(diff ** 2, sym_unc2)
                    else:
                        if diff > 0:
                            up_unc2 = max(diff ** 2, up_unc2)
                        else:
                            low_unc2 = max(diff ** 2, low_unc2)
                else:
                    print("combine strategy {} unknown".format(s["combine_strategy"]))
                    exit(1)

            if s["symmetrize"]:
                part_unc = part_unc_up = part_unc_low = math.sqrt(sym_unc2)
                partialRelSystUncUp[ivar].SetBinContent(ibin, part_unc / baseline.GetBinContent(ibin))
            else:
                part_unc_up2 = up_unc2
                part_unc_low2 = low_unc2
                part_unc_up = math.sqrt(part_unc_up2)
                part_unc_low = math.sqrt(part_unc_low2)
                partialRelSystUncUp[ivar].SetBinContent(ibin, part_unc_up / baseline.GetBinContent(ibin))
                partialRelSystUncLow[ivar].SetBinContent(ibin, part_unc_low / baseline.GetBinContent(ibin))
            ivar += 1
        for s in config["fixed_uncertainties"]:
            if not s["plot"]: continue
            partialRelSystUncUp[ivar].SetBinContent(ibin, s["uncertainty"])
            ivar += 1

    # smothening procedure
    for ibin in range(1, baseline.GetNbinsX() + 1):
        tot_syst_unc_up2 = 0
        tot_syst_unc_low2 = 0
        ivar = 0
        for s in config["sources"]:
            if not s["active"]: continue
            if partialRelSystUncUp[ivar]:
                n = 1
                part_unc_up = partialRelSystUncUp[ivar].GetBinContent(ibin)
                if ibin > 1:
                    part_unc_up += partialRelSystUncUp[ivar].GetBinContent(ibin - 1)
                    n += 1
                if ibin < baseline.GetNbinsX():
                    part_unc_up += partialRelSystUncUp[ivar].GetBinContent(ibin + 1)
                    n += 1
                part_unc_up /= n
                partialRelSystUncUp_smooth[ivar].SetBinContent(ibin, part_unc_up)
            if partialRelSystUncLow[ivar]:
                n = 1
                part_unc_low = partialRelSystUncLow[ivar].GetBinContent(ibin)
                if ibin > 1:
                    part_unc_low += partialRelSystUncLow[ivar].GetBinContent(ibin - 1)
                    n += 1
                if ibin < baseline.GetNbinsX():
                    part_unc_low += partialRelSystUncLow[ivar].GetBinContent(ibin + 1)
                    n += 1
                part_unc_low /= n
                partialRelSystUncLow_smooth[ivar].SetBinContent(ibin, part_unc_low)
            else:
                part_unc_low = part_unc_up
            ivar += 1
        for s in config["fixed_uncertainties"]:
            if not s["plot"]: continue
            partialRelSystUncUp_smooth[ivar].SetBinContent(ibin, partialRelSystUncUp[ivar].GetBinContent(ibin))
            ivar += 1

    fixed_unc2 = fixed_syst_unc ** 2
    for ibin in range(1, baseline.GetNbinsX() + 1):
        x.append(baseline.GetXaxis().GetBinCenter(ibin))
        xerrup.append(baseline.GetXaxis().GetBinWidth(ibin) / 2)
        xerrlow.append(baseline.GetXaxis().GetBinWidth(ibin) / 2)
        if baseline.GetBinContent(ibin) == 0:
            y.append(0)
            yerrup.append(0)
            yerrlow.append(0)
            continue
        tot_syst_unc_up2 = 0
        tot_syst_unc_low2 = 0
        ivar = 0
        for s in config["sources"]:
            if not s["active"]: continue
            if partialRelSystUncUp_smooth[ivar]:
                part_unc_up = partialRelSystUncUp_smooth[ivar].GetBinContent(ibin)
                if ibin < baseline.GetNbinsX():
                    if part_unc_up > partialRelSystUncUp_smooth[ivar].GetBinContent(ibin + 1): part_unc_up = math.floor(part_unc_up * 100) / 100
                    else: part_unc_up = math.floor(part_unc_up * 100 + 0.5) / 100
                else:
                    if part_unc_up < partialRelSystUncUp_smooth[ivar].GetBinContent(ibin - 1): part_unc_up = math.ceil(part_unc_up * 100) / 100
                    else: part_unc_up = math.floor(part_unc_up * 100 + 0.5) / 100
                partialRelSystUncUp_smooth[ivar].SetBinContent(ibin, part_unc_up)
            if partialRelSystUncLow_smooth[ivar]:
                part_unc_low = partialRelSystUncLow_smooth[ivar].GetBinContent(ibin)
                if ibin < baseline.GetNbinsX():
                    if part_unc_low > partialRelSystUncLow_smooth[ivar].GetBinContent(ibin + 1): part_unc_low = math.floor(part_unc_low * 100) / 100
                    else: part_unc_low = math.floor(part_unc_low * 100 + 0.5) / 100
                else:
                    if part_unc_low < partialRelSystUncLow_smooth[ivar].GetBinContent(ibin - 1): part_unc_low = math.ceil(part_unc_low * 100) / 100
                    else: part_unc_low = math.floor(part_unc_low * 100 + 0.5) / 100
                partialRelSystUncLow_smooth[ivar].SetBinContent(ibin, part_unc_low)
            else:
                part_unc_low = part_unc_up
            tot_syst_unc_up2 += part_unc_up ** 2
            tot_syst_unc_low2 += part_unc_low ** 2
            ivar += 1
        for s in config["fixed_uncertainties"]:
            if not s["plot"]: continue
            part_unc_up = partialRelSystUncUp_smooth[ivar].GetBinContent(ibin)
            part_unc_low = part_unc_up
            tot_syst_unc_up2 += part_unc_up ** 2
            tot_syst_unc_low2 += part_unc_low ** 2
            ivar += 1
        tot_syst_unc_up2 += fixed_unc2
        tot_syst_unc_low2 += fixed_unc2
        tot_syst_unc_up = math.sqrt(tot_syst_unc_up2)
        tot_syst_unc_low = math.sqrt(tot_syst_unc_low2)
        stat_unc = baseline.GetBinError(ibin) / baseline.GetBinContent(ibin)
        stat_unc2 = stat_unc ** 2
        tot_unc_up = math.sqrt(tot_syst_unc_up2 + stat_unc2)
        tot_unc_low = math.sqrt(tot_syst_unc_low2 + stat_unc2)
        totRelSystUncUp.SetBinContent(ibin, tot_syst_unc_up)
        totRelSystUncLow.SetBinContent(ibin, tot_syst_unc_low)
        statUnc.SetBinContent(ibin, stat_unc)
        totUncUp.SetBinContent(ibin, tot_unc_up)
        totUncLow.SetBinContent(ibin, tot_unc_low)

        y.append(baseline.GetBinContent(ibin))
        yerrup.append(tot_syst_unc_up * baseline.GetBinContent(ibin))
        yerrlow.append(tot_syst_unc_low * baseline.GetBinContent(ibin))

    centralSystUnc = ROOT.TGraphAsymmErrors(baseline.GetNbinsX(),
                                            numpy.array(x, dtype=float), numpy.array(y, dtype=float),
                                            numpy.array(xerrlow, dtype=float), numpy.array(xerrup, dtype=float),
                                            numpy.array(yerrlow, dtype=float), numpy.array(yerrup, dtype=float))
    centralSystUnc.SetName("central_syst_unc")
    centralSystUnc.SetTitle("{0} Systematics".format(baseline.GetTitle()))
    centralSystUnc.GetXaxis().SetTitle(xaxis_title)
    centralSystUnc.GetYaxis().SetTitle(yaxis_title)
    result = {"PartialSystematicUncertaintiesUp" : partialRelSystUncUp_smooth, "PartialSystematicUncertaintiesLow" : partialRelSystUncLow_smooth,
              "fixed_syst_unc" : fixed_syst_unc,
              totRelSystUncUp.GetName() : totRelSystUncUp, totRelSystUncLow.GetName() : totRelSystUncLow,
              statUnc.GetName() : statUnc, totUncUp.GetName() : totUncUp, totUncLow.GetName() : totUncLow,
              centralSystUnc.GetName() : centralSystUnc}
    return result

def PlotSystematicUncertaintySummary(name, results):
    sourcesUp = []
    sourcesLow = []
    colorsUp = []
    colorsLow = []
    h = results["Uncertainties"]["tot_unc_up"]
    tot_unc_up = h.Clone("{0}_copy".format(h.GetName()))
    sourcesUp.append(tot_unc_up)
    colorsUp.append(ROOT.kBlack)

    print("Source & \\multicolumn{{{0}}}{{c}}{{Uncertainty (\\%)}} \\\\ \\hline".format(tot_unc_up.GetNbinsX()))
    if "JetZDistr" in name or "JetZCrossSection" in name:
        print(" & ".join(["\\zpar\\"] + ["{0:.1f} - {1:.1f}".format(tot_unc_up.GetXaxis().GetBinLowEdge(ibin), tot_unc_up.GetXaxis().GetBinUpEdge(ibin)) for ibin in range(1, tot_unc_up.GetNbinsX() + 1)]) + "\\\\ \hline")
    elif "JetPtCrossSection" in name:
        print(" & ".join(["\\ptchjet\\ (\\GeVc)"] + ["{0:.0f} - {1:.0f}".format(tot_unc_up.GetXaxis().GetBinLowEdge(ibin), tot_unc_up.GetXaxis().GetBinUpEdge(ibin)) for ibin in range(1, tot_unc_up.GetNbinsX() + 1)]) + "\\\\ \hline")

    h = results["Uncertainties"]["tot_unc_low"]
    tot_unc_low = h.Clone("{0}_copy".format(h.GetName()))
    sourcesLow.append(tot_unc_low)
    colorsLow.append(ROOT.kBlack)

    h = results["Uncertainties"]["stat_unc"]
    stat_unc = h.Clone("{0}_copy".format(h.GetName()))
    sourcesUp.append(stat_unc)
    colorsUp.append(ROOT.kBlue + 2)

    h = results["Uncertainties"]["tot_rel_syst_unc_up"]
    tot_rel_syst_unc_up = h.Clone("{0}_copy".format(h.GetName()))
    sourcesUp.append(tot_rel_syst_unc_up)
    colorsUp.append(ROOT.kRed + 2)

    h = results["Uncertainties"]["tot_rel_syst_unc_low"]
    tot_rel_syst_unc_low = h.Clone("{0}_copy".format(h.GetName()))
    sourcesLow.append(tot_rel_syst_unc_low)
    colorsLow.append(ROOT.kRed + 2)

    colorsPart = [ROOT.kGreen + 2, ROOT.kOrange + 2, ROOT.kAzure + 2, ROOT.kMagenta + 2, ROOT.kCyan + 2, ROOT.kPink + 1, ROOT.kTeal + 2, ROOT.kYellow + 2]
    cont = 0
    for hUp, hLow in zip(results["Uncertainties"]["PartialSystematicUncertaintiesUp"], results["Uncertainties"]["PartialSystematicUncertaintiesLow"]):
        if hLow:
            print(" & ".join(["\multirow{{2}}{{*}}{{{}}}".format(hUp.GetTitle())] + ["+{0:.0f}".format(hUp.GetBinContent(ibin) * 100) for ibin in range(1, h.GetNbinsX() + 1)]) + "\\\\")
            print(" & ".join([" "] + ["-{0:.0f}".format(hLow.GetBinContent(ibin) * 100) for ibin in range(1, hLow.GetNbinsX() + 1)]) + "\\\\")
        else:
            print(" & ".join([hUp.GetTitle()] + ["{0:.0f}".format(hUp.GetBinContent(ibin) * 100) for ibin in range(1, h.GetNbinsX() + 1)]) + "\\\\")
        h_copy = hUp.Clone("{0}_copy".format(hUp.GetName()))
        sourcesUp.append(h_copy)
        colorsUp.append(colorsPart[cont])

        if hLow:
            h_copy = hLow.Clone("{0}_copy".format(hLow.GetName()))
            sourcesLow.append(h_copy)
            colorsLow.append(colorsPart[cont])

        cont += 1

    print("\\hline")
    print("Normalization (BR \& lumi) & \\multicolumn{{{0}}}{{c}}{{{1:.1f}}} \\\\".format(tot_unc_up.GetNbinsX(), results["Uncertainties"]["fixed_syst_unc"] * 100))
    print("\\hline")
    print(" & ".join(["\multirow{2}{*}{Total Systematic Uncertainty}"] + ["+{0:.1f}".format(tot_rel_syst_unc_up.GetBinContent(ibin) * 100) for ibin in range(1, tot_rel_syst_unc_up.GetNbinsX() + 1)]) + "\\\\")
    print(" & ".join([" "] + ["-{0:.1f}".format(tot_rel_syst_unc_low.GetBinContent(ibin) * 100) for ibin in range(1, tot_rel_syst_unc_low.GetNbinsX() + 1)]) + "\\\\")
    print("\\hline")
    print(" & ".join([stat_unc.GetTitle()] + ["{0:.1f}".format(stat_unc.GetBinContent(ibin) * 100) for ibin in range(1, stat_unc.GetNbinsX() + 1)]) + "\\\\")
    print("\\hline")
    print(" & ".join(["\multirow{2}{*}{Total Uncertainty}"] + ["+{0:.1f}".format(tot_unc_up.GetBinContent(ibin) * 100) for ibin in range(1, tot_unc_up.GetNbinsX() + 1)]) + "\\\\")
    print(" & ".join([" "] + ["-{0:.1f}".format(tot_unc_low.GetBinContent(ibin) * 100) for ibin in range(1, tot_unc_low.GetNbinsX() + 1)]) + "\\\\")

    globalList.extend(sourcesUp)
    globalList.extend(sourcesLow)
    comp = DMesonJetCompare.DMesonJetCompare("CompareUncertainties_{0}".format(name))
    comp.fOptSpectrum = "hist"
    comp.fOptSpectrumBaseline = "hist"
    comp.fX1LegSpectrum = 0.15
    comp.fX2LegSpectrum = 0.55
    comp.fLegLineHeight = 0.05
    comp.fLinUpperSpace = 1.5
    comp.fDoSpectraPlot = "lineary"
    comp.fDoRatioPlot = False
    comp.fColors = colorsUp
    comp.fLines = [1] * len(colorsUp)
    r = comp.CompareSpectra(sourcesUp[0], sourcesUp[1:])
    for obj in r:
        globalList.append(obj)
    comp.fOptSpectrumBaseline = "hist same"
    comp.fColors = colorsLow
    comp.fLines = [2] * len(colorsLow)
    comp.fDoSpectrumLegend = False
    r = comp.CompareSpectra(sourcesLow[0], sourcesLow[1:])
    for obj in r:
        globalList.append(obj)

def PlotSpectrumStatAndSyst(name, results):
    stat = results["Variations"]["default"]
    syst = results["Uncertainties"]["central_syst_unc"]
    canvas = ROOT.TCanvas(name, name)
    canvas.SetLeftMargin(0.13)
    canvas.cd()
    h = stat.DrawCopy("axis")
    h.GetYaxis().SetTitleOffset(1.5)

    h.GetXaxis().SetTitle(xaxis_title)
    h.GetYaxis().SetTitle(yaxis_title)
    if "JetPtCrossSection" in name:
        h.GetYaxis().SetRangeUser(5e-6, 4e-2)
        canvas.SetLogy()
    elif "JetZ" in name:
        h.GetXaxis().SetRangeUser(0.2, 1.0)
        if "JetPt_15_30" in name:
            if "CrossSection" in name:
                h.GetYaxis().SetRangeUser(-0.001, 0.016)
            elif "Distr" in name:
                h.GetYaxis().SetRangeUser(0, 5.0)
        else:
            if "CrossSection" in name:
                h.GetYaxis().SetRangeUser(-0.02, 0.36)
            elif "Distr" in name:
                h.GetYaxis().SetRangeUser(0, 3.5)

    syst_copy = syst.Clone("central_syst_unc_copy")
    syst_copy.Draw("e2")

    syst_copy.SetFillColor(ROOT.kGray)
    syst_copy.SetMarkerColor(ROOT.kGray)
    syst_copy.SetLineColor(ROOT.kGray)
    stat_copy = stat.DrawCopy("same p e0 x0")
    stat_copy.SetMarkerColor(ROOT.kRed + 2)
    stat_copy.SetLineColor(ROOT.kRed + 2)
    stat_copy.SetMarkerStyle(ROOT.kFullCircle)
    stat_copy.SetMarkerSize(0.9)
    leg = ROOT.TLegend(0.35, 0.71, 0.89, 0.89, "NB")
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(43)
    leg.SetTextSize(20)
    leg.AddEntry(stat_copy, "ALICE", "pe")
    leg.AddEntry(syst_copy, "Systematic Uncertainty", "f")
    leg.Draw()
    globalList.append(h)
    globalList.append(syst_copy)
    globalList.append(stat_copy)
    globalList.append(leg)
    globalList.append(canvas)
    result = dict()
    syst_copy_copy = syst_copy.Clone("CentralPointsSystematicUncertainty")
    stat_copy_copy = stat_copy.Clone("CentralPointsStatisticalUncertainty")
    result["CentralPointsSystematicUncertainty"] = syst_copy_copy
    result["CentralPointsStatisticalUncertainty"] = stat_copy_copy
    result["FinalSpectrumCanvas"] = canvas
    return result

def SaveCanvases(input_path):
    for obj in globalList:
        if isinstance(obj, ROOT.TCanvas):
            oname = obj.GetName().replace("/", "_")
            obj.SaveAs("{0}/{1}.pdf".format(input_path, oname))

def GenerateRootList(pdict, name):
    rlist = ROOT.TList()
    rlist.SetName(name)
    print("Loop in {0}".format(name))
    if isinstance(pdict, dict) or isinstance(pdict, OrderedDict):
        objects = pdict.iteritems()
    else:
        objects = zip([inner_obj.GetName() for inner_obj in pdict if inner_obj], pdict)
    for nobj, obj in objects:
        if isinstance(obj, ROOT.TObject):
            print("Adding in {0}".format(nobj))
            rlist.Add(obj)
        elif isinstance(obj, dict) or isinstance(obj, OrderedDict) or isinstance(obj, list):
            print("Recursion in {0}".format(nobj))
            rlist.Add(GenerateRootList(obj, nobj))
        else:
            print("Error: type of object {0} not recognized!".format(obj))
    return rlist

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Data systematics.')
    parser.add_argument('yaml', metavar='conf.yaml')

    args = parser.parse_args()

    f = open(args.yaml, 'r')
    __my_config__ = yaml.load(f)
    f.close()

    main(__my_config__)

    IPython.embed()
