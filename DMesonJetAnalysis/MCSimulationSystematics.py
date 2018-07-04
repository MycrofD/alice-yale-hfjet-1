#!/usr/bin/env python
# python script to calculate systematics on MC simulations

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
import DetectorResponseLoader
import copy
import os
from collections import OrderedDict
import HistogramNormalizator

globalList = []

def main(config, unfolding_debug):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gSystem.Load("libRooUnfold")

    if "b_response" in config:
        print("Opening the b->D0 response matrix file")
        bResponseFile = OpenResponseFile(config["input_path"], config["b_response"], False)

        print("Opening the b->D0 response matrix file (w/ efficiency)")
        bResponseFile_efficiency = OpenResponseFile(config["input_path"], config["b_response"], True)
    else:
        print("No b->D0 response provided.")
        bResponseFile = None
        bResponseFile_efficiency = None

    if "c_response" in config:
        if config["c_response"]["train"] == config["b_response"]["train"] and config["c_response"]["name"] == config["b_response"]["name"]:
            cResponseFile = bResponseFile
            cResponseFile_efficiency = bResponseFile_efficiency
        else:
            print("Opening the c->D0 response matrix file")
            cResponseFile = OpenResponseFile(config["input_path"], config["c_response"], False)

            print("Opening the c->D0 response matrix file (w/ efficiency)")
            cResponseFile_efficiency = OpenResponseFile(config["input_path"], config["c_response"], True)
    else:
        print("No c->D0 response provided.")
        cResponseFile = None
        cResponseFile_efficiency = None

    if bResponseFile_efficiency and cResponseFile_efficiency:
        do_unfolding = True
    else:
        do_unfolding = False

    results = OrderedDict()
    robjects = []
    for v in config["variations"]:
        if not v["active"]: continue
        if "scaling" in v:
            scaling_factor = v["scaling"]
        else:
            scaling_factor = 1
        name = v["name"]
        if "generator" in v:
            suffix = "_".join([v["generator"], str(v["ts"])])
        else:
            suffix = "_".join([config["generator"], str(v["ts"])])

        input_file_name = "{input_path}/FastSim_{suffix}/{name}_{suffix}.root".format(input_path=config["input_path"], suffix=suffix, name=config["analysis_name"])
        fd_histograms = LoadFDHistogram(input_file_name, config["spectra"], scaling_factor)
        if "inclusive" in config:
            input_inclusive_file_name = "{input_path}/FastSim_{generator}_{ts}/{name}.root".format(input_path=config["input_path"], generator=config["inclusive"]["generator"], ts=v["inclusive_ts"], name=config["inclusive"]["analysis_name"])
            inclusive_histograms = LoadInclusiveHistogram(input_inclusive_file_name, config["spectra"])
        else:
            inclusive_histograms = [None] * len(fd_histograms)
        results[name] = OrderedDict()
        for spectrum, hist_orig, hist_inclusive in zip(config["spectra"], fd_histograms, inclusive_histograms):
            results[name].update(PrepareFDhist(spectrum, v["ts"], hist_orig, hist_inclusive, bResponseFile, cResponseFile, bResponseFile_efficiency, cResponseFile_efficiency, unfolding_debug))
        robjects.append(GenerateRootList(results[name], name))

    spectrum_names = GenerateSpectrumNames(config["spectra"], do_unfolding)
    results["SystematicUncertainty"] = CompareVariations(config["variations"], spectrum_names, results)
    robjects.append(GenerateRootList(results["SystematicUncertainty"], "SystematicUncertainty"))
    PlotFSspectraAndSyst(spectrum_names, results)

    output_file_name = "{0}/{1}.root".format(config["input_path"], config["name"])
    outputFile = ROOT.TFile(output_file_name, "recreate")
    if not outputFile or outputFile.IsZombie():
        print("Could not open output file {0}".format(output_file_name))
        exit(1)
    for obj in robjects:
        obj.Write(obj.GetName(), ROOT.TObject.kSingleKey)
    outputFile.Close()

    SaveCanvases(config["name"], config["input_path"])

    print("Done: {0}.".format(output_file_name))

def GenerateSpectrumNames(spectra, do_unfolding):
    basic_spectra = ["{variable}Spectrum{kin_cuts}_{normalization}/GeneratorLevel_{variable}Spectrum"]
    efficiency_spectra = ["{variable}Spectrum{kin_cuts}_{normalization}/GeneratorLevel_{variable}Spectrum_bEfficiencyMultiply",
                            "{variable}Spectrum{kin_cuts}_{normalization}/GeneratorLevel_{variable}Spectrum_bEfficiencyMultiply_cEfficiencyDivide"]
    unfolding_spectra = ["{variable}Spectrum{kin_cuts}_{normalization}/DetectorLevel_{variable}Spectrum_bEfficiencyMultiply",
                 "{variable}Spectrum{kin_cuts}_{normalization}/DetectorLevel_{variable}Spectrum_bEfficiencyMultiply_cEfficiencyDivide",
                 "{variable}Spectrum{kin_cuts}_{normalization}/Unfolded_c_{variable}Spectrum_bEfficiencyMultiply_cEfficiencyDivide"]
    spectra_names = []
    for spectrum in spectra:
        if len(spectrum["d_pt_cuts"]) > 1:
            ptdmin = spectrum["d_pt_cuts"][0]
            d_cuts = "_DPt_{}_{}".format(spectrum["d_pt_cuts"][0] * 10, spectrum["d_pt_cuts"][1] * 10)
        elif len(spectrum["d_pt_cuts"]) == 1:
            ptdmin = spectrum["d_pt_cuts"][0]
            d_cuts = "_DPt_{}".format(int(ptdmin * 10))
        else:
            ptdmin = 0
            d_cuts = ""

        if len(spectrum["jet_pt_cuts"]) > 1:
            jet_cuts = "_JetPt_{}_{}".format(spectrum["jet_pt_cuts"][0], spectrum["jet_pt_cuts"][1])
        elif len(spectrum["jet_pt_cuts"]) == 1:
            jet_cuts = "_JetPt_{}".format(spectrum["jet_pt_cuts"][0])
        else:
            jet_cuts = ""

        kin_cuts = "{}{}".format(d_cuts, jet_cuts)

        format_spectra_names = []
        format_spectra_names.extend(basic_spectra)
        if do_unfolding and spectrum["normalization"] == "cross_section":
            format_spectra_names.extend(efficiency_spectra)
            if spectrum["variable_name"] != "DPt":
                format_spectra_names.extend(unfolding)
        spectra_names.extend([(sname.format(variable=spectrum["variable_name"], kin_cuts=kin_cuts, normalization=HistogramNormalizator.NormalizationType.from_string(spectrum["normalization"]).name), spectrum["normalization"]) for sname in format_spectra_names])

    return spectra_names

def PlotFSspectraAndSyst(spectra, results):
    for name, normalization in spectra:
        stat = GetSpectrum(results["default"], name)
        if not stat: continue
        hname = name[name.rfind("/") + 1:]
        syst = results["SystematicUncertainty"][name]["{0}_CentralAsymmSyst".format(hname)]
        canvas = ROOT.TCanvas("{0}_canvas".format(name), "{0}_canvas".format(name))
        if not "JetZ" in name: canvas.SetLogy()
        canvas.SetLeftMargin(0.13)
        canvas.cd()
        syst_copy = syst.Clone()
        syst_copy.SetFillColor(ROOT.kCyan + 1)
        syst_copy.SetMarkerColor(ROOT.kCyan + 1)
        syst_copy.GetYaxis().SetTitleOffset(1.5)
        syst_copy.Draw("a e2")
        stat_copy = stat.DrawCopy("same p e0 x0")
        stat_copy.SetMarkerColor(ROOT.kBlue + 1)
        stat_copy.SetMarkerStyle(ROOT.kFullCircle)
        stat_copy.SetMarkerSize(0.6)
        stat_copy.SetLineColor(ROOT.kBlue + 1)
        if "JetZ" in name:
            leg = ROOT.TLegend(0.15, 0.71, 0.70, 0.89, "NB")
        else:
            leg = ROOT.TLegend(0.35, 0.71, 0.89, 0.89, "NB")
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(43)
        leg.SetTextSize(20)
        leg.AddEntry(stat_copy, "Central values with stat. unc.", "pe")
        leg.AddEntry(syst_copy, "Systematic uncertainty", "f")
        leg.Draw()
        globalList.append(syst_copy)
        globalList.append(syst_copy.GetHistogram())
        globalList.append(stat_copy)
        globalList.append(leg)
        globalList.append(canvas)

def GetSpectrum(results, name):
    subdirs = name.split("/")
    if len(subdirs) < 1:
        print("{0} is not a valid spectrum name".format(name))
        exit(1)
    for subdir in subdirs:
        if not (isinstance(results, dict) or isinstance(results, OrderedDict)):
            print("ERROR: {0} is not a dictionary!".format(results))
            exit(1)
        if not subdir in results:
            print("Could not find key '{0}' in dictionary '{1}'".format(subdir, results))
            return None
        results = results[subdir]
    return results

def CompareVariationsForSpectrum(comp_template, variations, results, name, normalization):
    h = GetSpectrum(results["default"], name)
    if not h: return None
    comp = copy.deepcopy(comp_template)
    comp.fName = name
    baseline = h.Clone("{0}_copy".format(h.GetName()))
    baseline.SetTitle(variations[0]["title"])
    spectra = []  # all variations
    spectra_syst = []  # envelope is taken from this set of variations (usually all POWHEG paramter variations)
    spectra_syst_add = []  # added in quadrature (usually decay model, e.g. EvtGen vs PYTHIA)
    for v in variations:
        vname = v["name"]
        if not v["active"] or vname == "default": continue
        h = GetSpectrum(results[vname], name)
        h_copy = h.Clone("{0}_copy".format(h.GetName()))
        h_copy.SetTitle(v["title"])
        spectra.append(h_copy)
        if v["systematic"]:
            spectra_syst.append(h_copy)
        else:
            spectra_syst_add.append(h_copy)
    comp.CompareSpectra(baseline, spectra)
    globalList.append(baseline)
    globalList.extend(spectra)
    globalList.extend(comp.fResults)
    return GenerateSystematicUncertainty(baseline, spectra_syst, normalization, spectra_syst_add)

def GenerateSystematicUncertainty(baseline, spectra, normalization, spectra_add):
    hname = baseline.GetName().replace("_copy", "")
    upperLimitsHist = ROOT.TH1D("{0}_UpperSyst".format(hname), "{0}_UpperSyst".format(hname), baseline.GetNbinsX(), baseline.GetXaxis().GetXbins().GetArray())
    upperLimitsHist.Sumw2()
    lowerLimitsHist = ROOT.TH1D("{0}_LowerSyst".format(hname), "{0}_LowerSyst".format(hname), baseline.GetNbinsX(), baseline.GetXaxis().GetXbins().GetArray())
    lowerLimitsHist.Sumw2()
    symmetricLimitsHist = ROOT.TH1D("{0}_SymmSyst".format(hname), "{0}_SymmSyst".format(hname), baseline.GetNbinsX(), baseline.GetXaxis().GetXbins().GetArray())
    symmetricLimitsHist.Sumw2()
    relativeSystHist = ROOT.TH1D("{0}_RelSyst".format(hname), "{0}_RelSyst".format(hname), baseline.GetNbinsX(), baseline.GetXaxis().GetXbins().GetArray())

    result = OrderedDict()
    result[upperLimitsHist.GetName()] = upperLimitsHist
    result[lowerLimitsHist.GetName()] = lowerLimitsHist
    result[symmetricLimitsHist.GetName()] = symmetricLimitsHist
    result[relativeSystHist.GetName()] = relativeSystHist

    for ibin in range(1, baseline.GetNbinsX() + 1):
        centralValue = baseline.GetBinContent(ibin)
        centralValueErr2 = baseline.GetBinError(ibin) ** 2
        for var in spectra:
            diff = var.GetBinContent(ibin) - centralValue
            diffErr = math.sqrt(var.GetBinError(ibin) ** 2 + centralValueErr2)
            if diff > upperLimitsHist.GetBinContent(ibin):
                upperLimitsHist.SetBinContent(ibin, diff)
                upperLimitsHist.SetBinError(ibin, diffErr)
                print("Bin {0}, upper limit {1}".format(ibin, diff))
            if -diff > lowerLimitsHist.GetBinContent(ibin):
                lowerLimitsHist.SetBinContent(ibin, -diff)
                lowerLimitsHist.SetBinError(ibin, diffErr)
                print("Bin {0}, lower limit {1}".format(ibin, -diff))

        if upperLimitsHist.GetBinContent(ibin) > lowerLimitsHist.GetBinContent(ibin):
            symmetricLimitsHist.SetBinContent(ibin, upperLimitsHist.GetBinContent(ibin))
            symmetricLimitsHist.SetBinError(ibin, upperLimitsHist.GetBinError(ibin))
        else:
            symmetricLimitsHist.SetBinContent(ibin, lowerLimitsHist.GetBinContent(ibin))
            symmetricLimitsHist.SetBinError(ibin, lowerLimitsHist.GetBinError(ibin))

    tot_upperLimitsHist = upperLimitsHist.Clone("{0}_TotUpperSyst".format(hname))
    tot_lowerLimitsHist = lowerLimitsHist.Clone("{0}_TotLowerSyst".format(hname))
    tot_symmetricLimitsHist = symmetricLimitsHist.Clone("{0}_TotSymmSyst".format(hname))
    tot_relativeSystHist = relativeSystHist.Clone("{0}_TotRelSyst".format(hname))

    result[tot_upperLimitsHist.GetName()] = tot_upperLimitsHist
    result[tot_lowerLimitsHist.GetName()] = tot_lowerLimitsHist
    result[tot_symmetricLimitsHist.GetName()] = tot_symmetricLimitsHist
    result[tot_relativeSystHist.GetName()] = tot_relativeSystHist

    for ibin in range(1, baseline.GetNbinsX() + 1):
        centralValue = baseline.GetBinContent(ibin)
        tot_upperLimitsHist.SetBinError(ibin, 0)
        tot_lowerLimitsHist.SetBinError(ibin, 0)
        tot_symmetricLimitsHist.SetBinError(ibin, 0)
        tot_relativeSystHist.SetBinError(ibin, 0)
        for var in spectra_add:
            diff = var.GetBinContent(ibin) - centralValue
            if diff > 0:
                tot_upperLimitsHist.SetBinContent(ibin, math.sqrt(diff ** 2 + tot_upperLimitsHist.GetBinContent(ibin) ** 2))
            elif diff < 0:
                tot_lowerLimitsHist.SetBinContent(ibin, math.sqrt(diff ** 2 + tot_lowerLimitsHist.GetBinContent(ibin) ** 2))

            tot_symmetricLimitsHist.SetBinContent(ibin, math.sqrt(diff ** 2 + tot_symmetricLimitsHist.GetBinContent(ibin) ** 2))

        if baseline.GetBinContent(ibin) != 0:
            tot_relativeSystHist.SetBinContent(ibin, tot_symmetricLimitsHist.GetBinContent(ibin) / baseline.GetBinContent(ibin))

    xArray = numpy.array([baseline.GetXaxis().GetBinCenter(ibin) for ibin in range(1, baseline.GetNbinsX() + 1)], dtype=numpy.float32)
    xArrayErr = numpy.array([baseline.GetXaxis().GetBinWidth(ibin) / 2 for ibin in range(1, baseline.GetNbinsX() + 1)], dtype=numpy.float32)
    yArray = numpy.array([baseline.GetBinContent(ibin) for ibin in range(1, baseline.GetNbinsX() + 1)], dtype=numpy.float32)

    yArrayErrUp = numpy.array([upperLimitsHist.GetBinContent(ibin) for ibin in range(1, upperLimitsHist.GetNbinsX() + 1)], dtype=numpy.float32)
    yArrayErrLow = numpy.array([lowerLimitsHist.GetBinContent(ibin) for ibin in range(1, lowerLimitsHist.GetNbinsX() + 1)], dtype=numpy.float32)
    yArrayErrSym = numpy.array([symmetricLimitsHist.GetBinContent(ibin) for ibin in range(1, symmetricLimitsHist.GetNbinsX() + 1)], dtype=numpy.float32)

    yArrayTotErrUp = numpy.array([tot_upperLimitsHist.GetBinContent(ibin) for ibin in range(1, tot_upperLimitsHist.GetNbinsX() + 1)], dtype=numpy.float32)
    yArrayTotErrLow = numpy.array([tot_lowerLimitsHist.GetBinContent(ibin) for ibin in range(1, tot_lowerLimitsHist.GetNbinsX() + 1)], dtype=numpy.float32)
    yArrayTotErrSym = numpy.array([tot_symmetricLimitsHist.GetBinContent(ibin) for ibin in range(1, tot_symmetricLimitsHist.GetNbinsX() + 1)], dtype=numpy.float32)

    symmetricUncGraph = ROOT.TGraphErrors(baseline.GetNbinsX(), xArray, yArray, xArrayErr, yArrayErrSym)
    symmetricUncGraph.SetName("{0}_CentralSymmSyst".format(hname))
    symmetricUncGraph.GetXaxis().SetTitle(baseline.GetXaxis().GetTitle())
    asymmetricUncGraph = ROOT.TGraphAsymmErrors(baseline.GetNbinsX(), xArray, yArray, xArrayErr, xArrayErr, yArrayErrLow, yArrayErrUp)
    asymmetricUncGraph.SetName("{0}_CentralAsymmSyst".format(hname))
    asymmetricUncGraph.GetXaxis().SetTitle(baseline.GetXaxis().GetTitle())

    symmetricTotUncGraph = ROOT.TGraphErrors(baseline.GetNbinsX(), xArray, yArray, xArrayErr, yArrayTotErrSym)
    symmetricTotUncGraph.SetName("{0}_CentralTotSymmSyst".format(hname))
    symmetricTotUncGraph.GetXaxis().SetTitle(baseline.GetXaxis().GetTitle())
    asymmetricTotUncGraph = ROOT.TGraphAsymmErrors(baseline.GetNbinsX(), xArray, yArray, xArrayErr, xArrayErr, yArrayTotErrLow, yArrayTotErrUp)
    asymmetricTotUncGraph.SetName("{0}_CentralTotAsymmSyst".format(hname))
    asymmetricTotUncGraph.GetXaxis().SetTitle(baseline.GetXaxis().GetTitle())

    if "it{p}" in baseline.GetXaxis().GetTitle():
        if normalization == "cross_section":
            symmetricUncGraph.GetYaxis().SetTitle("#frac{d#sigma}{d#it{p}_{T}} [(mb) (GeV/#it{c})^{-1}]")
            asymmetricUncGraph.GetYaxis().SetTitle("#frac{d#sigma}{d#it{p}_{T}} [(mb) (GeV/#it{c})^{-1}]")
            symmetricTotUncGraph.GetYaxis().SetTitle("#frac{d#sigma}{d#it{p}_{T}} [(mb) (GeV/#it{c})^{-1}]")
            asymmetricTotUncGraph.GetYaxis().SetTitle("#frac{d#sigma}{d#it{p}_{T}} [(mb) (GeV/#it{c})^{-1}]")
        elif normalization == "distribution":
            symmetricUncGraph.GetYaxis().SetTitle("Probability Density (GeV/#it{c})^{-1}")
            asymmetricUncGraph.GetYaxis().SetTitle("Probability Density (GeV/#it{c})^{-1}")
            symmetricTotUncGraph.GetYaxis().SetTitle("Probability Density (GeV/#it{c})^{-1}")
            asymmetricTotUncGraph.GetYaxis().SetTitle("Probability Density (GeV/#it{c})^{-1}")
    elif "it{z}" in baseline.GetXaxis().GetTitle():
        if normalization == "cross_section":
            symmetricUncGraph.GetYaxis().SetTitle("#frac{d#sigma}{d#it{z}_{||}^{ch}} (mb)")
            asymmetricUncGraph.GetYaxis().SetTitle("#frac{d#sigma}{d#it{z}_{||}^{ch}} (mb)")
            symmetricTotUncGraph.GetYaxis().SetTitle("#frac{d#sigma}{d#it{z}_{||}^{ch}} (mb)")
            asymmetricTotUncGraph.GetYaxis().SetTitle("#frac{d#sigma}{d#it{z}_{||}^{ch}} (mb)")
        elif normalization == "distribution":
            symmetricUncGraph.GetYaxis().SetTitle("Probability Density")
            asymmetricUncGraph.GetYaxis().SetTitle("Probability Density")
            symmetricTotUncGraph.GetYaxis().SetTitle("Probability Density")
            asymmetricTotUncGraph.GetYaxis().SetTitle("Probability Density")

    result[symmetricUncGraph.GetName()] = symmetricUncGraph
    result[asymmetricUncGraph.GetName()] = asymmetricUncGraph
    result[symmetricTotUncGraph.GetName()] = symmetricTotUncGraph
    result[asymmetricTotUncGraph.GetName()] = asymmetricTotUncGraph
    return result

def CompareVariations(variations, spectra, results):
    comp_template = DMesonJetCompare.DMesonJetCompare("DMesonJetCompare")
    comp_template.fOptRatio = "hist"
    comp_template.fLogUpperSpace = 3
    comp_template.fNColsLegRatio = 2
    comp_template.fX1LegRatio = 0.2

    result = OrderedDict()

    for name, normalization in spectra:
        if "JetZ" in name:
            comp_template.fDoSpectraPlot = "lineary"
            comp_template.fX1LegSpectrum = 0.15
            comp_template.fX2LegSpectrum = 0.70
        else:
            comp_template.fDoSpectraPlot = "logy"
            comp_template.fX1LegSpectrum = 0.55
            comp_template.fX2LegSpectrum = 0.90
        r = CompareVariationsForSpectrum(comp_template, variations, results, name, normalization)
        if r: result[name] = r

    return result

def SaveCanvases(name, input_path):
    path = "{}/{}".format(input_path, name)
    if not os.path.exists(path):
        os.makedirs(path)
    for obj in globalList:
        if isinstance(obj, ROOT.TCanvas):
            oname = obj.GetName().replace("/", "_")
            obj.SaveAs("{}/{}.pdf".format(path, oname))

def GenerateRootList(pdict, name):
    rlist = ROOT.TList()
    rlist.SetName(name)
    print("Loop in {0}".format(name))
    for nobj, obj in pdict.iteritems():
        if isinstance(obj, ROOT.TObject):
            rlist.Add(obj)
            print("Adding {0}".format(nobj))
        elif isinstance(obj, dict) or isinstance(obj, OrderedDict):
            print("Recursion in {0}".format(nobj))
            rlist.Add(GenerateRootList(obj, nobj))
        else:
            print("Error: type of object {0} not recognized!".format(obj))
    return rlist

def PrepareFDhist(spectrum, ts, FDhistogram_old, inclusive_histogram, bResponseFile, cResponseFile, bResponseFile_efficiency, cResponseFile_efficiency, unfolding_debug):
    if spectrum["variable_name"] == "DPt":
        return PrepareFDhist_dpt(spectrum, ts, FDhistogram_old, bResponseFile, cResponseFile, bResponseFile_efficiency, cResponseFile_efficiency)
    else:
        return PrepareFDhist_jet(spectrum, ts, FDhistogram_old, inclusive_histogram, bResponseFile, cResponseFile, bResponseFile_efficiency, cResponseFile_efficiency, unfolding_debug)

def PrepareFDhist_dpt(spectrum, ts, FDhistogram_old, bResponseFile, cResponseFile, bResponseFile_efficiency, cResponseFile_efficiency):
    print("Preparing D pt FD histograms")

    if len(spectrum["jet_pt_cuts"]) > 1:
        jet_cuts = "_JetPt_{}_{}".format(spectrum["jet_pt_cuts"][0], spectrum["jet_pt_cuts"][1])
    elif len(spectrum["jet_pt_cuts"]) == 1:
        jet_cuts = "_JetPt_{}".format(spectrum["jet_pt_cuts"][0])
    else:
        jet_cuts = ""

    spectrum_name = "DPtSpectrum{}".format(jet_cuts)

    dmeson = spectrum["d_meson"]
    dmeson_prompt = "Prompt_{}".format(dmeson)
    dmeson_nonprompt = "NonPrompt_{}".format(dmeson)

    result = OrderedDict()

    dpt = OrderedDict()
    out_spectrum_name = "{}_{}".format(spectrum_name, HistogramNormalizator.NormalizationType.from_string(spectrum["normalization"]).name)
    result[out_spectrum_name] = dpt

    dptbins = spectrum["d_pt_bins"]

    print("Rebinning the FD original histogram")
    FDhistogram_orig = FDhistogram_old.Rebin(len(dptbins) - 1, FDhistogram_old.GetName(), array.array('d', dptbins))
    FDhistogram_orig.SetName("GeneratorLevel_DPtSpectrum")
    normalizator = HistogramNormalizator.MCSimulationNormalizator(FDhistogram_orig, spectrum["normalization"])
    if "normalization_configuration" in spectrum and ("xmin" in spectrum["normalization_configuration"] and "xmax" in spectrum["normalization_configuration"]):
        normalizator.fXmin = spectrum["normalization_configuration"]["xmin"]
        normalizator.fXmax = spectrum["normalization_configuration"]["xmax"]
    FDhistogram_orig = normalizator.NormalizeHistogram()
    dpt[FDhistogram_orig.GetName()] = FDhistogram_orig

    xaxis_title = "#it{p}_{T,D} (GeV/#it{c})"
    if spectrum["normalization"] == "cross_section":
        yaxis_title = "#frac{d^{2}#sigma}{d#it{p}_{T,D}d#it{#eta}} [mb (GeV/#it{c})^{-1}]"
    elif spectrum["normalization"] == "distribution":
        yaxis_title = "Probability Density (GeV/#it{c})^{-1}"
    
    FDhistogram_orig.GetXaxis().SetTitle(xaxis_title)
    FDhistogram_orig.GetYaxis().SetTitle(yaxis_title)

    if spectrum["normalization"] == "cross_section":
        if bResponseFile or cResponseFile:
            responseList = OrderedDict()
            dpt["DetectorResponse"] = responseList

        if bResponseFile:
            print("Loading the response matrix b->D0")
            temp, bEfficiency_dpt = LoadResponse(bResponseFile, dmeson_nonprompt, spectrum_name, "", "b", FDhistogram_orig.GetNbinsX(), FDhistogram_orig.GetXaxis().GetXbins().GetArray())
            bEfficiency_dpt.SetName("EfficiencyVsDPt_b")
            responseList[bEfficiency_dpt.GetName()] = bEfficiency_dpt
        else:
            bEfficiency_dpt = None

        if cResponseFile:
            print("Loading the response matrix c->D0")
            temp, cEfficiency_dpt = LoadResponse(cResponseFile, dmeson_prompt, spectrum_name, "", "c", FDhistogram_orig.GetNbinsX(), FDhistogram_orig.GetXaxis().GetXbins().GetArray())
            cEfficiency_dpt.SetName("EfficiencyVsDPt_c")
            responseList[cEfficiency_dpt.GetName()] = cEfficiency_dpt
        else:
            cEfficiency_dpt = None
    else:
        bEfficiency_dpt = None
        cEfficiency_dpt = None

    if bEfficiency_dpt:
        print("Applying the b->D0 reconstruction efficiency")
        FDhistogram_dpt = FDhistogram_orig.Clone("{0}_bEff".format(FDhistogram_orig.GetName()))
        FDhistogram_dpt.Multiply(bEfficiency_dpt)
        FDhistogram_dpt.SetName("GeneratorLevel_DPtSpectrum_bEfficiencyMultiply")
        dpt[FDhistogram_dpt.GetName()] = FDhistogram_dpt
    else:
        FDhistogram_dpt = None

    if FDhistogram_dpt and cEfficiency_dpt:
        print("Applying the correction for the c->D0 reconstruction efficiency")
        FDhistogram_efficiency_dpt = FDhistogram_dpt.Clone("{0}_efficiency".format(FDhistogram_dpt.GetName()))
        FDhistogram_efficiency_dpt.Divide(cEfficiency_dpt)
        FDhistogram_efficiency_dpt.SetName("GeneratorLevel_DPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide")
        dpt[FDhistogram_efficiency_dpt.GetName()] = FDhistogram_efficiency_dpt
    else:
        FDhistogram_efficiency_dpt = None

    return result

def PrepareFDhist_jet(spectrum, ts, FDhistogram_old, inclusive_histogram, bResponseFile, cResponseFile, bResponseFile_efficiency, cResponseFile_efficiency, unfolding_debug):
    dmeson = spectrum["d_meson"]
    dmeson_prompt = "Prompt_{}".format(dmeson)
    dmeson_nonprompt = "NonPrompt_{}".format(dmeson)
    jet_var_name = spectrum["variable_name"]
    if jet_var_name == "JetPt":
        jetxbins = spectrum["jet_pt_bins"]
        xaxis_title = "#it{p}_{T,jet}^{ch} (GeV/#it{c})"
        if spectrum["normalization"] == "cross_section":
            yaxis_title = "#frac{d^{2}#sigma}{d#it{p}_{T,jet}^{ch}d#it{#eta}} [mb (GeV/#it{c})^{-1}]"
        elif spectrum["normalization"] == "distribution":
            yaxis_title = "Probability Density (GeV/#it{c})^{-1}"
        elif spectrum["normalization"] == "ratio":
            yaxis_title = "#it{R}(#it{p}_{T,jet}^{ch})"
        else:
            print("Normalization option '{}' invalid".format(config["normalization"]))
            exit(1)
    elif jet_var_name == "JetZ":
        jetxbins = spectrum["d_z_bins"]
        xaxis_title = "#it{z}_{||}^{ch}"
        if spectrum["normalization"] == "cross_section":
            yaxis_title = "#frac{d^{2}#sigma}{d#it{z}_{||}^{ch}d#it{#eta}} (mb)"
        elif spectrum["normalization"] == "distribution":
            yaxis_title = "Probability Density"
        elif spectrum["normalization"] == "rate":
            yaxis_title = "#frac{#it{R}(#it{z}_{||}^{ch})}{#Delta#it{z}_{||}^{ch}}"
        else:
            print("Normalization option '{}' invalid".format(config["normalization"]))
            exit(1)
    else:
        print("Error: {} variable not known!".format(jet_var_name))
        exit(1)

    if len(spectrum["d_pt_cuts"]) > 1:
        ptdmin = spectrum["d_pt_cuts"][0]
        d_cuts = "_DPt_{}_{}".format(spectrum["d_pt_cuts"][0] * 10, spectrum["d_pt_cuts"][1] * 10)
    elif len(spectrum["d_pt_cuts"]) == 1:
        ptdmin = spectrum["d_pt_cuts"][0]
        d_cuts = "_DPt_{}".format(int(ptdmin * 10))
    else:
        ptdmin = 0
        d_cuts = ""

    if len(spectrum["jet_pt_cuts"]) > 1:
        jet_cuts = "_JetPt_{}_{}".format(spectrum["jet_pt_cuts"][0], spectrum["jet_pt_cuts"][1])
    elif len(spectrum["jet_pt_cuts"]) == 1:
        jet_cuts = "_JetPt_{}".format(spectrum["jet_pt_cuts"][0])
    else:
        jet_cuts = ""

    dptbins = spectrum["d_pt_bins"]

    print("Preparing {} FD histograms".format(jet_var_name))
    result = OrderedDict()

    def PrepareHistogramNormalizator(hist):
        normalizator = HistogramNormalizator.MCSimulationNormalizator(hist, spectrum["normalization"])
        if "normalization_configuration" in spectrum and ("xmin" in spectrum["normalization_configuration"] and "xmax" in spectrum["normalization_configuration"]):
            normalizator.fXmin = spectrum["normalization_configuration"]["xmin"]
            normalizator.fXmax = spectrum["normalization_configuration"]["xmax"]
        if spectrum["normalization"] == "ratio" or spectrum["normalization"] == "rate":
            if inclusive_histogram:
                normalizator.fNormalizationHistogram = inclusive_histogram
            else:
                print("Spectrum {}: the inclusive histogram was not provided. It is needed when normalization is '{}'".format(jet_var_name, spectrum["normalization"]))
                exit(1)
        return normalizator

    spectrum_2d_name = "{}DPtSpectrum".format(jet_var_name)

    jetxdpt = OrderedDict()
    result[spectrum_2d_name] = jetxdpt

    print("Rebinning the FD original histogram")
    FDhistogram_orig = DMesonJetUtils.Rebin2D_fromBins(FDhistogram_old, FDhistogram_old.GetName(), len(jetxbins) - 1, array.array('d', jetxbins), len(dptbins) - 1, array.array('d', dptbins))
    FDhistogram_orig.SetName("GeneratorLevel_{}".format(spectrum_2d_name))
    jetxdpt[FDhistogram_orig.GetName()] = FDhistogram_orig

    if spectrum["normalization"] == "cross_section":
        if bResponseFile or cResponseFile:
            responseList = OrderedDict()
            jetxdpt["DetectorResponse"] = responseList

        if bResponseFile:
            print("Loading the response matrix b->D0")
            temp, bEfficiency_dpt = LoadResponse(bResponseFile, dmeson_nonprompt, "DPtSpectrum{}".format(jet_cuts), "", "b", FDhistogram_orig.GetNbinsY(), FDhistogram_orig.GetYaxis().GetXbins().GetArray())
            bEfficiency_dpt.SetName("EfficiencyVsDPt_b")
            responseList[bEfficiency_dpt.GetName()] = bEfficiency_dpt
        else:
            bEfficiency_dpt = None

        if cResponseFile:
            print("Loading the response matrix c->D0")
            temp, cEfficiency_dpt = LoadResponse(cResponseFile, dmeson_prompt, "DPtSpectrum{}".format(jet_cuts), "", "c", FDhistogram_orig.GetNbinsY(), FDhistogram_orig.GetYaxis().GetXbins().GetArray())
            cEfficiency_dpt.SetName("EfficiencyVsDPt_c")
            responseList[cEfficiency_dpt.GetName()] = cEfficiency_dpt
        else:
            cEfficiency_dpt = None
    else:
        bEfficiency_dpt = None
        cEfficiency_dpt = None

    if bEfficiency_dpt:
        print("Applying the b->D0 reconstruction efficiency")
        FDhistogram = FDhistogram_orig.Clone("{0}_bEff".format(FDhistogram_orig.GetName()))
        ApplyEfficiency(FDhistogram, bEfficiency_dpt, False)
        FDhistogram.SetName("GeneratorLevel_{}_bEfficiencyMultiply".format(spectrum_2d_name))
        jetxdpt[FDhistogram.GetName()] = FDhistogram
    else:
       FDhistogram = None

    if FDhistogram and cEfficiency_dpt:
        print("Applying the correction for the c->D0 reconstruction efficiency")
        FDhistogram_efficiency = FDhistogram.Clone("{0}_efficiency".format(FDhistogram.GetName()))
        ApplyEfficiency(FDhistogram_efficiency, cEfficiency_dpt, True)
        FDhistogram_efficiency.SetName("GeneratorLevel_{}_bEfficiencyMultiply_cEfficiencyDivide".format(spectrum_2d_name))
        jetxdpt[FDhistogram_efficiency.GetName()] = FDhistogram_efficiency
    else:
        FDhistogram_efficiency = None

    if bEfficiency_dpt:
        # Finer bin versions used as priors in the folding/unfolding
        print("Applying the b->D0 reconstruction efficiency (fine bins)")
        FDhistogram_fineBins = DMesonJetUtils.Rebin2D_fromBins(FDhistogram_old, "{0}_fineBins_bEff".format(FDhistogram_old.GetName()), FDhistogram_old.GetNbinsX(), FDhistogram_old.GetXaxis().GetXbins().GetArray(), len(dptbins) - 1, array.array('d', dptbins))
        ApplyEfficiency(FDhistogram_fineBins, bEfficiency_dpt, False)
        FDhistogram_fineBins.SetName("GeneratorLevel_{}_FineBins_bEfficiencyMultiply".format(spectrum_2d_name))
        jetxdpt[FDhistogram_fineBins.GetName()] = FDhistogram_fineBins
    else:
        FDhistogram_fineBins = None

    if FDhistogram_fineBins and cEfficiency_dpt:
        print("Applying the correction for the c->D0 reconstruction efficiency (fine bins)")
        FDhistogram_fineBins_efficiency = FDhistogram_fineBins.Clone("{0}_efficiency".format(FDhistogram_fineBins.GetName()))
        ApplyEfficiency(FDhistogram_fineBins_efficiency, cEfficiency_dpt, True)
        FDhistogram_fineBins_efficiency.SetName("GeneratorLevel_{}_FineBins_bEfficiencyMultiply_cEfficiencyDivide".format(spectrum_2d_name))
        jetxdpt[FDhistogram_fineBins_efficiency.GetName()] = FDhistogram_fineBins_efficiency
    else:
        FDhistogram_fineBins_efficiency = None

    spectrumName = "{}Spectrum".format(jet_var_name)
    fullSpectrumName = "{}{}{}".format(spectrumName, d_cuts, jet_cuts)
    out_spectrum_name = "{}_{}".format(fullSpectrumName, HistogramNormalizator.NormalizationType.from_string(spectrum["normalization"]).name)

    ptdList = OrderedDict()
    result[out_spectrum_name] = ptdList

    print("Projecting FD histogram into the jet pt axis for D pt > {0} GeV/c (w/o b or c efficiency)".format(ptdmin))
    FDhistogram_jetx_orig = FDhistogram_orig.ProjectionX("{}_{}_DPt_{}".format(FDhistogram_orig.GetName(), jet_var_name, ptdmin * 10), FDhistogram_orig.GetYaxis().FindBin(ptdmin), FDhistogram_orig.GetNbinsY() + 1)
    FDhistogram_jetx_orig.SetName("GeneratorLevel_{}".format(spectrumName))
    normalizator = PrepareHistogramNormalizator(FDhistogram_jetx_orig)
    FDhistogram_jetx_orig = normalizator.NormalizeHistogram()
    FDhistogram_jetx_orig.GetXaxis().SetTitle(xaxis_title)
    FDhistogram_jetx_orig.GetYaxis().SetTitle(yaxis_title)
    ptdList[FDhistogram_jetx_orig.GetName()] = FDhistogram_jetx_orig

    if FDhistogram:
        print("Projecting FD histogram into the jet pt axis for D pt > {0} GeV/c (w/o efficiency)".format(ptdmin))
        FDhistogram_jetx = FDhistogram.ProjectionX("{}_{}_DPt_{}".format(FDhistogram.GetName(), jet_var_name, ptdmin * 10), FDhistogram.GetYaxis().FindBin(ptdmin), FDhistogram.GetNbinsY() + 1)
        FDhistogram_jetx.SetName("GeneratorLevel_{}_bEfficiencyMultiply".format(spectrumName))
        normalizator = PrepareHistogramNormalizator(FDhistogram_jetx)
        FDhistogram_jetx = normalizator.NormalizeHistogram()
        FDhistogram_jetx.GetXaxis().SetTitle(xaxis_title)
        FDhistogram_jetx.GetYaxis().SetTitle(yaxis_title)
        ptdList[FDhistogram_jetx.GetName()] = FDhistogram_jetx
    else:
        FDhistogram_jetx = None

    if FDhistogram_fineBins:
        print("Projecting FD histogram into the jet pt axis for D pt > {0} GeV/c (w/o efficiency, fine bins)".format(ptdmin))
        FDhistogram_fineBins_jetx = FDhistogram_fineBins.ProjectionX("{}_{}_DPt_{}".format(FDhistogram_fineBins.GetName(), jet_var_name, ptdmin * 10), FDhistogram_fineBins.GetYaxis().FindBin(ptdmin), FDhistogram_fineBins.GetNbinsY() + 1)
        FDhistogram_fineBins_jetx.SetName("GeneratorLevel_{}_FineBins_bEfficiencyMultiply".format(spectrumName))
        normalizator = PrepareHistogramNormalizator(FDhistogram_fineBins_jetx)
        FDhistogram_fineBins_jetx = normalizator.NormalizeHistogram()
        FDhistogram_fineBins_jetx.GetXaxis().SetTitle(xaxis_title)
        FDhistogram_fineBins_jetx.GetYaxis().SetTitle(yaxis_title)
        ptdList[FDhistogram_fineBins_jetx.GetName()] = FDhistogram_fineBins_jetx
    else:
        FDhistogram_fineBins_jetx = None

    if bResponseFile:
        print("Applying the jet b response matrix")
        FDhistogram_jetx_detector = ApplyResponse(FDhistogram_jetx, bResponseFile, FDhistogram_fineBins_jetx, fullSpectrumName, dmeson_nonprompt)
        FDhistogram_jetx_detector.SetName("DetectorLevel_{}_bEfficiencyMultiply".format(spectrumName))
        ptdList[FDhistogram_jetx_detector.GetName()] = FDhistogram_jetx_detector
    else:
        FDhistogram_jetx_detector = None

    if cResponseFile and FDhistogram_jetx_detector:
        print("Unfolding using the c response matrix")
        FDhistogram_jetx_unfolded = GetUnfoldedSpectrum("{0}_{1}_unfolded_c".format(ts, FDhistogram_jetx_detector.GetName()), FDhistogram_jetx_detector, cResponseFile, FDhistogram_fineBins_jetx, fullSpectrumName, dmeson_prompt, unfolding_debug)
        FDhistogram_jetx_unfolded.SetName("Unfolded_c_{}_bEfficiencyMultiply".format(spectrumName))
        ptdList[FDhistogram_jetx_unfolded.GetName()] = FDhistogram_jetx_unfolded
    else:
        FDhistogram_jetx_unfolded = None

    if FDhistogram_jetx_detector and bResponseFile:
        print("Unfolding using the b response matrix")
        FDhistogram_jetx_unfolded = GetUnfoldedSpectrum("{0}_{1}_unfolded_b".format(ts, FDhistogram_jetx_detector.GetName()), FDhistogram_jetx_detector, bResponseFile, FDhistogram_fineBins_jetx, fullSpectrumName, dmeson_nonprompt, unfolding_debug)
        FDhistogram_jetx_unfolded.SetName("Unfolded_b_{}_bEfficiencyMultiply".format(spectrumName))
        ptdList[FDhistogram_jetx_unfolded.GetName()] = FDhistogram_jetx_unfolded
    else:
        FDhistogram_jetx_unfolded = None

    if FDhistogram_efficiency:
        print("Projecting FD histogram into the jet axis for D pt > {0} GeV/c (w/ efficiency)".format(ptdmin))
        FDhistogram_efficiency_jetx = FDhistogram_efficiency.ProjectionX("{}_{}_DPt_{}".format(FDhistogram_efficiency.GetName(), jet_var_name, ptdmin * 10), FDhistogram_efficiency.GetYaxis().FindBin(ptdmin), FDhistogram_efficiency.GetNbinsY() + 1)
        FDhistogram_efficiency_jetx.SetName("GeneratorLevel_{}_bEfficiencyMultiply_cEfficiencyDivide".format(spectrumName))
        normalizator = PrepareHistogramNormalizator(FDhistogram_efficiency_jetx)
        FDhistogram_efficiency_jetx = normalizator.NormalizeHistogram()
        FDhistogram_efficiency_jetx.GetXaxis().SetTitle(xaxis_title)
        FDhistogram_efficiency_jetx.GetYaxis().SetTitle(yaxis_title)
        ptdList[FDhistogram_efficiency_jetx.GetName()] = FDhistogram_efficiency_jetx
    else:
        FDhistogram_efficiency_jetx = None

    if FDhistogram_fineBins_efficiency:
        print("Projecting FD histogram into the jet axis for D pt > {0} GeV/c (w/ efficiency, fine bins)".format(ptdmin))
        FDhistogram_fineBins_efficiency_jetx = FDhistogram_fineBins_efficiency.ProjectionX("{}_{}_DPt_{}".format(FDhistogram_fineBins_efficiency.GetName(), jet_var_name, ptdmin * 10), FDhistogram_fineBins_efficiency.GetYaxis().FindBin(ptdmin), FDhistogram_fineBins_efficiency.GetNbinsY() + 1)
        FDhistogram_fineBins_efficiency_jetx.SetName("GeneratorLevel_{}_FineBins_bEfficiencyMultiply_cEfficiencyDivide".format(spectrumName))
        normalizator = PrepareHistogramNormalizator(FDhistogram_fineBins_efficiency_jetx)
        FDhistogram_fineBins_efficiency_jetx = normalizator.NormalizeHistogram()
        FDhistogram_fineBins_efficiency_jetx.GetXaxis().SetTitle(xaxis_title)
        FDhistogram_fineBins_efficiency_jetx.GetYaxis().SetTitle(yaxis_title)
        ptdList[FDhistogram_fineBins_efficiency_jetx.GetName()] = FDhistogram_fineBins_efficiency_jetx
    else:
        FDhistogram_fineBins_efficiency_jetx = None

    if FDhistogram_efficiency_jetx and bResponseFile_efficiency:
        print("Applying the jet b response matrix")
        FDhistogram_efficiency_jetx_detector = ApplyResponse(FDhistogram_efficiency_jetx, bResponseFile_efficiency, FDhistogram_fineBins_efficiency_jetx, fullSpectrumName, dmeson_nonprompt)
        FDhistogram_efficiency_jetx_detector.SetName("DetectorLevel_{}_bEfficiencyMultiply_cEfficiencyDivide".format(spectrumName))
        ptdList[FDhistogram_efficiency_jetx_detector.GetName()] = FDhistogram_efficiency_jetx_detector
    else:
        FDhistogram_efficiency_jetx_detector = None

    if FDhistogram_efficiency_jetx_detector and cResponseFile_efficiency:
        print("Unfolding using the c response matrix")
        FDhistogram_efficiency_jetx_unfolded = GetUnfoldedSpectrum("{0}_{1}_unfolded_c".format(ts, FDhistogram_efficiency_jetx_detector.GetName()), FDhistogram_efficiency_jetx_detector, cResponseFile_efficiency, FDhistogram_fineBins_efficiency_jetx, fullSpectrumName, dmeson_prompt, unfolding_debug)
        FDhistogram_efficiency_jetx_unfolded.SetName("Unfolded_c_{}_bEfficiencyMultiply_cEfficiencyDivide".format(spectrumName))
        ptdList[FDhistogram_efficiency_jetx_unfolded.GetName()] = FDhistogram_efficiency_jetx_unfolded
    else:
        FDhistogram_efficiency_jetx_unfolded = None

    if FDhistogram_efficiency_jetx_detector and bResponseFile_efficiency:
        print("Unfolding using the b response matrix")
        FDhistogram_efficiency_jetx_unfolded_b = GetUnfoldedSpectrum("{0}_{1}_unfolded_b".format(ts, FDhistogram_efficiency_jetx_detector.GetName()), FDhistogram_efficiency_jetx_detector, bResponseFile_efficiency, FDhistogram_fineBins_efficiency_jetx, fullSpectrumName, dmeson_nonprompt, unfolding_debug)
        FDhistogram_efficiency_jetx_unfolded_b.SetName("Unfolded_b_{}_bEfficiencyMultiply_cEfficiencyDivide".format(spectrumName))
        ptdList[FDhistogram_efficiency_jetx_unfolded_b.GetName()] = FDhistogram_efficiency_jetx_unfolded_b
    else:
        FDhistogram_efficiency_jetx_unfolded_b = None

    return result

def GenerateUnfoldingEngine(name, responseFile, prior, spectrumName, nbins, dmeson):
    analysis = dict()
    analysis["name"] = name
    analysis["d_meson"] = "_d_meson_"
    analysis["variable"] = "_variable_"
    analysis["jet_type"] = "Charged"
    analysis["jet_radius"] = "R040"
    analysis["active"] = True
    analysis["kinematic_cuts"] = None
    analysis["raw_yield_method"] = None
    analysis["d_meson_response"] = dmeson
    analysis["spectrum_response"] = spectrumName
    analysis["priors"] = ["GeneratedSpectrum"]
    analysis["default_prior"] = "GeneratedSpectrum"
    analysis["default_method"] = "Svd"
    svd = dict()
    default_reg_svd = {"GeneratedSpectrum" : nbins}
    svd["default_reg"] = default_reg_svd
    bayes = dict()
    default_reg_bayes = {"GeneratedSpectrum" : 6}
    bayes["default_reg"] = default_reg_bayes
    bayes["iter_min"] = 1
    bayes["iter_max"] = 6
    bayes["iter_step"] = 1
    binBybin = dict()
    binBybin["default_reg"] = None
    analysis["methods"] = {"Svd" : svd, "Bayes" : bayes, "BinByBin" : binBybin}
    unf = DMesonJetUnfolding.DMesonJetUnfoldingEngine(None, None, "unfolding", analysis)
    unf.fDoErrorCrossChecks = False
    unf.SetCustomPrior("GeneratedSpectrum", prior)
    unf.fUseOverflow = True
    unf.fResponseFile = responseFile
    unf.LoadResponse(False)
    return unf

def GetUnfoldedSpectrum(name, histogram, responseFile, prior, spectrumName, dmeson, unfolding_debug):
    unf = GenerateUnfoldingEngine(name, responseFile, prior, spectrumName, histogram.GetXaxis().GetNbins(), dmeson)
    
    # Need to remove the bin width normalization
    h_copy = histogram.Clone()
    for ibin in range(h_copy.GetNbinsX() + 2):
        h_copy.SetBinContent(ibin, h_copy.GetBinContent(ibin) * h_copy.GetXaxis().GetBinWidth(ibin))
        h_copy.SetBinError(ibin, h_copy.GetBinError(ibin) * h_copy.GetXaxis().GetBinWidth(ibin))

    unf.fInputSpectrum = h_copy
    unf.fTruthSpectrum = None
    unf.Start(unfolding_debug)
    default_reg = unf.GetDefaultRegularization(unf.fDefaultMethod, unf.fDefaultPrior)
    default_unfolded = unf.fUnfoldedSpectra[unf.fDefaultMethod, default_reg, unf.fDefaultPrior]
    default_unfolded.SetName(histogram.GetName().replace("_detector", "_unfolded"))

    default_unfolded.Scale(1.0, "width")
    return default_unfolded

def ApplyEfficiency(hist, efficiency, reverse):
    print("Applying efficiency {0} to histogram {1}".format(efficiency.GetName(), hist.GetName()))
    for ybin in range(0, hist.GetNbinsY() + 2):
        if reverse and efficiency.GetBinContent(ybin) > 0:
            print("(bin {0:.2f} ({1:.2f}): dividing by {2}".format(hist.GetYaxis().GetBinUpEdge(ybin), efficiency.GetXaxis().GetBinUpEdge(ybin), efficiency.GetBinContent(ybin)))
            eff = 1. / efficiency.GetBinContent(ybin)
        else:
            print("(bin {0:.2f} ({1:.2f}): multiplying by {2}".format(hist.GetYaxis().GetBinUpEdge(ybin), efficiency.GetXaxis().GetBinUpEdge(ybin), efficiency.GetBinContent(ybin)))
            eff = efficiency.GetBinContent(ybin)
        for xbin in range(0, hist.GetNbinsX() + 2):
            hist.SetBinContent(xbin, ybin, hist.GetBinContent(xbin, ybin) * eff)
            hist.SetBinError(xbin, ybin, hist.GetBinError(xbin, ybin) * eff)

def ApplyResponse(truth, responseFile, prior, spectrumName, dmeson):
    # Need to remove the bin width normalization
    h_copy = truth.Clone()
    for ibin in range(0, h_copy.GetNbinsX() + 2):
        h_copy.SetBinContent(ibin, h_copy.GetBinContent(ibin) * h_copy.GetXaxis().GetBinWidth(ibin))
        h_copy.SetBinError(ibin, h_copy.GetBinError(ibin) * h_copy.GetXaxis().GetBinWidth(ibin))

    unf = GenerateUnfoldingEngine(truth.GetName(), responseFile, prior, spectrumName, truth.GetXaxis().GetNbins(), dmeson)

    unf.fInputSpectrum = h_copy
    unf.fTruthSpectrum = None
    unf.GenerateResponse()
    resp = unf.fResponseMatrices["GeneratedSpectrum"].fNormResponse.Clone("temp")
    for ybin in range(0, resp.GetYaxis().GetNbins() + 2):
        inty = resp.Integral(0, resp.GetXaxis().GetNbins() + 1, ybin, ybin)
        if inty == 0:
            continue
        print("Integral is {0:.3f}".format(inty))
        scaling = h_copy.GetBinContent(ybin) / inty
        for xbin in range(0, resp.GetXaxis().GetNbins() + 2):
            binValue = resp.GetBinContent(xbin, ybin) * scaling
            binErr = resp.GetBinError(xbin, ybin) * scaling
            resp.SetBinContent(xbin, ybin, binValue)
            resp.SetBinError(xbin, ybin, binErr)
    result = resp.ProjectionX("{0}_detector".format(truth.GetName()), 0, -1)
    result.GetYaxis().SetTitle(truth.GetYaxis().GetTitle())

    result.Scale(1.0, "width")

    return result

def OpenResponseFile(input_path, response, efficiency):
    if efficiency:
        file_name = "{0}/{1}/{2}_efficiency.root".format(input_path, response["train"], response["name"])
    else:
        file_name = "{0}/{1}/{2}.root".format(input_path, response["train"], response["name"])
    file = ROOT.TFile(file_name)
    print("Response file: {0}".format(file_name))
    if not file or file.IsZombie():
        print("Could not open input file {0}".format(file_name))
        exit(1)
    return file

def LoadResponse(responseFile, dmeson, spectrumName, suffix_in, suffix_out, detAxisNbins, detAxisBin):
    loader = DetectorResponseLoader.DetectorResponseLoader(responseFile, dmeson, "Jet_AKTChargedR040_pt_scheme", spectrumName, suffix_in, detAxisNbins, detAxisBin)
    resp_coarse = loader.GetResponseMatrixObject()
    resp_coarse.SetName("{}_coarse_{}".format(resp_coarse.GetName(), suffix_out))
    eff_coarse = loader.GetEfficiencyObject()
    eff_coarse.SetName("{}_coarse_{}".format(eff_coarse.GetName(), suffix_out))
    return resp_coarse, eff_coarse

def LoadInclusiveHistogram(file_name, spectra):
    result = []
    for spectrum in spectra:
        if "normalization_configuration" in spectrum and (spectrum["normalization"] == "rate" or spectrum["normalization"] == "ratio"):
            spectrum_name = spectrum["normalization_configuration"]["spectrum_name"]
            file = ROOT.TFile(file_name)
            if not file or file.IsZombie():
                print("Could not open file {0}".format(file_name))
                exit(1)
            else:
                print("File {0} open".format(file_name))
            objname = "Charged_R040/{spectrum_name}".format(spectrum_name=spectrum_name)
            h = DMesonJetUtils.GetObject(file, objname)
            if not h:
                print("Could not find histogram '{}'!".format(objname))
                exit(1)
            result.append(h)
        else:
            result.append(None)
    return result

def LoadFDHistogram(file_name, spectra, scaling_factor):
    result = []
    for spectrum in spectra:
        jet_var_name = spectrum["variable_name"]
        if len(spectrum["jet_pt_cuts"]) > 1:
            jet_cuts = "_JetPt_{}_{}".format(spectrum["jet_pt_cuts"][0], spectrum["jet_pt_cuts"][1])
        elif len(spectrum["jet_pt_cuts"]) == 1:
            jet_cuts = "_JetPt_{}".format(spectrum["jet_pt_cuts"][0])
        else:
            jet_cuts = ""
        label = spectrum["variable_title"]
        if not jet_var_name == "DPt":
            spectrum_name = "{}DPtSpectrum{}".format(jet_var_name, jet_cuts)
        else:
            spectrum_name = "{}Spectrum{}".format(jet_var_name, jet_cuts)
        file = ROOT.TFile(file_name)
        if not file or file.IsZombie():
            print("Could not open input file {0}".format(file_name))
            exit(1)
        else:
            print("File {0} open".format(file_name))
        dlist = file.Get("D0_MCTruth")
        jlist = dlist.FindObject("Charged_R040")
        slist = jlist.FindObject("D0_MCTruth_Charged_R040_{}".format(spectrum_name))
        if not slist:
            print("Could not find list '{}' ({})!".format(spectrum_name, jet_var_name))
            jlist.Print()
            exit(1)
        jet_hist = slist.FindObject("D0_MCTruth_Charged_R040_{}".format(spectrum_name))
        if not jet_hist:
            print("Could not find FD histogram '{}' ({})!".format(spectrum_name, jet_var_name))
            slist.Print()
            exit(1)
        jet_hist.Scale(scaling_factor)
        jet_hist.GetZaxis().SetTitle("#frac{{d#sigma}}{{d{lab}}} #times #Delta{lab} (mb)".format(lab=label))
        print("Histogram {} loaded".format(jet_hist.GetName()))
        result.append(jet_hist)

    return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Prepare B feed-down correction file.')
    parser.add_argument('yaml', metavar='conf.yaml')
    parser.add_argument("--unfolding-debug", action='store_const',
                        default=False, const=True,
                        help='Unfolding debug plots.')

    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.unfolding_debug)

    IPython.embed()
