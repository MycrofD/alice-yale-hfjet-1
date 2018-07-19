#!/usr/bin/env python

from collections import OrderedDict
import ROOT
import DMesonJetUtils
import HistogramNormalizator

def GetD0JetTheoryCrossSectionAll(config, axis):
    if "scale" in config:
        scale_base = config["scale"]
    else:
        scale_base = 1.0

    if "data_minx" in config and "data_maxx" in config:
        data_minx = config["data_minx"]
        data_maxx = config["data_maxx"]
    else:
        data_minx = 0
        data_maxx = -1

    if "max_jet_pt" in config and "min_jet_pt" in config:
        max_jet_pt = config["max_jet_pt"]
        min_jet_pt = config["min_jet_pt"]
    else:
        max_jet_pt = 0
        min_jet_pt = -1

    for t in config["theory"]:
        if not t["active"]:
            continue

        if not "spectrum_name" in t:
            t["spectrum_name"] = config["theory_spectrum"]

        if not "normalization" in config:
            config["normalization"] = "CrossSection"

        if not "jet_type" in t:
            t["jet_type"] = config["jet_type"]

        if not "type" in t:
            t["type"] = "stat-only"
        
        scale = scale_base
        if "scale" in t and config["normalization"] != "Distribution":
            scale *= t["scale"]

        if t["type"] == "stat-only":
            if config["normalization"] == "Distribution" and data_maxx < data_minx:
                print("LoadTheoryCrossSection.GetD0JetTheoryCrossSectionAll: requsted normalization, but did not provide minx and maxx via the fields 'data_minx' and 'data_maxx'. Aborting.")
                exit(1)
            
            if config["normalization"] == "Ratio" or config["normalization"] == "Rate":
                if not "inclusive" in t or not "histogram" in t["inclusive"]:
                    print("LoadTheoryCrossSection.GetD0JetTheoryCrossSectionAll: requested {} normalization but did not load inclusive jet spectrum to normalize.".format(config["normalization"]))
                    exit(1)
                if max_jet_pt < min_jet_pt:
                    print("LoadTheoryCrossSection.GetD0JetTheoryCrossSectionAll: requested {} normalization but did not load jet pt limits.".format(config["normalization"]))
                    exit(1)

            h = GetD0JetTheoryCrossSection(config, t, scale, axis, data_minx, data_maxx, min_jet_pt, max_jet_pt)
            t["histogram"] = h

        elif t["type"] == "stat+syst":
            file_name = t["file_name"]
            hStat, hSyst = GetD0JetTheoryCrossSectionStatSyst(config["input_path"], file_name, scale, t["spectrum_name"], config["normalization"])
            t["histogram"] = hStat
            t["systematics"] = hSyst
            if "ratio" in t and t["ratio"]:
                hStat, hSyst = GetD0JetTheoryCrossSectionStatSyst(config["input_path"], file_name, scale, t["spectrum_name"], "Ratio")
                t["ratio_histogram"] = hStat
                t["ratio_systematics"] = hSyst

def GetD0JetTheoryCrossSectionStatSyst(input_path, file_name, scale, spectrum, normalization):
    fname = "{}/{}".format(input_path, file_name)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    underscore = spectrum.find("_")
    spectrum_name_prefix = spectrum[:underscore]
    spectrum_name_suffix = spectrum[underscore + 1:]
    if "DPt_3" in spectrum_name_suffix and not "DPt_30" in spectrum_name_suffix:
        spectrum_name_suffix = spectrum_name_suffix.replace("DPt_3", "DPt_30")
    if "DPt_2" in spectrum_name_suffix and not "DPt_20" in spectrum_name_suffix:
        spectrum_name_suffix = spectrum_name_suffix.replace("DPt_2", "DPt_20")
    if "DPt_6" in spectrum_name_suffix and not "DPt_60" in spectrum_name_suffix:
        spectrum_name_suffix = spectrum_name_suffix.replace("DPt_6", "DPt_60")
    hStat = DMesonJetUtils.GetObject(file, "default/{prefix}_{suffix}_{normalization}/GeneratorLevel_{prefix}".format(prefix=spectrum_name_prefix, suffix=spectrum_name_suffix, normalization=normalization))
    if not hStat:
        print("Cannot get theory cross section with statistical uncertainty!")
        exit(1)
    hSyst = DMesonJetUtils.GetObject(file, "SystematicUncertainty/{prefix}_{suffix}_{normalization}//GeneratorLevel_{prefix}/GeneratorLevel_{prefix}_CentralAsymmSyst".format(prefix=spectrum_name_prefix, suffix=spectrum_name_suffix, normalization=normalization))
    if not hSyst:
        print("Cannot get theory cross section lower systematic uncertainty!")
        exit(1)

    hStat.Scale(scale)
    for ipoint in range(0, hSyst.GetN()):
        hSyst.SetPointEYlow(ipoint, hSyst.GetErrorYlow(ipoint) * scale)
        hSyst.SetPointEYhigh(ipoint, hSyst.GetErrorYhigh(ipoint) * scale)
        hSyst.SetPoint(ipoint, hSyst.GetX()[ipoint], hSyst.GetY()[ipoint] * scale)

    return hStat, hSyst

def GetD0JetTheoryCrossSection(config, t, scale, axis, data_minx, data_maxx, min_jet_pt, max_jet_pt):
    input_path = config["input_path"]
    gen = t["gen"]
    proc = t["proc"]
    ts = t["ts"]
    spectrum = t["spectrum_name"]
    jet_type = t["jet_type"]
    normalization = config["normalization"]
    fname = "{input_path}/FastSim_{gen}_{proc}_{ts}/FastSimAnalysis_ccbar_{gen}_{proc}_{ts}.root".format(input_path=input_path, gen=gen, proc=proc, ts=ts)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    h_orig = DMesonJetUtils.GetObject(file, "D0_MCTruth/{jet_type}/D0_MCTruth_{jet_type}_{spectrum}/D0_MCTruth_{jet_type}_{spectrum}".format(spectrum=spectrum, jet_type=jet_type))
    if not h_orig:
        print("Cannot get theory cross section with statistical uncertainty!")
        exit(1)

    h = DMesonJetUtils.Rebin1D(h_orig, axis)
    normalizator = HistogramNormalizator.MCSimulationNormalizator(h, normalization)
    normalizator.fScale = scale
    normalizator.fXmin = data_minx
    normalizator.fXmax = data_maxx
    normalizator.fNormalizationXmin = min_jet_pt
    normalizator.fNormalizationXmax = max_jet_pt
    if normalizator.fNormalizationType == HistogramNormalizator.NormalizationType.Rate or normalizator.fNormalizationType == HistogramNormalizator.NormalizationType.Ratio:
        normalizator.fNormalizationHistogram = t["inclusive"]["histogram"]
    h_normalized = normalizator.NormalizeHistogram()

    return h_normalized

def GetInclusiveJetTheoryCrossSectionAll(config):
    inclusive_jet_cross_sections = OrderedDict()
    if "inclusive_theory_spectrum" in config:
        spectrum_name = config["inclusive_theory_spectrum"]
    else:
        spectrum_name = "JetPt"
    for t in config["theory"]:
        if not t["active"]:
            continue
        if not "inclusive" in t or not t["inclusive"]:
            continue
        if "jet_type" in t["inclusive"]:
            jet_type = t["inclusive"]["jet_type"]
        elif "jet_type" in t:
            jet_type = t["jet_type"]
        else:
            jet_type = config["jet_type"]
        if "scale" in t["inclusive"]:
            scale = t["inclusive"]["scale"]
        else:
            scale = None
        if not "color" in t["inclusive"] and "color" in t:
            t["inclusive"]["color"] = t["color"]
        if not "title" in t["inclusive"] and "title" in t:
            t["inclusive"]["title"] = t["title"]
        if not "line" in t["inclusive"] and "line" in t:
            t["inclusive"]["line"] = t["line"]

        jet_key = (t["inclusive"]["ts"], jet_type)

        if jet_key in inclusive_jet_cross_sections:
            t["inclusive"] = inclusive_jet_cross_sections[jet_key]
            continue

        h = GetInclusiveJetTheoryCrossSection(config["input_path"], t["inclusive"]["gen"], t["inclusive"]["proc"], t["inclusive"]["ts"], scale, t["inclusive"]["file_name"], jet_type, spectrum_name)
        t["inclusive"]["histogram"] = h
        inclusive_jet_cross_sections[jet_key] = t["inclusive"]
        if h:
            print("Histogram '{}' loaded".format(jet_key))
        else:
            print("Could not load histogram '{}'!".format(jet_key))
    return inclusive_jet_cross_sections

def GetInclusiveJetTheoryCrossSection(input_path, gen, proc, ts, scale, file_name, jet_type, spectrum_name):
    fname = "{input_path}/FastSim_{gen}_{proc}_{ts}/{file_name}".format(input_path=input_path, gen=gen, proc=proc, ts=ts, file_name=file_name)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    if jet_type:
        objname = "{jet_type}/{spectrum_name}".format(jet_type=jet_type, spectrum_name=spectrum_name)
    else:
        objname = spectrum_name
    h_orig = DMesonJetUtils.GetObject(file, objname)
    if not h_orig:
        print("Cannot get theory cross section with statistical uncertainty!")
        exit(1)

    h = h_orig.Clone()

    if scale:
        h.Scale(scale)

    return h
