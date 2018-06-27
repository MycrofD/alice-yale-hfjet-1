#!/usr/bin/env python

import ROOT
import DMesonJetUtils
import HistogramNormalizator

def GetD0JetTheoryCrossSectionAll(config, axis):
    if "scale" in config:
        scale = config["scale"]
    else:
        scale = 1.0
    for t in config["theory"]:
        if not t["active"]: continue
        if "spectrum_name" in t:
            spectrum_name = t["spectrum_name"]
        else:
            spectrum_name = config["theory_spectrum"]
        if "jet_type" in t:
            jet_type = t["jet_type"]
        else:
            jet_type = config["jet_type"]
        if "normalization" in config:
            normalization = config["normalization"]
        else:
            normalization = "CrossSection"
        if "scale" in t and normalization != "Distribution":
            scale *= t["scale"]

        if "data_minx" in config and "data_maxx" in config:
            data_minx = config["data_minx"]
            data_maxx = config["data_maxx"]
        else:
            data_minx = 0
            data_maxx = -1
            
        if normalization == "Distribution" and data_maxx < data_minx:
            print("Error: requsted normalization, but did not provide minx and maxx via the fields 'data_minx' and 'data_maxx'. Aborting.")
            exit(1)

        h = GetD0JetTheoryCrossSection(config["input_path"], t["gen"], t["proc"], t["ts"], scale, spectrum_name, jet_type, axis, normalization, data_minx, data_maxx)
        t["histogram"] = h

def GetD0JetTheoryCrossSection(input_path, gen, proc, ts, scale, spectrum, jet_type, axis, normalization, data_minx, data_maxx):
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
    h_normalized = normalizator.NormalizeHistogram()

    return h_normalized

def GetInclusiveJetTheoryCrossSectionAll(config):
    inclusive_jet_cross_sections = dict()
    if "inclusive_theory_spectrum" in config:
        spectrum_name = config["inclusive_theory_spectrum"]
    else:
        spectrum_name = "JetPt"
    for t in config["theory"]:
        if not t["active"]: continue
        if not t["inclusive"]: continue
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
