#!/usr/bin/env python
# python script to compare inclusive jet spectra
# us it with
# InclusiveJetCompareTheory.yaml
# InclusiveJetComparePowhegKTCutSuppFact20.yaml
# InclusiveJetComparePowhegKTCutSuppFact60.yaml
# InclusiveJetComparePowhegKTCutSuppFact90.yaml
# InclusiveJetComparePowhegSuppFactKT1.yaml
# InclusiveJetComparePowhegSuppFactKT10.yaml
# InclusiveJetComparePowhegSuppFactKT20.yaml

import argparse
import numpy
import yaml
import IPython
import ROOT
import DMesonJetCompare
import DMesonJetUtils
import LoadInclusiveJetSpectrum

globalList = []

def LoadHistograms(config, ts, file_name, prefix, title, jet_type, scale):
    result = dict()
    full_file_name = "{}/{}_{}/{}".format(config["input_path"], prefix, ts, file_name)
    myfile = ROOT.TFile(full_file_name, "read")
    for hdef in config["histograms"]:
        if jet_type:
            hname = "{}/{}".format(jet_type, hdef["name"])
        else:
            hname = hdef["name"]
            jet_type = "Charged_R040"
        h = DMesonJetUtils.GetObject(myfile, hname)
        if not h: 
            continue
        if "bins" in hdef:
            h_orig = h
            h = DMesonJetUtils.Rebin1D_fromBins(h_orig, "{}_rebinned".format(h_orig.GetName()), len(hdef["bins"])-1, numpy.array(hdef["bins"], dtype=numpy.float))
        if "max" in hdef and "min" in hdef:
            h.GetXaxis().SetRangeUser(hdef["min"], hdef["max"])
        if scale:
            h.Scale(scale)
        h.SetTitle(title)
        hname_result = "{}/{}".format(jet_type, hdef["name"])
        result[hname_result] = h.Clone(hname_result)
        globalList.append(result[hname_result])
    return result

def main(config):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    measured_inclusive_cross_section, _ = LoadInclusiveJetSpectrum.GetCrossSection("original")

    histograms = []
    for element in config["list"]:
        if not element["active"]: 
            continue
        if "file_name" in element:
            file_name = element["file_name"]
        else:
            file_name = config["file_name"]
        if "prefix" in element:
            prefix = element["prefix"]
        else:
            prefix = config["prefix"]
        if "jet_type" in element:
            jet_type = element["jet_type"]
        else:
            jet_type = config["jet_type"]
        if "scale" in element:
            scale = element["scale"]
        else:
            scale = 1.0
        histograms.append(LoadHistograms(config, element["ts"], file_name, prefix, element["title"], jet_type, scale))
    for hdef in config["histograms"]:
        if "jet_type" in hdef:
            jet_type = hdef["jet_type"]
        else:
            jet_type = config["jet_type"]
        hname = "{}/{}".format(jet_type, hdef["name"])
        cname = "{}_{}_{}".format(config["name"], jet_type, hdef["name"])
        if "suffix" in hdef:
            cname += hdef["suffix"]
        histo_to_compare = [hlist[hname] for hlist in histograms if hname in hlist]
        globalList.extend(histo_to_compare)
        comp = DMesonJetCompare.DMesonJetCompare(cname)
        comp.fX1LegRatio = 0.40
        comp.fX2LegRatio = 0.90
        comp.fX1LegSpectrum = 0.40
        comp.fX2LegSpectrum = 0.90
        comp.fLogUpperSpace = 50
        comp.fLinUpperSpace = 0.4

        if measured_inclusive_cross_section and hdef["measured"]:
            if not histo_to_compare or len(histo_to_compare) < 1:
                print("No histograms to compare!")
                continue
            if "bins" in hdef:
                measured_inclusive_cross_section_copy = DMesonJetUtils.Rebin1D_fromBins(measured_inclusive_cross_section, "measured_inclusive_cross_section_copy", len(hdef["bins"])-1, numpy.array(hdef["bins"], dtype=numpy.float))
            else:
                measured_inclusive_cross_section_copy = measured_inclusive_cross_section.Clone("measured_inclusive_cross_section_copy")
            measured_inclusive_cross_section_copy.GetXaxis().SetTitle(histo_to_compare[0].GetXaxis().GetTitle())
            measured_inclusive_cross_section_copy.GetYaxis().SetTitle(histo_to_compare[0].GetYaxis().GetTitle())
            globalList.append(measured_inclusive_cross_section_copy)
            r = comp.CompareSpectra(measured_inclusive_cross_section_copy, histo_to_compare)
        else:
            if len(histo_to_compare) <= 1:
                print("No histograms to compare!")
                continue
            r = comp.CompareSpectra(histo_to_compare[0], histo_to_compare[1:])
        if "min_ratio" in hdef:
            comp.fMainRatioHistogram.SetMinimum(hdef["min_ratio"])
        if "max_ratio" in hdef:
            comp.fMainRatioHistogram.SetMaximum(hdef["max_ratio"])

        for obj in r:
            globalList.append(obj)
            if isinstance(obj, ROOT.TCanvas):
                obj.SaveAs("{}/{}.pdf".format(config["input_path"], obj.GetName()))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compare inclusive jet spectra.')
    parser.add_argument('config',
                        help='YAML file')
    args = parser.parse_args()

    f = open(args.config, 'r')
    yconfig = yaml.load(f)
    f.close()

    main(yconfig)

    IPython.embed()