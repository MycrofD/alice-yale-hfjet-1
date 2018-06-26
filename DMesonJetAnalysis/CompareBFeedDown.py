#!/usr/bin/env python

import argparse
import yaml
import IPython
import ROOT
import DMesonJetCompare
import DMesonJetUtils
import MCSimulationSystematics

globalList = []

def main(configs, name):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    CompareFeedDown(configs, name)

def OpenFiles(configs):
    for c in configs:
        fname = "{}/{}.root".format(c["input_path"], c["name"])
        c["file"] = ROOT.TFile(fname, "read")
        c["default_list"] = c["file"].Get("default")

def CompareFeedDown(configs, name):
    spectrum_names = MCSimulationSystematics.GenerateSpectrumNames(configs[0]["spectra"], True)
    OpenFiles(configs)
    if not name: name = "Comparison"
    for sname,_ in spectrum_names:
        cname = ""
        histos = []
        print("Spectrum {}...".format(sname))
        for c in configs:
            print("Loading spectrum for '{}'".format(c["name"]))
            h = DMesonJetUtils.GetObject(c["default_list"], sname)
            if not h: continue
            if cname: cname = "{}_{}".format(cname, c["name"])
            else: cname = c["name"]
            h.SetTitle(c["name"])
            globalList.append(h)
            histos.append(h)

        cname = "{}_{}_{}".format(name, sname.replace("/", "_"), cname)
        if len(histos) < 2: continue
        comp = DMesonJetCompare.DMesonJetCompare(cname)
        comp.fOptRatio = "hist"
        comp.fX1LegRatio = 0.15
        comp.fX1LegSpectrum = 0.25
        comp.fLogUpperSpace = 5  # this factor will be used to adjust the y axis in log scale
        comp.fLogLowerSpace = 1.5  # this factor will be used to adjust the y axis in log scale
        comp.fLinUpperSpace = 0.45  # this factor will be used to adjust the y axis in linear scale
        comp.fLinLowerSpace = 0.15  # this factor will be used to adjust the y axis in linear scale
        comp.fGridyRatio = True
        if "JetZ" in sname:
            comp.fDoSpectraPlot = "lineary"
        else:
            comp.fDoSpectraPlot = "logy"
        r = comp.CompareSpectra(histos[0], histos[1:])
        for obj in r:
            if not obj in globalList:
                globalList.append(obj)

        comp.fCanvasSpectra.SaveAs("{0}/{1}.pdf".format(configs[0]["input_path"], comp.fCanvasSpectra.GetName()))
        comp.fCanvasRatio.SaveAs("{0}/{1}.pdf".format(configs[0]["input_path"], comp.fCanvasRatio.GetName()))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='B feed-down comparison')
    parser.add_argument('yaml', nargs='*',
                        help='List of YAML configuration files')
    parser.add_argument('--name', metavar='NAME',
                        default=None)
    args = parser.parse_args()

    yconfigs = []

    for fname in args.yaml:
        f = open(fname, 'r')
        c = yaml.load(f)
        f.close()
        yconfigs.append(c)

    main(yconfigs, args.name)

    IPython.embed()
