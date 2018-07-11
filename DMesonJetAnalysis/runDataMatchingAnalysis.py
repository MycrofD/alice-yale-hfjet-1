#!/usr/bin/env python
# python script to run the D meson jet analysis on 2010 pp data

import argparse
import subprocess
import yaml
import IPython
import DMesonJetProjectors
import DMesonJetMatchingAnalysis
import ROOT
print("ROOT {} imported".format(ROOT.gROOT.GetVersionInt()))
ROOT.PyConfig.IgnoreCommandLineOptions = True

globalList = []

# To mimic ROOT5 behavior
if ROOT.gROOT.GetVersionInt() >= 60000: ROOT.ROOT.Math.IntegratorOneDimOptions.SetDefaultIntegrator("Gauss")

def main(config, maxEvents, fmt, gen, proc, ts, stage, ask, bg):
    if bg:
        ROOT.gROOT.SetBatch(True)

    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    subprocess.call("make")
    ROOT.gSystem.Load("MassFitter.so")

    if config["train"] == "FastSim" or config["train"] == "FastSimOld":
        suffix = "{0}_{1}_{2}".format(gen, proc, ts)

        name = "{0}_{1}".format(config["name"], suffix)

        input_path = "{0}/FastSim_{1}/".format(config["input_path"], suffix)

        if stage >= 0:
            input_path += "stage_{0}/output".format(stage)
            if config["train"] == "FastSimOld":
                file_name = "AnalysisResults_FastSim_{0}.root".format(suffix)
            else:
                file_name = "AnalysisResults_FastSim_{0}_{1}.root".format(gen, proc)
        else:
            input_path += "output"
            file_name = "AnalysisResults_FastSim_{0}_{1}.root".format(gen, proc)
        train = ""
        output_path = "{0}/FastSim_{1}/".format(config["input_path"], suffix)
    else:
        name = config["name"]
        train = config["train"]
        input_path = config["input_path"]
        file_name = config["file_name"]
        output_path = input_path

    if "normalization_factor" in config:
        norm_factor = config["normalization_factor"]
    else:
        norm_factor = 1

    if "tree_type" in config:
        tree_type = config["tree_type"]
    else:
        tree_type = "simple"

    ana = DMesonJetMatchingAnalysis.DMesonJetMatchingAnalysis(name)
    projector = DMesonJetProjectors.DMesonJetProjector(input_path, train, file_name, config["task_name"], tree_type, config["merging_type"], norm_factor, maxEvents)
    projector.fDoNotAsk = not ask
    ana.SetProjector(projector)

    ana.StartAnalysis(config)

    ana.SaveRootFile("{0}/{1}".format(output_path, train))
    ana.SavePlots("{0}/{1}".format(output_path, train), fmt)
    ana.SavePlots("{0}/{1}".format(output_path, train), "C")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='D meson jet analysis for 2010 pp data.')
    parser.add_argument('yaml', metavar='config.yaml',
                        help='YAML configuration file')
    parser.add_argument('--events', metavar='N',
                        default=-1, type=int)
    parser.add_argument('--format', metavar='pdf',
                        default="pdf")
    parser.add_argument('--gen', metavar='GEN',
                        default=None)
    parser.add_argument('--proc', metavar='PROC',
                        default=None)
    parser.add_argument('--ts', metavar='TS',
                        default=None)
    parser.add_argument('--stage', metavar='N',
                        default=-1, type=int)
    parser.add_argument('-a', metavar='a', help='Ask before starting each analysis cycle',
                        action='store_const', default=False, const=True)
    parser.add_argument('-b', action='store_const',
                        default=False, const=True)
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    yconfig = yaml.load(f)
    f.close()

    main(yconfig, args.events, args.format, args.gen, args.proc, args.ts, args.stage, args.a, args.b)

    IPython.embed()
