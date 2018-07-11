#!/usr/bin/env python
# python script to calculate average jet pt

import numpy
import IPython
import ROOT
import DMesonJetUtils

globalList = []

def GetMeasuredCrossSection(input_path, file_name):
    fname = "{}/{}.root".format(input_path, file_name)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    hStat = DMesonJetUtils.GetObject(file, "FinalSpectrum/CentralPointsStatisticalUncertainty")
    if not hStat:
        print("Cannot get measured cross section with statistical uncertainty!")
        exit(1)
    return hStat

def Calculate(dataStat):
    bins = [5, 6, 8, 10, 15, 20, 30]
    new_axis = ROOT.TAxis(len(bins) - 1, numpy.array(bins, numpy.float))
    h_new, results = DMesonJetUtils.FitAndRebin(dataStat, new_axis)
    globalList.extend(results)
    h = h_new.Clone()
    for ibin in range(1, h.GetNbinsX() + 1):
        h.SetBinContent(ibin, h.GetBinContent(ibin) * h.GetXaxis().GetBinWidth(ibin))
    h.GetXaxis().SetRangeUser(5, 30)
    print("The mean in the range 5,30 is {} +/- {}".format(h.GetMean(), h.GetMeanError()))
    h.GetXaxis().SetRangeUser(5, 15)
    print("The mean in the range 5,15 is {} +/- {}".format(h.GetMean(), h.GetMeanError()))
    h.GetXaxis().SetRangeUser(15, 30)
    print("The mean in the range 15,30 is {} +/- {}".format(h.GetMean(), h.GetMeanError()))

def main():
    dataStat = GetMeasuredCrossSection("/Volumes/DATA/ALICE/JetResults", "JetPtCrossSection_DPt_30_Systematics")
    Calculate(dataStat)

if __name__ == '__main__':
    main()

    IPython.embed()
