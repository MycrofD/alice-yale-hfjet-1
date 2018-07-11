# Load inclusive jet spectrum

import math
import numpy
import ROOT
import DMesonJetUtils

globalList = []

def GetCrossSection(jetptbins=[5, 6, 8, 10, 14, 20, 30]):
    sigmapp = 73.2  # published x-section in mbarn
    triggerEff = 62.2 / 73.2  # eff of MB OR trigger to be applied to data
    bin0CorrData = 0.91  # ratio event after vertexNcontributors cut to evts after phys selection

    fname = "../obusch/outData_spec_Bayes_combPtH_rebinned.root"
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {}".format(fname))
        exit(1)
    hStat_old = file.Get("hSpecComb")
    if not hStat_old:
        print("Cannot get inclusive cross section with statistical uncertainty!")
        exit(1)
    file.Close()

    hStat_old.Scale(sigmapp * triggerEff)
    hStat_old.Scale(bin0CorrData)

    fname = "../obusch/outData_spec_Bayes_combPtH.root"
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {}".format(fname))
        exit(1)
    hStat_orig = file.Get("hSpecComb")
    if not hStat_orig:
        print("Cannot get inclusive cross section with statistical uncertainty (not rebinned)!")
        exit(1)
    file.Close()

    hStat_orig.Scale(sigmapp * triggerEff)
    hStat_orig.Scale(bin0CorrData)

    fname = "../obusch/outSys_tot_spec.root"
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {}".format(fname))
        exit(1)
    hSyst_orig = file.Get("grSysErrTotNoUE")
    if not hSyst_orig:
        print("Cannot get inclusive cross section with systematic uncertainty!")
        exit(1)
    file.Close()

    if isinstance(jetptbins, list):
        new_axis = ROOT.TAxis(len(jetptbins) - 1, numpy.array(jetptbins, dtype=numpy.float))
        axis_comparison_old = DMesonJetUtils.AxisCompare.CheckConsistency(new_axis, hStat_old.GetXaxis())
        axis_comparison_orig = DMesonJetUtils.AxisCompare.CheckConsistency(new_axis, hStat_orig.GetXaxis())

        print("Binning of the original spectrum")
        binning = "["
        for ibin in range(1, hStat_orig.GetXaxis().GetNbins() + 1):
            binning += "{}, ".format(hStat_orig.GetXaxis().GetBinLowEdge(ibin))
        binning += "{}]".format(hStat_orig.GetXaxis().GetBinUpEdge(hStat_orig.GetXaxis().GetNbins()))
        print(binning)

        print("Binning of the rebinned spectrum")
        binning = "["
        for ibin in range(1, hStat_old.GetXaxis().GetNbins() + 1):
            binning += "{}, ".format(hStat_old.GetXaxis().GetBinLowEdge(ibin))
        binning += "{}]".format(hStat_old.GetXaxis().GetBinUpEdge(hStat_old.GetXaxis().GetNbins()))
        print(binning)

        print("Requested binning")
        binning = "["
        for ibin in range(1, new_axis.GetNbins() + 1):
            binning += "{}, ".format(new_axis.GetBinLowEdge(ibin))
        binning += "{}]".format(new_axis.GetBinUpEdge(new_axis.GetNbins()))
        print(binning)

        print("The relation between the original and requested binning is '{}'".format(axis_comparison_orig))
        print("The relation between the rebinned and requested binning is '{}'".format(axis_comparison_old))

        if axis_comparison_orig == DMesonJetUtils.AxisCompare.Identical:
            return_original = True

        elif axis_comparison_old == DMesonJetUtils.AxisCompare.Identical:
            return_original = False
            hStat = hStat_old.Clone("{}_copy".format(hStat_old.GetName()))
            
        elif axis_comparison_orig == DMesonJetUtils.AxisCompare.IsContainedSameBinning:
            return_original = True
            hStat_orig = DMesonJetUtils.Rebin1D_fromBins(hStat_orig, "{}_rebinned".format(hStat_orig.GetName()), len(jetptbins) - 1, numpy.array(jetptbins, dtype=numpy.float))
            points_to_be_removed = set()
            for ipoint in range(0, hSyst_orig.GetN()):
                if hSyst_orig.GetX()[ipoint] < jetptbins[0] or hSyst_orig.GetX()[ipoint] > jetptbins[-1]: points_to_be_removed.add(ipoint)
            for ipoint in points_to_be_removed:
                hSyst_orig.RemovePoint(ipoint)

        elif axis_comparison_old == DMesonJetUtils.AxisCompare.IsContainedSameBinning:
            return_original = False
            hStat = DMesonJetUtils.Rebin1D_fromBins(hStat_old, "{}_rebinned".format(hStat_old.GetName()), len(jetptbins) - 1, numpy.array(jetptbins, dtype=numpy.float))

        elif axis_comparison_orig == DMesonJetUtils.AxisCompare.IsContained or axis_comparison_orig == DMesonJetUtils.AxisCompare.SameLimits:
            return_original = False
            print("I will attempt to fit the original histogram and use the requested binning.")
            hStat, fit_results = DMesonJetUtils.FitAndRebin(hStat_orig, new_axis)
            globalList.extend(fit_results)
        
        else:
            print("I am unable to convert binning of either original or rebinned spectrum to the requested binning.")
            return None, None
 
    elif jetptbins == "rebinned":
        return_original = False
        hStat = hStat_old.Clone("{}_copy".format(hStat_old.GetName()))
    elif jetptbins == "original":
        return_original = True
    else:
        print("Option '{}' unknown!".format(jetptbins))
        return None, None

    if return_original:
        return hStat_orig.Clone(), hSyst_orig.Clone()

    if not hStat:
        return None, None

    hSyst = ROOT.TGraphErrors(hStat.GetNbinsX())

    for i in range(1, hStat.GetNbinsX() + 1):
        x = hStat.GetXaxis().GetBinCenter(i)
        y = hStat.GetBinContent(i)
        xerr = hStat.GetXaxis().GetBinWidth(i) / 2
        print("hStat Bin {}, x in [{}, {}], y = {}".format(i, x - xerr, x + xerr, y))
        yerr_rel = 0.
        y_check = 0.
        for iorig in range(0, hSyst_orig.GetN()):
            if hSyst_orig.GetX()[iorig] < hStat.GetXaxis().GetBinLowEdge(i): continue
            if hSyst_orig.GetX()[iorig] > hStat.GetXaxis().GetBinUpEdge(i): break
            yerr_rel += hSyst_orig.GetErrorY(iorig) * hStat_orig.GetBinContent(iorig + 1) * hStat_orig.GetXaxis().GetBinWidth(iorig + 1)
            y_check += hStat_orig.GetBinContent(iorig + 1) * hStat_orig.GetXaxis().GetBinWidth(iorig + 1)
            print("hSyst_orig Bin {}, x = {}, err_y_rel = {}".format(iorig, hSyst_orig.GetX()[iorig], hSyst_orig.GetErrorY(iorig)))
            print("hStat_orig Bin {}, x = {}, y = {}, err_y_rel = {}".format(iorig, hStat_orig.GetXaxis().GetBinCenter(iorig + 1), hStat_orig.GetBinContent(iorig + 1), hStat_orig.GetBinError(iorig + 1)))
        yerr_rel /= y_check
        y_check /= hStat.GetXaxis().GetBinWidth(i)
        yerr = yerr_rel * y
        hSyst.SetPoint(i - 1, x, y)
        hSyst.SetPointError(i - 1, xerr, yerr)
        print("Bin {} Tot Error rel {}".format(i, yerr_rel))
        print("Bin {}, tot y {}, check {}\n".format(i, y, y_check))

    return hStat, hSyst
