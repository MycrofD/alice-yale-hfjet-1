#!/usr/bin/env python
# python utilities for D Meson jet analysis

import os
import ROOT
import math
import array
from collections import OrderedDict
from enum import Enum

class AxisCompare(Enum):
    Identical = 0
    SameLimits = 1
    IsContainedSameBinning = 2
    ContainsSameBinning = 3
    OverlapsSameBinning = 4
    IsContained = 5
    Contains = 6
    Overlaps = 7
    NoOverlap = 8

    @classmethod
    def CheckConsistency(cls, axis1, axis2):
        isContained = False
        contains = False
        overlaps = False
        if axis1.GetBinLowEdge(1) <= axis2.GetBinLowEdge(1):
            if axis1.GetBinUpEdge(axis1.GetNbins()) >= axis2.GetBinUpEdge(axis2.GetNbins()):
                contains = True
            else:
                overlaps = True
        if axis2.GetBinLowEdge(1) <= axis1.GetBinLowEdge(1):
            if axis2.GetBinUpEdge(axis2.GetNbins()) >= axis1.GetBinUpEdge(axis1.GetNbins()):
                isContained = True
            else:
                overlaps = True

        if not contains and not isContained and not overlaps:
            return cls.NoOverlap

        sameBinning = True
        for ibin1 in range(1, axis1.GetNbins() + 1):
            if axis1.GetBinLowEdge(ibin1) >= axis2.GetBinLowEdge(1): break
            ibin1 += 1

        for ibin2 in range(1, axis2.GetNbins() + 1):
            if axis2.GetBinLowEdge(ibin2) >= axis1.GetBinLowEdge(1): break
            ibin2 += 1

        while(ibin1 <= axis1.GetNbins() and ibin2 <= axis2.GetNbins()):
            if axis1.GetBinLowEdge(ibin1) != axis2.GetBinLowEdge(ibin2):
                sameBinning = False
                break
            if axis1.GetBinUpEdge(ibin1) != axis2.GetBinUpEdge(ibin2):
                sameBinning = False
                break
            ibin1 += 1
            ibin2 += 1

        if contains and isContained:
            if sameBinning:
                return cls.Identical
            else:
                return cls.SameLimits
        elif contains:
            if sameBinning:
                return cls.ContainsSameBinning
            else:
                return cls.Contains
        elif isContained:
            if sameBinning:
                return cls.IsContainedSameBinning
            else:
                return cls.IsContained
        else:
            if sameBinning:
                return cls.OverlapsSameBinning
            else:
                return cls.Overlaps


def soft_clone(origin, name, title=None, yaxisTitle=None):
    if not title: title = name
    if not yaxisTitle: yaxisTitle = origin.GetYaxis().GetTitle()
    h = ROOT.TH1D(name, title, origin.GetNbinsX(), origin.GetXaxis().GetXbins().GetArray())
    h.GetXaxis().SetTitle(origin.GetXaxis().GetTitle())
    h.GetYaxis().SetTitle(yaxisTitle)
    return h

def GetRelativeUncertaintyHistogram(h):
    h_unc = soft_clone(h, "{0}_unc".format(h.GetName()), "{0} Rel. Unc.".format(h.GetTitle()), "rel. unc.")
    for ibin in xrange(0, h.GetNbinsX() + 2):
        if h.GetBinContent(ibin) == 0: continue
        h_unc.SetBinContent(ibin, h.GetBinError(ibin) / h.GetBinContent(ibin))
    return h_unc

def ConvertDMesonName(dmeson):
    if "D0" in dmeson:
        return "D^{0} #rightarrow K^{-}#pi^{+} and c.c."
    else:
        return dmeson

def binom(n, k):
    return math.factorial(n) / math.factorial(n - k) / math.factorial(k)

def GetObject(obj, name):
    slash = 0
    while(slash >= 0):
        slash = name.find("/", slash)
        if name[slash + 1] == '/':
            slash += 2
        else:
            break

    if slash < 0:
        name_lookup = name.replace("//", "/")
        name = None
    else:
        name_lookup = name[:slash].replace("//", "/")
        name = name[slash + 1:]
    res = None
    name_obj = None
    if isinstance(obj, ROOT.TCollection):
        res = obj.FindObject(name_lookup)
        name_obj = obj.GetName()
    elif isinstance(obj, dict) or isinstance(obj, OrderedDict):
        res = obj[name_lookup]
        name_obj = obj
    elif isinstance(obj, ROOT.TDirectory):
        res = obj.Get(name_lookup)
        name_obj = obj.GetName()
    if not res:
        if name_obj:
            print("Could not find object {0} in collection '{1}'".format(name_lookup, name_obj))
            if isinstance(obj, ROOT.TObject): obj.ls()
        else:
            print("Could not find object {} in following collection".format(name_lookup))
            print(obj)
        return None
    if name:
        return GetObject(res, name)
    else:
        return res

def GetObjectAndMerge(fileList, name):
    res = None
    for f in fileList:
        obj = GetObject(f, name)
        if res:
            res.Add(obj)
        else:
            res = obj.Clone()
    return res

def GenerateMultiCanvas(cname, n):
    rows = int(math.floor(math.sqrt(n)))
    cols = int(math.ceil(float(n) / rows))
    c = ROOT.TCanvas(cname, cname, cols * 400, rows * 400)
    c.Divide(cols, rows)
    return c

def find_file(path, file_name):
    for root, dirs, files in os.walk(path):
        for file in files:
            if file == file_name:
                yield os.path.join(root, file)

def FindMinimum(histogram, limit=0., errors=True):
    m = None
    if histogram.GetDimension() == 1:
        for ibin in xrange(1, histogram.GetNbinsX() + 1):
            if errors:
                cont = histogram.GetBinContent(ibin) - histogram.GetBinError(ibin)
            else:
                cont = histogram.GetBinContent(ibin)
            if cont <= limit: continue
            if m is None or cont < m: m = cont
    elif histogram.GetDimension() == 2:
        for xbin in xrange(1, histogram.GetNbinsX() + 1):
            for ybin in xrange(1, histogram.GetNbinsY() + 1):
                if errors:
                    cont = histogram.GetBinContent(xbin, ybin) - histogram.GetBinError(xbin, ybin)
                else:
                    cont = histogram.GetBinContent(xbin, ybin)
                if cont <= limit: continue
                if m is None or cont < m: m = cont
    return m

def FindMaximum(histogram, limit=0., errors=True):
    m = None
    if histogram.GetDimension() == 1:
        for ibin in xrange(1, histogram.GetNbinsX() + 1):
            if errors:
                cont = histogram.GetBinContent(ibin) + histogram.GetBinError(ibin)
            else:
                cont = histogram.GetBinContent(ibin)
            if cont <= limit: continue
            if m is None or cont > m: m = cont
    elif histogram.GetDimension() == 2:
        for xbin in xrange(1, histogram.GetNbinsX() + 1):
            for ybin in xrange(1, histogram.GetNbinsY() + 1):
                if errors:
                    cont = histogram.GetBinContent(xbin, ybin) + histogram.GetBinError(xbin, ybin)
                else:
                    cont = histogram.GetBinContent(xbin, ybin)
                if cont <= limit: continue
                if m is None or cont > m: m = cont
    return m

def DivideNoErrors(ratio, den):
    if not ratio.GetNbinsX() == den.GetNbinsX():
        print("DMesonJetUtils.DivideNoErrors: histograms have different number of bins!")
        return False

    for ibin in xrange(0, ratio.GetNbinsX() + 2):
        if den.GetBinContent(ibin) == 0:
            ratio.SetBinContent(ibin, 0)
        else:
            ratio.SetBinContent(ibin, ratio.GetBinContent(ibin) / den.GetBinContent(ibin))

    return True

def V2TH1(vect):
    result = ROOT.TH1D("vect", "vect", len(vect) - 2, 1, len(vect) - 2)
    for ibin in xrange(0, result.GetNbinsX() + 2):
        result.SetBinContent(ibin, vect[ibin])
    return result

def BuildHistogram(axis, name, yaxis):
    if len(axis) == 1:
        hist = ROOT.TH1D(name, name, len(axis[0].fBins) - 1, array.array('d', axis[0].fBins))
        hist.GetXaxis().SetTitle(axis[0].GetTitle())
        hist.GetYaxis().SetTitle(yaxis)
        hist.Sumw2()
    elif len(axis) == 2:
        hist = ROOT.TH2D(name, name, len(axis[0].fBins) - 1, array.array('d', axis[0].fBins), len(axis[1].fBins) - 1, array.array('d', axis[1].fBins))
        hist.GetXaxis().SetTitle(axis[0].GetTitle())
        hist.GetYaxis().SetTitle(axis[1].GetTitle())
        hist.GetZaxis().SetTitle(yaxis)
        hist.Sumw2()
    else:
        hist = ROOT.TH3D(name, name, len(axis[0].fBins) - 1, array.array('d', axis[0].fBins), len(axis[1].fBins) - 1, array.array('d', axis[1].fBins), len(axis[2].fBins) - 1, array.array('d', axis[2].fBins))
        hist.GetXaxis().SetTitle(axis[0].GetTitle())
        hist.GetYaxis().SetTitle(axis[1].GetTitle())
        hist.GetZaxis().SetTitle(axis[2].GetTitle())
        hist.Sumw2()
    return hist

def Rebin1D(hist, xaxis, warnings=False):
    return Rebin1D_fromBins(hist, hist.GetName(), xaxis.GetNbins(), xaxis.GetXbins().GetArray(), warnings)

def Rebin1D_fromBins(hist, name, nbinsX, binsX, warnings=False):
    axis = ROOT.TAxis(nbinsX, binsX)
    compareAxis = AxisCompare.CheckConsistency(hist.GetXaxis(), axis)
    if compareAxis == AxisCompare.Identical:
        r = hist.Clone(name)
        return r
    r = ROOT.TH1D(name, name, nbinsX, binsX)
    r.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
    r.GetYaxis().SetTitle(hist.GetYaxis().GetTitle())
    for xbin in xrange(0, hist.GetXaxis().GetNbins() + 2):
        xbinCenter = hist.GetXaxis().GetBinCenter(xbin)
        rxbin = r.GetXaxis().FindBin(xbinCenter)
        binValue = hist.GetBinContent(xbin) + r.GetBinContent(rxbin)
        binError = math.sqrt(hist.GetBinError(xbin) ** 2 + r.GetBinError(rxbin) ** 2)
        r.SetBinContent(rxbin, binValue)
        r.SetBinError(rxbin, binError)
        if binValue > 0:
            relErr = binError / binValue
            if relErr > 0.9 and warnings:
                print("Bin ({0}) has rel stat err = {1}. This is VERY dangerous!".format(xbin, relErr))
    return r

def Rebin2D(hist, xaxis, yaxis, warnings=False):
    return Rebin2D_fromBins(hist, hist.GetName(), xaxis.GetNbins(), xaxis.GetXbins().GetArray(), yaxis.GetNbins(), yaxis.GetXbins().GetArray(), warnings)

def Rebin2D_fromBins(hist, name, nbinsX, binsX, nbinsY, binsY, warnings=False):
    xaxis = ROOT.TAxis(nbinsX, binsX)
    xcompareAxis = AxisCompare.CheckConsistency(hist.GetXaxis(), xaxis)
    yaxis = ROOT.TAxis(nbinsY, binsY)
    ycompareAxis = AxisCompare.CheckConsistency(hist.GetYaxis(), yaxis)
    if xcompareAxis == AxisCompare.Identical and ycompareAxis == AxisCompare.Identical:
        r = hist.Clone(name)
        return r
    r = ROOT.TH2D(name, name, nbinsX, binsX, nbinsY, binsY)
    r.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
    r.GetYaxis().SetTitle(hist.GetYaxis().GetTitle())
    r.GetZaxis().SetTitle(hist.GetZaxis().GetTitle())
    r.Sumw2()
    hist_GetXaxis_GetBinCenter = hist.GetXaxis().GetBinCenter
    hist_GetYaxis_GetBinCenter = hist.GetYaxis().GetBinCenter
    r_GetXaxis_GetBinUpEdge = r.GetXaxis().GetBinUpEdge
    r_GetYaxis_GetBinUpEdge = r.GetYaxis().GetBinUpEdge
    hist_At = hist.At
    r_At = r.At
    r_SetAt = r.SetAt
    if hist.GetSumw2().GetSize() > 0:
        hist_GetSumw2_At = hist.GetSumw2().At
    else:
        hist_GetSumw2_At = hist.At
    r_GetSumw2_At = r.GetSumw2().At
    r_GetSumw2_SetAt = r.GetSumw2().SetAt
    hist_GetXaxis_GetNbins_2 = hist.GetXaxis().GetNbins() + 2
    hist_GetYaxis_GetNbins_2 = hist.GetYaxis().GetNbins() + 2
    hist_bin = 0
    rybin = 0
    r_AddBinContent = r.AddBinContent
    for ybin in xrange(0, hist_GetYaxis_GetNbins_2):
        ybinCenter = hist_GetYaxis_GetBinCenter(ybin)
        while (ybinCenter > r_GetYaxis_GetBinUpEdge(rybin) and rybin < nbinsY + 1): rybin += 1
        rxbin = 0
        for xbin in xrange(0, hist_GetXaxis_GetNbins_2):
            xbinCenter = hist_GetXaxis_GetBinCenter(xbin)
            while (xbinCenter > r_GetXaxis_GetBinUpEdge(rxbin) and rxbin < nbinsX + 1): rxbin += 1
            r_bin = rxbin + (nbinsX + 2) * rybin
            r_SetAt(r_At(r_bin) + hist_At(hist_bin), r_bin)
            r_GetSumw2_SetAt(r_GetSumw2_At(r_bin) + hist_GetSumw2_At(hist_bin), r_bin)
            hist_bin += 1

    r.SetEntries(hist.GetEntries())

    if warnings:
        for xbin in xrange(0, r.GetXaxis().GetNbins() + 2):
            for ybin in xrange(0, r.GetYaxis().GetNbins() + 2):
                binValue = r.GetBinContent(xbin, ybin)
                binError = r.GetBinError(xbin, ybin)
                if binValue > 0:
                    relErr = binError / binValue
                    if relErr > 0.9:
                        print("Bin ({0},{1}) has rel stat err = {2}. This is VERY dangerous!".format(xbin, ybin, relErr))
    return r

def frange(start, stop, step, closed=False):
    i = start
    if closed: stop += step
    while i < stop:
        yield i
        i += step

def FitAndRebin(hstat, new_axis):
    results = []
    xmin = new_axis.GetBinLowEdge(1) - 5
    if xmin < 5:
        xmin = 5
    xmax = new_axis.GetBinUpEdge(new_axis.GetNbins()) + 5
    fit_func = ROOT.TF1("myfit", "expo(0)+expo(2)", xmin, xmax)
    results.append(fit_func)
    cname = "Fit_{}".format(hstat.GetName())
    canvas = ROOT.TCanvas(cname, cname)
    canvas.SetLogy()
    results.append(canvas)
    hcopy = hstat.DrawCopy()
    hcopy.GetXaxis().SetRangeUser(new_axis.GetBinLowEdge(1), new_axis.GetBinUpEdge(new_axis.GetNbins()))
    results.append(hcopy)
    fit_func.SetParameter(0, 1)
    fit_func.SetParameter(1, -1)
    fit_func.SetParameter(2, -99999999)
    fit_func.SetParameter(3, 0)
    fit_func.SetParameter(0, 1 + math.log(hstat.GetBinContent(hstat.GetXaxis().FindBin(xmin)) / fit_func.Eval(xmin)))
    fit_func.SetParameter(2, 1)
    fit_func.SetParameter(3, -0.5)
    fit_func.SetParameter(2, 1 + math.log(hstat.GetBinContent(hstat.GetXaxis().FindBin((xmax - xmin) / 2)) / fit_func.Eval((xmax - xmin) / 2)))
    fitR = hcopy.Fit(fit_func, "S", "", xmin, xmax)
    fitOk = int(fitR)
    if fitOk != 0:
        print("The fit was unsuccessfull!")
        return None
    h_fit = ROOT.TH1D("{}_rebinned".format(hstat.GetName()), hstat.GetTitle(), new_axis.GetNbins(), new_axis.GetXbins().GetArray())
    results.append(h_fit)
    for ibin in range(1, h_fit.GetNbinsX() + 1):
        valErr = fit_func.IntegralError(h_fit.GetXaxis().GetBinLowEdge(ibin), h_fit.GetXaxis().GetBinUpEdge(ibin)) / (h_fit.GetXaxis().GetBinUpEdge(ibin) - h_fit.GetXaxis().GetBinLowEdge(ibin))
        val = fit_func.Integral(h_fit.GetXaxis().GetBinLowEdge(ibin), h_fit.GetXaxis().GetBinUpEdge(ibin)) / (h_fit.GetXaxis().GetBinUpEdge(ibin) - h_fit.GetXaxis().GetBinLowEdge(ibin))
        print("integral = {0:.5f}, central = {1:.5f}".format(val, fit_func.Eval((h_fit.GetXaxis().GetBinCenter(ibin)))))
        h_fit.SetBinContent(ibin, val)
        h_fit.SetBinError(ibin, valErr)
    h_fit.Draw("same")

    ratio = hcopy.Clone("ratio")
    results.append(ratio)
    ratio.Divide(fit_func)
    cname = "FitRatio_{}".format(hstat.GetName())
    canvas_ratio = ROOT.TCanvas(cname, cname)
    results.append(canvas_ratio)
    ratio.Draw()
    return h_fit, results