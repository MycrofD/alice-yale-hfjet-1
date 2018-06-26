#!/usr/bin/env python

import yaml
import IPython
import ROOT
import DMesonJetUtils
import DMesonJetCompare

globalList = []

def ResolutionComparison(config):
    fname = "{0}/{1}/{2}.root".format(config["input_path"], config["train"], config["name"])
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)

    spectrumName = "JetPtSpectrum_DPt_30"
    jetName = "Jet_AKTChargedR040_pt_scheme"
    dmesonName = "D0_D0toKpiCuts"
    prefix = "Prompt_{}_{}_{}".format(dmesonName, jetName, spectrumName)

    pt_lim = [(5, 6), (8, 10), (20, 30)]
    histos = []
    for (minJetPt, maxJetPt) in pt_lim:
        resolutionName = "{0}/DetectorResponse/{0}_DetectorResponse_{1}_{2}".format(prefix, minJetPt * 10, maxJetPt * 10)
        h = DMesonJetUtils.GetObject(file, resolutionName)
        h.SetTitle("{} < #it{{p}}_{{T,gen jet}}^{{ch}} < {} GeV/#it{{c}}".format(minJetPt, maxJetPt))
        globalList.append(h)
        histos.append(h)

    cname = "ResolutionVsJetPt_Paper"
    comp = DMesonJetCompare.DMesonJetCompare(cname)
    comp.fDoSpectraPlot = "lineary"
    comp.fDoRatioPlot = None
    comp.fMarkerSize = 1.5
    comp.fX1LegSpectrum = 0.14
    comp.fX2LegSpectrum = 0.41
    comp.fY1LegSpectrum = 0.53
    comp.fLinUpperSpace = 0.50
    comp.fLegLineHeight = 0.075
    comp.fColors = [ROOT.kRed + 2, ROOT.kBlue + 2, ROOT.kGreen + 2]
    comp.fMarkers = [ROOT.kOpenCircle, ROOT.kOpenSquare, ROOT.kOpenDiamond]
    r = comp.CompareSpectra(histos[0], histos[1:])
    for obj in r:
        globalList.append(obj)

    canvas = comp.fCanvasSpectra
    canvas.SetTicks(1, 1)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.05)
    canvas.SetTopMargin(0.05)
    canvas.SetBottomMargin(0.15)
    canvas.cd()

    h = comp.fMainHistogram

    h.GetYaxis().SetTitle("Probability Density")
    h.GetXaxis().SetTitle("(#it{p}_{T,det jet}^{ch} #font[122]{-} #it{p}_{T,gen jet}^{ch}) / #it{p}_{T,gen jet}^{ch}")
    h.GetXaxis().SetTitleFont(43)
    h.GetXaxis().SetTitleSize(26)
    h.GetXaxis().SetLabelFont(43)
    h.GetXaxis().SetLabelSize(22)
    h.GetXaxis().SetLabelOffset(0.008)
    h.GetXaxis().SetRangeUser(-0.6, 0.6)
    h.GetYaxis().SetTitleFont(43)
    h.GetYaxis().SetTitleSize(26)
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelSize(22)
    h.GetYaxis().SetTitleOffset(0.9)
    h.GetYaxis().SetRangeUser(0, 14)

    paveALICE = ROOT.TPaveText(0.14, 0.62, 0.53, 0.95, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(21)
    paveALICE.SetTextAlign(13)
    paveALICE.AddText("ALICE PYTHIA 6")
    paveALICE.AddText("pp, #sqrt{#it{s}} = 7 TeV")
    paveALICE.AddText("Prompt D^{0} #rightarrow K^{#font[122]{-}}#pi^{+}")
    paveALICE.AddText("and charge conj.")
    paveALICE.AddText("#it{p}_{T,D} > 3 GeV/#it{c}")
    paveALICE.Draw()
    
    paveALICE = ROOT.TPaveText(0.65, 0.75, 0.90, 0.95, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(21)
    paveALICE.SetTextAlign(13)
    paveALICE.AddText("Charged Jets")
    paveALICE.AddText("Anti-#it{k}_{T}, #it{R} = 0.4")
    paveALICE.AddText("|#it{#eta}_{jet}| < 0.5")
    paveALICE.Draw()

    return canvas

def main():
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    f = open("LHC15i2response_Train1399_efficiency.yaml", 'r')
    config = yaml.load(f)
    f.close()

    canvas = ResolutionComparison(config)
    canvas.SaveAs("{0}/ResolutionVsJetPt_Paper.pdf".format(config["input_path"]))
    canvas.SaveAs("{0}/ResolutionVsJetPt_Paper.C".format(config["input_path"]))


if __name__ == '__main__':

    main()

    IPython.embed()
