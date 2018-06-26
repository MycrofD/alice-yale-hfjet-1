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

    spectrumName = "JetZSpectrum"
    jetName = "Jet_AKTChargedR040_pt_scheme"
    dmesonName = "D0_D0toKpiCuts"
    prefix = "Prompt_{}_{}_{}".format(dmesonName, jetName, spectrumName)

    pt_lim = [(0.2, 0.4, 5, 15, "DPt_20_JetPt_5_15"), (0.6, 0.8, 5, 15, "DPt_20_JetPt_5_15"), (0.8, 1.0, 5, 15, "DPt_20_JetPt_5_15"),
              (0.6, 0.8, 15, 30, "DPt_60_JetPt_15_30"), (0.8, 1.0, 15, 30, "DPt_60_JetPt_15_30")]
    histos = []
    for (minZ, maxZ, minJetPt, maxJetPt, kincuts) in pt_lim:
        resolutionName = "{0}_{1}/DetectorResponse/{0}_{1}_DetectorResponse_{2}_{3}".format(prefix, kincuts, int(minZ * 10), int(maxZ * 10))
        h = DMesonJetUtils.GetObject(file, resolutionName)
        h.SetTitle("{} < #it{{z}}_{{||,D}}^{{part}} < {} and {} < #it{{p}}_{{T,ch jet}}^{{part}} < {}".format(minZ, maxZ, minJetPt, maxJetPt))
        globalList.append(h)
        histos.append(h)

    cname = "ResolutionVsZ_Paper"
    comp = DMesonJetCompare.DMesonJetCompare(cname)
    comp.fDoSpectraPlot = "lineary"
    comp.fDoRatioPlot = None
    comp.fDoSpectrumLegend = False
    comp.fLinUpperSpace = 0.50
    comp.fMarkerSize = 1.5
    comp.fColors = [ROOT.kBlue + 2, ROOT.kRed + 2, ROOT.kGreen + 2,
                    ROOT.kRed + 2, ROOT.kGreen + 2]
    comp.fMarkers = [ROOT.kFullCircle, ROOT.kFullSquare, ROOT.kFullDiamond,
                     ROOT.kOpenSquare, ROOT.kOpenDiamond]
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
    h.GetXaxis().SetTitleFont(43)
    h.GetXaxis().SetTitleSize(26)
    h.GetXaxis().SetLabelFont(43)
    h.GetXaxis().SetLabelSize(22)
    h.GetXaxis().SetLabelOffset(0.008)
    h.GetXaxis().SetRangeUser(-0.6, 0.6)
    h.GetXaxis().SetTitle("(#it{z}_{||, det}^{ch} #font[122]{-} #it{z}_{||, gen}^{ch}) / #it{z}_{||, gen}^{ch}")
    h.GetYaxis().SetTitleFont(43)
    h.GetYaxis().SetTitleSize(26)
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelSize(22)
    h.GetYaxis().SetTitleOffset(0.9)
    h.GetYaxis().SetRangeUser(0, 17)

    paveALICE = ROOT.TPaveText(0.14, 0.68, 0.53, 0.95, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(21)
    paveALICE.SetTextAlign(13)
    paveALICE.AddText("ALICE PYTHIA 6")
    paveALICE.AddText("pp, #sqrt{#it{s}} = 7 TeV")
    paveALICE.AddText("with D^{0} #rightarrow K^{#font[122]{-}}#pi^{+}")
    paveALICE.AddText("and charge conj.")
    paveALICE.Draw()

    paveALICE = ROOT.TPaveText(0.61, 0.75, 0.90, 0.95, "NB NDC")
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

    leg = ROOT.TLegend(0.61, 0.29, 0.86, 0.66)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(43)
    leg.SetTextSize(19)
    leg.SetMargin(0.15)
    entry = leg.AddEntry(None, "5 < #it{p}_{T,gen jet}^{ch} < 15 GeV/#it{c}", "")
    entry = leg.AddEntry(None, "#it{p}_{T,D} > 2 GeV/#it{c}", "")
    entry = leg.AddEntry(None, "0.2 < #it{z}_{||, gen}^{ch} < 0.4", "p")
    entry.SetLineColor(ROOT.kBlue + 2)
    entry.SetMarkerColor(ROOT.kBlue + 2)
    entry.SetMarkerStyle(ROOT.kFullCircle)
    entry = leg.AddEntry(None, "0.6 < #it{z}_{||, gen}^{ch} < 0.8", "p")
    entry.SetLineColor(ROOT.kRed + 2)
    entry.SetMarkerColor(ROOT.kRed + 2)
    entry.SetMarkerStyle(ROOT.kFullSquare)
    entry = leg.AddEntry(None, "0.8 < #it{z}_{||, gen}^{ch} < 1.0", "p")
    entry.SetLineColor(ROOT.kGreen + 2)
    entry.SetMarkerColor(ROOT.kGreen + 2)
    entry.SetMarkerStyle(ROOT.kFullDiamond)
    leg.Draw()
    globalList.append(leg)

    leg = ROOT.TLegend(0.15, 0.29, 0.43, 0.57)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(43)
    leg.SetTextSize(19)
    leg.SetMargin(0.15)
    entry = leg.AddEntry(None, "15 < #it{p}_{T,gen jet}^{ch} < 30 GeV/#it{c}", "")
    entry = leg.AddEntry(None, "#it{p}_{T,D} > 6 GeV/#it{c}", "")
    entry = leg.AddEntry(None, "0.6 < #it{z}_{||, gen}^{ch} < 0.8", "p")
    entry.SetLineColor(ROOT.kRed + 2)
    entry.SetMarkerColor(ROOT.kRed + 2)
    entry.SetMarkerStyle(ROOT.kOpenSquare)
    entry = leg.AddEntry(None, "0.8 < #it{z}_{||, gen}^{ch} < 1.0", "p")
    entry.SetLineColor(ROOT.kGreen + 2)
    entry.SetMarkerColor(ROOT.kGreen + 2)
    entry.SetMarkerStyle(ROOT.kOpenDiamond)
    leg.Draw()
    globalList.append(leg)

    return canvas


def main():
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    f = open("LHC15i2response_Train1399_efficiency.yaml", 'r')
    config = yaml.load(f)
    f.close()

    canvas = ResolutionComparison(config)
    canvas.SaveAs("{0}/ResolutionVsZ_Paper.pdf".format(config["input_path"]))
    canvas.SaveAs("{0}/ResolutionVsZ_Paper.C".format(config["input_path"]))
    canvas.SaveAs("{0}/ResolutionVsZ_Paper.eps".format(config["input_path"]))

if __name__ == '__main__':

    main()

    IPython.embed()
