#!/usr/bin/env python
# python program to do extract B feed down correction factors

import ROOT

class DMesonJetFDCorrection:

    def __init__(self, config, spectrumName, inputPath=None, dmeson=None, jetType=None, jetRadius=None):
        self.fInputPath = inputPath
        self.fSpectrumName = spectrumName
        self.fFDHistogramVariations = dict()
        # For the moment these parameters are not used
        #####
        self.fDMeson = dmeson
        self.fJetType = jetType
        self.fJetRadius = jetRadius
        #####
        if config:
            self.fFileName = config["file_name"]
            self.fCentralPoints = config["central_points"]
            self.fFDSpectrumName = config["spectrum"]
            self.LoadFD()
        else:
            self.fFileName = None
            self.fCentralPoints = None
            self.fFDSpectrumName = None
            self.fFDHistogram = None
            self.fFDUpSystUncHistogram = None
            self.fFDLowSystUncHistogram = None
            self.fFDTotUpSystUncHistogram = None
            self.fFDTotLowSystUncHistogram = None
            self.fFDSystUncGraph = None
            self.fFDTotSystUncGraph = None

    def GetFDHistogram(self, variation_name):
        if variation_name == self.fCentralPoints: return self.fFDHistogram
        if variation_name in self.fFDHistogramVariations: return self.fFDHistogramVariations[variation_name]
        h = self.LoadVariation(variation_name)
        self.fFDHistogramVariations[variation_name] = h
        return h

    def LoadVariation(self, variation_name):
        filename = "{0}/{1}".format(self.fInputPath, self.fFileName)
        file = ROOT.TFile(filename)
        if not file or file.IsZombie():
            print("Could not open file '{0}'".format(filename))
            return
        else:
            print("File '{0}' open for FD correction".format(filename))
        rlist = file.Get(variation_name)
        if not rlist:
            print("Could not find list '{0}'".format(variation_name))
            file.ls()
            return
        else:
            print("List '{0}' loaded".format(variation_name))
        rlist2 = rlist.FindObject("{}_CrossSection".format(self.fSpectrumName))
        if not rlist2:
            print("Could not find list '{0}'".format(self.fSpectrumName))
            rlist.Print()
            return
        else:
            print("List '{0}' loaded".format(self.fSpectrumName))
        h = rlist2.FindObject(self.fFDSpectrumName)
        if not h:
            print("Could not find histogram '{0}'".format(self.fFDSpectrumName))
            rlist2.Print()
            return
        else:
            print("Histogram '{0}' loaded".format(self.fFDSpectrumName))
        h_copy = h.Clone("FD")
        return h_copy

    def LoadFD(self):
        self.fFDHistogram = self.LoadVariation(self.fCentralPoints)

        filename = "{0}/{1}".format(self.fInputPath, self.fFileName)
        file = ROOT.TFile(filename)
        if not file or file.IsZombie():
            print("Could not open file {0}".format(filename))
            return
        else:
            print("File '{0}' open for FD correction".format(filename))

        # Loading systematic uncertainty
        systListName = "SystematicUncertainty"
        systList = file.Get(systListName)
        if not systList:
            print("Could not find list '{0}'".format(systListName))
            file.ls()
            return
        else:
            print("List '{0}' loaded".format(systListName))

        detLevListName = "{0}_CrossSection/{1}".format(self.fSpectrumName, self.fFDSpectrumName)
        detLevList = systList.FindObject(detLevListName)
        if not detLevList:
            print("Could not find list '{0}'".format(detLevListName))
            systList.Print()
            return
        else:
            print("List '{0}' loaded".format(detLevListName))

        detLevUpSystName = "{0}_UpperSyst".format(self.fFDSpectrumName)
        h = detLevList.FindObject(detLevUpSystName)
        if not h:
            print("Could not find hist '{0}'".format(detLevUpSystName))
            detLevList.Print()
            return
        else:
            print("Hist '{0}' loaded".format(detLevUpSystName))
        self.fFDUpSystUncHistogram = h.Clone("FD_DetectorLevelUpSyst")

        detLevLowSystName = "{0}_LowerSyst".format(self.fFDSpectrumName)
        h = detLevList.FindObject(detLevLowSystName)
        if not h:
            print("Could not find hist '{0}'".format(detLevLowSystName))
            detLevList.Print()
            return
        else:
            print("Hist '{0}' loaded".format(detLevLowSystName))
        self.fFDLowSystUncHistogram = h.Clone("FD_DetectorLevelLowSyst")

        tot_detLevUpSystName = "{0}_TotUpperSyst".format(self.fFDSpectrumName)
        h = detLevList.FindObject(tot_detLevUpSystName)
        if not h:
            print("Could not find hist '{0}'".format(tot_detLevUpSystName))
            detLevList.Print()
            return
        else:
            print("Hist '{0}' loaded".format(tot_detLevUpSystName))
        self.fFDTotUpSystUncHistogram = h.Clone("FD_DetectorLevelTotUpSyst")

        tot_detLevLowSystName = "{0}_TotLowerSyst".format(self.fFDSpectrumName)
        h = detLevList.FindObject(tot_detLevLowSystName)
        if not h:
            print("Could not find hist '{0}'".format(tot_detLevLowSystName))
            detLevList.Print()
            return
        else:
            print("Hist '{0}' loaded".format(tot_detLevLowSystName))
        self.fFDTotLowSystUncHistogram = h.Clone("FD_DetectorLevelTotLowSyst")

        detLevSystName = "{0}_CentralAsymmSyst".format(self.fFDSpectrumName)
        h = detLevList.FindObject(detLevSystName)
        if not h:
            print("Could not find hist '{0}'".format(detLevSystName))
            detLevList.Print()
            return
        else:
            print("Hist '{0}' loaded".format(detLevSystName))
        self.fFDSystUncGraph = h.Clone("FD_DetectorLevelSyst")

        tot_detLevSystName = "{0}_CentralTotAsymmSyst".format(self.fFDSpectrumName)
        h = detLevList.FindObject(tot_detLevSystName)
        if not h:
            print("Could not find hist '{0}'".format(tot_detLevSystName))
            detLevList.Print()
            return
        else:
            print("Hist '{0}' loaded".format(tot_detLevSystName))
        self.fFDTotSystUncGraph = h.Clone("FD_DetectorLevelTotSyst")

        file.Close()
