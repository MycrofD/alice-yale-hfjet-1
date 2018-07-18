#!/usr/bin/env python
# python program to project D meson jet trees into invariant mass histograms

import math
import os
import time

import ROOT

import DMesonJetUtils
import DetectorResponse
import DetectorResponseLoader

class DMesonJetProjector:

    def __init__(self, inputPath, train, fileName, taskName, tree_type, merging_type, norm_factor, maxEvents):
        self.fInputPath = inputPath
        self.fTrain = train
        self.fFileName = fileName
        self.fTaskName = taskName
        self.fMaxEvents = maxEvents
        self.fChain = None
        self.fCurrentFileName = None
        self.fCurrentRootFile = None
        self.fCurrentTreeNumber = -1
        self.fPeriod = None
        self.fPtHardBin = None
        self.fWeights = dict()
        self.fWeight = 1
        self.fNFiles = 0
        self.fTreeType = tree_type
        self.fMergingType = merging_type
        self.fHistEvents = None
        self.fNormFactor = norm_factor
        self.fDoNotAsk = True
        self.fFilePaths = None
        if self.fNormFactor > 1 or self.fNormFactor < 1:
            print("An artificial normalization factor of {} is applied.".format(self.fNormFactor))

    def GenerateChain(self, treeName):
        self.fChain = ROOT.TChain(treeName)
        path = "{0}/{1}".format(self.fInputPath, self.fTrain)

        print("Looking for file {0} in path {1}".format(self.fFileName, path))
        self.fFilePaths = list(DMesonJetUtils.find_file(path, self.fFileName))
        self.fNFiles = len(self.fFilePaths)

        for myfile in self.fFilePaths:
            print("Adding file {0}...".format(myfile))
            self.fChain.Add(myfile)

    def ExtractCurrentFileInfo(self):
        self.fCurrentRootFile = self.fChain.GetCurrentFile()
        self.fCurrentFileName = self.fCurrentRootFile.GetName()
        self.GetInfoFromFileName(self.fCurrentFileName)

    def GetInfoFromFileName(self, fname):
        lastSlash = fname.rfind('/')
        secondLastSlash = fname.rfind('/', 0, lastSlash - 1)
        thirdLastSlash = fname.rfind('/', 0, secondLastSlash - 1)
        ptHardBinStr = fname[secondLastSlash + 1:lastSlash]
        if ptHardBinStr.isdigit():
            self.fPtHardBin = int(ptHardBinStr)
        else:
            self.fPtHardBin = -1
        self.fPeriod = fname[thirdLastSlash + 1:secondLastSlash]

    def InitAnalysis(self, trigger, DMesonDef, DMesonDefSuffix):
        if trigger:
            listName = "{0}_{1}_histos".format(self.fTaskName, trigger)
        else:
            listName = "{0}_histos".format(self.fTaskName)

        for myfilename in self.fFilePaths:
            print("Extracting weights from file: '{}'.".format(myfilename))
            myfile = ROOT.TFile(myfilename, "read")

            if not myfile:
                print("Could not open file '{}'".format(myfilename))
                return 1

            self.CalculateEvents(myfile, trigger, DMesonDef, DMesonDefSuffix)

            if self.fMergingType != "explicit_weight": 
                weight = self.CalculateWeight(myfile, trigger)
                self.fWeights[myfilename] = weight
                print("Weight: {} (merging type = {})".format(weight, self.fMergingType))

        self.CalculateNormalizedEvents()

    def CalculateEvents(self, myfile, trigger, DMesonDef, DMesonDefSuffix):
        if trigger:
            listName = "{0}_{1}_histos".format(self.fTaskName, trigger)
            subListName = "histos{0}_{1}".format(self.fTaskName, trigger)
        else:
            listName = "{0}_histos".format(self.fTaskName)
            subListName = "histos{0}".format(self.fTaskName)

        if DMesonDefSuffix:
            subSubListName = "{}_{}".format(DMesonDef, DMesonDefSuffix)
        else:
            subSubListName = DMesonDef

        hlist = myfile.Get(listName)

        if not hlist:
            print("Could not get list '{0}' from file '{1}'".format(listName, myfile.GetName()))
            myfile.ls()
            return

        hSubList = hlist.FindObject(subListName)

        if not hSubList:
            print("Could not get list '{0}' from list '{1}'".format(subListName, listName))
            hlist.Print()
            return

        hSubSubList = hSubList.FindObject(subSubListName)

        if not hSubSubList:
            print("Could not get list '{0}' from list '{1}'".format(subSubListName, subListName))
            hSubList.Print()
            return

        print("Calculating events for file: '{}'".format(myfile.GetName()))
        events = self.ExtractEventsFromHistogramList(hSubSubList)

        if self.fHistEvents: self.fHistEvents.Add(events)
        else: self.fHistEvents = events

    def ExtractEventsFromHistogramList(self, hlist):
        accEventsHist = hlist.FindObject("fHistNEvents")
        rejEventsHist = hlist.FindObject("fHistEventRejectionReasons")

        if not (accEventsHist and rejEventsHist):
            print("Could not find event book-keeping histogram!")
            hlist.Print()
            return 0

        accepted = accEventsHist.GetBinContent(accEventsHist.GetXaxis().FindBin("Accepted"));
        nRecoVertVz_GT10 = rejEventsHist.GetBinContent(rejEventsHist.GetXaxis().FindBin("ZVtxOutFid"))
        nNoVert = rejEventsHist.GetBinContent(rejEventsHist.GetXaxis().FindBin("NoVertex")) + \
        rejEventsHist.GetBinContent(rejEventsHist.GetXaxis().FindBin("TooFewVtxContrib"))

        events = ROOT.TH1F("histEvents", "histEvents", 4, 0, 4)
        events.Fill("Accepted", accepted)
        events.Fill("ZVtxOutFid", nRecoVertVz_GT10)
        events.Fill("NoVert", nNoVert)

        print("Total number of accepted events with reconstructed vertex Vz < 10 cm: {0:e}".format(accepted))
        print("Total number of rejected events with reconstructed vertex Vz > 10 cm: {0:e}".format(nRecoVertVz_GT10))
        print("Total number of rejected events without reconstructed vertex: {0:e}".format(nNoVert))

        return events

    def CalculateNormalizedEvents(self):
        if not self.fHistEvents: return

        self.fAcceptedEvents = self.fHistEvents.GetBinContent(self.fHistEvents.GetXaxis().FindBin("Accepted"))
        self.fZVtxOutFidEvents = self.fHistEvents.GetBinContent(self.fHistEvents.GetXaxis().FindBin("ZVtxOutFid"))
        self.fNoVertEvents = self.fHistEvents.GetBinContent(self.fHistEvents.GetXaxis().FindBin("NoVert"))

        nRecoVert = self.fAcceptedEvents + self.fZVtxOutFidEvents
        self.fNormalizedEvents = self.fAcceptedEvents + self.fNoVertEvents * (1 - self.fZVtxOutFidEvents / nRecoVert)
        self.fHistEvents.Fill("Normalized", self.fNormalizedEvents)

        print("Calculation of the normalized events. Summary.")
        print("Total number of accepted events with reconstructed vertex Vz < 10 cm: {0:e}".format(self.fAcceptedEvents))
        print("Total number of rejected events with reconstructed vertex Vz > 10 cm: {0:e}".format(self.fZVtxOutFidEvents))
        print("Total number of rejected events without reconstructed vertex: {0:e}".format(self.fNoVertEvents))
        print("Correction factor: {}".format(self.fNormalizedEvents / self.fAcceptedEvents))
        print("Total number of normalized events: {0:e}".format(self.fNormalizedEvents))

    def CalculateWeight(self, myfile, trigger):
        if self.fMergingType == "simple_sum":
            weight = 1
        else:
            if trigger:
                listName = "{0}_{1}_histos".format(self.fTaskName, trigger)
            else:
                listName = "{0}_histos".format(self.fTaskName)
            hlist = myfile.Get(listName)

            if not hlist:
                print("Could not get list '{0}' from file '{1}'".format(listName, myfile.GetName()))
                weight = 1

            if "weighted_sum_old" in self.fMergingType:
                weight = self.ExtractWeightFromHistogramListOld(hlist)
            else:
                weight = self.ExtractWeightFromHistogramList(hlist)

            if self.fMergingType == "average":
                averageFactor = 1. / self.fNFiles

                if self.fMaxEvents > 0:
                    totEvents = self.fChain.GetEntries()
                    if self.fMaxEvents < totEvents:
                        averageFactor *= totEvents / self.fMaxEvents
                weight *= averageFactor
            else:
                weight *= self.fNormFactor
        return weight

    def ExtractWeightFromHistogramListOld(self, hlist):
        xsection = hlist.FindObject("fHistXsection")
        trials = hlist.FindObject("fHistTrials")

        if not trials or not xsection:
            print("Could not find trial and x-section information!")
            hlist.Print()
            return 1

        valNTRIALS = trials.GetBinContent(self.fPtHardBin + 1)
        valXSEC = xsection.GetBinContent(self.fPtHardBin + 1)
        if valNTRIALS > 0:
            return valXSEC / valNTRIALS

    def ExtractWeightFromHistogramList(self, hlist):
        xsection = hlist.FindObject("fHistXsectionNoSel")
        trials = hlist.FindObject("fHistTrialsVsPtHardNoSel")

        if not trials or not xsection:
            print("Falling back to secondary method for x-section and trials...")
            xsection = hlist.FindObject("fHistXsectionVsPtHardNoSel")
            trials = hlist.FindObject("fHistTrialsVsPtHardNoSel")

        if not trials or not xsection:
            print("Falling back to tertiary method for x-section and trials...")
            xsection = hlist.FindObject("fHistXsectionAfterSel")
            trials = hlist.FindObject("fHistTrialsAfterSel")

        if not trials or not xsection:
            print("Could not find trial and x-section information (not necessarily a bad thing)!")
            hlist.Print()
            return 1

        valNTRIALS = trials.Integral()
        valXSEC = xsection.GetMean(2)
        if valNTRIALS > 0:
            return valXSEC / valNTRIALS

    def GetExplicitWeight(self, entry):
        event = entry.Event
        self.fWeight = event.fWeight / self.fNormalizedEvents
        # print("w = {:e}".format(self.fWeight))

    def StartProjection(self, trigger, DMesonDef, DMesonDefSuffix, output, norm=1.0):
        treeName = "_".join([obj for obj in [self.fTaskName, trigger, DMesonDef, DMesonDefSuffix] if obj])

        print("Running analysis on tree {}.".format(treeName))
        if self.fDoNotAsk:
            print("Waiting 3 seconds...")
            time.sleep(3)
        else:
            raw_input("Press <ENTER> to continue...")

        self.fHistEvents = None
        self.GenerateChain(treeName)
        self.InitAnalysis(trigger, DMesonDef, DMesonDefSuffix)

        print("Total number of entries is {}".format(self.fChain.GetEntries()))
        if self.fMaxEvents > 0:
            print("The analysis will stop at the {0} entry.".format(self.fMaxEvents))

        if self.fTreeType == "simple":
            self.DoProjectionSimple(trigger, DMesonDef, DMesonDefSuffix, output, norm)
        elif self.fTreeType == "extended":
            self.DoProjectionExtended(trigger, DMesonDef, DMesonDefSuffix, output, norm)
        else:
            print("Tree type '{}' not recognized.".format(self.fTreeType))
            exit(1)

    def DoProjectionExtended(self, trigger, DMesonDef, DMesonDefSuffix, output, norm):
        class DMesonEvent:
            pass
        jet_branch_names = output.GetJetBranches()
        dmesonEvent = DMesonEvent()
        idmeson = 0
        for i, event in enumerate(self.fChain):
            if i % 10000 == 0:
                print("Event n. {0}".format(i))
                if self.fMaxEvents > 0 and i >= self.fMaxEvents:
                    print("Stopping the analysis.")
                    break
            self.OnFileChange()
            if self.fMergingType == "explicit_weight": 
                self.GetExplicitWeight(event)
            jet_branches = [getattr(event, jet_branch_name) for jet_branch_name in jet_branch_names]
            for dmeson_and_jets in zip(event.Dmesons, *jet_branches):
                if idmeson % 10000 == 0:
                    print("D meson candidate n. {0}".format(idmeson))
                dmeson = dmeson_and_jets[0]
                jets = dmeson_and_jets[1:]
                dmesonEvent.DmesonJet = dmeson
                for jet, jet_branch_name in zip(jets, jet_branch_names):
                    setattr(dmesonEvent, jet_branch_name, jet)
                output.Fill(dmesonEvent, self.fWeight * norm)
                idmeson += 1

    def DoProjectionSimple(self, trigger, DMesonDef, DMesonDefSuffix, output, norm):
        for i, dmesonEvent in enumerate(self.fChain):
            if i % 10000 == 0:
                print("D meson candidate n. {0}".format(i))
                if self.fMaxEvents > 0 and i >= self.fMaxEvents:
                    print("Stopping the analysis.")
                    break
            self.OnFileChange()
            output.Fill(dmesonEvent, self.fWeight * norm)

    def OnFileChange(self):
        if self.fCurrentTreeNumber == self.fChain.GetTreeNumber():
            return

        self.fCurrentTreeNumber = self.fChain.GetTreeNumber()
        self.ExtractCurrentFileInfo()

        if self.fMergingType != "explicit_weight":
            self.fWeight = self.fWeights[self.fCurrentFileName]
            print("File: {}\nPeriod: {}\nPt hard bin: {}\nWeight: {} (merging type = {})".format(self.fCurrentFileName, self.fPeriod, self.fPtHardBin, self.fWeight, self.fMergingType))
