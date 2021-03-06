#!/usr/bin/env python
# python script to test Jet D meson analysis

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import argparse
import helperFunctions
import yaml


def AddDMesonJetTask(mgr, config, doRecLevel, doSignalOnly, doMCTruth, doWrongPID, doResponse, cmin, cmax):
    nOutputTrees = 0

    nDmesonTypes = 0;
    if config["dzero"]: nDmesonTypes += 2
    if config["dstar"]: nDmesonTypes += 1

    if doRecLevel:
        nOutputTrees += nDmesonTypes

    if config["MC"] and doSignalOnly:
        nOutputTrees += nDmesonTypes

    if config["MC"] and doMCTruth:
        nOutputTrees += nDmesonTypes

    if config["MC"] and doWrongPID and config["dzero"]:
        nOutputTrees += 1

    if not cmin is None and not cmax is None:
        suffix = "Cent_{0}_{1}".format(int(cmin), int(cmax))
    else:
        suffix = ""

    if config["ue_sub"]:
        rhoName = "Rho"
        rhoGenName = "RhoGen"
    else:
        rhoName = ""
        rhoGenName = ""

    if nOutputTrees > 0:
        if config["MC"]:
            if doResponse:
                if config["full_jets"]:
                    pDMesonJetsTask = ROOT.AliAnalysisTaskDmesonJetsDetectorResponse.AddTaskDmesonJetsDetectorResponse("usedefault", "usedefault", "usedefault", 3, suffix)
                    pDMesonJetsTask.SetNeedEmcalGeom(True)
                else:
                    pDMesonJetsTask = ROOT.AliAnalysisTaskDmesonJetsDetectorResponse.AddTaskDmesonJetsDetectorResponse("usedefault", "", "usedefault", 3, suffix)
            else:
                if config["full_jets"]:
                    pDMesonJetsTask = ROOT.AliAnalysisTaskDmesonJets.AddTaskDmesonJets("usedefault", "usedefault", "usedefault", nOutputTrees, suffix)
                    pDMesonJetsTask.SetNeedEmcalGeom(True)
                else:
                    pDMesonJetsTask = ROOT.AliAnalysisTaskDmesonJets.AddTaskDmesonJets("usedefault", "", "usedefault", nOutputTrees, suffix)
                pDMesonJetsTask.SetOutputType(ROOT.AliAnalysisTaskDmesonJets.kTreeExtendedOutput)
        else:
            if config["full_jets"]:
                pDMesonJetsTask = ROOT.AliAnalysisTaskDmesonJets.AddTaskDmesonJets("usedefault", "usedefault", "", nOutputTrees, suffix)
                pDMesonJetsTask.SetNeedEmcalGeom(True)
            else:
                pDMesonJetsTask = ROOT.AliAnalysisTaskDmesonJets.AddTaskDmesonJets("usedefault", "", "", nOutputTrees, suffix)
            pDMesonJetsTask.SetOutputType(ROOT.AliAnalysisTaskDmesonJets.kTreeExtendedOutput)

        pDMesonJetsTask.SetApplyKinematicCuts(True)
        if not cmin is None and not cmax is None:
            print("Setting centrality range ({}, {}) for task '{}'".format(cmin, cmax, pDMesonJetsTask.GetName()))
            pDMesonJetsTask.SetCentRange(cmin, cmax)

        if doRecLevel and not doResponse:
            # D0
            if config["dzero"]:
                if config["charged_jets"]:
                    eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "", ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kChargedJet, 0.4, rhoName)
                    # eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "", ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kChargedJet, 0.6, rhoName)
                    eng.SetJetPtRange(5, 500)
                    eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "loosest_nopid", ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kChargedJet, 0.4, rhoName)
                    # eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "loosest_nopid", ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kChargedJet, 0.6, rhoName)
                    eng.SetJetPtRange(5, 500)
                    eng.SetD0Extended(True)

                if config["full_jets"]:
                    eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "", ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kFullJet, 0.2, rhoName)
                    # eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "", ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kFullJet, 0.4, rhoName)
                    if not doResponse: eng.SetJetPtRange(5, 500)
                    eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "loosest_nopid", ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kFullJet, 0.2, rhoName)
                    # eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "loosest_nopid", ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kFullJet, 0.4, rhoName)
                    if not doResponse: eng.SetJetPtRange(5, 500)

            # D*
            if config["dstar"]:
                if config["charged_jets"]:
                    eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, config["rdhf_cuts_dstar"], "", ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kChargedJet, 0.4, rhoName)
                    # eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, config["rdhf_cuts_dstar"], "", ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kChargedJet, 0.6, rhoName)
                    if not doResponse: eng.SetJetPtRange(5, 500)

                if config["full_jets"]:
                    eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, config["rdhf_cuts_dstar"], "", ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kFullJet, 0.2, rhoName)
                    # eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, config["rdhf_cuts_dstar"], "", ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kFullJet, 0.4, rhoName)
                    if not doResponse: eng.SetJetPtRange(5, 500)

        if config["MC"] and doSignalOnly:
            # D0
            if config["dzero"]:
                if config["charged_jets"]:
                    eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "", ROOT.AliAnalysisTaskDmesonJets.kSignalOnly, ROOT.AliJetContainer.kChargedJet, 0.4, rhoName)
                    # eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "", ROOT.AliAnalysisTaskDmesonJets.kSignalOnly, ROOT.AliJetContainer.kChargedJet, 0.6, rhoName)
                    if not doResponse: eng.SetJetPtRange(5, 500)
                    eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "loosest_nopid", ROOT.AliAnalysisTaskDmesonJets.kSignalOnly, ROOT.AliJetContainer.kChargedJet, 0.4, rhoName)
                    # eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "loosest_nopid", ROOT.AliAnalysisTaskDmesonJets.kSignalOnly, ROOT.AliJetContainer.kChargedJet, 0.6, rhoName)
                    if not doResponse: eng.SetJetPtRange(5, 500)

                if config["full_jets"]:
                    eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "", ROOT.AliAnalysisTaskDmesonJets.kSignalOnly, ROOT.AliJetContainer.kFullJet, 0.2, rhoName)
                    # eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "", ROOT.AliAnalysisTaskDmesonJets.kSignalOnly, ROOT.AliJetContainer.kFullJet, 0.4, rhoName)
                    if not doResponse: eng.SetJetPtRange(5, 500)
                    eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "loosest_nopid", ROOT.AliAnalysisTaskDmesonJets.kSignalOnly, ROOT.AliJetContainer.kFullJet, 0.2, rhoName)
                    # eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "loosest_nopid", ROOT.AliAnalysisTaskDmesonJets.kSignalOnly, ROOT.AliJetContainer.kFullJet, 0.4, rhoName)
                    if not doResponse: eng.SetJetPtRange(5, 500)

            # D*
            if config["dstar"]:
                if config["charged_jets"]:
                    eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, config["rdhf_cuts_dstar"], "", ROOT.AliAnalysisTaskDmesonJets.kSignalOnly, ROOT.AliJetContainer.kChargedJet, 0.4, rhoName)
                    if not doResponse: eng.SetJetPtRange(5, 500)
                    # eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, config["rdhf_cuts_dstar"], "", ROOT.AliAnalysisTaskDmesonJets.kSignalOnly, ROOT.AliJetContainer.kChargedJet, 0.6, rhoName)

                if config["full_jets"]:
                    eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, config["rdhf_cuts_dstar"], "", ROOT.AliAnalysisTaskDmesonJets.kSignalOnly, ROOT.AliJetContainer.kFullJet, 0.2, rhoName)
                    if not doResponse: eng.SetJetPtRange(5, 500)
                    # eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, config["rdhf_cuts_dstar"], "", ROOT.AliAnalysisTaskDmesonJets.kSignalOnly, ROOT.AliJetContainer.kFullJet, 0.4, rhoName)

        if config["MC"] and doWrongPID and not doResponse and config["dzero"]:
            # D0
            if config["charged_jets"]:
                eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "", ROOT.AliAnalysisTaskDmesonJets.kD0Reflection, ROOT.AliJetContainer.kChargedJet, 0.4, rhoName)
                eng.SetJetPtRange(5, 500)
                # eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "", ROOT.AliAnalysisTaskDmesonJets.kD0Reflection, ROOT.AliJetContainer.kChargedJet, 0.6, rhoName)
                eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "loosest_nopid", ROOT.AliAnalysisTaskDmesonJets.kD0Reflection, ROOT.AliJetContainer.kChargedJet, 0.4, rhoName)
                eng.SetJetPtRange(5, 500)
                # eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "loosest_nopid", ROOT.AliAnalysisTaskDmesonJets.kD0Reflection, ROOT.AliJetContainer.kChargedJet, 0.6, rhoName)

            if config["full_jets"]:
                eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "", ROOT.AliAnalysisTaskDmesonJets.kD0Reflection, ROOT.AliJetContainer.kFullJet, 0.2, rhoName)
                # eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "", ROOT.AliAnalysisTaskDmesonJets.kD0Reflection, ROOT.AliJetContainer.kFullJet, 0.4, rhoName)
                if not doResponse: eng.SetJetPtRange(5, 500)
                eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "loosest_nopid", ROOT.AliAnalysisTaskDmesonJets.kD0Reflection, ROOT.AliJetContainer.kFullJet, 0.2, rhoName)
                # eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "loosest_nopid", ROOT.AliAnalysisTaskDmesonJets.kD0Reflection, ROOT.AliJetContainer.kFullJet, 0.4, rhoName)
                if not doResponse: eng.SetJetPtRange(5, 500)

        if config["MC"] and doMCTruth:
            # D0
            if config["dzero"]:
                if config["charged_jets"]:
                    eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "", ROOT.AliAnalysisTaskDmesonJets.kMCTruth, ROOT.AliJetContainer.kChargedJet, 0.4, rhoGenName)
                    # eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "", ROOT.AliAnalysisTaskDmesonJets.kMCTruth, ROOT.AliJetContainer.kChargedJet, 0.6, rhoGenName)
                    if not doResponse: eng.SetJetPtRange(5, 500)

                if config["full_jets"]:
                    eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "", ROOT.AliAnalysisTaskDmesonJets.kMCTruth, ROOT.AliJetContainer.kFullJet, 0.2, rhoGenName)
                    # eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, config["rdhf_cuts_dzero"], "", ROOT.AliAnalysisTaskDmesonJets.kMCTruth, ROOT.AliJetContainer.kFullJet, 0.4, rhoGenName)
                    if not doResponse: eng.SetJetPtRange(5, 500)

            # D*
            if config["dstar"]:
                pass
                if config["charged_jets"]:
                    eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, config["rdhf_cuts_dstar"], "", ROOT.AliAnalysisTaskDmesonJets.kMCTruth, ROOT.AliJetContainer.kChargedJet, 0.4, rhoGenName)
                    # eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, config["rdhf_cuts_dstar"], "", ROOT.AliAnalysisTaskDmesonJets.kMCTruth, ROOT.AliJetContainer.kChargedJet, 0.6, rhoGenName)
                    if not doResponse: eng.SetJetPtRange(5, 500)

                if config["full_jets"]:
                    eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, config["rdhf_cuts_dstar"], "", ROOT.AliAnalysisTaskDmesonJets.kMCTruth, ROOT.AliJetContainer.kFullJet, 0.2, rhoGenName)
                    # eng = pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, config["rdhf_cuts_dstar"], "", ROOT.AliAnalysisTaskDmesonJets.kMCTruth, ROOT.AliJetContainer.kFullJet, 0.4, rhoGenName)
                    if not doResponse: eng.SetJetPtRange(5, 500)

    return pDMesonJetsTask


def ExtractTriggerSelection(triggerList):
    r = 0
    for t in triggerList:
        if hasattr(ROOT.AliVEvent, t):
            r = r | getattr(ROOT.AliVEvent, t)
            print("Trigger '{}' added".format(t))
        else:
            print("Error: could not parse trigger name '{}'".format(t))
            exit(1)
    return r


def main(configFileName, nFiles, nEvents, d2h, doRecLevel, doSignalOnly, doMCTruth, doWrongPID, doResponse, noInclusiveJets, efficiency, selectGen, taskName="JetDmesonAna", debugLevel=0):

    f = open(configFileName, 'r')
    config = yaml.load(f)
    f.close()

    physSel = ExtractTriggerSelection(config["trigger"])
    print("Trigger selection is {}".format(physSel))

    ROOT.gSystem.Load("libCGAL")

    ROOT.AliTrackContainer.SetDefTrackCutsPeriod(config["run_period"])

    mode = None
    if config["mode"] == "AOD":
        mode = helperFunctions.AnaMode.AOD
    elif config["mode"] == "ESD":
        mode = helperFunctions.AnaMode.ESD
    else:
        print "Error: mode has to be either ESD or AOD!"
        exit(1)

    print("{0} analysis chosen.".format(config["mode"]))
    print("Setting local analysis for {0} files from list {1} max events = {2}".format(nFiles, config["file_list"], nEvents))

    # Analysis manager
    mgr = ROOT.AliAnalysisManager(taskName)

    helperFunctions.LoadMacros()

    if mode is helperFunctions.AnaMode.AOD:
        helperFunctions.AddAODHandler()
    elif mode is helperFunctions.AnaMode.ESD:
        helperFunctions.AddESDHandler()

    # task = helperFunctions.AddTaskCDBConnect()
    # task.SetFallBackToRaw(True)

    # Physics selection task
    if not config["MC"]:
        ROOT.AddTaskPhysicsSelection()

    if config["cent_type"] == "new":
        ROOT.AddTaskMultSelection(False)
    elif config["cent_type"] == "old":
        print("Not configured for old centrality framework!")
        return

    if config["full_jets"]:
        helperFunctions.PrepareEMCAL("userQAconfiguration.yaml")

    if config["MC"]:
        print "Running on a MC production"
    else:
        print "Running on data"

    # PID response
    helperFunctions.AddTaskPIDResponse(config["MC"], True, True, config["reco_pass"])

    if config["full_jets"]:
        pSpectraTask = ROOT.AddTaskEmcalJetQA("usedefault", "usedefault", "usedefault")
        pSpectraTask.SetNeedEmcalGeom(True)
    else:
        pSpectraTask = ROOT.AddTaskEmcalJetQA("usedefault", "", "")
        pSpectraTask.SetNeedEmcalGeom(False)

    pSpectraTask.SelectCollisionCandidates(physSel)
    pSpectraTask.SetPtBin(1, 150)
    if config["beam_type"] == "PbPb": pSpectraTask.SetCentRange(0, 90)

    if config["ue_sub"]:
        pKtChJetTask = ROOT.AddTaskEmcalJet("usedefault", "", ROOT.AliJetContainer.kt_algorithm, 0.4, ROOT.AliJetContainer.kChargedJet, 0.15, 0., 0.005, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
        pKtChJetTask.SelectCollisionCandidates(physSel)

        pChJetTask = ROOT.AddTaskEmcalJet("usedefault", "", ROOT.AliJetContainer.antikt_algorithm, 0.4, ROOT.AliJetContainer.kChargedJet, 0.15, 0., 0.005, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
        pChJetTask.SelectCollisionCandidates(physSel)

        pRhoTask = ROOT.AliAnalysisTaskRhoDev.AddTaskRhoDev("usedefault", 0.15, "", 0.30, "Rho", 0.4, ROOT.AliEmcalJet.kTPCfid, ROOT.AliJetContainer.kChargedJet, True)
        pRhoTask.SelectCollisionCandidates(physSel)
        pRhoTask.SetVzRange(-10, 10)
        pRhoTask.SetEventSelectionAfterRun(True)
        jetCont = ROOT.AliJetContainer(ROOT.AliJetContainer.kChargedJet, ROOT.AliJetContainer.kt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.4, pRhoTask.GetParticleContainer("tracks"), 0)
        jetCont.SetJetPtCut(1)
        jetCont.SetJetAcceptanceType(ROOT.AliEmcalJet.kTPCfid)
        jetCont.SetName("Signal")
        pRhoTask.AdoptJetContainer(jetCont)
        if config["beam_type"] == "pp":
            pRhoTask.SetRhoSparse(True)
        elif config["beam_type"] == "PbPb":
            pRhoTask.SetRhoSparse(False)
        elif config["beam_type"] == "pPb":
            pRhoTask.SetRhoSparse(True)

        pRhoExclLeadJetsTask = ROOT.AliAnalysisTaskRhoDev.AddTaskRhoDev("usedefault", 0.15, "", 0.30, "RhoExclLeadJets", 0.4, ROOT.AliEmcalJet.kTPCfid, ROOT.AliJetContainer.kChargedJet, True)
        pRhoExclLeadJetsTask.SelectCollisionCandidates(physSel)
        pRhoExclLeadJetsTask.SetVzRange(-10, 10)
        pRhoExclLeadJetsTask.SetEventSelectionAfterRun(True)
        pRhoExclLeadJetsTask.SetExcludeLeadJets(2)
        jetCont = ROOT.AliJetContainer(ROOT.AliJetContainer.kChargedJet, ROOT.AliJetContainer.kt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.4, pRhoExclLeadJetsTask.GetParticleContainer("tracks"), 0)
        jetCont.SetJetPtCut(1)
        jetCont.SetJetAcceptanceType(ROOT.AliEmcalJet.kTPCfid)
        jetCont.SetName("Signal")
        pRhoExclLeadJetsTask.AdoptJetContainer(jetCont)
        if config["beam_type"] == "pp":
            pRhoExclLeadJetsTask.SetRhoSparse(True)
        elif config["beam_type"] == "PbPb":
            pRhoExclLeadJetsTask.SetRhoSparse(False)
        elif config["beam_type"] == "pPb":
            pRhoExclLeadJetsTask.SetRhoSparse(True)

        pRhoTransTask = ROOT.AliAnalysisTaskRhoTransDev.AddTaskRhoTransDev("usedefault", 0.15, "", 0.30, "RhoTrans", 0.4, ROOT.AliEmcalJet.kTPCfid, ROOT.AliJetContainer.kChargedJet, True)
        pRhoTransTask.SelectCollisionCandidates(physSel)
        pRhoTransTask.SetVzRange(-10, 10)
        pRhoTransTask.SetEventSelectionAfterRun(True)

        pUEstudies = ROOT.AliAnalysisTaskJetUEStudies.AddTaskJetUEStudies("usedefault", "", 0.15, 0.30)
        jetCont = pUEstudies.AddJetContainer(ROOT.AliJetContainer.kChargedJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.4, ROOT.AliJetContainer.kTPCfid, "tracks", "")
        jetCont.SetRhoName("Rho")
        pUEstudies.AddAltRho("RhoTrans")
        pUEstudies.AddAltRho("RhoExclLeadJets")

        if config["MC"]:
            pGenKtChJetTask = ROOT.AddTaskEmcalJet("mcparticles", "", ROOT.AliJetContainer.kt_algorithm, 0.4, ROOT.AliJetContainer.kChargedJet, 0., 0., 0.005, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
            pGenKtChJetTask.SelectCollisionCandidates(physSel)

            pGenChJetTask = ROOT.AddTaskEmcalJet("mcparticles", "", ROOT.AliJetContainer.antikt_algorithm, 0.4, ROOT.AliJetContainer.kChargedJet, 0., 0., 0.005, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
            pGenChJetTask.SelectCollisionCandidates(physSel)

            pGenRhoTask = ROOT.AliAnalysisTaskRhoDev.AddTaskRhoDev("mcparticles", 0, "", 0, "RhoGen", 0.4, ROOT.AliEmcalJet.kTPCfid, ROOT.AliJetContainer.kChargedJet, True)
            pGenRhoTask.SelectCollisionCandidates(physSel)
            pGenRhoTask.SetVzRange(-10, 10)
            pGenRhoTask.SetEventSelectionAfterRun(True)
            jetCont = ROOT.AliJetContainer(ROOT.AliJetContainer.kChargedJet, ROOT.AliJetContainer.kt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.4, pGenRhoTask.GetParticleContainer("mcparticles"), 0)
            jetCont.SetJetPtCut(1)
            jetCont.SetJetAcceptanceType(ROOT.AliEmcalJet.kTPCfid)
            jetCont.SetName("Signal")
            pGenRhoTask.AdoptJetContainer(jetCont)
            if config["beam_type"] == "pp":
                pGenRhoTask.SetRhoSparse(True)
            elif config["beam_type"] == "PbPb":
                pGenRhoTask.SetRhoSparse(False)
            elif config["beam_type"] == "pPb":
                pGenRhoTask.SetRhoSparse(True)

            pGenRhoExclLeadJetsTask = ROOT.AliAnalysisTaskRhoDev.AddTaskRhoDev("mcparticles", 0, "", 0, "RhoExclLeadJetsGen", 0.4, ROOT.AliEmcalJet.kTPCfid, ROOT.AliJetContainer.kChargedJet, True)
            pGenRhoExclLeadJetsTask.SelectCollisionCandidates(physSel)
            pGenRhoExclLeadJetsTask.SetVzRange(-10, 10)
            pGenRhoExclLeadJetsTask.SetEventSelectionAfterRun(True)
            pGenRhoExclLeadJetsTask.SetExcludeLeadJets(2)
            jetCont = ROOT.AliJetContainer(ROOT.AliJetContainer.kChargedJet, ROOT.AliJetContainer.kt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.4, pGenRhoExclLeadJetsTask.GetParticleContainer("mcparticles"), 0)
            jetCont.SetJetPtCut(1)
            jetCont.SetJetAcceptanceType(ROOT.AliEmcalJet.kTPCfid)
            jetCont.SetName("Signal")
            pGenRhoExclLeadJetsTask.AdoptJetContainer(jetCont)
            if config["beam_type"] == "pp":
                pGenRhoExclLeadJetsTask.SetRhoSparse(True)
            elif config["beam_type"] == "PbPb":
                pGenRhoExclLeadJetsTask.SetRhoSparse(False)
            elif config["beam_type"] == "pPb":
                pGenRhoExclLeadJetsTask.SetRhoSparse(True)

            pRhoGenTransTask = ROOT.AliAnalysisTaskRhoTransDev.AddTaskRhoTransDev("mcparticles", 0, "", 0, "RhoTransGen", 0.4, ROOT.AliEmcalJet.kTPCfid, ROOT.AliJetContainer.kChargedJet, True)
            pRhoGenTransTask.SelectCollisionCandidates(physSel)
            pRhoGenTransTask.SetVzRange(-10, 10)
            pRhoGenTransTask.SetEventSelectionAfterRun(True)

            pUEstudies = ROOT.AliAnalysisTaskJetUEStudies.AddTaskJetUEStudies("mcparticles", "", 0., 0., "gen")
            jetCont = pUEstudies.AddJetContainer(ROOT.AliJetContainer.kChargedJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.4, ROOT.AliJetContainer.kTPCfid, "mcparticles", "")
            jetCont.SetRhoName("RhoGen")
            pUEstudies.AddAltRho("RhoTransGen")
            pUEstudies.AddAltRho("RhoExclLeadJetsGen")
            pUEstudies.SetMaxMomentumThridJet(15)
            pUEstudies.SetBackToBackJetPtFraction(0.8)

    if not noInclusiveJets:
        # Charged jet analysis
        if config["charged_jets"]:
            pChJetTask = ROOT.AddTaskEmcalJet("usedefault", "", ROOT.AliJetContainer.antikt_algorithm, 0.4, ROOT.AliJetContainer.kChargedJet, 0.15, 0., 0.005, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
            pChJetTask.SelectCollisionCandidates(physSel)

            # pChJetTask = ROOT.AddTaskEmcalJet("usedefault", "", 1, 0.6, ROOT.AliJetContainer.kChargedJet, 0.15, 0., 0.005, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
            # pChJetTask.SelectCollisionCandidates(physSel)

            if config["MC"]:
                pGenChJetTask = ROOT.AddTaskEmcalJet("mcparticles", "", ROOT.AliJetContainer.antikt_algorithm, 0.4, ROOT.AliJetContainer.kChargedJet, 0., 0., 0.005, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
                pGenChJetTask.SelectCollisionCandidates(physSel)

        # Full jet analysis
        if config["full_jets"]:
            pJetTask = ROOT.AddTaskEmcalJet("usedefault", "usedefault", ROOT.AliJetContainer.antikt_algorithm, 0.2, ROOT.AliJetContainer.kFullJet, 0.15, 0.30, 0.005, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
            pJetTask.SelectCollisionCandidates(physSel)

            # pJetTask = ROOT.AddTaskEmcalJet("usedefault", "usedefault", 1, 0.4, ROOT.AliJetContainer.kFullJet, 0.15, 0.30, 0.005, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
            # pJetTask.SelectCollisionCandidates(physSel)

            if config["MC"]:
                pGenJetTask = ROOT.AddTaskEmcalJet("mcparticles", "", ROOT.AliJetContainer.antikt_algorithm, 0.4, ROOT.AliJetContainer.kFullJet, 0., 0., 0.005, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
                pGenJetTask.SelectCollisionCandidates(physSel)

        ROOT.AliAnalysisManager.SetCommonFileName("AnalysisResults_jets.root")
        if config["full_jets"]:
            if config["ue_sub"]:
                anaType = ROOT.AliAnalysisTaskEmcalJetTreeBase.kJetPbPb
                rhoName = "Rho"
                rhoGenName = "RhoGen"
            else:
                anaType = ROOT.AliAnalysisTaskEmcalJetTreeBase.kJetPP
                rhoName = ""
                rhoGenName = ""
            pJetSpectraTask = ROOT.AliAnalysisTaskEmcalJetTreeBase.AddTaskEmcalJetTree("usedefault", "usedefault", 0.15, 0.30, anaType)
            pJetSpectraTask.SetNeedEmcalGeom(True)
        else:
            if config["ue_sub"]:
                anaType = ROOT.AliAnalysisTaskEmcalJetTreeBase.kJetPbPbCharged
                rhoName = "Rho"
                rhoGenName = "RhoGen"
            else:
                anaType = ROOT.AliAnalysisTaskEmcalJetTreeBase.kJetPPCharged
                rhoName = ""
                rhoGenName = ""
            pJetSpectraTask = ROOT.AliAnalysisTaskEmcalJetTreeBase.AddTaskEmcalJetTree("usedefault", "", 0.15, 0.30, anaType)
            pJetSpectraTask.SetNeedEmcalGeom(False)
        ROOT.AliAnalysisManager.SetCommonFileName("AnalysisResults.root")

        if config["MC"]:
            partCont = pJetSpectraTask.AddParticleContainer("mcparticles")
            partCont.SetMinPt(0)
        pJetSpectraTask.SelectCollisionCandidates(physSel)
        if config["beam_type"] == "PbPb": pJetSpectraTask.SetCentRange(0, 90)

        if config["charged_jets"]:
            jetCont = pJetSpectraTask.AddJetContainer(ROOT.AliJetContainer.kChargedJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.4, ROOT.AliJetContainer.kTPCfid, "tracks", "")
            jetCont.SetRhoName(rhoName)
            # jetCont = pJetSpectraTask.AddJetContainer(ROOT.AliJetContainer.kChargedJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.6, ROOT.AliJetContainer.kTPCfid)
            # jetCont.SetRhoName(rhoName)
            if config["MC"]:
                jetCont = pJetSpectraTask.AddJetContainer(ROOT.AliJetContainer.kChargedJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.4, ROOT.AliJetContainer.kTPCfid, "mcparticles", "")
                jetCont.SetRhoName(rhoGenName)

        if config["full_jets"]:
            jetCont = pJetSpectraTask.AddJetContainer(ROOT.AliJetContainer.kFullJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.2, ROOT.AliJetContainer.kEMCALfid, "tracks", "caloClusters")
            jetCont.SetRhoName(rhoName)
            # jetCont = pJetSpectraTask.AddJetContainer(ROOT.AliJetContainer.kFullJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.4, ROOT.AliJetContainer.kEMCALfid)
            # jetCont.SetRhoName(rhoName)
            if config["MC"]:
                jetCont = pJetSpectraTask.AddJetContainer(ROOT.AliJetContainer.kFullJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.4, ROOT.AliJetContainer.kTPCfid, "mcparticles", "")
                jetCont.SetRhoName(rhoGenName)

    # ROOT.AddTaskCleanupVertexingHF()

    if d2h: ROOT.AddTaskD0Mass(0, config["MC"], False, False, 0, 0, 0, 0, "Loose", config["rdhf_cuts_dzero"])

    if config["cent_bins"]:
        for cmin, cnmax in zip(config["cent_bins"][:-1], config["cent_bins"][1:]):
            pDMesonJetsTask = AddDMesonJetTask(mgr, config, doRecLevel, doSignalOnly, doMCTruth, doWrongPID, doResponse, cmin, cnmax)
            pDMesonJetsTask.SelectCollisionCandidates(physSel)
            pDMesonJetsTask.SetTrackEfficiency(efficiency)
    else:
        pDMesonJetsTask = AddDMesonJetTask(mgr, config, doRecLevel, doSignalOnly, doMCTruth, doWrongPID, doResponse, None, None)
        pDMesonJetsTask.SelectCollisionCandidates(physSel)
        pDMesonJetsTask.SetTrackEfficiency(efficiency)

    tasks = mgr.GetTasks()
    for task in tasks:
        if isinstance(task, ROOT.AliAnalysisTaskEmcalLight):
            task.SetIsPythia(config["is_pythia"])
            if selectGen: pDMesonJetsTask.SelectGeneratorName(selectGen)
            if config["is_pythia"]: task.SetMaxMinimumBiasPtHard(5)
            if config["run_period"] == "LHC15o": task.SetSwitchOffLHC15oFaultyBranches(True)
            if config["beam_type"] == "pp":
                task.SetForceBeamType(ROOT.AliAnalysisTaskEmcalLight.kpp)
            elif config["beam_type"] == "PbPb":
                task.SetForceBeamType(ROOT.AliAnalysisTaskEmcalLight.kAA)
            elif config["beam_type"] == "pPb":
                task.SetForceBeamType(ROOT.AliAnalysisTaskEmcalLight.kpA)
            if config["cent_type"] == "new":
                task.SetCentralityEstimation(ROOT.AliAnalysisTaskEmcalLight.kNewCentrality)
            elif config["cent_type"] == "old":
                task.SetCentralityEstimation(ROOT.AliAnalysisTaskEmcalLight.kOldCentrality)
            else:
                task.SetCentralityEstimation(ROOT.AliAnalysisTaskEmcalLight.kNoCentrality)

        if isinstance(task, ROOT.AliAnalysisTaskEmcal):
            if config["beam_type"] == "pp":
                task.SetForceBeamType(ROOT.AliAnalysisTaskEmcal.kpp)
            elif config["beam_type"] == "PbPb":
                task.SetForceBeamType(ROOT.AliAnalysisTaskEmcal.kAA)
                task.SetUseNewCentralityEstimation(True)

    res = mgr.InitAnalysis()

    if not res:
        print "Error initializing the analysis!"
        exit(1)

    mgr.PrintStatus()

    outFile = ROOT.TFile("train.root", "RECREATE")
    outFile.cd()
    mgr.Write()
    outFile.Close()

    chain = None
    if mode is helperFunctions.AnaMode.AOD:
        ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateAODChain.C")
        chain = ROOT.CreateAODChain(config["file_list"], nFiles, 0, False , "AliAOD.VertexingHF.root")
    elif mode is helperFunctions.AnaMode.ESD:
        ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
        chain = ROOT.CreateESDChain(config["file_list"], nFiles, 0, False)

    if debugLevel == 0:
        mgr.SetUseProgressBar(1, 250)
    else:
        # To have more debug info
        mgr.AddClassDebug("AliAnalysisTaskDmesonJets", ROOT.AliLog.kDebug + 100)
        mgr.AddClassDebug("AliAnalysisTaskDmesonJets::AnalysisEngine", ROOT.AliLog.kDebug + 100)

    mgr.SetDebugLevel(debugLevel)

    # start analysis
    print "Starting Analysis..."
    mgr.StartAnalysis("local", chain, nEvents)


if __name__ == '__main__':
    # runJetDmeson.py executed as script

    parser = argparse.ArgumentParser(description='Jet D meson analysis.')
    parser.add_argument('config',
                        help='YAML configuration file')
    parser.add_argument('-f', '--n-files',
                        default=100,
                        type=int,
                        help='Number of files to be analyzed')
    parser.add_argument('-e', '--n-events',
                        type=int,
                        default=2000,
                        help='Number of events to be analyzed')
    parser.add_argument('--no-signal-only', action='store_const',
                        default=False, const=True,
                        help='Signal only analysis')
    parser.add_argument('--no-rec-level', action='store_const',
                        default=False, const=True,
                        help='Reconstructed level analysis')
    parser.add_argument('--no-mc-truth', action='store_const',
                        default=False, const=True,
                        help='MC truth analysis')
    parser.add_argument('--no-wrong-pid', action='store_const',
                        default=False, const=True,
                        help='Wrong PID analysis (reflections)')
    parser.add_argument('--no-incl-jets', action='store_const',
                        default=False, const=True,
                        help='No inclusive jets')
    parser.add_argument('--d2h', action='store_const',
                        default=False, const=True,
                        help='D2H task')
    parser.add_argument('--response', action='store_const',
                        default=False, const=True,
                        help='Run response matrix')
    parser.add_argument('--efficiency',
                        default=0.0, type=float,
                        help='Artificial tracking inefficiency')
    parser.add_argument('--select-gen',
                        default=None,
                        help='Select events with a generator name')
    parser.add_argument('--task-name',
                        default="JetDmesonAna",
                        help='Task name')
    parser.add_argument('-d', '--debug-level',
                        default=0,
                        type=int,
                        help='Debug level')
    args = parser.parse_args()

    main(args.config, args.n_files, args.n_events, args.d2h, not args.no_rec_level, not args.no_signal_only, not args.no_mc_truth, not args.no_wrong_pid,
         args.response, args.no_incl_jets, args.efficiency, args.select_gen,
         args.task_name, args.debug_level)
