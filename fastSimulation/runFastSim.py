#!/usr/bin/env python

import time
import datetime
import platform
import os
import shutil
import subprocess
import argparse
import random
import sys
import yaml
import GeneratePowhegInput


class PowhegResult:

    def __init__(self, events_generated, lhe_file):
        self.lhe_file = lhe_file
        self.events_generated = events_generated


def RunPowhegParallel(powhegExe, pythiaEvents, powheg_stage, job_number):
    print("Running POWHEG simulation at stage {}!".format(powheg_stage))

    with open("powheg.input", 'r') as fin:
        powheg_input = fin.read().splitlines()
    for line in powheg_input:
        print(line)

    print("Running POWHEG...")
    LogFileName = "Powheg_Stage_{}_Job_{:03d}.log".format(powheg_stage, job_number)
    with open(LogFileName, "w") as myfile:
        print([powhegExe, str(job_number)])
        p = subprocess.Popen([powhegExe], stdout=myfile, stderr=myfile, stdin=subprocess.PIPE)
        p.communicate(input=str(job_number))

    if powheg_stage == 4:
        result = PowhegResult(True, "pwgevents.lhe")
    else:
        result = PowhegResult(False, "")

    return result


def RunPowhegSingle(powhegExe, pythiaEvents, yamlConfigFile):
    print("Running POWHEG simulation!")

    GeneratePowhegInput.main("./", pythiaEvents, 0, yamlConfigFile)

    with open("powheg.input", 'r') as fin:
        powheg_input = fin.read().splitlines()
    for line in powheg_input:
        print(line)

    print("Running POWHEG...")
    with open("powheg.log", "w") as myfile:
        subprocess.call([powhegExe], stdout=myfile, stderr=myfile)

    result = PowhegResult(True, "pwgevents.lhe")

    return result


def main(pythiaEvents, powheg_stage, job_number, yamlConfigFile, batch_job, LHEfile, minpthard, maxpthard, debug_level):
    print("------------------ job starts ---------------------")
    dateNow = datetime.datetime.now()
    print(dateNow)

    try:
        rootPath = subprocess.check_output(["which", "root"]).rstrip()
        alirootPath = subprocess.check_output(["which", "aliroot"]).rstrip()
    except subprocess.CalledProcessError:
        print "Environment is not configured correctly!"
        exit()

    print "Root: " + rootPath
    print "AliRoot: " + alirootPath

    f = open(args.config, 'r')
    config = yaml.load(f)
    f.close()

    gen = config["gen"]
    proc = config["proc"]
    beamType = config["beam_type"]
    ebeam1 = config["ebeam1"]
    ebeam2 = config["ebeam2"]
    nPDFset = config["nPDFset"]
    nPDFerrSet = config["nPDFerrSet"]
    if "rejectISR" in config:
        rejectISR = config["rejectISR"]
    else:
        rejectISR = False

    if batch_job == "grid":
        fname = "{0}_{1}".format(gen, proc)
    elif batch_job == "lbnl3":
        fname = "{}_{}_{:03d}".format(gen, proc, job_number)
    else:
        unixTS = int(time.time())
        fname = "{0}_{1}_{2}".format(gen, proc, unixTS)

    print("Running {0} MC production on: {1}".format(proc, " ".join(platform.uname())))

    if "powheg" in gen:
        if os.path.isfile("pwgevents.lhe") or os.path.isfile("powheg.input") or os.path.isfile("powheg.log"):
            print("Before running POWHEG again you must delete or move the following files: pwgevents.lhe, powheg.input, powheg.log")
            exit(1)

        if LHEfile:
            print("Using previously generated POWHEG events from file {0}!".format(LHEfile))
            runPOWHEG = False
        else:
            runPOWHEG = True
    else:
        runPOWHEG = False

    if runPOWHEG:
        if batch_job == "local":
            powhegExe = "./POWHEG_bins/{0}".format(powhegExe)

        if powheg_stage > 0 and powheg_stage <= 4:
            powheg_result = RunPowhegParallel(powhegExe, pythiaEvents, powheg_stage, job_number)
        else:
            powheg_result = RunPowhegSingle(powhegExe, pythiaEvents, yamlConfigFile)

        if not powheg_result.events_generated:
            if powheg_stage > 0 and powheg_stage <= 3:
                print("POWHEG stage {} completed. Exiting.".format(powheg_stage))
                os.unlink("powheg.input")
                os.rename("powheg.log", "powheg-stage{}-{:03d}.log".format(powheg_stage, job_number))
                exit(0)
            else:
                print("POWHEG at stage {} did not produce any event!!!".format(powheg_stage))
                exit(1)

        if not os.path.isfile(powheg_result.lhe_file):
            print("Could not find POWHEG output pwgevents.lhe. Something went wrong, aborting...")
            if os.path.isfile("powheg.log"):
                print("Check log file below.")
                with open("powheg.log", "r") as myfile:
                    powheg_log = myfile.read().splitlines()
                for line in powheg_log:
                    print(line)
            else:
                print("No log file was found.")
            exit(1)

        if batch_job == "grid" or batch_job == "lbnl3":
            LHEfile = powheg_result.lhe_file
        else:
            LHEfile = "pwgevents_{0}.lhe".format(fname)
            os.rename("powheg.input", "{0}.input".format(fname))
            print("POWHEG configuration backed up in {0}.input".format(fname))
            os.rename(powheg_result.lhe_file, LHEfile)
            print("POWHEG events backed up in {0}".format(LHEfile))
            os.rename("powheg.log", "{0}.log".format(fname))
            print("POWHEG log backed up in {0}.log".format(fname))
            # cleaning the working directory and archiving POWHEG files
            subprocess.call(["./clean_powheg.sh", "{0}.tar".format(fname)])

    rnd = random.randint(0, 1073741824)  # 2^30
    print("Setting PYTHIA seed to {0}".format(rnd))

    print("Compiling analysis code...")
    subprocess.call(["make"])

    print("Running PYTHIA simulation...")
    with open("sim_{0}.log".format(fname), "w") as myfile:
        subprocess.call(["aliroot", "-b", "-l", "-q", "start_simulation.C(\"{0}\", {1}, \"{2}\", \"{3}\", {4}, \"{5}\", \"{6}\", {7}, {8}, {9}, {10}, {11}, {12})".format(fname, pythiaEvents, proc, gen, rnd, LHEfile, beamType, ebeam1, ebeam2, int(rejectISR), minpthard, maxpthard, debug_level)], stdout=myfile, stderr=myfile)

    print("Done")
    print("...see results in the log files")

    subprocess.call(["ls", "-l"])

    print("############### SCRATCH SIZE ####################")
    subprocess.call(["du", "-sh"])

    print("------------------ job ends ----------------------")
    dateNow = datetime.datetime.now()
    print(dateNow)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run fast simulation.')
    parser.add_argument('config', metavar='config.yaml',
                        default="default.yaml", help='YAML configuration file')
    parser.add_argument('--numevents', metavar='NEVT',
                        default=50000, type=int)
    parser.add_argument('--lhe', metavar='LHE',
                        default='')
    parser.add_argument('--minpthard', metavar="MINPTHARD",
                        default=-1, type=float)
    parser.add_argument('--maxpthard', metavar='MAXPTHARD',
                        default=-1, type=float)
    parser.add_argument('--batch-job', metavar='job',
                        default='local')
    parser.add_argument('--powheg-stage', metavar='STAGE',
                        default=0, type=int)
    parser.add_argument('--job-number', metavar='STAGE',
                        default=0, type=int)
    parser.add_argument('-d', metavar='debug_level',
                        default=0, type=int)
    args = parser.parse_args()

    main(args.numevents, args.powheg_stage, args.job_number, args.config, args.lhe, args.minpthard, args.maxpthard, args.batch_job, args.d)
