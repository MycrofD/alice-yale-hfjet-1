#!/usr/bin/env python

from ROOT import TFileMerger

def MergeFiles(output, fileList, n=20):
    merger = TFileMerger(False)
    merger.OutputFile(output);
    merger.SetMaxOpenedFiles(n);

    print "Total number of files is {0}".format(len(fileList))
    
    for fileName in fileList:
        print "Adding file {0}".format(fileName)
        merger.AddFile(fileName)

    merger.AddObjectNames("AnalysisCutsDzero");
    merger.AddObjectNames("AnalysisCutsDStar");
  
    merger.PrintFiles("");
    r = merger.PartialMerge(TFileMerger.kAllIncremental | TFileMerger.kSkipListed);

    if not r:
        print "Merge error!"

    return r