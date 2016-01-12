// Base class for D jet correlation analysis
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

#include <algorithm>

#include <TList.h>
#include <TFile.h>
#include <TString.h>
#include <THnSparse.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <Riostream.h>
#include <TParameter.h>
#include <TObjString.h>
#include <TString.h>
#include <TDirectoryFile.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLine.h>

#include "DJetCorrAnalysisParams.h"

#include "DJetCorrBase.h"

const Double_t DJetCorrBase::fgkEpsilon = 1E-6;

ClassImp(DJetCorrBase);

//____________________________________________________________________________________
DJetCorrBase::DJetCorrBase() :
  TNamed(),
  fTrainName(),
  fInputPath(),
  fInputFileName(),
  fInputDirFileName(),
  fOutputPath(),
  fOutputFileName(),
  fOverwrite(kFALSE),
  fAnalysisParams(0),
  fPlotFormat("pdf"),
  fSavePlots(kFALSE),
  fAddTrainToCanvasName(kFALSE),
  fAnaType(DJetCorrAnalysisParams::KUndefinedAna),
  fTHnSparseAxisMaps(),
  fInputFile(0),
  fInputDirectoryFile(0),
  fInputList(0),
  fOutputList(0),
  fCanvases(0),
  fEvents(0.)
{
  // Default ctr.
  
}

//____________________________________________________________________________________
DJetCorrBase::DJetCorrBase(const char* train, const char* path) :
  TNamed(train, train),
  fTrainName(train),
  fInputPath(path),
  fInputFileName(),
  fInputDirFileName(),
  fOutputPath("../data/"),
  fOutputFileName("DJetCorr.root"),
  fOverwrite(kFALSE),
  fAnalysisParams(new TList()),
  fPlotFormat("pdf"),
  fSavePlots(kFALSE),
  fAddTrainToCanvasName(kFALSE),
  fAnaType(DJetCorrAnalysisParams::KUndefinedAna),
  fTHnSparseAxisMaps(),
  fInputFile(0),
  fInputDirectoryFile(0),
  fInputList(0),
  fOutputList(0),
  fCanvases(0),
  fEvents(0.)
{
  // Standard ctr.

}

//____________________________________________________________________________________
Bool_t DJetCorrBase::ClearInputData()
{
  // Clear the input data.

  fTHnSparseAxisMaps.Clear();

  return kTRUE;
}

//____________________________________________________________________________________
DJetCorrAnalysisParams* DJetCorrBase::AddAnalysisParams(const char* dmeson, const char* jetType, const char* jetRadius, const char* tracksName, Bool_t isMC, Bool_t isBkgSub)
{
  // Add the analysis params for the D meson type, jet type and jet radius provided.

  if (!fAnalysisParams) fAnalysisParams = new TList();

  DJetCorrAnalysisParams* params = new DJetCorrAnalysisParams(dmeson, jetType, jetRadius, tracksName, fAnaType, isMC, isBkgSub);

  fAnalysisParams->Add(params);

  return params;
}

//____________________________________________________________________________________
DJetCorrAnalysisParams* DJetCorrBase::AddAnalysisParams(DJetCorrAnalysisParams* params)
{
  // Add the analysis params.

  if (!fAnalysisParams) fAnalysisParams = new TList();

  DJetCorrAnalysisParams* paramsCopy = new DJetCorrAnalysisParams(*params);
  
  fAnalysisParams->Add(paramsCopy);

  return paramsCopy;
}

//____________________________________________________________________________________
TMap* DJetCorrBase::GenerateAxisMap(THnSparse* hn)
{
  if (!hn) return 0;
  
  TMap* axisMap = new TMap();
  axisMap->SetOwnerKeyValue();
  
  for (Int_t i = 0; i < hn->GetNdimensions(); i++) {
    TObjString* key = new TObjString(hn->GetAxis(i)->GetTitle());
    TParameter<Int_t>* value = new TParameter<Int_t>(hn->GetAxis(i)->GetTitle(), i);
    axisMap->Add(key, value);
  }

  fTHnSparseAxisMaps.Add(axisMap);
  
  return axisMap;
}

//____________________________________________________________________________________
Bool_t DJetCorrBase::Init()
{
  // Init class.

  Printf("Info-DJetCorrBase::Init : Now initializing.");
  
  if (fOutputList) {
    delete fOutputList;
    fOutputList = 0;
  }
  fOutputList = new THashList();

  if (fCanvases) {
    delete fCanvases;
    fCanvases = 0;
  }
  fCanvases = new TList();  
    
  Printf("Info-DJetCorrBase::Init : Initialization done.");

  return kTRUE;
}

//____________________________________________________________________________________
TFile* DJetCorrBase::OpenOutputFile()
{
  TString fname(fOutputPath);
  
  if (fTrainName.IsNull()) {
    fname += "/test/";
  }
  else {
    fname += "/";
    fname += fTrainName;
    fname += "/";
  }

  fname += fOutputFileName;
  
  TFile* outputFile = TFile::Open(fname);

  if (!outputFile || outputFile->IsZombie()) {
    Printf("Error-DJetCorrBase::LoadOutputHistograms : Could not open file '%s' to read.", fname.Data()); 
    outputFile = 0;
  }

  return outputFile;
}

//____________________________________________________________________________________
Bool_t DJetCorrBase::LoadOutputHistograms()
{
  TFile* outputFile = OpenOutputFile();

  if (!outputFile) return kFALSE;

  fOutputList = static_cast<THashList*>(outputFile->Get("fOutputList"));

  outputFile->Close();
  delete outputFile;
  outputFile = 0;

  return kTRUE;
}

//____________________________________________________________________________________
void DJetCorrBase::SavePlot(TCanvas* canvas)
{
  if (!canvas) return;
  
  TString fname(fOutputPath);

  if (fTrainName.IsNull()) {
    fname += "/test/";
  }
  else {
    fname += "/";
    fname += fTrainName;
    fname += "/";
  }

  fname += canvas->GetTitle();
  fname += ".";
  fname += fPlotFormat;

  canvas->SaveAs(fname);
}

//____________________________________________________________________________________
TVirtualPad* DJetCorrBase::SetUpPad(TVirtualPad* pad,
                                    const char* xTitle, Double_t minX, Double_t maxX, Bool_t logX,
                                    const char* yTitle, Double_t minY, Double_t maxY, Bool_t logY,
                                    Double_t lmar, Double_t rmar, Double_t bmar, Double_t tmar)
{
  if (!pad) return 0;

  pad->SetLeftMargin(lmar);
  pad->SetRightMargin(rmar);
  pad->SetBottomMargin(bmar);
  pad->SetTopMargin(tmar);
  
  pad->SetLogx(logX);
  pad->SetLogy(logY);

  pad->SetTicks(1, 1);

  pad->cd();

  TString blankHistName(Form("%s_blankHist", pad->GetName()));
  TH1* blankHist = new TH1D(blankHistName, blankHistName, 1000, minX, maxX);

  if (logY) {
    if (maxY <= 0 && minY > 0) {
      maxY = minY * 100;
    }
    else if (minY <= 0 && maxY > 0) {
      minY = maxY / 100;
    }
    else {
      minY = 1;
      maxY = 1;
    }
  }

  blankHist->SetMinimum(minY);
  blankHist->SetMaximum(maxY);
  
  blankHist->GetXaxis()->SetTitle(xTitle);
  blankHist->GetYaxis()->SetTitle(yTitle);

  blankHist->GetXaxis()->SetTitleOffset(1.12);
  blankHist->GetYaxis()->SetTitleOffset(1.14);

  blankHist->GetXaxis()->SetTitleFont(43);
  blankHist->GetYaxis()->SetTitleFont(43);

  blankHist->GetXaxis()->SetTitleSize(20);
  blankHist->GetYaxis()->SetTitleSize(20);

  blankHist->GetXaxis()->SetLabelOffset(0.010);
  blankHist->GetYaxis()->SetLabelOffset(0.010);

  blankHist->Draw("AXIS");

  return pad;
}

//____________________________________________________________________________________
TCanvas* DJetCorrBase::SetUpCanvas(TH1* histo, Bool_t logX, Bool_t logY,
                                   Double_t w, Double_t h, Int_t rows, Int_t cols,
                                   Double_t lmar, Double_t rmar, Double_t bmar, Double_t tmar)
{
  return SetUpCanvas(histo->GetName(),
                     histo->GetXaxis()->GetTitle(), histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax(), logX,
                     histo->GetYaxis()->GetTitle(), histo->GetYaxis()->GetXmin(), histo->GetYaxis()->GetXmax(), logY,
                     w, h, rows, cols, lmar, rmar, bmar, tmar);
}

//____________________________________________________________________________________
TCanvas* DJetCorrBase::SetUpCanvas(const char* name,
                                   const char* xTitle, Double_t minX, Double_t maxX, Bool_t logX,
                                   const char* yTitle, Double_t minY, Double_t maxY, Bool_t logY,
                                   Double_t w, Double_t h, Int_t rows, Int_t cols,
                                   Double_t lmar, Double_t rmar, Double_t bmar, Double_t tmar)
{
  Printf("Info-DJetCorrBase::SetUpCanvas : Setting up canvas '%s'", name);

  TString cname;

  if (fAddTrainToCanvasName) {
    cname = Form("%s_%s", GetName(), name);
  }
  else {
    cname = name;
  }

  TCanvas* canvas = static_cast<TCanvas*>(fCanvases->FindObject(cname));
  if (canvas) {
    Printf("Warning-DJetCorrBase::SetUpCanvas : Canvas '%s' already exists.", cname.Data());
    return canvas;
  }
  
  canvas = new TCanvas(cname, cname, w, h);
  
  if (rows == 1 && cols == 1) {
    SetUpPad(canvas, xTitle, minX, maxX, logX, yTitle, minY, maxY, logY);
    canvas->SetLeftMargin(lmar);
    canvas->SetRightMargin(rmar);
    canvas->SetBottomMargin(bmar);
    canvas->SetTopMargin(tmar);
  }
  else {
    canvas->Divide(cols, rows);
    Int_t n = rows * cols;
    for (Int_t i = 1; i <= n; i++) {
      TVirtualPad* pad = canvas->cd(i);
      SetUpPad(pad, xTitle, minX, maxX, logX, yTitle, minY, maxY, logY,
               lmar, rmar, bmar, tmar); 
    }
  }

  fCanvases->Add(canvas);
  
  return canvas;
}

//____________________________________________________________________________________
TLegend* DJetCorrBase::SetUpLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Int_t textSize)
{
  TLegend* leg = new TLegend(x1, y1, x2, y2);
  leg->SetTextFont(43);
  leg->SetTextSize(textSize);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  return leg;
}

//____________________________________________________________________________________
TPaveText* DJetCorrBase::SetUpPaveText(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Int_t textSize, const char* text)
{
  TPaveText* pave = new TPaveText(x1, y1, x2, y2, "brNDC");
  pave->SetTextFont(43);
  pave->SetTextAlign(11);
  pave->SetTextSize(textSize);
  pave->SetFillStyle(0);
  pave->SetBorderSize(0);
  
  if (text) {
    TString stext(text);
    TObjArray* lines = stext.Tokenize("\n");
    TIter next(lines);
    TObject *obj = 0;
    while ((obj = next())) {
      pave->AddText(obj->GetName());
    }
    delete lines;
  }

  return pave;
}

//____________________________________________________________________________________
Bool_t DJetCorrBase::PlotObservable(DJetCorrAnalysisParams* params, TString obsName, Double_t xmin, Double_t xmax,
                                        Double_t minDPt, Double_t maxDPt, Double_t minJetPt, Double_t maxJetPt, Double_t minZ, Double_t maxZ,
                                        Int_t step, Int_t rebin, Int_t norm, Int_t plotStats)
{
  if (!fOutputList) return kFALSE;

  TString cname(Form("fig_%s_%s",obsName.Data(), params->GetName()));
  
  TString jetCuts;
  TString dCuts;
  TString zCuts;

  Int_t nj = params->GetNJetPtBins();
  Int_t nd = params->GetNDPtBins();
  Int_t nz = params->GetNzBins();

  if (maxJetPt > 0) {
    jetCuts = Form("JetPt_%03.0f_%03.0f", minJetPt, maxJetPt);
    jetCuts.ReplaceAll(".", "");
    nj = 1;
  }

  if (maxDPt > 0) {
    dCuts = Form("DPt_%02.0f_%02.0f", minDPt, maxDPt);
    dCuts.ReplaceAll(".", "");
    nd = 1;
  }

  if (maxZ > 0) {
    zCuts = Form("z_%.1f_%.1f", minZ, maxZ);
    zCuts.ReplaceAll(".", "");
    nz = 1;
  }
  
  if (!jetCuts.IsNull()) {
    cname += "_";
    cname += jetCuts;
  }

  if (!dCuts.IsNull()) {
    cname += "_";
    cname += dCuts;
  }

  if (!zCuts.IsNull()) {
    cname += "_";
    cname += zCuts;
  }

  TString hname(obsName);
  TString prefix(params->GetName());

  TObjArray histos;
  Int_t ih = 0;
  TString jetTitle;
  TString dTitle;
  TString zTitle;
  
  for (Int_t ij = 0; ij < nj; ij+=step) {
    if (maxJetPt < 0) {
      jetCuts = Form("JetPt_%03.0f_%03.0f", params->GetJetPtBin(ij), params->GetJetPtBin(ij+1));
      jetCuts.ReplaceAll(".", "");
      jetTitle = Form("%.1f < #it{p}_{T,jet} < %.1f GeV/#it{c}", params->GetJetPtBin(ij), params->GetJetPtBin(ij+1));
    }
    for (Int_t id = 0; id < nd; id+=step) {
      if (maxDPt < 0) {
        dCuts = Form("DPt_%02.0f_%02.0f", params->GetDPtBin(id), params->GetDPtBin(id+1));
        dCuts.ReplaceAll(".", "");
        dTitle = Form("%.1f < #it{p}_{T,D} < %.1f GeV/#it{c}", params->GetDPtBin(id), params->GetDPtBin(id+1));
      }
      for (Int_t iz = 0; iz < nz; iz+=step) {
        if (maxZ < 0) {
          zCuts = Form("z_%.1f_%.1f", params->GetzBin(iz), params->GetzBin(iz+1));
          zCuts.ReplaceAll(".", "");
          zTitle = Form("%.1f < z < %.1f", params->GetzBin(iz), params->GetzBin(iz+1));
        }
            
        TString cuts(Form("%s_%s_%s", dCuts.Data(), jetCuts.Data(), zCuts.Data()));
        TString htitle;

        htitle = jetTitle;
        if (!dTitle.IsNull()) {
          if (!htitle.IsNull()) htitle += ", ";
          htitle += dTitle;
        }
        if (!zTitle.IsNull()) {
          if (!htitle.IsNull()) htitle += ", ";
          htitle += zTitle;
        }
    
        TString objname(Form("h%s_%s_%s_Matched", prefix.Data(), hname.Data(), cuts.Data()));
        //Printf("Info-DJetCorrAnalysis::PlotObservable : Retrieving histogram '%s'", objname.Data());
        TH1* h = dynamic_cast<TH1*>(fOutputList->FindObject(objname));
        if (!h) {
          Printf("Error-DJetCorrAnalysis::PlotObservable : Histogram '%s' not found!", objname.Data());
          continue;
        }
        TString newName(objname);
        newName += "_copy";
        //Printf("Info-DJetCorrAnalysis::PlotObservable : Cloning histogram '%s'", objname.Data());
        TH1* hist = static_cast<TH1*>(h->Clone(newName));
        histos.Add(hist);
        if (norm == 1) {
          hist->Scale(1. / hist->Integral(), "width");
          hist->GetYaxis()->SetTitle("Probability density");
        }
        else if (norm == 2) {
          hist->Scale(1. / hist->Integral(), "");
          TF1 f1("f1", "TMath::TwoPi() * x", 0, hist->GetXaxis()->GetXmax()*1.5);
          for (Int_t ibin = 1; ibin <= hist->GetNbinsX(); ibin++) {
            Double_t integ = f1.Integral(hist->GetBinLowEdge(ibin), hist->GetBinLowEdge(ibin+1));
            hist->SetBinContent(ibin, hist->GetBinContent(ibin) / integ);
            hist->SetBinError(ibin, hist->GetBinError(ibin) / integ);
          }
          hist->GetYaxis()->SetTitle("Prob. density / 2#pi#Delta R");
        }
        if (rebin > 1) hist->Rebin(4);
        hist->SetTitle(htitle);
        ih++;
      }
    }
  }

  if (ih == 0) return kFALSE;
  
  return Plot1DHistos(cname, histos, xmin, xmax, plotStats);
}

//____________________________________________________________________________________
Bool_t DJetCorrBase::Plot1DHistos(TString cname, TObjArray& histos, Double_t xmin, Double_t xmax, Int_t plotStats)
{
  Color_t colors[12] = {kRed+1, kBlue+1, kMagenta, kGreen+2, kOrange+1, kCyan+2, kPink+1, kTeal+1, kViolet, kAzure+1, kYellow+2, kSpring+3};

  Int_t n = histos.GetEntries();
  if (n == 0) return kFALSE;

  TH1* first = static_cast<TH1*>(histos.At(0));
  
  TString xTitle(first->GetXaxis()->GetTitle());
  TString yTitle(first->GetYaxis()->GetTitle());
  if (xmax - xmin < fgkEpsilon) {
    xmin = first->GetXaxis()->GetXmin();
    xmax = first->GetXaxis()->GetXmax();
  }
  TCanvas *canvas = SetUpCanvas(cname, xTitle, xmin, xmax, kFALSE, yTitle, 0, 1, kFALSE);

  Printf("Getting blank histogram");
  TH1* hblank = dynamic_cast<TH1*>(canvas->GetListOfPrimitives()->At(0));

  Printf("Setting up legend");
  TLegend* leg = 0;
  Int_t neff = n;
  if (plotStats) neff *= 2;
  if (neff > 8) leg = SetUpLegend(0.38, 0.45, 0.88, 0.88, 15);
  else leg = SetUpLegend(0.38, 0.69, 0.88, 0.88, 15);
  leg->SetNColumns(2);
  Double_t maxY = 0;

  Printf("Starting loop on histograms");
  TIter next(&histos);
  TH1* hist = 0;
  Int_t i = 0;
  while ((hist = static_cast<TH1*>(next()))) {
    Printf("Working on histogram %s", hist->GetName());
    if (maxY < hist->GetMaximum()) maxY = hist->GetMaximum();
    hist->SetLineColor(colors[i]);
    hist->SetMarkerColor(colors[i]);
    hist->SetMarkerStyle(kFullCircle);
    hist->SetMarkerSize(0.8);
    hist->Draw("same");
    TLegendEntry* legEntry = leg->AddEntry(hist, hist->GetTitle(), "pe");
    legEntry->SetLineColor(colors[i]);
    if (plotStats) {
      TString statsStr(Form("#mu = %.3f, #sigma = %.3f", hist->GetMean(), hist->GetRMS()));
      leg->AddEntry((TObject*)0, statsStr, "");
    }
    i++;
  }

  if (hblank && maxY > 0) {
    hblank->GetYaxis()->SetRangeUser(0, maxY*1.6);
  }

  leg->Draw();

  if (fSavePlots) {
    SavePlot(canvas);
  }
  
  return kTRUE;
}

//____________________________________________________________________________________
Int_t DJetCorrBase::GetAxisIndex(TString title, THnSparse* hn, Bool_t messageOnFail)
{
  TMap* axisMap = static_cast<TMap*>(fTHnSparseAxisMaps.FindObject(hn->GetName()));
  if (!axisMap) {
    axisMap = GenerateAxisMap(hn);
    if (!axisMap) return -1;
  }

  TParameter<Int_t>* par = static_cast<TParameter<Int_t>*>(axisMap->GetValue(title));
  if (!par) {
    if (messageOnFail) Printf("Warning-DJetCorrBase::GetAxisIndex : Could not find axis with title '%s' in histogram '%s'.", title.Data(), hn->GetName());
    return -1;
  }
  
  return par->GetVal();
}

//____________________________________________________________________________________
Bool_t DJetCorrBase::GenerateRatios(const char* nname, const char* dname)
{
  TIter next(fOutputList);
  TObject* obj = 0;

  while ((obj = next())) {
    TH1* hist1D = dynamic_cast<TH1*>(obj);
    if (!hist1D) continue;

    TString hdname(hist1D->GetName());
    if (!hdname.Contains(dname)) continue;

    TString hnname(hdname);
    hnname.ReplaceAll(dname, nname);

    TH1* num = dynamic_cast<TH1*>(fOutputList->FindObject(hnname));

    if (!num) continue;

    TString rname(Form("Ratio_%s_%s", hnname.Data(), hdname.Data()));

    if (fOutputList->Contains(rname)) continue;

    Printf("Info-DJetCorrBase::GenerateRatios : Now calcutaling ratio '%s' over '%s'", hnname.Data(), hdname.Data());

    if (hist1D->GetSumw2N() == 0) hist1D->Sumw2();
    if (num->GetSumw2N() == 0) num->Sumw2();

    TH2* hist2D = dynamic_cast<TH2*>(hist1D);
    if (hist2D) {
      TH2* num2D = dynamic_cast<TH2*>(num);
      if (!num2D) continue;
      TH2* ratio2D = static_cast<TH2*>(num2D->Clone(rname));
      ratio2D->Divide(hist2D);
      fOutputList->Add(ratio2D);
    }
    else {
      TGraphAsymmErrors* graph = new TGraphAsymmErrors(num, hist1D);
      graph->SetName(rname);
      fOutputList->Add(graph);
    }
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrBase::OpenInputFile()
{
  // Open the input file and set the input list.
  
  TString fname;
  
  fname = fInputPath;
  fname += "/";
  fname += fTrainName;
  fname += "/";
  fname += fInputFileName;

  if (fInputDirectoryFile == fInputFile) {
    fInputDirectoryFile = 0; // to avoid double delete
  }

  if (!fInputFile) {
    Printf("Info-DJetCorrBase::OpenInputFile : Opening file '%s'", fname.Data()); 
    fInputFile = TFile::Open(fname);
  }
  else {
    Printf("Info-DJetCorrBase::OpenInputFile : File '%s' already open", fInputFile->GetName()); 
  }

  if (!fInputFile || fInputFile->IsZombie()) {
    Printf("Error-DJetCorrBase::OpenInputFile : Could not open file '%s'", fname.Data()); 
    fInputFile = 0;
    return kFALSE;
  }

  if (!fInputDirectoryFile || fInputDirFileName != fInputDirectoryFile->GetName()) {
    if (fInputDirFileName.IsNull()) {
      fInputDirectoryFile = fInputFile;
    }
    else {
      ClearInputData();
    
      delete fInputDirectoryFile;
      Printf("Info-DJetCorrBase::OpenInputFile : Getting directory '%s' from file '%s'", fInputDirFileName.Data(), fInputFile->GetName()); 
      fInputDirectoryFile = dynamic_cast<TDirectoryFile*>(fInputFile->Get(fInputDirFileName));

      delete fInputList;
      fInputList = 0;

      if (!fInputDirectoryFile) {
        Printf("Error-DJetCorrBase::OpenInputFile : Could not get directory '%s' from file '%s'", fInputDirFileName.Data(), fInputFile->GetName()); 
        return kFALSE;
      }
    }
  }

  Printf("Info-DJetCorrBase::OpenInputFile : Success.");

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrBase::LoadInputList(const char* inputListName)
{
  // Load the input list.

  if (!OpenInputFile()) return kFALSE;
  
  if (!fInputList || strcmp(inputListName, fInputList->GetName()) != 0) {
    ClearInputData();
    
    delete fInputList;
    Printf("Info-DJetCorrBase::OpenInputFile : Getting list '%s' from directory '%s' of file '%s'", inputListName, fInputDirectoryFile->GetName(), fInputFile->GetName()); 
    fInputList = dynamic_cast<TList*>(fInputDirectoryFile->Get(inputListName));
  }
  
  if (!fInputList) {
    Printf("Error-DJetCorrBase::OpenInputFile : Could not get list '%s' from directory '%s' of file '%s'", inputListName, fInputDirectoryFile->GetName(), fInputFile->GetName());
    return kFALSE;
  }

  Printf("Info-DJetCorrBase::LoadInputList : Success.");

  GetEvents();

  Printf("Info-DJetCorrBase::LoadInputList : Total number of events: %.0f.", fEvents);

  return kTRUE;
}

//____________________________________________________________________________________
void DJetCorrBase::CloseInputFile()
{
  // Close the input file.

  ClearInputData();
  
  if (fInputFile) {
    fInputFile->Close();
    delete fInputFile;
    fInputFile = 0;
  }

  if (fInputDirectoryFile) {
    delete fInputDirectoryFile;
    fInputDirectoryFile = 0;
  }

  if (fInputList) {
    delete fInputList;
    fInputList = 0;
  }
}

//____________________________________________________________________________________
Bool_t DJetCorrBase::SaveOutputFile()
{
  // Save the output in a file.

  TObjArray arr;
  arr.SetOwner(kFALSE);

  fOutputList->SetName("fOutputList");
    
  arr.Add(fOutputList);

  return SaveOutputFile(arr);
}

//____________________________________________________________________________________
TString DJetCorrBase::GetOutpuFileName() const
{
  TString fname(fOutputPath);

  if (fTrainName.IsNull()) {
    fname += "/test/";
  }
  else {
    fname += "/";
    fname += fTrainName;
    fname += "/";
  }

  fname += fOutputFileName;

  return fname;
}

//____________________________________________________________________________________
Bool_t DJetCorrBase::SaveOutputFile(TObjArray& arr)
{
  // Save the content arr in a file.

  TString opt("create");
  if (fOverwrite) opt = "recreate";

  TString fname = GetOutpuFileName();

  TString path(fname);
  path.Remove(path.Last('/'));
  if (gSystem->AccessPathName(path)) {
    Printf("Info-DJetCorrBase::SaveOutputFile : creating directory '%s'.", path.Data()); 
    gSystem->mkdir(path.Data(), kTRUE);
  }
  
  TFile* outputFile = TFile::Open(fname, opt);

  if (!outputFile || outputFile->IsZombie()) {
    Printf("Error-DJetCorrBase::SaveOutputFile : could not open file '%s' to write (if file exists, you may need to set the overwrite option).", fname.Data()); 
    outputFile = 0;
    return kFALSE;
  }

  outputFile->cd();

  Printf("Info-DJetCorrBase::SaveOutputFile : Now streaming results.");
  TIter next(&arr);
  TObject* obj = 0;
  while ((obj = next())) {
    TCollection* coll = dynamic_cast<TCollection*>(obj);
    if (coll) {
      coll->Write(coll->GetName(), TObject::kSingleKey);
    }
    else {
      obj->Write();
    }
  }

  Printf("Info-DJetCorrBase::SaveOutputFile : Closing the output file."); 
  outputFile->Close();
  delete outputFile;
  outputFile = 0;
  
  return kTRUE;
}

//____________________________________________________________________________________
Double_t DJetCorrBase::GetEvents(Bool_t recalculate)
{
  if (fEvents == 0 || recalculate) {
    fEvents = 0;

    if (fOutputList) {

      TH1* hevents = dynamic_cast<TH1*>(fOutputList->FindObject("hEvents"));

      if (!hevents && fInputList) {
        TH1* hevents_temp = static_cast<TH1*>(fInputList->FindObject("fHistEventCount"));
        if (hevents_temp) {
          hevents = static_cast<TH1*>(hevents_temp->Clone("hEvents"));
          hevents->SetTitle("hEvents");
          fOutputList->Add(hevents);
        }
      }
      
      if (hevents) {
        fEvents = hevents->GetBinContent(1);
      }
    }
  }

  return fEvents;
}

//____________________________________________________________________________________
void DJetCorrBase::FitObjectInPad(TObject* obj, TVirtualPad* pad, Option_t* opt, Bool_t copyAxisTitle, Double_t extraFactor)
{
  if (obj->InheritsFrom("TGraph")) {
    FitGraphInPad(static_cast<TGraph*>(obj), pad, opt, copyAxisTitle, extraFactor);
  }
  else if (obj->InheritsFrom("TH1")) {
    FitHistogramInPad(static_cast<TH1*>(obj), pad, opt, copyAxisTitle, extraFactor);
  }
  else {
    Printf("Don't know how to draw this object!");
  }
}

//____________________________________________________________________________________
void DJetCorrBase::FitGraphInPad(TGraph* graph, TVirtualPad* pad, Option_t* opt, Bool_t copyAxisTitle, Double_t extraFactor)
{
  TString strOpt(opt);
  if (!strOpt.Contains("same")) strOpt += "same";

  TH1* blankHist = dynamic_cast<TH1*>(pad->GetListOfPrimitives()->At(0));
  if (blankHist) {
    if (copyAxisTitle) {
      blankHist->GetXaxis()->SetTitle(graph->GetXaxis()->GetTitle());
      blankHist->GetYaxis()->SetTitle(graph->GetYaxis()->GetTitle());
    }

    Double_t miny = blankHist->GetMinimum();
    if (pad->GetLogy()) miny *= extraFactor;
    Double_t maxy = blankHist->GetMaximum() / extraFactor;
  
    GetMinMax(graph, miny, maxy);
    
    if (pad->GetLogy()) {
      if (miny <= 0) miny = blankHist->GetMinimum() / extraFactor;
      if (maxy <= 0) maxy = blankHist->GetMaximum() * extraFactor;

      miny /= extraFactor;
    }
    else {
      miny = 0;
    }

    maxy *= extraFactor;

    blankHist->SetMinimum(miny);
    blankHist->SetMaximum(maxy);

    pad->cd();
    graph->Draw(strOpt);
  }
  else {
    Printf("Error-DJetCorrBase::FitGraphInPad : Could not find blank histogram!");
  }
}

//____________________________________________________________________________________
void DJetCorrBase::FitHistogramInPad(TH1* hist, TVirtualPad* pad, Option_t* opt, Bool_t copyAxisTitle, Double_t extraFactor)
{
  TString strOpt(opt);
  if (!strOpt.Contains("same")) strOpt += "same";

  TH1* blankHist = dynamic_cast<TH1*>(pad->GetListOfPrimitives()->At(0));
  if (blankHist) {
    if (copyAxisTitle) {
      blankHist->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
      blankHist->GetYaxis()->SetTitle(hist->GetYaxis()->GetTitle());
    }

    Double_t miny = blankHist->GetMinimum();
    if (pad->GetLogy()) miny *= extraFactor;
    Double_t maxy = blankHist->GetMaximum() / extraFactor;

    GetMinMax(hist, miny, maxy);

    if (pad->GetLogy()) {
      if (miny <= 0) miny = blankHist->GetMinimum() / extraFactor;
      if (maxy <= 0) maxy = blankHist->GetMaximum() * extraFactor;

      miny /= extraFactor;
    }
    else {
      miny = 0;
    }

    maxy *= extraFactor;

    blankHist->SetMinimum(miny);
    blankHist->SetMaximum(maxy);

    pad->cd();
    hist->Draw(strOpt);
  }
  else {
    Printf("Error-DJetCorrBase::FitHistogramInPad : Could not find blank histogram!");
  }
}

//____________________________________________________________________________________
void DJetCorrBase::GetMinMax(TGraph* graph, Double_t& miny, Double_t& maxy)
{
  Double_t* array = graph->GetY();

  for (Int_t i = 0; i < graph->GetN(); i++) {
    if (miny > graph->GetY()[i] - graph->GetEYlow()[i]) miny = graph->GetY()[i] - graph->GetEYlow()[i];
    if (maxy < graph->GetY()[i] + graph->GetEYhigh()[i]) maxy = graph->GetY()[i] + graph->GetEYhigh()[i];
  } 
}

//____________________________________________________________________________________
void DJetCorrBase::GetMinMax(TH1* hist, Double_t& miny, Double_t& maxy)
{
  Int_t minBin = hist->GetMinimumBin();
  miny = TMath::Min(hist->GetBinContent(minBin) - hist->GetBinError(minBin), miny);

  Int_t maxBin = hist->GetMaximumBin();
  maxy = TMath::Max(hist->GetBinContent(maxBin) + hist->GetBinError(maxBin), maxy);
}

//____________________________________________________________________________________
TLegend* DJetCorrBase::GetLegend(TPad* pad)
{
  TLegend* leg = 0;

  for (Int_t i = 0; i < pad->GetListOfPrimitives()->GetEntries(); i++) {
    leg = dynamic_cast<TLegend*>(pad->GetListOfPrimitives()->At(i));
    if (leg) break;
  }

  return leg;
}

//____________________________________________________________________________________
Bool_t DJetCorrBase::CheckExactRebin(TAxis* orig, TAxis* dest)
{
  for (Int_t i = 1; i <= orig->GetNbins(); i++) {
    Double_t xlow    = orig->GetBinLowEdge(i);
    Double_t xup     = orig->GetBinUpEdge(i);
    Int_t    xlowBin = dest->FindBin(xlow);
    Int_t    xupBin  = dest->FindBin(xup);
    if (TMath::Abs(xup - dest->GetBinLowEdge(xupBin)) < fgkEpsilon && xupBin > 0) {
      xupBin--;
    }
    if (xlowBin != xupBin) {
      Printf("Bin %d = [%.2f, %.2f] -> %d = [%.2f, %.2f], %d = [%.2f, %.2f]", i, xlow, xup, xlowBin, dest->GetBinLowEdge(xlowBin), dest->GetBinUpEdge(xlowBin), xupBin, dest->GetBinLowEdge(xupBin), dest->GetBinUpEdge(xupBin));
      return kFALSE;
    }
  }

  return kTRUE;
}

//____________________________________________________________________________________
void DJetCorrBase::GetBinCenter(THnSparse* hn, Int_t* coord_ind, Double_t* coord)
{
  for (Int_t idim = 0; idim < hn->GetNdimensions(); idim++) {
    coord[idim] = hn->GetAxis(idim)->GetBinCenter(coord_ind[idim]);
  }
}

//____________________________________________________________________________________
THnSparse* DJetCorrBase::Rebin(THnSparse* orig, const char* name, const Int_t* nbins, const Double_t** bins)
{
  THnSparse* result = new THnSparseD(name, orig->GetTitle(), orig->GetNdimensions(), nbins);
  for (Int_t idim = 0; idim < orig->GetNdimensions(); idim++) {
    result->GetAxis(idim)->SetTitle(orig->GetAxis(idim)->GetTitle());
    result->GetAxis(idim)->Set(nbins[idim], bins[idim]);
    if (!CheckExactRebin(orig->GetAxis(idim), result->GetAxis(idim))) {
      Printf("WARNING: unable to exact rebin axis %s of histogram %s", orig->GetAxis(idim)->GetTitle(), name);
    }
  }

  Int_t* coord_ind = new Int_t[orig->GetNdimensions()];
  Double_t* coord = new Double_t[orig->GetNdimensions()];
    
  for (Int_t ibin = 0; ibin < orig->GetNbins(); ibin++) {
    Double_t content = orig->GetBinContent(ibin, coord_ind);
    GetBinCenter(orig, coord_ind, coord);
    Int_t newbin = result->GetBin(coord);
    result->AddBinContent(newbin, content);
    result->AddBinError2(newbin, orig->GetBinError2(ibin));
  }

  delete[] coord_ind;
  delete[] coord;

  return result;
}

//____________________________________________________________________________________
TH1* DJetCorrBase::Rebin(TH1* orig, const char* name, Int_t nbins, const Double_t* bins)
{
  TH1* dest = new TH1D(name, orig->GetTitle(), nbins, bins);

  if (!CheckExactRebin(orig->GetXaxis(), dest->GetXaxis())) Printf("WARNING: unable to exact rebin axis %s of histogram %s", orig->GetXaxis()->GetTitle(), name);
  dest->GetXaxis()->SetTitle(orig->GetXaxis()->GetTitle());
  
  for (Int_t xBin = 0; xBin <= orig->GetNbinsX()+1; xBin++) {
    Double_t xVal     = orig->GetXaxis()->GetBinCenter(xBin);
    Int_t    xDestBin = dest->GetXaxis()->FindBin(xVal);
    
    Double_t content = orig->GetBinContent(xBin);
    Double_t err     = orig->GetBinError(xBin);

    dest->SetBinContent(xDestBin, dest->GetBinContent(xDestBin)+content);
    dest->SetBinError(xDestBin, TMath::Sqrt(dest->GetBinError(xDestBin)*dest->GetBinError(xDestBin)+err*err));
  }

  return dest;
}

//____________________________________________________________________________________
TH2* DJetCorrBase::Rebin(TH2* orig, const char* name, Int_t nbinsx, const Double_t* binsx, Int_t nbinsy, const Double_t* binsy)
{
  TH2* dest = new TH2D(name, orig->GetTitle(), nbinsx, binsx, nbinsy, binsy);

  if (!CheckExactRebin(orig->GetXaxis(), dest->GetXaxis())) Printf("WARNING: unable to exact rebin axis %s of histogram %s", orig->GetXaxis()->GetTitle(), name);
  if (!CheckExactRebin(orig->GetYaxis(), dest->GetYaxis())) Printf("WARNING: unable to exact rebin axis %s of histogram %s", orig->GetYaxis()->GetTitle(), name);

  dest->GetXaxis()->SetTitle(orig->GetXaxis()->GetTitle());
  dest->GetYaxis()->SetTitle(orig->GetYaxis()->GetTitle());
  
  for (Int_t xBin = 0; xBin <= orig->GetNbinsX(); xBin++) {
    Double_t xVal     = orig->GetXaxis()->GetBinCenter(xBin);
    Int_t    xDestBin = dest->GetXaxis()->FindBin(xVal);

    for (Int_t yBin = 0; yBin <= orig->GetNbinsY(); yBin++) {
      Double_t yVal     = orig->GetYaxis()->GetBinCenter(yBin);
      Int_t    yDestBin = dest->GetYaxis()->FindBin(yVal);
    
      Double_t content = orig->GetBinContent(xBin, yBin);
      Double_t err     = orig->GetBinError(xBin, yBin);

      dest->SetBinContent(xDestBin, yDestBin, dest->GetBinContent(xDestBin, yDestBin)+content);
      dest->SetBinError(xDestBin, yDestBin, TMath::Sqrt(dest->GetBinError(xDestBin, yDestBin)*dest->GetBinError(xDestBin, yDestBin)+err*err));
    }
  }

  return dest;
}

//____________________________________________________________________________________
TH1* DJetCorrBase::Rebin(TH1* orig, const char* name, Int_t nbins, Double_t min, Double_t max)
{
  Double_t* bins = GenerateFixedArray(nbins, min, max);
  TH1* res = Rebin(orig, name, nbins, bins);
  delete[] bins;
  return res;
}

//____________________________________________________________________________________
TH2* DJetCorrBase::Rebin(TH2* orig, const char* name, Int_t nbinsx, Double_t minX, Double_t maxX, Int_t nbinsy, Double_t minY, Double_t maxY)
{
  Double_t* binsX = GenerateFixedArray(nbinsx, minX, maxX);
  Double_t* binsY = GenerateFixedArray(nbinsy, minY, maxY);

  TH2* res = Rebin(orig, name, nbinsx, binsX, nbinsy, binsY);
  delete[] binsX;
  delete[] binsY;
  return res;
}

//____________________________________________________________________________________
Double_t* DJetCorrBase::GenerateFixedArray(Int_t n, Double_t min, Double_t max)
{
  Double_t* res = new Double_t[n+1];
  Double_t binWidth = (max-min) / n;
  Printf("Max = %.1f, min = %.1f, bins = %d, binWidth = %.1f", max, min, n, binWidth);
  res[0] = min;
  for (Int_t i = 1; i <= n; i++) {
    res[i] = res[i-1] + binWidth;
  }
  return res;
}

//____________________________________________________________________________________
TH2* DJetCorrBase::GetTruth(Int_t p, Bool_t copy)
{
  TString hname = GetTruthName(p);

  if (hname.IsNull()) return 0;

  TH2* hist = dynamic_cast<TH2*>(GetOutputHistogram(hname));

  if (copy && hist) {
    TString hname = hist->GetName();
    hname += "_copy";

    hist = static_cast<TH2*>(hist->Clone(hname));
  }

  if (!hist) Printf("Histogram %s not found!", hname.Data());

  return hist;
}

//____________________________________________________________________________________
TH1* DJetCorrBase::GetDzTruth(Int_t p, Bool_t copy)
{
  TString hname = GetDzTruthName(p);

  if (hname.IsNull()) return 0;

  TH2* hist = dynamic_cast<TH2*>(GetOutputHistogram(hname));

  if (copy && hist) {
    TString hname = hist->GetName();
    hname += "_copy";

    hist = static_cast<TH2*>(hist->Clone(hname));
  }

  if (!hist) Printf("Histogram %s not found!", hname.Data());

  return hist;
}

//____________________________________________________________________________________
TH1* DJetCorrBase::GetDzMeasured(Int_t p, Bool_t copy)
{
  TString hname = GetDzMeasuredName(p);

  if (hname.IsNull()) return 0;

  TH2* hist = dynamic_cast<TH2*>(GetOutputHistogram(hname));

  if (copy && hist) {
    TString hname = hist->GetName();
    hname += "_copy";

    hist = static_cast<TH2*>(hist->Clone(hname));
  }

  if (!hist) Printf("Histogram %s not found!", hname.Data());

  return hist;
}

//____________________________________________________________________________________
TH2* DJetCorrBase::GetMeasured(Int_t p, Bool_t copy)
{
  TString hname = GetMeasuredName(p);

  if (hname.IsNull()) return 0;

  TH2* hist = dynamic_cast<TH2*>(GetOutputHistogram(hname));

  if (copy && hist) {
    TString hname = hist->GetName();
    hname += "_copy";

    hist = static_cast<TH2*>(hist->Clone(hname));
  }

  if (!hist) Printf("Histogram %s not found!", hname.Data());

  return hist;
}

//____________________________________________________________________________________
TH1* DJetCorrBase::GetDPtTruth(Int_t p, Bool_t copy, const char* matching)
{
  TString hname = GetDPtTruthName(p, matching);

  if (hname.IsNull()) return 0;

  TH1* hist = dynamic_cast<TH1*>(GetOutputHistogram(hname));

  if (copy && hist) {
    TString hname = hist->GetName();
    hname += "_copy";

    hist = static_cast<TH1*>(hist->Clone(hname));
  }

  if (!hist) Printf("Histogram %s not found!", hname.Data());

  return hist;
}

//____________________________________________________________________________________
TH1* DJetCorrBase::GetDPtMeasured(Int_t p, Bool_t copy, const char* matching)
{
  TString hname = GetDPtMeasuredName(p, matching);

  if (hname.IsNull()) return 0;

  TH1* hist = dynamic_cast<TH1*>(GetOutputHistogram(hname));

  if (copy && hist) {
    TString hname = hist->GetName();
    hname += "_copy";

    hist = static_cast<TH1*>(hist->Clone(hname));
  }

  if (!hist) Printf("Histogram %s not found!", hname.Data());

  return hist;
}

//____________________________________________________________________________________
TH1* DJetCorrBase::GetDEtaTruth(Int_t p, Bool_t copy, const char* matching)
{
  TString hname = GetDEtaTruthName(p, matching);

  if (hname.IsNull()) return 0;

  TH1* hist = dynamic_cast<TH1*>(GetOutputHistogram(hname));

  if (copy && hist) {
    TString hname = hist->GetName();
    hname += "_copy";

    hist = static_cast<TH1*>(hist->Clone(hname));
  }

  if (!hist) Printf("Histogram %s not found!", hname.Data());

  return hist;
}

//____________________________________________________________________________________
TH1* DJetCorrBase::GetDEtaMeasured(Int_t p, Bool_t copy, const char* matching)
{
  TString hname = GetDEtaMeasuredName(p, matching);

  if (hname.IsNull()) return 0;

  TH1* hist = dynamic_cast<TH1*>(GetOutputHistogram(hname));

  if (copy && hist) {
    TString hname = hist->GetName();
    hname += "_copy";

    hist = static_cast<TH1*>(hist->Clone(hname));
  }

  if (!hist) Printf("Histogram %s not found!", hname.Data());

  return hist;
}

//____________________________________________________________________________________
TH1* DJetCorrBase::GetJetPtTruth(Int_t p, Bool_t copy)
{
  TString hname = GetJetPtTruthName(p);

  if (hname.IsNull()) return 0;

  TH1* hist = dynamic_cast<TH1*>(GetOutputHistogram(hname));

  if (copy && hist) {
    TString hname = hist->GetName();
    hname += "_copy";

    hist = static_cast<TH1*>(hist->Clone(hname));
  }

  if (!hist) Printf("Histogram %s not found!", hname.Data());

  return hist;
}

//____________________________________________________________________________________
TH1* DJetCorrBase::GetJetPtMeasured(Int_t p, Bool_t copy)
{
  TString hname = GetJetPtMeasuredName(p);

  if (hname.IsNull()) return 0;

  TH1* hist = dynamic_cast<TH1*>(GetOutputHistogram(hname));

  if (copy && hist) {
    TString hname = hist->GetName();
    hname += "_copy";

    hist = static_cast<TH1*>(hist->Clone(hname));
  }

  if (!hist) Printf("Histogram %s not found!", hname.Data());

  return hist;
}

//____________________________________________________________________________________
void DJetCorrBase::MakeBinomialConsistent(TH1* pass, TH1* total)
{
  Int_t nbinsx, nbinsy, nbinsz, nbins;

  nbinsx = pass->GetNbinsX();
  nbinsy = pass->GetNbinsY();
  nbinsz = pass->GetNbinsZ();

  switch(pass->GetDimension()) {
    case 1: nbins = nbinsx + 2; break;
    case 2: nbins = (nbinsx + 2) * (nbinsy + 2); break;
    case 3: nbins = (nbinsx + 2) * (nbinsy + 2) * (nbinsz + 2); break;
    default: nbins = 0;
  }

  for(Int_t i = 0; i < nbins; ++i) {
    if (pass->GetBinContent(i) == 0 || total->GetBinContent(i) == 0) {
      pass->SetBinContent(i, 0);
      total->SetBinContent(i, 0);
    }

    if(pass->GetBinContent(i) > total->GetBinContent(i)) {
      pass->SetBinContent(i, total->GetBinContent(i));
    }
  }
}

//____________________________________________________________________________________
const char* DJetCorrBase::GetParamName(Int_t i) const
{
  TObject* obj = fAnalysisParams->At(i);
  if (obj) {
    return obj->GetName();
  }
  else {
    return 0;
  }
}

//____________________________________________________________________________________
TH2* DJetCorrBase::Normalize(TH2* orig, const char* name)
{
  TH2* res = static_cast<TH2*>(orig->Clone(name));

  for (Int_t y = 0; y <= res->GetNbinsY(); y++) {
    Double_t integral = res->Integral(0, -1, y, y);
    if (integral == 0.) continue;
    for (Int_t x = 0; x <= res->GetNbinsX(); x++) {
      res->SetBinContent(x, y, res->GetBinContent(x, y) / integral);
      res->SetBinError(x, y, res->GetBinError(x, y) / integral);
    }
  }

  return res;
}

//____________________________________________________________________________________
void DJetCorrBase::UpdateAllCanvases()
{
  TIter next(fCanvases);
  TVirtualPad* pad = 0;

  while ((pad = static_cast<TVirtualPad*>(next()))) {
    pad->Update();
    gSystem->ProcessEvents();
  }
}
