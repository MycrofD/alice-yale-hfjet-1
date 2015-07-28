// Helper class to set the style of histograms.
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

#include <TAttLine.h>
#include <TAttFill.h>
#include <TAttMarker.h>

class TH1;

class HistoStyler : public TAttLine, public TAttFill, public TAttMarker {
  
 public:
  enum EMarkerColor { kMarkerColorFixed, kMarkerColorVariable };
  enum EMarkerStyle { kMarkerStyleFixed, kMarkerStyleVariable , kMarkerStyleOpen, kMarkerStyleFull };
  enum ELineColor   { kLineColorFixed  , kLineColorVariable   };
  enum ELineStyle   { kLineStyleFixed  , kLineStyleVariable   };
  enum EFillColor   { kFillColorFixed  , kFillColorVariable   };
  enum EFillStyle   { kFillStyleFixed  , kFillStyleVariable   };
  
  HistoStyler();

  void       SetVariableMarkerColor(EMarkerColor v=kMarkerColorVariable)   { fVariableMarkerColor  = v; }
  void       SetVariableMarkerStyle(EMarkerStyle v=kMarkerStyleOpen)       { fVariableMarkerStyle  = v; }

  void       SetVariableLineColor(ELineColor v=kLineColorVariable)         { fVariableLineColor    = v; }
  void       SetVariableLineStyle(ELineStyle v=kLineStyleVariable)         { fVariableLineStyle    = v; }
  
  void       SetVariableFillColor(EFillColor v=kFillColorVariable)         { fVariableFillColor    = v; }
  void       SetVariableFillStyle(EFillStyle v=kFillStyleVariable)         { fVariableFillStyle    = v; }

  using TAttMarker::GetMarkerColor;
  Color_t    GetMarkerColor(UShort_t i) const;

  using TAttMarker::GetMarkerStyle;
  Style_t    GetMarkerStyle(UShort_t i) const;

  using TAttLine::GetLineColor;
  Color_t    GetLineColor(UShort_t i) const;
  
  using TAttLine::GetLineStyle;
  Style_t    GetLineStyle(UShort_t i) const;

  using TAttFill::GetFillColor;
  Color_t    GetFillColor(UShort_t i) const;

  using TAttFill::GetFillStyle;
  Style_t    GetFillStyle(UShort_t i) const;

  void       Apply(Int_t n, TH1** histo) const;

  static const Int_t fgkNColors = 10;
  static const Int_t fgkNOpenMarkerStyles = 6;
  static const Int_t fgkNFullMarkerStyles = 6;
  
 protected:
  EMarkerColor     fVariableMarkerColor;
  EMarkerStyle     fVariableMarkerStyle;
  ELineColor       fVariableLineColor;
  ELineStyle       fVariableLineStyle;
  EFillColor       fVariableFillColor;
  EFillStyle       fVariableFillStyle;

  static const Color_t fgkColors[fgkNColors];
  static const Color_t fgkOpenMarkerStyles[fgkNOpenMarkerStyles];
  static const Color_t fgkFullMarkerStyles[fgkNFullMarkerStyles];
  
 private:
   
  HistoStyler(const HistoStyler &source);
  HistoStyler& operator=(const HistoStyler& source); 

  ClassDef(HistoStyler, 1);
};
