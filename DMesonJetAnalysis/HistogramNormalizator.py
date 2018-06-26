from enum import Enum
import math
import DMesonJetUtils

CROSS_SECTION = 62.2  # mb CINT1
BRANCHING_RATIO = 0.0393  # D0->Kpi

class NormalizationType(Enum):
    Unknown = 0
    CrossSection = 1
    Distribution = 2

    @classmethod
    def from_string(cls, norm_type_str):
        if norm_type_str == "cross_section":
            return cls.CrossSection
        elif norm_type_str == "distribution":
            return cls.Distribution
        else:
            print("Unknown normalization option '{}'".format(norm_type_str))
            return cls.Unknown

class Normalizator:
    def __init__(self, h, normalization):
        self.fHistogram = h
        self.fNormalizationType = NormalizationType.from_string(normalization)
        self.fNormalizationOption = ""
        self.fNormalizedHistogram = None
        print("Normalizator '{}' for histogram '{}'".format(self.fNormalizationType, self.fHistogram.GetName()))

    def GetCrossSectionNormalization(self):
        return 1.0

    def CalculateIntegralAndError(self, exclude_bin):
        integral = 0
        error2 = 0
        for ibin in range(1, self.fHistogram.GetNbinsX() + 1):
            if ibin == exclude_bin:
                continue
            integral += self.fHistogram.GetBinContent(ibin)
            error2 += self.fHistogram.GetBinError(ibin) ** 2
        return integral, error2

    def NormalizeHistogram(self):
        print("Performing normalization '{}' for histogram '{}'".format(self.fNormalizationType, self.fHistogram.GetName()))
        if self.fNormalizationType == NormalizationType.CrossSection:
            self.fNormalizedHistogram = self.fHistogram.Clone(self.fHistogram.GetName())
            print("Scaling by {} with option '{}'".format(self.GetCrossSectionNormalization(), self.fNormalizationOption))
            self.fNormalizedHistogram.Scale(self.GetCrossSectionNormalization(), self.fNormalizationOption)
        elif self.fNormalizationType == NormalizationType.Distribution:
            self.fNormalizedHistogram = DMesonJetUtils.soft_clone(self.fHistogram, self.fHistogram.GetName())
            self.fNormalizedHistogram.Sumw2()

            tot_integral, tot_integral_err2 = self.CalculateIntegralAndError(-1)
            tot_integral_err = math.sqrt(tot_integral_err2)
            print("\n\n")
            print("Total integral is {} +/- {} (rel err {})".format(tot_integral, tot_integral_err, tot_integral_err / tot_integral))
            for ibin in range(1, self.fHistogram.GetNbinsX() + 1):
                print("\n")
                print("Bin {}: {} +/- {} (rel err {})".format(ibin, self.fHistogram.GetBinContent(ibin), self.fHistogram.GetBinError(ibin), self.fHistogram.GetBinError(ibin) / self.fHistogram.GetBinContent(ibin)))
                if self.fHistogram.GetBinContent(ibin) == 0:
                    continue
                I_minus_A, I_minus_A_err2 = self.CalculateIntegralAndError(ibin)
                if I_minus_A == 0:
                    break
                I_minus_A_over_A = I_minus_A / self.fHistogram.GetBinContent(ibin)
                I_minus_A_over_A_err2 = (I_minus_A_err2 / (I_minus_A ** 2) + (self.fHistogram.GetBinError(ibin) ** 2) / (self.fHistogram.GetBinContent(ibin) ** 2)) * (I_minus_A_over_A ** 2)
                I_minus_A_over_A_err = math.sqrt(I_minus_A_over_A_err2)
                I_over_A = 1 + I_minus_A_over_A
                I_over_A_err2 = I_minus_A_over_A_err2
                I_over_A_err = math.sqrt(I_over_A_err2)
                A_over_I = 1 / I_over_A
                A_over_I_err = I_over_A_err / I_over_A * A_over_I
                self.fNormalizedHistogram.SetBinContent(ibin, A_over_I)
                self.fNormalizedHistogram.SetBinError(ibin, A_over_I_err)

                print("I_minus_A_over_A = {} +/- {} (rel err {})".format(I_minus_A_over_A, I_minus_A_over_A_err, I_minus_A_over_A_err / I_minus_A_over_A))
                print("I_over_A = {} +/- {} ({})".format(I_over_A, I_over_A_err, tot_integral / self.fHistogram.GetBinContent(ibin)))
                print("A_over_I = {} +/- {} (rel err {})".format(A_over_I, A_over_I_err, A_over_I_err / A_over_I))
                print("Old calculation: A_over_I = {} +/- {} (rel err {})".format(self.fHistogram.GetBinContent(ibin) / tot_integral, self.fHistogram.GetBinError(ibin) / tot_integral, self.fHistogram.GetBinError(ibin) / self.fHistogram.GetBinContent(ibin)))
            
            if self.fNormalizationOption:
                self.fNormalizedHistogram.Scale(1.0, self.fNormalizationOption)
        else:
            print("Unknown normalization option '{}'".format(self.fNormalizationType))
            return None

        return self.fNormalizedHistogram

class DataNormalizator(Normalizator):
    def __init__(self, h, normalization, events):
        Normalizator.__init__(self, h, normalization)
        self.fEvents = events
        self.fNormalizationOption = "width"
        self.GetCrossSectionNormalization = self.GetCrossSectionNormalizationData

    def GetCrossSectionNormalizationData(self):
        return CROSS_SECTION / (self.fEvents * BRANCHING_RATIO)

class MCSimulationNormalizator(Normalizator):
    def __init__(self, h, normalization):
        Normalizator.__init__(self, h, normalization)
        self.fNormalizationOption = "width"
        self.GetCrossSectionNormalization = self.GetCrossSectionNormalizationMCSimulation

    def GetCrossSectionNormalizationMCSimulation(self):
        return 1.0
