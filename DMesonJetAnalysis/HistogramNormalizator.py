from enum import Enum
import math
import DMesonJetUtils

CROSS_SECTION = 62.2  # mb CINT1
BRANCHING_RATIO = 0.0389  # D0->Kpi

class NormalizationType(Enum):
    Unknown = 0
    CrossSection = 1
    Distribution = 2
    Ratio = 3
    Rate = 4

    @classmethod
    def from_string(cls, norm_type_str):
        if norm_type_str in cls.__members__:
            return cls[norm_type_str]
        elif norm_type_str == "cross_section":
            return cls.CrossSection
        elif norm_type_str == "distribution":
            return cls.Distribution
        elif norm_type_str == "ratio":
            return cls.Ratio
        elif norm_type_str == "rate":
            return cls.Rate
        else:
            print("Unknown normalization option '{}'".format(norm_type_str))
            return cls.Unknown

class Normalizator:
    def __init__(self, h, normalization):
        self.fHistogram = h
        self.fGraph = None
        self.fNormalizationType = NormalizationType.from_string(normalization)
        self.fNormalizationOption = ""
        self.fNormalizedHistogram = None
        self.fNormalizedGraph = None
        self.fScale = 1.0
        self.fXmin = 0
        self.fXmax = -1
        self.GetCrossSectionNormalization = self._GetCrossSectionNormalization
        
        self.fNormalizationHistogram = None
        self.fNormalizationGraph = None
        self.fNormalizationXmin = 0
        self.fNormalizationXmax = -1
        self.fTotalCrossSection = None
        self.fTotalCrossSectionStatUnc = 0
        self.fTotalCrossSectionSystUnc = 0

        print("Normalizator '{}' for histogram '{}'".format(self.fNormalizationType, self.fHistogram.GetName()))

    def _GetCrossSectionNormalization(self):
        return self.fScale, 0, 0

    def CalculateIntegralAndError(self, exclude_bin):
        integral = 0
        error2 = 0
        for ibin in range(1, self.fHistogram.GetNbinsX() + 1):
            if ibin == exclude_bin:
                continue
            if self.fXmin < self.fXmax:
                if self.fHistogram.GetXaxis().GetBinCenter(ibin) < self.fXmin or self.fHistogram.GetXaxis().GetBinCenter(ibin) > self.fXmax:
                    continue
            integral += self.fHistogram.GetBinContent(ibin)
            error2 += self.fHistogram.GetBinError(ibin) ** 2
        return integral, error2

    def NormalizeHistogramCrossSection(self):
        scaling, scaling_stat, scaling_syst = self.GetCrossSectionNormalization()
        print("Scaling by {} +/- {} (stat) +/- {} (syst) with option '{}'".format(scaling, scaling_stat*scaling, scaling_syst*scaling, self.fNormalizationOption))
        self.fNormalizedHistogram = self.fHistogram.Clone(self.fHistogram.GetName())
        if self.fXmin < self.fXmax:
            for ibin in range(1, self.fHistogram.GetNbinsX() + 1):
                if self.fNormalizedHistogram.GetXaxis().GetBinCenter(ibin) < self.fXmin or self.fNormalizedHistogram.GetXaxis().GetBinCenter(ibin) > self.fXmax:
                    self.fNormalizedHistogram.SetBinContent(ibin, 0)
                    self.fNormalizedHistogram.SetBinError(ibin, 0)
        if scaling_stat:
            for ibin in range(1, self.fNormalizedHistogram.GetNbinsX() + 1):
                if self.fNormalizedHistogram.GetBinContent(ibin) == 0:
                    continue
                stat_err = math.sqrt((self.fNormalizedHistogram.GetBinError(ibin) / self.fNormalizedHistogram.GetBinContent(ibin)) ** 2 + scaling_stat ** 2) *  self.fNormalizedHistogram.GetBinContent(ibin)
                self.fNormalizedHistogram.SetBinError(ibin, stat_err)
        self.fNormalizedHistogram.Scale(scaling, self.fNormalizationOption)
        if self.fGraph:
            self.fNormalizedGraph =  self.fGraph.Clone(self.fGraph.GetName())
            for ipoint in range(0, self.fGraph.GetN()):
                if self.fGraph.GetX()[ipoint] < self.fGraph.GetX()[ipoint] > self.fXmax:
                    self.fNormalizedGraph.RemovePoint(ipoint)
            for ipoint in range(0, self.fNormalizedGraph.GetN()):
                if self.fNormalizationOption == "width":
                    scaling *= self.fNormalizedGraph.GetErrorX(ipoint) * 2
                syst_err = math.sqrt((self.fNormalizedGraph.GetErrorY(ipoint) / self.fNormalizedGraph.GetY()[ipoint]) ** 2 + scaling_syst ** 2) * self.fNormalizedGraph.GetY()[ipoint]
                self.fNormalizedGraph.SetPointError(ipoint, self.fNormalizedGraph.GetErrorX(ipoint), self.fNormalizedGraph.GetErrorX(ipoint), syst_err * scaling, syst_err * scaling)
                self.fNormalizedGraph.SetPoint(ipoint, self.fNormalizedGraph.GetX()[ipoint],  self.fNormalizedGraph.GetY()[ipoint] * scaling)

    def NormalizeHistogramDistribution(self):
        if self.fGraph:
            print("NormalizeHistogramDistribution not implemented for systematic uncertainties.")
            exit(1)

        self.fNormalizedHistogram = DMesonJetUtils.soft_clone(self.fHistogram, self.fHistogram.GetName())
        self.fNormalizedHistogram.Sumw2()

        tot_integral, tot_integral_err2 = self.CalculateIntegralAndError(-1)
        tot_integral_err = math.sqrt(tot_integral_err2)
        print("\n")
        print("Total integral is {} +/- {} (rel err {})".format(tot_integral, tot_integral_err, tot_integral_err / tot_integral))
        for ibin in range(1, self.fHistogram.GetNbinsX() + 1):
            print("\n")
            print("Bin {}: {} +/- {} (rel err {})".format(ibin, self.fHistogram.GetBinContent(ibin), self.fHistogram.GetBinError(ibin), self.fHistogram.GetBinError(ibin) / self.fHistogram.GetBinContent(ibin)))
            if self.fXmin < self.fXmax:
                if self.fHistogram.GetXaxis().GetBinCenter(ibin) < self.fXmin or self.fHistogram.GetXaxis().GetBinCenter(ibin) > self.fXmax:
                    print("Skipping this bin because outside of requested range [{}, {}]".format(self.fXmin, self.fXmax))
                    continue
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

    def NormalizeHistogramRatio(self):
        if not self.fNormalizationHistogram:
            print("For normalization type '{}' the field 'fNormalizationHistogram' must be set!".format(self.fNormalizationType))
            exit(1)
        
        if self.fHistogram.GetNbinsX() != self.fNormalizationHistogram.GetNbinsX():
            print("Error: graph '{}' has {} points, while histogram '{}' has {} bins.".format(self.fNormalizationHistogram.GetName(), self.fNormalizationHistogram.GetNbinsX(), self.fHistogram.GetName(), self.fHistogram.GetNbinsX()))
            exit(1)

        self.fNormalizedHistogram = self.fHistogram.Clone(self.fHistogram.GetName())
        if self.fGraph:
            if not self.fNormalizationGraph:
                print("Did not provide normalization graph (syst unc) for histogram '{}'".format(self.fHistogram.GetName()))
                exit(1)
            self.fNormalizedGraph = self.fGraph.Clone(self.fGraph.GetName())
        else:
            self.fNormalizedGraph = None

        xsec_incl_tot = 0.0
        stat_xsec_den_tot2 = 0.0
        syst_xsec_den_tot = 0.0

        xsec_d0_tot = 0.0
        stat_xsec_num_tot2 = 0.0
        syst_xsec_num_tot = 0.0

        avg_tot = 0.0
        stat_avg_tot = 0.0
        syst_avg_tot = 0.0

        avg_flat = 0.0
        sum_of_w_avg_flat = 0.0
        stat_avg_flat2 = 0.0
        syst_avg_flat2 = 0.0

        for ibin in range(0, self.fNormalizedHistogram.GetNbinsX()):
            if self.fNormalizedGraph:
                self.fNormalizedGraph.SetPoint(ibin, self.fNormalizedGraph.GetX()[ibin], self.fGraph.GetY()[ibin] / self.fNormalizationGraph.GetY()[ibin])
                syst_erry2 = ((self.fGraph.GetErrorY(ibin) / self.fGraph.GetY()[ibin]) ** 2 + (self.fNormalizationGraph.GetErrorY(ibin) / self.fNormalizationGraph.GetY()[ibin]) ** 2) * self.fNormalizedGraph.GetY()[ibin] ** 2
                syst_erry = math.sqrt(syst_erry2)
                syst_xsec_num_tot += self.fGraph.GetErrorY(ibin)  # take the weithed average of the rel unc
                self.fNormalizedGraph.SetPointError(ibin, self.fNormalizedGraph.GetErrorX(ibin), self.fNormalizedGraph.GetErrorX(ibin), syst_erry, syst_erry)
                syst_xsec_den_tot += self.fNormalizationGraph.GetErrorY(ibin)  # take the weighted average of the rel unc
            else:
                syst_erry2 = 0
                syst_erry = 0
            xsec_incl_tot += self.fNormalizationHistogram.GetBinContent(ibin + 1)
            stat_xsec_den_tot2 += self.fNormalizationHistogram.GetBinError(ibin + 1) ** 2

            self.fNormalizedHistogram.SetBinContent(ibin + 1, self.fHistogram.GetBinContent(ibin + 1) / self.fNormalizationHistogram.GetBinContent(ibin + 1))
            stat_err_y2 = ((self.fHistogram.GetBinError(ibin + 1) / self.fHistogram.GetBinContent(ibin + 1)) ** 2 + (self.fNormalizationHistogram.GetBinError(ibin + 1) / self.fNormalizationHistogram.GetBinContent(ibin + 1)) ** 2) * self.fNormalizedHistogram.GetBinContent(ibin + 1) ** 2
            stat_err_y = math.sqrt(stat_err_y2)
            self.fNormalizedHistogram.SetBinError(ibin + 1, stat_err_y)

            tot_err_y2 = stat_err_y2 + syst_erry2
            tot_err_y = math.sqrt(tot_err_y2)

            xsec_d0_tot += self.fHistogram.GetBinContent(ibin + 1)
            stat_xsec_num_tot2 += self.fHistogram.GetBinError(ibin + 1) ** 2
            

            if self.fNormalizedHistogram.GetXaxis().GetBinCenter(ibin+1) > 8:
                avg_flat += self.fNormalizedHistogram.GetBinContent(ibin + 1) / tot_err_y
                sum_of_w_avg_flat += 1.0 / tot_err_y
                stat_avg_flat2 += stat_err_y2 / tot_err_y2
                syst_avg_flat2 += syst_erry2 / tot_err_y2


            print("Bin {}, x = {}, ratio = {} +/- {} (stat) +/- {} (syst)".format(ibin, self.fNormalizedHistogram.GetXaxis().GetBinCenter(ibin+1), self.fNormalizedHistogram.GetBinContent(ibin + 1), stat_err_y, syst_erry))

        stat_xsec_incl_tot = math.sqrt(stat_xsec_den_tot2)
        stat_xsec_d0_tot = math.sqrt(stat_xsec_num_tot2)

        avg_tot = xsec_d0_tot / xsec_incl_tot
        stat_avg_tot = math.sqrt(stat_xsec_den_tot2 / xsec_incl_tot ** 2 + stat_xsec_num_tot2 / xsec_d0_tot ** 2) * avg_tot

        syst_xsec_incl_tot2 = syst_xsec_den_tot ** 2
        syst_xsec_num_tot2 = syst_xsec_num_tot ** 2

        syst_avg_tot = math.sqrt(syst_xsec_incl_tot2 / xsec_incl_tot ** 2 + syst_xsec_num_tot2 / xsec_d0_tot ** 2) * avg_tot

        avg_flat /= sum_of_w_avg_flat
        stat_avg_flat = math.sqrt(stat_avg_flat2) / sum_of_w_avg_flat
        syst_avg_flat = math.sqrt(syst_avg_flat2) / sum_of_w_avg_flat

        print("Average ratio in full range: {} +/- {} (stat) +/- {} (syst)".format(avg_tot, stat_avg_tot, syst_avg_tot))
        print("Average ratio in flat region (pt > 8 GeV/c): {} +/- {} (stat) +/- {} (syst)".format(avg_flat, stat_avg_flat, syst_avg_flat))

    
    def GetRateNormalization(self):
        if not self.fNormalizationHistogram:
            print("For normalization type '{}' the field 'fNormalizationHistogram' must be set!".format(self.fNormalizationType))
            exit(1)
        
        if not self.fTotalCrossSection:
            xsec_tot = 0
            stat_xsec_tot2 = 0
            syst_xsec_tot = 0
            for ibin in range(1, self.fNormalizationHistogram.GetNbinsX() + 1):
                if self.fNormalizationXmin < self.fNormalizationXmax:
                    if self.fNormalizationHistogram.GetXaxis().GetBinCenter(ibin) < self.fNormalizationXmin:
                        continue
                    if self.fNormalizationHistogram.GetXaxis().GetBinCenter(ibin) > self.fNormalizationXmax:
                        break
                binw = self.fNormalizationHistogram.GetXaxis().GetBinWidth(ibin)
                xsec = self.fNormalizationHistogram.GetBinContent(ibin)
                xsec_tot += xsec * binw
                stat_xsec_tot2 += (self.fNormalizationHistogram.GetBinError(ibin) * binw) ** 2
                if self.fNormalizationGraph:
                    syst_xsec_tot += self.fNormalizationGraph.GetErrorY(ibin - 1) * binw  # take the weighted average of the rel unc
                print("Cross section in bin [{}, {}] is {} +/- {}".format(self.fNormalizationHistogram.GetXaxis().GetBinLowEdge(ibin), self.fNormalizationHistogram.GetXaxis().GetBinUpEdge(ibin), xsec * binw, self.fNormalizationHistogram.GetBinError(ibin) * binw))

            stat_xsec_tot = math.sqrt(stat_xsec_tot2)

            print("The total cross section for '{}' is {:.4f} +/- {:.4f} (stat) +/- {:.4f} (syst)".format(self.fNormalizationHistogram.GetName(), xsec_tot, stat_xsec_tot, syst_xsec_tot))
            self.fTotalCrossSection = xsec_tot
            self.fTotalCrossSectionStatUnc = stat_xsec_tot
            self.fTotalCrossSectionSystUnc = syst_xsec_tot

        return 1.0 / self.fTotalCrossSection, self.fTotalCrossSectionStatUnc / self.fTotalCrossSection, self.fTotalCrossSectionSystUnc  / self.fTotalCrossSection

    def NormalizeHistogram(self, graph=None):
        print("Performing normalization '{}' for histogram '{}'".format(self.fNormalizationType, self.fHistogram.GetName()))
        self.fGraph = graph
        if self.fGraph:
            if self.fGraph.GetN() != self.fHistogram.GetNbinsX():
                print("Error: graph '{}' has {} points, while histogram '{}' has {} bins.".format(self.fGraph.GetName(), self.fGraph.GetN(), self.fHistogram.GetName(), self.fHistogram.GetNbinsX()))
                exit(1)
        if self.fNormalizationType == NormalizationType.CrossSection:
             self.NormalizeHistogramCrossSection()
        elif self.fNormalizationType == NormalizationType.Distribution:
            self.NormalizeHistogramDistribution()
        elif self.fNormalizationType == NormalizationType.Rate:
            self.GetCrossSectionNormalization = self.GetRateNormalization
            self.NormalizeHistogramCrossSection()
        elif self.fNormalizationType == NormalizationType.Ratio:
            self.NormalizeHistogramRatio()
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
        return CROSS_SECTION / (self.fEvents * BRANCHING_RATIO), 0, 0

class MCSimulationNormalizator(Normalizator):
    def __init__(self, h, normalization):
        Normalizator.__init__(self, h, normalization)
        self.fNormalizationOption = "width"
        self.GetCrossSectionNormalization = self._GetCrossSectionNormalization
