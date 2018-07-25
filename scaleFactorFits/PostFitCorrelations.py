import ROOT
from ROOT import TColor
import math
from copy import deepcopy
from array import array

def PlotPostFitCorrelations(model, dictIn, outputName):
  
    par_names = []

    for o in model.get_observables():
        for p in model.get_processes(o):
            parameters = model.get_histogram_function(o, p).parameters
            for par in parameters:
                if par not in par_names:
                    par_names.append(par)
        break

    for o in model.get_observables():
        for p in model.get_processes(o):
            coeffs = model.get_coeff(o,p)
            for c in coeffs.factors:
                if c not in par_names:
                    par_names.append(c)
    
    best_fit = 0
    if '__chi2' in dictIn:
        best_chi2 =9999999999999 
        for i, chi2_val in enumerate(dictIn['__chi2']):
            if chi2_val < best_chi2:
                best_chi2 = chi2_val
                best_fit = i
 
    Nbins = len(par_names);
    h_corr = ROOT.TH2D("","", Nbins, 0, Nbins, Nbins, 0, Nbins)

    cov = deepcopy(dictIn['__cov'][best_fit])

    lables = []
    for i in range(0,len(par_names)):
        for j in range(0,len(par_names)):
            if math.fabs(cov[i][i] - dictIn[par_names[j]][best_fit][1]**2) < 10**(-13):
                lables.append(par_names[j])

    for i in range(0,Nbins):
        for j in range(0,Nbins):
            h_corr.SetBinContent(i+1,j+1, cov[i][j]/math.sqrt(cov[i][i])/math.sqrt(cov[j][j]))

    #histogram settings
    Red = array('d', [0.00, 1.00, 0.90])
    Green = array('d', [0.90, 1.00, 0.20])
    Blue = array('d', [ 0.00, 1.00, 0.20])
    Length = array('d', [0.00, 0.50, 1.00])

    ROOT.TColor.CreateGradientColorTable(3,Length,Red,Green,Blue,50);

    ROOT.gStyle.SetPaintTextFormat("4.2f")
    ROOT.gStyle.SetOptStat(0)
    #h_corr.SetPaintTextFormat("4.2f");
    xAxis = h_corr.GetXaxis()
    yAxis = h_corr.GetYaxis()
    for i in range(0,Nbins):
        xAxis.SetBinLabel(i+1,lables[i])
        yAxis.SetBinLabel(i+1,lables[i])
    h_corr.SetMaximum(1.)
    h_corr.SetMinimum(-1.)
    h_corr.SetMarkerSize(0.6)
    h_corr.GetXaxis().SetLabelSize(0.015)
    h_corr.GetYaxis().SetLabelSize(0.02)
    can = ROOT.TCanvas()
    xAxis.SetTicks("")
    yAxis.SetTicks("")
    can.SetTickx(0)
    can.SetTicky(0)
    #can.SetRightMargin(0.2)
    can.SetLeftMargin(0.2)

    
    h_corr.Draw("COLZ")                          
    h_corr.Draw("TEXT SAME") 
   
    can.SaveAs(outputName+".eps")
