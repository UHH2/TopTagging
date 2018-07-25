# -*- coding: utf-8 -*-

import numpy
import ROOT
from array import array

import os
import sys

sys.path.append('./')
from PostFitUncertainties import *
from PostFitCorrelations import *
from SFcalculation import *
from ThetaPostFitPlot import *


def run(fname_stat, fname_sys, outname_stat, outname_sys, calc_sfs = True, write_report = False):
    if write_report:report.reopen_file()

    inputpath = "thetaFilesNewSYS/mass_sub/"
    print 'use files in', inputpath

    model_stat = build_model_from_rootfile(inputpath+fname_stat)
    model_stat.fill_histogram_zerobins()

    model_sys = build_model_from_rootfile(inputpath+fname_sys)
    model_sys.fill_histogram_zerobins()

    print '================================================='
    print fname_sys
    print '================================================='

    rate_unc_high = math.log(2.0)
    rate_unc_low = math.log(2.0)

    model_stat.add_lognormal_uncertainty('TTbar_mergedTop_Pass_rate', rate_unc_low , procname='TTbar_mergedTop',obsname='Mass_pass')
    model_stat.add_lognormal_uncertainty('TTbar_semimerged_Pass_rate',  rate_unc_high, procname='TTbar_semimerged',obsname='Mass_pass')
    model_stat.add_lognormal_uncertainty('TTbar_notmerged_Pass_rate', rate_unc_high , procname='TTbar_notmerged',obsname='Mass_pass')

    model_sys.add_lognormal_uncertainty('TTbar_mergedTop_Pass_rate', rate_unc_low , procname='TTbar_mergedTop',obsname='Mass_pass')
    model_sys.add_lognormal_uncertainty('TTbar_semimerged_Pass_rate',  rate_unc_high, procname='TTbar_semimerged',obsname='Mass_pass')
    model_sys.add_lognormal_uncertainty('TTbar_notmerged_Pass_rate', rate_unc_high , procname='TTbar_notmerged',obsname='Mass_pass')


    model_stat.add_lognormal_uncertainty('TTbar_mergedTop_Fail_rate', rate_unc_high, procname='TTbar_mergedTop',obsname='Mass_fail')
    model_stat.add_lognormal_uncertainty('TTbar_semimerged_Fail_rate', rate_unc_low, procname='TTbar_semimerged',obsname='Mass_fail')
    model_stat.add_lognormal_uncertainty('TTbar_notmerged_Fail_rate', rate_unc_low, procname='TTbar_notmerged',obsname='Mass_fail') 

    model_sys.add_lognormal_uncertainty('TTbar_mergedTop_Fail_rate', rate_unc_high, procname='TTbar_mergedTop',obsname='Mass_fail')
    model_sys.add_lognormal_uncertainty('TTbar_semimerged_Fail_rate', rate_unc_low, procname='TTbar_semimerged',obsname='Mass_fail')
    model_sys.add_lognormal_uncertainty('TTbar_notmerged_Fail_rate', rate_unc_low, procname='TTbar_notmerged',obsname='Mass_fail') 


    options = Options()
    options.set('minimizer', 'strategy', 'robust')

    n_name_stat = deepcopy(outname_stat)
    n_name_stat = n_name_stat.replace(".root","_")
    n_name_sys = deepcopy(outname_sys)
    n_name_sys = n_name_sys.replace(".root","_")
 
    print '============'
    print 'stat'
    print '============'
    mle_output_stat = mle(model_stat, input='data', n=100, with_covariance=True, chi2=True, options=options, signal_process_groups ={'background_only':[]})
    print '============'
    print 'sys'
    print '============'
    mle_output_sys = mle(model_sys, input='data', n=100, with_covariance=True, chi2=True, options=options, signal_process_groups ={'background_only':[]})

    PlotPostFitCorrelations(model_stat, mle_output_stat['background_only'], "fitResultsNewSYS/nuissance/Corr_"+n_name_stat)
    PlotPostFitCorrelations(model_sys, mle_output_sys['background_only'], "fitResultsNewSYS/nuissance/Corr_"+n_name_sys)
   
    writeOutputFile(inputpath+fname_stat,  "fitResultsNewSYS/mass_sub/"+outname_stat, mle_output_stat['background_only'], model_stat)
    writeOutputFile(inputpath+fname_sys,  "fitResultsNewSYS/mass_sub/"+outname_sys, mle_output_sys['background_only'], model_sys)

    #writeOutputFile("thetaFilesNew+fname,  "fitResultsNewSYS/pt/"+outname, mle_output['background_only'], model, False)
    #writeOutputFile("thetaFilesNewSYS/tau32/"+fname,  "fitResultsNewSYS/tau32/"+outname, mle_output['background_only'], model, False)
    
    if calc_sfs:
        sfs = SFcalculation(inputpath+fname_stat, inputpath+fname_sys, mle_output_stat['background_only'], mle_output_sys['background_only'], model_stat, model_sys)

        if 'CHS' in fname_stat:
            sfs.setMassWindow(105,220)
        elif 'PUPPI' in fname_stat:
            sfs.setMassWindow(105,210)
        else:
            sfs.setMassWindow(140,220)
        if not 'HOTVR' in fname_stat:
            print 'set extra mass files'
            sfs.setMassEffHists("thetaFilesNewSYS/mass_sub/fine/"+fname_stat);
            sfs.setMassEffHistsSys("thetaFilesNewSYS/mass_sub/fine/"+fname_sys);

        dictOut_stat = sfs.calcEfficiencies('stat')
        dictOut_sys = sfs.calcEfficiencies('sys')
    
#        for pf_vals in mle_output_stat.itervalues():
#            del pf_vals['__nll']
#            del pf_vals['__cov']
#        for pf_vals in mle_output_sys.itervalues():
#            del pf_vals['__nll']
#            del pf_vals['__cov']
     
#        postfit_stat = ThetaPostFitPlot(mle_output_stat)
#        postfit_stat.make_plots("fitResultsNewSYS/nuissance/",n_name_stat)

#       postfit_sys = ThetaPostFitPlot(mle_output_sys)
#       postfit_sys.make_plots("fitResultsNewSYS/nuissance/",n_name_sys)

        return [dictOut_stat, dictOut_sys]


###loop over taggers and wps


bins = array('d', [300, 400, 480, 600, 1100])
calculate_scaleFactors = True

#run("thetaFile_400_PUPPI_sys_.root", "Hists_400_PUPPI_sys.root")

wps = ["_wp1", "_wp2", "_wp3", "_wp4", "_wp5", "_wp1_btag", "_wp2_btag", "_wp3_btag", "_wp4_btag", "_wp5_btag"]

for wp in wps: 

    d_300to400 = run("thetaFile_300to400_PUPPI_stat"+wp+".root", "thetaFile_300to400_PUPPI_sys"+wp+".root", "Hists_300to400_PUPPI_stat"+wp+".root", "Hists_300to400_PUPPI_sys"+wp+".root", calculate_scaleFactors)
    d_400to480 = run("thetaFile_400to480_PUPPI_stat"+wp+".root", "thetaFile_400to480_PUPPI_sys"+wp+".root", "Hists_400to480_PUPPI_stat"+wp+".root", "Hists_400to480_PUPPI_sys"+wp+".root", calculate_scaleFactors)
    d_480to600 = run("thetaFile_480to600_PUPPI_stat"+wp+".root", "thetaFile_480to600_PUPPI_sys"+wp+".root", "Hists_480to600_PUPPI_stat"+wp+".root", "Hists_480to600_PUPPI_sys"+wp+".root", calculate_scaleFactors)
    d_600 = run("thetaFile_600_PUPPI_stat"+wp+".root", "thetaFile_600_PUPPI_sys"+wp+".root", "Hists_600_PUPPI_stat"+wp+".root", "Hists_600_PUPPI_sys"+wp+".root", calculate_scaleFactors)

    if calculate_scaleFactors:
        dicts_stat = [d_300to400[0], d_400to480[0], d_480to600[0], d_600[0]]
        dicts_sys = [d_300to400[1], d_400to480[1], d_480to600[1], d_600[1]]
        WriteEffGraphs_separate("ScaleFactors_NoIso/eff_hists_PUPPI"+wp+".root", dicts_stat, dicts_sys, bins)

d_300to400_HOTVR = run("thetaFile_300to400_HOTVR_stat__modelUnc.root", "thetaFile_300to400_HOTVR_sys__modelUnc.root", "Hists_300to400_HOTVR_stat__modelUnc.root", "Hists_300to400_HOTVR_sys__modelUnc.root", calculate_scaleFactors)
d_400to480_HOTVR = run("thetaFile_400to480_HOTVR_stat__modelUnc.root", "thetaFile_400to480_HOTVR_sys__modelUnc.root", "Hists_400to480_HOTVR_stat__modelUnc.root", "Hists_400to480_HOTVR_sys__modelUnc.root", calculate_scaleFactors)
d_480to600_HOTVR = run("thetaFile_480to600_HOTVR_stat__modelUnc.root", "thetaFile_480to600_HOTVR_sys__modelUnc.root", "Hists_480to600_HOTVR_stat__modelUnc.root", "Hists_480to600_HOTVR_sys__modelUnc.root", calculate_scaleFactors)
d_600_HOTVR = run("thetaFile_600_HOTVR_stat__modelUnc.root", "thetaFile_600_HOTVR_sys__modelUnc.root", "Hists_600_HOTVR_stat__modelUnc.root", "Hists_600_HOTVR_sys__modelUnc.root", calculate_scaleFactors)

if calculate_scaleFactors:
    dicts_stat = [d_300to400_HOTVR[0], d_400to480_HOTVR[0], d_480to600_HOTVR[0], d_600_HOTVR[0]]
    dicts_sys = [d_300to400_HOTVR[1], d_400to480_HOTVR[1], d_480to600_HOTVR[1], d_600_HOTVR[1]]
    WriteEffGraphs_separate("ScaleFactors_NoIso/eff_hists_HOTVR.root", dicts_stat, dicts_sys, bins)


