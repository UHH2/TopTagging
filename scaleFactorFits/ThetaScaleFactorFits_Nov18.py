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

    inputpath = "thetaFiles_FSR/mass_sub/"
    print 'use files in', inputpath

    model_stat = build_model_from_rootfile(inputpath+fname_stat, include_mc_uncertainties=True)
    model_stat.fill_histogram_zerobins()

    model_sys = build_model_from_rootfile(inputpath+fname_sys,  histogram_filter=(lambda s:( not s.count('Btag') and not s.count('MounID') and not s.count('Trigger')) ) , include_mc_uncertainties=True)
    model_sys.fill_histogram_zerobins()

    print '================================================='
    print fname_sys
    print '================================================='

    rate_unc_TTbar = math.log(2.0)
    rate_unc_TTbar_notmerged_fail = math.log(2.0)
    rate_unc_TTbar_semimerged_pass = math.log(2.0)

    if '400to480' in fname_sys and 'PUPPI' in fname_sys and 'wp2' in fname_sys and not '_btag' in fname_sys:
        rate_unc_TTbar = math.log(1.8)
        rate_unc_TTbar_notmerged_fail = math.log(1.8)
        rate_unc_TTbar_semimerged_pass = math.log(1.8)

    if '480to600' in fname_sys and 'PUPPI' in fname_sys:
        rate_unc_TTbar = math.log(2.0)
        rate_unc_TTbar_notmerged_fail = math.log(2.0)
        rate_unc_TTbar_semimerged_pass = math.log(3.0) 

    if '_600' in fname_sys and 'CHS' in fname_sys:
        rate_unc_TTbar = math.log(2.0)
        rate_unc_TTbar_notmerged_fail = math.log(2.0)
        rate_unc_TTbar_semimerged_pass = math.log(3.0)

    if '300to400' in fname_sys and 'PUPPI' in fname_sys:
        rate_unc_TTbar = math.log(1.9)
        rate_unc_TTbar_notmerged_fail = math.log(1.9)
        rate_unc_TTbar_semimerged_pass = math.log(1.9)
 
        if 'wp5' in fname_sys and not '_btag' in fname_sys:
            rate_unc_TTbar = math.log(1.7)
            rate_unc_TTbar_notmerged_fail = math.log(1.7)
            rate_unc_TTbar_semimerged_pass = math.log(1.7)
            
        if 'wp2_btag' in fname_sys:
            rate_unc_TTbar = math.log(1.8)
            rate_unc_TTbar_notmerged_fail = math.log(1.8)
            rate_unc_TTbar_semimerged_pass = math.log(1.8)
            
        if 'wp3_btag' in fname_sys:
            rate_unc_TTbar = math.log(1.9)
            rate_unc_TTbar_notmerged_fail = math.log(1.15)
            rate_unc_TTbar_semimerged_pass = math.log(1.9)

        if 'wp4_btag' in fname_sys:
            rate_unc_TTbar = math.log(1.8)
            rate_unc_TTbar_notmerged_fail = math.log(1.8)
            rate_unc_TTbar_semimerged_pass = math.log(1.8)

        if 'wp5_btag' in fname_sys:
            rate_unc_TTbar = math.log(1.6)
            rate_unc_TTbar_notmerged_fail = math.log(1.6)
            rate_unc_TTbar_semimerged_pass = math.log(1.6)

    if '300to400' in fname_sys and 'CHS' in fname_sys:
        rate_unc_TTbar = math.log(1.9)
        rate_unc_TTbar_notmerged_fail = math.log(1.9)
        rate_unc_TTbar_semimerged_pass = math.log(1.9)
#        rate_unc_TTbar_notmerged_pass = math.log(1.9)

        if 'wp2_btag' in fname_sys:
            rate_unc_TTbar = math.log(1.7)
            rate_unc_TTbar_notmerged_fail = math.log(1.7)
            rate_unc_TTbar_semimerged_pass = math.log(1.7)
 
        if 'wp3_btag' in fname_sys:
            rate_unc_TTbar = math.log(1.8)
            rate_unc_TTbar_notmerged_fail = math.log(1.8)
            rate_unc_TTbar_semimerged_pass = math.log(1.8)


    print rate_unc_TTbar

    model_stat.add_lognormal_uncertainty('TTbar_mergedTop_Pass_rate', rate_unc_TTbar, procname='TTbar_mergedTop',obsname='Mass_pass')
    model_stat.add_lognormal_uncertainty('TTbar_semimerged_Pass_rate', rate_unc_TTbar, procname='TTbar_semimerged',obsname='Mass_pass')
    model_stat.add_lognormal_uncertainty('TTbar_notmerged_Pass_rate', rate_unc_TTbar_semimerged_pass, procname='TTbar_notmerged',obsname='Mass_pass')
 
    model_sys.add_lognormal_uncertainty('TTbar_mergedTop_Pass_rate', rate_unc_TTbar , procname='TTbar_mergedTop',obsname='Mass_pass')
    model_sys.add_lognormal_uncertainty('TTbar_semimerged_Pass_rate', rate_unc_TTbar, procname='TTbar_semimerged',obsname='Mass_pass')
    model_sys.add_lognormal_uncertainty('TTbar_notmerged_Pass_rate', rate_unc_TTbar_semimerged_pass, procname='TTbar_notmerged',obsname='Mass_pass')

    model_stat.add_lognormal_uncertainty('TTbar_mergedTop_Fail_rate', rate_unc_TTbar, procname='TTbar_mergedTop',obsname='Mass_fail')
    model_stat.add_lognormal_uncertainty('TTbar_semimerged_Fail_rate', rate_unc_TTbar, procname='TTbar_semimerged',obsname='Mass_fail')
    model_stat.add_lognormal_uncertainty('TTbar_notmerged_Fail_rate', rate_unc_TTbar_notmerged_fail, procname='TTbar_notmerged',obsname='Mass_fail') 

    model_sys.add_lognormal_uncertainty('TTbar_mergedTop_Fail_rate',  rate_unc_TTbar , procname='TTbar_mergedTop',obsname='Mass_fail')
    model_sys.add_lognormal_uncertainty('TTbar_semimerged_Fail_rate', rate_unc_TTbar, procname='TTbar_semimerged',obsname='Mass_fail')
    model_sys.add_lognormal_uncertainty('TTbar_notmerged_Fail_rate', rate_unc_TTbar_notmerged_fail, procname='TTbar_notmerged',obsname='Mass_fail') 

 
    model_sys.add_lognormal_uncertainty('WJets_rate', math.log(1.5) , procname='WJets',obsname='*')
    model_sys.add_lognormal_uncertainty('SingleTop_rate', math.log(1.2) , procname='ST',obsname='*')
    model_sys.add_lognormal_uncertainty('DY_rate', math.log(1.5) , procname='DYJets',obsname='*')
    model_sys.add_lognormal_uncertainty('QCD_rate', math.log(2.0), procname='QCD',obsname='*')

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

    PlotPostFitCorrelations(model_stat, mle_output_stat['background_only'], "fitResults_FSR_f/nuissance/Corr_"+n_name_stat)
    PlotPostFitCorrelations(model_sys, mle_output_sys['background_only'], "fitResults_FSR_f/nuissance/Corr_"+n_name_sys)
   
    writeOutputFile(inputpath+fname_stat,  "fitResults_FSR_f/mass_sub/"+outname_stat, mle_output_stat['background_only'], model_stat)
    writeOutputFile(inputpath+fname_sys,  "fitResults_FSR_f/mass_sub/"+outname_sys, mle_output_sys['background_only'], model_sys)

    #writeOutputFile("thetaFiles_FSRNew+fname,  "fitResults_FSR_fNewSYS/pt/"+outname, mle_output['background_only'], model, False)
    #writeOutputFile("thetaFiles_FSRNewSYS/tau32/"+fname,  "fitResults_FSR_fNewSYS/tau32/"+outname, mle_output['background_only'], model, False)
    
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
            sfs.setMassEffHists("thetaFiles_FSR/mass_sub/fine/"+fname_stat);
            sfs.setMassEffHistsSys("thetaFiles_FSR/mass_sub/fine/"+fname_sys);

        dictOut_stat = sfs.calcEfficiencies('stat')
        dictOut_sys = sfs.calcEfficiencies('sys')
    
    for pf_vals in mle_output_stat.itervalues():
        del pf_vals['__nll']
        del pf_vals['__cov']
    for pf_vals in mle_output_sys.itervalues():
        del pf_vals['__nll']
        del pf_vals['__cov']
     
    postfit_stat = ThetaPostFitPlot(mle_output_stat)
    postfit_stat.make_plots("fitResults_FSR_f/nuissance/",n_name_stat)

    postfit_sys = ThetaPostFitPlot(mle_output_sys)
    postfit_sys.make_plots("fitResults_FSR_f/nuissance/",n_name_sys)

    if calc_sfs:
        return [dictOut_stat, dictOut_sys]


###loop over taggers and wps


bins = array('d', [300, 400, 480, 600, 1100])
calculate_scaleFactors = False

#run("thetaFile_400_PUPPI_sys_.root", "Hists_400_PUPPI_sys.root")

wps = ["_wp1", "_wp2", "_wp3", "_wp4", "_wp5", "_wp1_btag", "_wp2_btag", "_wp3_btag", "_wp4_btag", "_wp5_btag"]
#wps = ["_wp1", "_wp2", "_wp3", "_wp1_btag", "_wp2_btag"]
#wps = [  "_wp3_btag"]

#wps = ["_wp3_btag"]


#wps = []


for wp in wps: 

#    d_300to400 = run("thetaFile_300to400_PUPPI_stat"+wp+".root", "thetaFile_300to400_PUPPI_sys"+wp+".root", "Hists_300to400_PUPPI_stat"+wp+".root", "Hists_300to400_PUPPI_sys"+wp+".root", calculate_scaleFactors)
#    d_400to480 = run("thetaFile_400to480_PUPPI_stat"+wp+".root", "thetaFile_400to480_PUPPI_sys"+wp+".root", "Hists_400to480_PUPPI_stat"+wp+".root", "Hists_400to480_PUPPI_sys"+wp+".root", calculate_scaleFactors)
    d_480to600 = run("thetaFile_480to600_PUPPI_stat"+wp+".root", "thetaFile_480to600_PUPPI_sys"+wp+".root", "Hists_480to600_PUPPI_stat"+wp+".root", "Hists_480to600_PUPPI_sys"+wp+".root", calculate_scaleFactors)
#    d_600 = run("thetaFile_600_PUPPI_stat"+wp+".root", "thetaFile_600_PUPPI_sys"+wp+".root", "Hists_600_PUPPI_stat"+wp+".root", "Hists_600_PUPPI_sys"+wp+".root", calculate_scaleFactors)

    if calculate_scaleFactors:
        dicts_stat = [d_300to400[0], d_400to480[0], d_480to600[0], d_600[0]]
        dicts_sys = [d_300to400[1], d_400to480[1], d_480to600[1], d_600[1]]
        WriteEffGraphs_separate("ScaleFactors_FSR_f/eff_hists_PUPPI"+wp+".root", dicts_stat, dicts_sys, bins)


wpsCHS = ["_wp2", "_wp3", "_wp4", "_wp5", "_wp2_btag", "_wp3_btag", "_wp4_btag", "_wp5_btag"]
#wpsCHS = ["_wp2", "_wp3", "_wp2_btag", "_wp3_btag"]
#wpsCHS = []

for wp in wpsCHS: 

 #   d_300to400_CHS = run("thetaFile_300to400_CHS_stat"+wp+".root", "thetaFile_300to400_CHS_sys"+wp+".root", "Hists_300to400_CHS_stat"+wp+".root", "Hists_300to400_CHS_sys"+wp+".root", calculate_scaleFactors)
 #   d_400to480_CHS = run("thetaFile_400to480_CHS_stat"+wp+".root", "thetaFile_400to480_CHS_sys"+wp+".root", "Hists_400to480_CHS_stat"+wp+".root", "Hists_400to480_CHS_sys"+wp+".root", calculate_scaleFactors)
    d_480to600_CHS = run("thetaFile_480to600_CHS_stat"+wp+".root", "thetaFile_480to600_CHS_sys"+wp+".root", "Hists_480to600_CHS_stat"+wp+".root", "Hists_480to600_CHS_sys"+wp+".root", calculate_scaleFactors)
 #   d_600_CHS = run("thetaFile_600_CHS_stat"+wp+".root", "thetaFile_600_CHS_sys"+wp+".root", "Hists_600_CHS_stat"+wp+".root", "Hists_600_CHS_sys"+wp+".root", calculate_scaleFactors)

    if calculate_scaleFactors:
        dicts_stat_CHS = [d_300to400_CHS[0], d_400to480_CHS[0], d_480to600_CHS[0], d_600_CHS[0]]
        dicts_sys_CHS = [d_300to400_CHS[1], d_400to480_CHS[1], d_480to600_CHS[1], d_600_CHS[1]]
        WriteEffGraphs_separate("ScaleFactors_FSR_f/eff_hists_CHS"+wp+".root", dicts_stat_CHS, dicts_sys_CHS, bins)

