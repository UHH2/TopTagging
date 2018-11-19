import ROOT
import numpy
from copy import deepcopy
from array import array

from efficiencyClasses import *
from distribution import *


    
class SFcalculation(object):
    def __init__(self,fileName_stat, fileName_sys, dictIn_stat, dictIn_sys, model_stat, model_sys): 
    
        self.Infile_stat = ROOT.TFile(fileName_stat, "READ")
        fileName_stat_no = deepcopy(fileName_stat)
        fileName_stat_no = fileName_stat_no.replace(".root", "_NotScaled.root")
        self.InfileNotScaled_stat = ROOT.TFile(fileName_stat_no, "READ")

        self.Infile_sys = ROOT.TFile(fileName_sys, "READ")
        fileName_sys_no = deepcopy(fileName_sys)
        fileName_sys_no = fileName_sys_no.replace(".root", "_NotScaled.root")
        self.InfileNotScaled_sys = ROOT.TFile(fileName_sys_no, "READ")
 
        self.dictIn_stat = dictIn_stat
        self.dictIn_sys = dictIn_sys
        self.model_stat = model_stat
        self.model_sys = model_sys

        self.extra_mass_files = False
        self.extra_mass_files_sys = False

    def __del__(self):
        self.Infile_stat.Close()
        self.InfileNotScaled_stat.Close()

        self.Infile_sys.Close()
        self.InfileNotScaled_sys.Close()

        if(self.extra_mass_files):
            self.Infile_Mass_stat.Close()
        if(self.extra_mass_files_sys):
            self.Infile_Mass_sys.Close()
            

    def setMassEffHists(self,fileName):
        self.Infile_Mass_stat = ROOT.TFile(fileName, "READ")
        self.extra_mass_files = True

    def setMassEffHistsSys(self,fileName):
        self.Infile_Mass_sys = ROOT.TFile(fileName, "READ")
        self.extra_mass_files_sys = True    
        
    def setMassWindow(self, mass_min, mass_max):
        self.mass_min = mass_min
        self.mass_max = mass_max

    def getDistribution(self,observable, process, stat_sys = 'stat', name_postfix='', fine = False):

        fileIn = self.Infile_stat
        model = self.model_stat
        dictIn = self.dictIn_stat
        if stat_sys == 'sys':
            fileIn = self.Infile_sys
            model = self.model_sys
            dictIn = self.dictIn_sys 
        
        if fine:
            if stat_sys == 'sys':
                fileIn = self.Infile_Mass_sys   
            else:
                fileIn = self.Infile_Mass_stat

        #set the input histograms 
        hist = ROOT.TH1F()
        fileIn.GetObject(observable+'__'+process, hist)

        distr = distribution(name_postfix)
        distr.set(hist)  

        best_fit = 0

        if '__chi2' in dictIn:
            best_chi2 =9999999999999 
            for i, chi2_val in enumerate(dictIn['__chi2']):
                if chi2_val < best_chi2:
                    best_chi2 = chi2_val
                    best_fit = i

 
        #apply shape morphing and uncertainties
        parameters = model.get_histogram_function(observable, process).parameters
        for par in parameters:
            shape = dictIn[par][best_fit][0]
            shape_err = dictIn[par][best_fit][1]

            h_up = ROOT.TH1F()
            h_down = ROOT.TH1F()
            fileIn.GetObject(observable+'__'+process+'__'+par+'__plus',h_up)
            fileIn.GetObject(observable+'__'+process+'__'+par+'__minus',h_down)
                
            distr.addShapeVar(par, h_up, h_down, shape, shape_err)

              
        #apply rate scaling and uncertainties
        coeffs = model.get_coeff(observable,process)
        for c in coeffs.factors:
            scn = dictIn[c][best_fit][0]
            scn_err = dictIn[c][best_fit][1]
            lambda_scale = coeffs.factors[c]['lambda_plus']
        
            scale = numpy.exp(scn*lambda_scale)
            scale_err = scale*lambda_scale*scn_err

            distr.scale(c,scale,scale_err)

        return distr


    def calcMassEff(self, observable, pre_post, stat_sys):

        fileIn = self.Infile_stat
        model = self.model_stat
        if self.extra_mass_files: 
            fileIn = self.Infile_Mass_stat
        if self.extra_mass_files_sys and stat_sys == 'sys': 
            fileIn = self.Infile_Mass_sys
            model = self.model_sys
     
        if pre_post is 'post':
        #get the data histogram
            h_data = ROOT.TH1F()
            fileIn.GetObject(observable+'_pass__DATA', h_data)
            distr_DATA = distribution()
            distr_DATA.set(h_data)

            eff_mass_post = efficiencyMass(distr_DATA, self.mass_min, self.mass_max)
            return eff_mass_post


        if pre_post is 'pre':
            distr_MC = distribution()

            MCset = False
     
            for p in model.get_processes(observable+'_pass'):
                print p
                distr_pass = self.getDistribution(observable+'_pass' ,p, stat_sys, '_mass',  self.extra_mass_files)

                if not MCset:
                    distr_MC = deepcopy(distr_pass)
                    MCset = True
                else:
                    distr_MC.add(distr_pass.name, distr_pass)

            eff_mass_pre = efficiencyMass(distr_MC, self.mass_min, self.mass_max)
            return eff_mass_pre

    
    def calcEfficiencies(self, stat_sys):
        #output dictionary
        d_out = dict()

        model = self.model_stat
        InfileNotScaled  = self.InfileNotScaled_stat
        if stat_sys == 'sys':
            model = self.model_sys
            InfileNotScaled  = self.InfileNotScaled_sys

        observables = []
        for o in model.get_observables():
            if 'pass' in o:
                observables.append(o)


        for o in observables:
            d_obs = dict()

            o_mass = o.replace('_pass','')
            eff_mass_pre = self.calcMassEff(o_mass, 'pre', stat_sys)
            eff_mass_post = self.calcMassEff(o_mass, 'post', stat_sys)

            for p in model.get_processes(o):

                #read histograms
                hist_fail = ROOT.TH1F()
                hist_pass = ROOT.TH1F()

                o = o.replace('_pass','')
 
                InfileNotScaled.GetObject(o+'_fail__'+p, hist_fail)
                InfileNotScaled.GetObject(o+'_pass__'+p, hist_pass)

                #define pre fit and post fit distributions
                d_pass = distribution("_pre")
                d_fail = distribution("_pre")
    
                d_pass.set(hist_pass)
                d_fail.set(hist_fail) 

                d_pass_post = self.getDistribution(o+'_pass',p, stat_sys)
                d_fail_post = self.getDistribution(o+'_fail',p, stat_sys)

                #calculate efficiencies
                eff_pre = efficiency(d_pass, d_fail) 
                eff_post = efficiency(d_pass_post, d_fail_post)
   
 #               eff_mass_pre = self.calcMassEff(o, 'pre', stat_sys)
 #               eff_mass_post = self.calcMassEff(o, 'post', stat_sys)
 
                eff_tot_pre = deepcopy(eff_pre)
                eff_tot_pre.add(eff_mass_pre)

                eff_tot_post = deepcopy(eff_post)
                eff_tot_post.add(eff_mass_post)
 
                #calculate scale facotrs
                sf = scaleFactor(eff_pre, eff_post)
                sf_tot = scaleFactor(eff_tot_pre, eff_tot_post)
               
                d_eff = {p: { 'effPreFit_NoMass':[eff_pre.eff, eff_pre.unc_tot] ,  'effPostFit_NoMass': [eff_post.eff, eff_post.unc_tot],  'effPreFit_Mass':[eff_mass_pre.eff, eff_mass_pre.unc_tot] ,  'effPostFit_Mass': [eff_mass_post.eff, eff_mass_post.unc_tot], 'effPreFit':[eff_tot_pre.eff, eff_tot_pre.unc_tot] ,  'effPostFit': [eff_tot_post.eff, eff_tot_post.unc_tot], 'sf_NoMass':[sf.sf, sf.unc_tot], 'sf' : [sf_tot.sf, sf_tot.unc_tot] }, }

                d_obs.update(d_eff)

            d_out.update( {o: d_obs})
          
        return d_out


def WriteEffGraphs_separate(name, dicts, dicts_sys, bins):

    file_out = ROOT.TFile(name, "RECREATE")

    observables = []
    for i in range(0, len(dicts)):
        for obs in dicts[i]:
            observables.append(obs)
        
    processes = []
    for proc in dicts[0][observables[0]]:
        processes.append(proc)

    variables = []
    for var in dicts[0][observables[0]][ processes[0]]:

        variables.append(var)

    for proc in processes:
        for var in variables:

            title = ''
            if 'eff' in var:
                title = ";p_{T};Efficieny"
            if 'ttSF' in var:
                title = ";p_{T};Scale factor"
          
            eff_stat = array('d', [dicts[bin-1][observables[bin-1]][proc][var][0] for bin in range(1, len(bins))])
            eff = array('d', [dicts_sys[bin-1][observables[bin-1]][proc][var][0] for bin in range(1, len(bins))])        
            e_high = array('d', [dicts_sys[bin-1][observables[bin-1]][proc][var][1] for bin in range(1, len(bins))])
            e_low = array('d', [dicts_sys[bin-1][observables[bin-1]][proc][var][1] for bin in range(1, len(bins))])
            e_stat_high = array('d', [dicts[bin-1][observables[bin-1]][proc][var][1] for bin in range(1, len(bins))])
            e_stat_low = array('d', [dicts[bin-1][observables[bin-1]][proc][var][1] for bin in range(1, len(bins))])
            e_x = array('d', [(bins[bin-1] + 0.5*(bins[bin]-bins[bin-1])) for bin in range(1, len(bins))])
            e_x_low = array('d', [(0.5*(bins[bin]-bins[bin-1])) for  bin in range(1, len(bins))])
            e_x_high = array('d', [(0.5*(bins[bin]-bins[bin-1])) for  bin in range(1, len(bins))])

            graph_stat = ROOT.TGraphAsymmErrors(len(bins)-1, e_x, eff_stat, e_x_low, e_x_high, e_stat_low, e_stat_high)
            graph_sys = ROOT.TGraphAsymmErrors(len(bins)-1, e_x, eff, e_x_low, e_x_high, e_low, e_high)
            
            graph_stat.Write(var+'_'+proc+'_stat')
            graph_sys.Write(var+'_'+proc+'_tot')

    file_out.Close()

