import ROOT
from copy import deepcopy
from distribution import *


class efficiency(object):
    
    def __init__(self, d_pass, d_fail):
 
        int_min = 1
        int_max = d_pass.distr.GetNbinsX()+1
        overflowBin =  d_pass.distr.GetNbinsX()+1

        Npass = d_pass.distr.Integral(1,overflowBin )
        Nfail = d_fail.distr.Integral(1,overflowBin )

        e_pass_tot2 = d_pass.cov_tot().Integral(1, overflowBin, 1, overflowBin)
        e_fail_tot2 = d_fail.cov_tot().Integral(1, overflowBin, 1, overflowBin)

        #calculate derivatives 
        d_pass = (Nfail)/( (Npass+Nfail)**2)
        d_fail = (-1)*Npass/( (Npass+Nfail)**2)

        self.eff = Npass/(Npass+Nfail)
        self.unc_tot = math.sqrt(d_pass**2 * e_pass_tot2 + d_fail**2 * e_fail_tot2 - 2*d_pass*d_fail*math.sqrt(e_fail_tot2)*math.sqrt(e_pass_tot2))  

         
    def add(self, eff_add):
        self.unc_tot = math.sqrt(eff_add.eff**2 * self.unc_tot**2 + self.eff**2 * eff_add.unc_tot**2)
        self.eff = self.eff * eff_add.eff  

class efficiencyMass(object):

    def __init__(self, d_pass, mass_min = -1, mass_max = -1):

        int_min = 1
        int_max = d_pass.distr.GetNbinsX()+1
        overflowBin =  d_pass.distr.GetNbinsX()+1
        if mass_min >= 0:
            int_min = d_pass.distr.FindBin(mass_min+0.001)
        if mass_max >= 0:
            int_max = d_pass.distr.FindBin(mass_max-0.001)

        #get the Numbers for passed and failed events
        Npass = d_pass.distr.Integral(int_min, int_max)
        Nfail = 0
        if int_min > 1: 
            Nfail = Nfail + d_pass.distr.Integral(1, int_min-1)
        if int_max < d_pass.distr.GetNbinsX()+1:
            Nfail = Nfail + d_pass.distr.Integral(int_max+1, overflowBin )

        Nall  = d_pass.distr.Integral(1,overflowBin )
 
        #get the uncertainties 
        e_pass_tot2 = d_pass.cov_tot_normalized().Integral(int_min, int_max, int_min, int_max)
        e_fail_tot2 = 0
        if int_min > 1: 
            e_fail_tot2 = e_fail_tot2 + d_pass.cov_tot_normalized().Integral(1, int_min-1, 1, int_min-1)
        if int_max < d_pass.distr.GetNbinsX()+1:
            e_fail_tot2 = e_fail_tot2 + d_pass.cov_tot_normalized().Integral(int_max+1, overflowBin, int_max+1, overflowBin)
   
        #calculate derivatives 
        d_pass = (Nfail)/( (Npass+Nfail)**2)
        d_fail = (-1)*Npass/( (Npass+Nfail)**2)

        self.eff = Npass/(Npass+Nfail)
        self.unc_tot = math.sqrt(d_pass**2 * e_pass_tot2 + d_fail**2 * e_fail_tot2 - 2*d_pass*d_fail*math.sqrt(e_fail_tot2)*math.sqrt(e_pass_tot2))  
 
    def add(self, eff_add):
        self.unc_tot = math.sqrt(eff_add.eff**2 * self.unc_tot**2 + self.eff**2 * eff_add.unc_tot**2) 
        self.eff = self.eff * eff_add.eff  


class scaleFactor:
    
    def __init__(self, e_pre, e_post):
        if e_pre.eff > 0:
            self.sf = e_post.eff/e_pre.eff
            self.unc_tot = math.sqrt(1/e_pre.eff**2 * e_post.unc_tot**2 + (e_post.eff/(e_pre.eff**2))**2 * e_pre.unc_tot**2)
        else: 
            self.sf = 1.
            self.unc_tot = 0.

    def add(self, sf_add):
        self.unc_tot = math.sqrt(sf_add.sf**2 * self.unc_tot**2 + self.sf**2 * sf_add.unc_tot**2) 
        self.sf = self.sf * sf_add.sf
