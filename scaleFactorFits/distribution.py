import ROOT
import math
from copy import deepcopy

class distribution(object):

    def __init__(self, _name_postfix = ""):
        self.distr = ROOT.TH1F()
        self.d_cov_sys = dict()
        self.d_shift_sys = dict()
        self.name_postfix = _name_postfix

    def set(self, hist):
        self.name = hist.GetName()
        self.distr = hist.Clone("distr")
        self.distr_nom =  hist.Clone("distr_nominal")
        self.cov_stat = ROOT.TH2D("cov_stat"+hist.GetName()+self.name_postfix,";Bin;Bin",self.distr.GetNbinsX(),0.5,self.distr.GetNbinsX()+0.5,self.distr.GetNbinsX(), 0.5,self.distr.GetNbinsX()+0.5 )
        for bin1 in range(1, self.distr.GetNbinsX()+1):
            for bin2 in range(1, self.distr.GetNbinsX()+1):
                if bin1 is bin2:
                    self.cov_stat.SetBinContent(bin1,bin2,self.distr.GetBinError(bin1)*self.distr.GetBinError(bin2))

    def set_correlations(self,hist):
        self.correlations = hist
  
    def scale(self, name, scale, err):

        h_shift = self.distr.Clone("h_shift")
        h_shift.Reset()
        cov = self.cov_stat.Clone("cov_scale")
        cov.Reset()
        for bin1 in range(1, self.distr.GetNbinsX()+1):
            h_shift.SetBinContent(bin1, self.distr.GetBinContent(bin1)*err)
            for bin2 in range(1, self.distr.GetNbinsX()+1):
                cov.SetBinContent(bin1,bin2,err**2*self.distr.GetBinContent(bin1)*self.distr.GetBinContent(bin2))

        self.distr.Scale(scale)
        self.cov_stat.Scale(scale**2)
        for key in self.d_cov_sys:
            self.d_cov_sys[key].Scale(scale**2)

        for key in self.d_shift_sys:
            self.d_shift_sys[key].Scale(scale)

        d_shift = {name: h_shift,}
        self.d_shift_sys.update(d_shift)
        d_cov = {name: cov,}
        self.d_cov_sys.update(d_cov)

    def addShapeVar(self, name, h_up, h_down, scale, err):

        h_diff = self.distr.Clone("h_diff")
        h_diff.Reset()
        h_shift = self.distr.Clone("h_shift")
        h_shift.Reset()
        cov = self.cov_stat.Clone("cov_shape")
        cov.Reset()

        #calculate diff
        for bin in range(1, h_diff.GetNbinsX()+1):

            up = h_up.GetBinContent(bin)
            down = h_down.GetBinContent(bin)
            nom = self.distr_nom.GetBinContent(bin)
            diff = 0.

      #      if math.fabs(scale) > 1:
      #          if scale > 0:
      #              diff = scale*(up-nom);
      #          else:
      #              diff = scale*(nom-down);                             
      #      else:
            diff = 0.5*scale*(up - down) + (scale**2 - 0.5*math.fabs(scale**3))*(up + down - 2*nom)
                
            h_diff.SetBinContent(bin, diff)
        
        
        #calculate covariance
        for bin1 in range(1, h_diff.GetNbinsX()+1):
            for bin2 in range(1, h_diff.GetNbinsX()+1):

                up1 = h_up.GetBinContent(bin1)
                down1 = h_down.GetBinContent(bin1)
                nom1 = self.distr_nom.GetBinContent(bin1)

                up2 = h_up.GetBinContent(bin2)
                down2 = h_down.GetBinContent(bin2)
                nom2 = self.distr_nom.GetBinContent(bin2)
                d_diff1 = 0.
                d_diff2 = 0.

             #   if math.fabs(scale) > 1:
             #       if scale > 0:
             #           d_diff1 = up1-nom;
             #           d_diff2 = up2-nom;
             #       else:
             #           d_diff1 = nom-down1;
             #           d_diff2 = nom-down2;
             #   else:
                d_diff1 = 0.5*(up1 - down1)+(2*scale-1.5*scale**2)*(up1+down1-2*nom1)
                d_diff2 = 0.5*(up2 - down2)+(2*scale-1.5*scale**2)*(up2+down2-2*nom2)

                vij = d_diff1*d_diff2*err**2
                cov.SetBinContent(bin1,bin2, vij)
                h_shift.SetBinContent(bin1, d_diff1*err)
        
        self.distr.Add(h_diff)
        d_cov = {name: cov,}
        d_shift = {name: h_shift,}
        self.d_cov_sys.update(d_cov)
        self.d_shift_sys.update(d_shift)
        
    def cov_statistical(self):
        cov_stati = self.cov_stat.Clone("cov_stati")
        return cov_stati
    
    def cov_sys(self):
        cov_sys = self.cov_stat.Clone("cov_sys")
        cov_sys.Reset()
        for key in self.d_cov_sys:
            cov_sys.Add(d_cov_sys[key])
        return cov_sys

    def cov_tot(self):
        cov_total = self.cov_stat.Clone("cov_tot")
        for key in self.d_cov_sys:
            cov_total.Add(self.d_cov_sys[key])
        return cov_total

    def cov_tot_normalized(self):
        cov =  self.cov_stat.Clone("cov_tot_norm")
        for key in self.d_shift_sys:
            cov_sys =  self.cov_stat.Clone("cov_sys_norm")
            h_shift_norm = self.d_shift_sys[key].Clone("shift_sys_norm")

            integr_nom = self.distr.Integral(0, self.distr.GetNbinsX()+1)
            intrgr_shift =  self.d_shift_sys[key].Integral(0, self.d_shift_sys[key].GetNbinsX()+1)

            scale = integr_nom /(integr_nom+intrgr_shift)

            for bin1 in range(1, h_shift_norm.GetNbinsX()+1):
                norm_shift = self.d_shift_sys[key].GetBinContent(bin1)*scale + self.distr.GetBinContent(bin1)*(scale-1)
                h_shift_norm.SetBinContent(bin1, norm_shift)

            for bin1 in range(1, h_shift_norm.GetNbinsX()+1):
                for bin2 in range(1, h_shift_norm.GetNbinsX()+1):
                    cov_sys.SetBinContent(bin1,bin2,h_shift_norm.GetBinContent(bin1)*h_shift_norm.GetBinContent(bin2))
            
            cov.Add(cov_sys)
        
        return cov 
                                         
        
    def add(self, name, d_in, scale = 1.0):       
        self.distr.Add(d_in.distr, scale)
        self.cov_stat.Add(d_in.cov_stat, scale**2)
        for key in d_in.d_cov_sys:
            cov_scaled = d_in.d_cov_sys[key].Clone("cov_in")
            cov_scaled.Scale(scale**2)
            d_cov = {key+'_'+name: cov_scaled,}
            self.d_cov_sys.update(d_cov)

        
    def writeToFile(self,OutFile):
        OutFile.cd()
        self.distr.Write(self.name)
         
        for key in self.d_cov_sys:
            
            h_shift = self.distr.Clone("h_shift")
            h_shift.Reset()  
             
            for bin in range(1,h_shift.GetNbinsX()+1):
                h_shift.SetBinContent(bin, math.sqrt(self.d_cov_sys[key].GetBinContent(bin,bin)))
                
            h_up = self.distr.Clone("h_up")
            h_up.Add(h_shift)
            
            h_down = self.distr.Clone("h_down")
            h_down.Add(h_shift,-1)
        
            h_up.Write(self.name+'__'+key+'__plus')
            h_down.Write(self.name+'__'+key+'__minus')
  
