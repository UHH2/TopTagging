import ROOT 

class ThetaPostFitPlot():
    def __init__(self,input_dicts):
        self.input_dicts = input_dicts

    @staticmethod
    def prepare_post_fit_items(post_fit_dict):
	#print post_fit_dict
	#print " "


        best_fit = 0

        if '__chi2' in post_fit_dict:
            best_chi2 =9999999999999 
        for i, chi2_val in enumerate(post_fit_dict['__chi2']):
            if chi2_val < best_chi2:
                best_chi2 = chi2_val
                best_fit = i
        print '+-+-+-+', best_fit

	result = []
	for entry in post_fit_dict:
	  if entry in "__nll": continue
	  if entry in "__chi2": continue
	  if entry in "__cov": continue
	  #print entry, post_fit_dict[entry][0]
          #print '---------',entry, post_fit_dict[entry]
	  result.append((entry,post_fit_dict[entry][best_fit]))
	return result

        mylist =list((name, val_err)
                     for name, (val_err,) in sorted(post_fit_dict.iteritems())
                     if name not in ('__nll'))
        print mylist
        #for k, v in sorted(post_fit_dict.iteritems()):
            
        return list(
            (name, val_err)
            for name, (val_err,) in sorted(post_fit_dict.iteritems())
            if name not in ('__nll')
        )

    @staticmethod
    def prepare_pull_graph(n_items, post_fit_items):
        g = ROOT.TGraphAsymmErrors(n_items)
        for i, (_, (val, err)) in enumerate(post_fit_items):
            x, y = val, i + 1.5
            g.SetPoint(i, x, y)
            g.SetPointEXlow(i, err)
            g.SetPointEXhigh(i, err)

        g.SetLineStyle(1)
        g.SetLineWidth(1)
        g.SetLineColor(1)
        g.SetMarkerStyle(21)
        g.SetMarkerSize(1.25)
        return g

    @staticmethod
    def prepare_band_graphs(n_items):
        g68 = ROOT.TGraph(2*n_items+7)
        g95 = ROOT.TGraph(2*n_items+7)
        for a in xrange(0, n_items+3):
            g68.SetPoint(a, -1, a)
            g95.SetPoint(a, -2, a)
            g68.SetPoint(a+1+n_items+2, 1, n_items+2-a)
            g95.SetPoint(a+1+n_items+2, 2, n_items+2-a)
        g68.SetFillColor(ROOT.kGreen)
        g95.SetFillColor(ROOT.kYellow)
        return g68, g95

    @staticmethod
    def prepare_canvas(name):
        c_name = 'cnv_post_fit_' + name
        c = ROOT.TCanvas(c_name, c_name, 600, 700)
        c.SetTopMargin(0.06)
        c.SetRightMargin(0.02)
        c.SetBottomMargin(0.12)
        c.SetLeftMargin(0.35*700/650)
        c.SetTickx()
        c.SetTicky()
        return c

    @staticmethod
    def put_axis_foo(n_items, prim_graph, post_fit_items):
        prim_hist = prim_graph.GetHistogram()
        ax_1 = prim_hist.GetYaxis()
        ax_2 = prim_hist.GetXaxis()

        prim_graph.SetTitle('')
        ax_2.SetTitle('post-fit nuisance parameters values')
        #ax_2.SetTitle('deviation in units of #sigma')
        ax_1.SetTitleSize(0.050)
        ax_2.SetTitleSize(0.050)
        ax_1.SetTitleOffset(1.4)
        ax_2.SetTitleOffset(1.0)
        ax_1.SetLabelSize(0.03)
        #ax_2.SetLabelSize(0.05)
        ax_1.SetRangeUser(0, n_items+2)
        ax_2.SetRangeUser(-2.4, 2.4)

        ax_1.Set(n_items+2, 0, n_items+2)
        ax_1.SetNdivisions(-414) 
        for i, (uncert_name, _) in enumerate(post_fit_items):
            ax_1.SetBinLabel(i+2, uncert_name)

    def mk_canvas(self, sig_name, post_fit_dict,saveDir, name):
        n = len(post_fit_dict)-1
        #print post_fit_dict
        items = self.prepare_post_fit_items(post_fit_dict)

        g = self.prepare_pull_graph(n, items)
        g68, g95 = self.prepare_band_graphs(n)
        cnv = self.prepare_canvas(sig_name)
        
        cnv.Print(saveDir+'/'+sig_name+"_"+name+"_PostFitNuissancePar.pdf[")

        cnv.cd()
        g95.Draw('AF')
        g68.Draw('F')
        g.Draw('P')
        
        self.put_axis_foo(n, g95, items)
        g95.GetHistogram().Draw('axis,same')
        cnv.Modified()
        cnv.Update()

        cnv.Print(saveDir+'/'+sig_name+"_"+name+"_PostFitNuissancePar.pdf")
        cnv.Print(saveDir+'/'+sig_name+"_"+name+"_PostFitNuissancePar.pdf]")

        return cnv

    def make_plots(self,saveDir,name):
        #print self.input_dicts
        for k,v in self.input_dicts.iteritems():
            canvas = self.mk_canvas(k,self.input_dicts[k],saveDir,name)
            
        
        
