import ROOT
import numpy
from copy import deepcopy
from distribution import *



def writeOutputFile(InputFileName, OutputFileName, InputDict, model, pdf = True):  

    InputFile = ROOT.TFile(InputFileName, "READ")
    OutputFile = ROOT.TFile(OutputFileName, "RECREATE")

    observables = model.get_observables()
    for o in observables:
        d_o = dict()

        #write data hist 
        h_data = ROOT.TH1F()
        InputFile.GetObject(o+'__DATA', h_data)
        h_data.Write()


        best_fit = 0

        if '__chi2' in InputDict:
            best_chi2 =9999999999999 
            for i, chi2_val in enumerate(InputDict['__chi2']):
                if chi2_val < best_chi2:
                    best_chi2 = chi2_val
                    best_fit = i

        for p in model.get_processes(o):

            #set the input histogram 
            h_in = ROOT.TH1F()
            InputFile.GetObject(o+'__'+p, h_in)

            distr = distribution()
            distr.set(h_in)  

            #apply shape morphing and uncertainties
            parameters = model.get_histogram_function(o, p).parameters
            for par in parameters:
                if not pdf and 'PDF' in par:
                    continue
                shape = InputDict[par][best_fit][0]
                shape_err = InputDict[par][best_fit][1]

                h_up = ROOT.TH1F()
                h_down = ROOT.TH1F()
                InputFile.GetObject(o+'__'+p+'__'+par+'__plus',h_up)
                InputFile.GetObject(o+'__'+p+'__'+par+'__minus',h_down)
                
                distr.addShapeVar(par, h_up, h_down, shape, shape_err)

            #apply rate scaling and uncertainties
            coeffs = model.get_coeff(o,p)
            for c in coeffs.factors:
                scn = InputDict[c][best_fit][0]
                scn_err = InputDict[c][best_fit][1]
                lambda_scale = coeffs.factors[c]['lambda_plus']
        
                scale = numpy.exp(scn*lambda_scale)
                scale_err = scale*lambda_scale*scn_err

                distr.scale(c,scale,scale_err) 
        
            #write to file
            distr.writeToFile(OutputFile)    
            


          
   

          

