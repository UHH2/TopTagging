
void norm(TH1F* h_sys, TFile* f_nom, TFile* f_sys, TString name);
void scaleShape( TH1F* h_shape, TH1F* h_nom, double scale = 1.);
double getTTbarScale(TFile *dataFile, std::vector<std::vector<TFile*>> MCFiles, TString dirName, TString histName); 
vector<TH1F*> GetPDFVariations(TFile* f_PDF, TFile* f_nominal, TString in_name, TString out_name, bool rebin);
TH1F* calcPDFunc(std::vector<TH1F*> PDFhists, TH1F *h_nominal);
void write_symmetric_uncertainty(TH1F *hist, TH1F* h_nominal, TFile* outputFile);

void write_theta_files(){

  //===========================
  //steering options
  //===========================
 
  TString InputPathPUPPI = "/nfs/dust/cms/user/dreyert/RunII_HeavyResonance/PostSelection/CMSTopTaggerPuppi_NoIso/";
  //TString InputPathPUPPI = "/nfs/dust/cms/user/dreyert/RunII_HeavyResonance/PostSelection/CMSTopTaggerPuppi/";
  TString InputPathCHS = "/nfs/dust/cms/user/dreyert/RunII_80X_Moriond17/PostSelection/AK8CHS/new_data2/";
  // TString InputPathHOTVR = "/nfs/dust/cms/user/dreyert/RunII_80X_Moriond17/PostSelection/HOTVR/newMatching/";
  TString InputPathHOTVR = "/nfs/dust/cms/user/dreyert/RunII_HeavyResonance/PostSelection/HOTVR_PUPPI_NoIso/";
  TString OutputPath = "./";

  TString PDF_dir = "fill_PDF_TRUE";

  //std::vector<TString> MCNames {"TTbar_mergedTop", "TTbar_semimerged", "TTbar_notmerged", "QCD", "DYJets", "ST", "WJets"};
  std::vector<TString> MCNames {"TTbar_mergedTop", "TTbar_semimerged", "TTbar_notmerged", "QCD", "ST", "WJets"};
 
  std::vector< vector<TString> > systematics {
    {"Btag_bc__plus", "BTag_variation_up_bc"}, {"Btag_udsg__plus", "BTag_variation_up_udsg"}, 
    {"Btag_bc__minus", "BTag_variation_down_bc"}, {"Btag_udsg__minus", "BTag_variation_down_udsg"},
    {"MounID__plus", "MuonID_variation_up"}, {"MounID__minus", "MuonID_variation_down"},
    {"Trigger__plus", "MuonTrigger_variation_up"}, {"Trigger__minus", "MuonTrigger_variation_down"}, 
    {"PU__plus", "PU_variation_up"}, {"PU__minus", "PU_variation_down"},
    {"ScaleVariationMuF__plus", "ScaleVariationMuF_up"},  {"ScaleVariationMuF__minus", "ScaleVariationMuF_down"},
    {"ScaleMuR__plus", "ScaleVariationMuR_up"}, {"ScaleMuR__minus", "ScaleVariationMuR_down"}, 
    {"JEC__plus", "jecsmear_direction_up"}, {"JEC__minus", "jecsmear_direction_down"}, 
    {"JER__plus", "jersmear_direction_up"}, {"JER__minus", "jersmear_direction_down"}};

  std::vector< vector<TString> > model_systematics {//};
   {"generator", "aMCatNLO"},
   {"shower_model","Herwig"}};

  std::vector< vector<TString> > observables {
    //{"ProbeJet_pt400_WPXXX_PASSFAIL", "Mass_PASSFAIL__", "400_PASSFAIL"},
    {"ProbeJet_pt300to400_WPXXX_PASSFAIL", "Mass_PASSFAIL__", "300to400_PASSFAIL"},
    {"ProbeJet_pt400to480_WPXXX_PASSFAIL", "Mass_PASSFAIL__", "400to480_PASSFAIL"}, 
    {"ProbeJet_pt480to600_WPXXX_PASSFAIL", "Mass_PASSFAIL__", "480to600_PASSFAIL"}, 
    {"ProbeJet_pt600_WPXXX_PASSFAIL", "Mass_PASSFAIL__",  "600_PASSFAIL"}};
 
  // std::vector<TString> wps {"", "wp1","wp2", "wp3", "wp4", "wp5", "wp1_btag","wp2_btag", "wp3_btag", "wp4_btag", "wp5_btag"};
  //std::vector<TString> wps {"wp1"};
  // std::vector<TString> wps {"wp2"};
  /// std::vector<TString> wps {"wp1_mass_btag","wp2_mass_btag", "wp3_mass_btag", "wp4_mass_btag", "wp5_mass_btag"};
  std::vector<TString> wps { ""};

  //std::vector<TString> variables = { "pt", "tau32"};

  std::vector<TString> variables = {"mass_sub"};
  //std::vector<TString> variables = {"tau32"};
  
  std::vector<bool> rebinning = {true, false};

  std::vector<TString> statsys = {"stat", "sys"};
  //std::vector<TString> statsys = {"stat"};
  
  std::vector<bool> ttbarscaling  = {true, false};

  //std::vector<TString> JetCollections = {"PUPPI", "CHS"};

  //std::vector<TString> JetCollections = {"PUPPI"};
  std::vector<TString> JetCollections = {"HOTVR"};

  std::vector<TString> vPassFail = {"pass", "fail"};

  double bins_pt[] = {300, 400, 440, 480, 520, 560, 600,640,680,720,760,800,840,880,920,960,1000,1060,1120,1180,1260, 1340, 1420, 1500, 1600, 1800};
  int Nbins_pt = sizeof(bins_pt)/sizeof(*bins_pt)-1;

  bool separate = false;

  bool normalizeSYS = true;


  //================
  //start the loop
  //================

  for ( const auto & JetCollection : JetCollections){
    
    TString Path = "";
    if(JetCollection == "PUPPI") Path = InputPathPUPPI;
    if(JetCollection == "CHS") Path = InputPathCHS;
    if(JetCollection == "HOTVR") Path = InputPathHOTVR;

    //=======================
    //load input files
    //=======================
    TFile* dataFile = new TFile(Path+"uhh2.AnalysisModuleRunner.DATA.DATA.root","READ");

    std::vector<std::vector<TFile*>> MCFiles;
    std::vector<std::vector<TFile*>> modelFiles;
    std::vector<TFile*> PDFFiles;

    for( const auto & MCName: MCNames){
      std::vector<TFile*> vFiles; 
      TFile* nominalFile = new TFile(Path+"uhh2.AnalysisModuleRunner.MC."+MCName+".root","READ");
      nominalFile->SetName(MCName);
      vFiles.emplace_back(nominalFile);     
      for(const auto & systematic: systematics){
	TFile* file = new TFile(Path+systematic.at(1)+"/uhh2.AnalysisModuleRunner.MC."+MCName+".root","READ");
	file->SetName(MCName+"__"+systematic.at(0));
	vFiles.emplace_back(file);
      }
      if(MCName.Contains("TTbar")){
	for(const auto & model_systematic: model_systematics){
	  TFile* file = new TFile(Path+model_systematic.at(1)+"/uhh2.AnalysisModuleRunner.MC."+MCName+".root","READ");
	  cout <<Path+model_systematic.at(1)+"/uhh2.AnalysisModuleRunner.MC."+MCName+".root" << endl;
	  file->SetName(MCName+"__"+model_systematic.at(0));
	  //file->SetName(MCName+"_"+model_systematic.at(0));
	  vFiles.emplace_back(file);
	}
      }
      MCFiles.emplace_back(vFiles);

      /*     std::vector<TFile*> vFiles_model;
	     for(const auto & model_systematic: model_systematics){
	     TFile* file = new TFile(Path+model_systematic.at(1)+"/uhh2.AnalysisModuleRunner.MC."+MCName+".root","READ");
	     file->SetName(MCName+"_"+model_systematic.at(0));
	     }
	     modelFiles.emplace_back(vFiles_model);
      */
      if(MCName.Contains("TTbar")){
	TFile* PDFFile = new TFile(Path+PDF_dir+"/uhh2.AnalysisModuleRunner.MC."+MCName+".root" , "READ");
	PDFFile->SetName(MCName+"_PDF");
	PDFFiles.emplace_back(PDFFile);
      }
    }
	

    //=========================================================
    //loop over different variables, working points and pt bins
    //=========================================================

    for(const auto & variable: variables){
      for(const auto & Statsys: statsys){
	for(const auto & rebin: rebinning){
	  for(const auto & scaleTTbar: ttbarscaling){
	    for( const auto & wp: wps){
	      for( unsigned int bin = 0; bin < observables.size(); ++bin){

		//==========================================
		// naming conventions (backward compatible)
		//=========================================
		TString fine = "";
		if(!rebin) fine = "fine/";
		TString scaled = "";
		if(!scaleTTbar) scaled = "_NotScaled";
		
		//set output file in case of combined pass and fail region files
		TFile *outputFile;
		if(!separate){
		  TString binName =  observables.at(bin).at(2);
		  binName.ReplaceAll("_PASSFAIL", "");
		  
		  //TString outName = "thetaFilesTest/"+variable+"/"+fine+"thetaFile_"+binName+"_"+JetCollection+"_"+Statsys+"_"+wp+scaled+".root";
		  TString outName = "thetaFilesNewSYS/"+variable+"/"+fine+"thetaFile_"+binName+"_"+JetCollection+"_"+Statsys+"_"+wp+"_modelUnc"+scaled+".root";
		  outputFile = new TFile(outName,"RECREATE");
		}

		for( const auto & passFail: vPassFail){

		  //====================================
		  // more naming
		  //===================================
		  if(JetCollection != "HOTVR" && wp == "" && passFail == "fail") continue;

		  cout << passFail << endl;
		
		  TString dirName = observables.at(bin).at(0);
		  if(wp.Contains("btag")) dirName.ReplaceAll("WPXXX",wp);
		  else if(wp =="")dirName.ReplaceAll("WPXXX",wp+"all");
		  else dirName.ReplaceAll("WPXXX",wp+"_all");
		  dirName.ReplaceAll("PASSFAIL", passFail);
		
		  TString categoryName =  observables.at(bin).at(1);
		  categoryName.ReplaceAll("PASSFAIL", passFail);

		  TString histName = variable;
		  /* if(JetCollection == "HOTVR" && variable == "mass_sub") {
		    histName = "mass";
		    }*/
		
		  //=======================
		  //get data histogram
		  //=======================
		  TH1F* h_data; 
		  dataFile->GetObject(dirName+"/"+histName, h_data);
		  h_data->SetName(categoryName+"DATA");
 		  if(rebin){
		    if(variable == "pt")  h_data = (TH1F*)h_data->Rebin(Nbins_pt,h_data->GetName(),bins_pt);
		    else h_data->Rebin(2);
		  }
		  outputFile->cd();
		  h_data->Write();

		  //======================
		  //get ttbar scale
		  //======================
		  double scale = 1.;
		  if (scaleTTbar){
		    // if(JetCollection == "HOTVR") scale = getTTbarScale(dataFile, MCFiles, dirName, "mass");
		    // else 
		    scale = getTTbarScale(dataFile, MCFiles, dirName, "mass_sub");
		  }

		  //=======================
		  //write MC
		  //=======================
		  for( unsigned int i = 0; i < MCFiles.size(); ++i){
		    TH1F* h_nominal;
		    for( unsigned int j = 0; j < MCFiles.at(i).size(); ++j){

		      //j == 0 : nominal
		      //j >  0 : systematic

		      if(Statsys == "stat" && j != 0) continue;

		      TString MCName = MCFiles.at(i).at(j)->GetName();

		      cout <<  MCName << endl;

		      TH1F* hist_file;
		      MCFiles.at(i).at(j)->GetObject(dirName+"/"+histName, hist_file);
		      TH1F* hist = (TH1F*)hist_file->Clone(categoryName+MCName);

		      if(rebin){
			if(variable == "pt")  hist = (TH1F*)hist->Rebin(Nbins_pt,hist->GetName(),bins_pt);
			else hist->Rebin(2);
		      }

		      //normalize the systematic variations
		      //(scale with respecto to rate differences between nominal and sys for pass and fail combined)
		      if(j != 0 && normalizeSYS && MCName.Contains("TTbar")){ 
			norm( hist, MCFiles.at(i).at(0), MCFiles.at(i).at(j), dirName+"/"+histName);
		      }

		      // scale TTbar to data -> better starting values for the fit
		      if( MCName.Contains("TTbar") && scaleTTbar) hist->Scale(scale);
		    
		      hist->SetName(categoryName+MCName);

		      if(j == 0) h_nominal = (TH1F*)hist->Clone("h_nominal");

		      //write empty histogram and break if nominal histogram is empty
		      if(j == 0 && variable == "mass_sub" && (hist->Integral() <= 1) ){
			outputFile->cd();
			hist->Reset();
			hist->Write();
			break;
		      }

		      if( j != 0 && !MCName.Contains("__plus") && !MCName.Contains("__minus") ){
			if(MCName.Contains("TTbar")) write_symmetric_uncertainty(hist, h_nominal, outputFile);
		      }else{		    
			outputFile->cd();		
			hist->Write();
		      }
		    }
	
		    //=====
		    //PDF 
		    //=====
		    if (variable == "mass_sub" && ((TString) MCFiles.at(i).at(0)->GetName()).Contains("TTbar") && Statsys == "sys"){

		      TString TTbarName = MCFiles.at(i).at(0)->GetName();
		      if( !((TString)PDFFiles.at(i)->GetName()).Contains(TTbarName) )
			cout << "WARING: wrong PDF uncertainty for "+TTbarName << endl;
		    
		      //get PDF histogram
		      TString inNamePDF = dirName+"/"+histName;
		      std::vector<TH1F*> PDFhists = GetPDFVariations(PDFFiles.at(i), MCFiles.at(i).at(0), inNamePDF, categoryName+TTbarName, rebin);
		      
		      TString inNamePDF_extra = inNamePDF.Copy();
		      bool extraPDFFile = false; 
		      if(inNamePDF_extra.Contains("pass")){
			inNamePDF_extra.ReplaceAll("pass", "fail");
			extraPDFFile = true;
		      }
		      else if(inNamePDF_extra.Contains("fail")){
			inNamePDF_extra.ReplaceAll("fail", "pass");
			extraPDFFile = true;
		      }
		      std::vector<TH1F*> PDFhists_extra = GetPDFVariations(PDFFiles.at(i), MCFiles.at(i).at(0), inNamePDF_extra, categoryName+TTbarName, rebin);
		      
		      outputFile->cd(); 
		      for(unsigned int var = 0; var < 2; ++var){
			TH1F* h_pdf = (TH1F*)PDFhists.at(var)->Clone(PDFhists.at(var)->GetName());

			if(normalizeSYS){
			  TH1F* h_nom = (TH1F*)PDFhists.at(2)->Clone("h_nom");
			  h_nom->Add(PDFhists_extra.at(2));
     
			  TH1F* h_pdf_total = (TH1F*)PDFhists.at(var)->Clone("h_pdf_total");
			  h_pdf_total->Add(PDFhists_extra.at(var));
		   
			  h_pdf->Scale(h_nom->Integral(0,h_nom->GetNbinsX()+1)/h_pdf_total->Integral(0, h_pdf_total->GetNbinsX()+1));
			}
			if(scaleTTbar) h_pdf->Scale(scale);
			
			h_pdf->Write();
		   
		      }
		    }
		  }
		}
	    
		outputFile->Close();
		   
	      }
	    }
	  }
	}
      }
    }

    //======================
    //close input files
    //=====================
    dataFile->Close();
    for( const auto & vFiles : MCFiles){
      for( const auto & file : vFiles){
	file->Close();
      }
    }
    for(const auto & file : PDFFiles){
      file->Close();
    }
    
  }  
}




void norm(TH1F* h_sys, TFile* f_nom, TFile* f_sys, TString name){

  TString name_pass = name;
  TString name_fail = name;

  bool passfail = false;

  if(name.Contains("pass")){
    name_fail.ReplaceAll("pass", "fail");
    passfail = true;
  }
  else if(name.Contains("fail")){
    name_pass.ReplaceAll("fail", "pass");
    passfail = true;
  }
  else cout << "WARNING: no pass or fail region specified. Systematic template will not be normalized!" << endl;
    
  if(passfail){
    TH1F* h_nom_pass, *h_nom_fail, * h_sys_pass, * h_sys_fail;
    f_nom->GetObject(name_pass, h_nom_pass); 
    f_nom->GetObject(name_fail, h_nom_fail); 
    f_sys->GetObject(name_pass, h_sys_pass); 
    f_sys->GetObject(name_fail, h_sys_fail); 

    TH1F* h_nom_total = (TH1F*)h_nom_pass->Clone("h_nom_total");
    h_nom_total->Add(h_nom_fail);

    TH1F* h_sys_total = (TH1F*)h_sys_pass->Clone("h_sys_total");
    h_sys_total->Add(h_sys_fail);

    double N_nom = h_nom_total->Integral(0,h_nom_total->GetNbinsX()+1);
    double N_sys = h_sys_total->Integral(0,h_sys_total->GetNbinsX()+1);
    double scale = N_nom/N_sys;
    
    h_sys->Scale(scale);
  }
}

void scaleShape( TH1F* h_shape, TH1F* h_nom, double scale = 1.){

  for(unsigned int bin = 0; bin < h_shape->GetNbinsX(); ++bin){
    
    double binContent = h_nom->GetBinContent(bin+1) + (scale * (h_shape->GetBinContent(bin+1) - h_nom->GetBinContent(bin+1)));
    h_shape->SetBinContent(bin+1, binContent);

  }

}

double getTTbarScale(TFile* dataFile, std::vector<std::vector<TFile*>> MCFiles, TString dirName, TString histName){
  
  double scale = 1.;

  TH1F *hDataScale = (TH1F*)dataFile->Get(dirName+"/"+histName);
  double int_data = hDataScale->Integral();

  double int_bkg = 0.;
  double int_tt = 0.;
  for(const auto & mc: MCFiles){
    TH1F * hist = (TH1F*)mc.at(0)->Get(dirName+"/"+histName);
    if( ((TString)mc.at(0)->GetName()).Contains("TTbar")){
      int_tt += hist->Integral();
    }else{
      int_bkg += hist->Integral();
    }
  }
  scale = (int_data-int_bkg)/int_tt;

  return scale;
}

vector<TH1F*> GetPDFVariations(TFile* f_PDF, TFile* f_nominal, TString in_name, TString out_name, bool rebin){
 
  std::vector<TH1F*> PDFhists;

  TH1F* h_nominal_FILE;
  f_nominal->GetObject(in_name, h_nominal_FILE);
  TH1F* h_nominal = (TH1F*) h_nominal_FILE->Clone("h_nominalPDF");
  if(rebin) h_nominal->Rebin(2); 
    
  for(unsigned int k = 0; k < 100; ++k){
    TH1F* h_PDF1;
    f_PDF->GetObject(in_name+"_PDF_"+std::to_string(k), h_PDF1);
    TH1F* h_PDF = (TH1F*)h_PDF1->Clone("h_PDF");
    if(rebin) h_PDF->Rebin(2);
    PDFhists.emplace_back(h_PDF);
  }
		      
  TH1F *h_PDF_scale = calcPDFunc(PDFhists, h_nominal);
  TH1F* h_PDF_up = (TH1F*)h_nominal->Clone(out_name+"__PDF__plus");
  TH1F* h_PDF_down = (TH1F*)h_nominal->Clone(out_name+"__PDF__minus");
  h_PDF_up->Add(h_PDF_scale);
  h_PDF_down->Add(h_PDF_scale, -1);

  std::vector<TH1F*> outhists {h_PDF_up, h_PDF_down, h_nominal};

  return outhists;

}
 
TH1F* calcPDFunc(std::vector<TH1F*> PDFhists, TH1F* h_nominal){
  
  TH1F* h_PDF = PDFhists.at(0);
  h_PDF->Add(h_nominal,-1);
  h_PDF->Multiply(h_PDF);
  
  for( unsigned int i = 1; i < PDFhists.size(); ++i){
    TH1F* h_PDF2 = PDFhists.at(i);
    //  cout << "add nominal" << h_PDF2->GetNbinsX()  << h_nominal->GetNbinsX() << endl;
    h_PDF2->Add(h_nominal,-1);
    h_PDF2->Multiply(h_PDF2);
    h_PDF->Add(h_PDF2);
  }
  TH1F* h_PDFscale = (TH1F*)h_PDF->Clone("h_PDF_scale");
  h_PDFscale->Reset();
  for(unsigned int bin = 1; bin <= h_PDFscale->GetNbinsX(); ++bin){
    double PDFscale = sqrt(h_PDF->GetBinContent(bin)/100);
    h_PDFscale->SetBinContent(bin,PDFscale);
  }

  return h_PDFscale;
}

void write_symmetric_uncertainty(TH1F *hist, TH1F* h_nominal, TFile* outputFile){

  cout << "start" << endl;
  
  TString histName = hist->GetName();
  cout << h_nominal->GetName() << endl;


  cout << "clone hists" << endl;
  TH1F* hist_up = (TH1F*)hist->Clone(histName+"__plus");
  // TH1F* hist_up = (TH1F*)hist->Clone(histName);
  TH1F* hist_down = (TH1F*)h_nominal->Clone(histName+"__minus");

  cout << "calc" << endl;

  TH1F* hist_diff = (TH1F*)hist->Clone(histName+"__diff"); 
  hist_diff->Add(h_nominal, -1.);
  hist_down->Add(hist_diff, -1.);

  cout << "write hists" << endl;
 
  outputFile->cd();
  hist_up->Write();
  hist_down->Write();

  cout << "DONE" << endl;
 
}
