#include "UHH2/TopTagging/include/TopTaggingUtils.h"
#include "UHH2/core/include/LorentzVector.h"

bool TopJetLeptonDeltaRCleaner::process(uhh2::Event & event) {

  assert(event.topjets);
  std::vector<TopJet> cleaned_topjets;

  for(const auto & tjet : *event.topjets){
    bool skip_tjet(false);

    if(event.muons){
      for(const auto & muo : *event.muons)
        if(uhh2::deltaR(tjet, muo) < minDR_) skip_tjet = true;
    }

    if(skip_tjet) continue;

    if(event.electrons){
      for(const auto & ele : *event.electrons)
        if(uhh2::deltaR(tjet, ele) < minDR_) skip_tjet = true;
    }

    if(!skip_tjet) cleaned_topjets.push_back(tjet);
  }

  event.topjets->clear();
  event.topjets->reserve(cleaned_topjets.size());
  for(auto & j : cleaned_topjets) event.topjets->push_back(j);

  return true;
}

TopJetCorrectionModules::TopJetCorrectionModules(uhh2::Context & ctx, jet_type type) {

  _type = type;
  is_mc = (ctx.get("dataset_type") == "MC");

  std::vector<std::string> JEC_TOP_MC, JEC_TOP_B, JEC_TOP_C, JEC_TOP_D, JEC_TOP_E, JEC_TOP_F;
  std::vector<std::string> JEC_SUB_MC, JEC_SUB_B, JEC_SUB_C, JEC_SUB_D, JEC_SUB_E, JEC_SUB_F;
  
  if(_type == AK8_PUPPI || _type == CA15_PUPPI){
    if(is_mc){
      JEC_TOP_MC = JERFiles::Fall17_17Nov2017_V8_L123_AK8PFPuppi_MC;
      JEC_SUB_MC = JERFiles::Fall17_17Nov2017_V8_L123_AK4PFPuppi_MC;
    }else{
      // JEC_TOP_MC = JERFiles::Fall17_17Nov2017_V4_L123_AK8PFPuppi_MC;
      // JEC_SUB_MC = JERFiles::Fall17_17Nov2017_V4_L123_AK4PFPuppi_MC;
     
      JEC_TOP_B = JERFiles::Fall17_17Nov2017_V8_B_L123_AK8PFPuppi_DATA;
      JEC_TOP_C = JERFiles::Fall17_17Nov2017_V8_C_L123_AK8PFPuppi_DATA;
      JEC_TOP_D = JERFiles::Fall17_17Nov2017_V8_D_L123_AK8PFPuppi_DATA;
      JEC_TOP_E = JERFiles::Fall17_17Nov2017_V8_E_L123_AK8PFPuppi_DATA;
      JEC_TOP_F = JERFiles::Fall17_17Nov2017_V8_F_L123_AK8PFPuppi_DATA;
     
      JEC_SUB_B = JERFiles::Fall17_17Nov2017_V8_B_L123_AK4PFPuppi_DATA;
      JEC_SUB_C = JERFiles::Fall17_17Nov2017_V8_C_L123_AK4PFPuppi_DATA;
      JEC_SUB_D = JERFiles::Fall17_17Nov2017_V8_D_L123_AK4PFPuppi_DATA;
      JEC_SUB_E = JERFiles::Fall17_17Nov2017_V8_E_L123_AK4PFPuppi_DATA;
      JEC_SUB_F = JERFiles::Fall17_17Nov2017_V8_F_L123_AK4PFPuppi_DATA;
    }
  }
  
  if(_type == AK8_CHS || _type == CA15_CHS){
    if(is_mc){
      JEC_TOP_MC = JERFiles::Fall17_17Nov2017_V8_L123_AK8PFchs_MC;
      JEC_SUB_MC = JERFiles::Fall17_17Nov2017_V8_L123_AK4PFchs_MC;
    }else{
      //JEC_TOP_MC = JERFiles::Fall17_17Nov2017_V4_L123_AK8PFchs_MC;
      //JEC_SUB_MC = JERFiles::Fall17_17Nov2017_V4_L123_AK4PFchs_MC;
 
      JEC_TOP_B = JERFiles::Fall17_17Nov2017_V8_B_L123_AK8PFchs_DATA;
      JEC_TOP_C = JERFiles::Fall17_17Nov2017_V8_C_L123_AK8PFchs_DATA;
      JEC_TOP_D = JERFiles::Fall17_17Nov2017_V8_D_L123_AK8PFchs_DATA;
      JEC_TOP_E = JERFiles::Fall17_17Nov2017_V8_E_L123_AK8PFchs_DATA;
      JEC_TOP_F = JERFiles::Fall17_17Nov2017_V8_F_L123_AK8PFchs_DATA;
     
      JEC_SUB_B = JERFiles::Fall17_17Nov2017_V8_B_L123_AK4PFchs_DATA;
      JEC_SUB_C = JERFiles::Fall17_17Nov2017_V8_C_L123_AK4PFchs_DATA;
      JEC_SUB_D = JERFiles::Fall17_17Nov2017_V8_D_L123_AK4PFchs_DATA;
      JEC_SUB_E = JERFiles::Fall17_17Nov2017_V8_E_L123_AK4PFchs_DATA;
      JEC_SUB_F = JERFiles::Fall17_17Nov2017_V8_F_L123_AK4PFchs_DATA;
    }
  }

  if(_type == HOTVR_CHS){
    if(is_mc){
      JEC_TOP_MC = JERFiles::Fall17_17Nov2017_V8_L123_AK8PFchs_MC;
      JEC_SUB_MC = JERFiles::Fall17_17Nov2017_V8_L123_AK4PFchs_MC;
    }else{
      JEC_TOP_MC = JERFiles::Fall17_17Nov2017_V8_L123_AK8PFchs_MC;
      JEC_SUB_MC = JERFiles::Fall17_17Nov2017_V8_L123_AK4PFchs_MC;
      // JEC_TOP_BCD = JERFiles::Summer16_23Sep2016_V4_BCD_L23_AK8PFchs_DATA;
      // JEC_TOP_EF = JERFiles::Summer16_23Sep2016_V4_EF_L23_AK8PFchs_DATA;
      // JEC_TOP_FG = JERFiles::Summer16_23Sep2016_V4_G_L23_AK8PFchs_DATA;
      // JEC_TOP_H = JERFiles::Summer16_23Sep2016_V4_H_L23_AK8PFchs_DATA;
     
      // JEC_SUB_BCD = JERFiles::Summer16_23Sep2016_V4_BCD_L23_AK4PFchs_DATA;
      // JEC_SUB_EF = JERFiles::Summer16_23Sep2016_V4_EF_L23_AK4PFchs_DATA;
      // JEC_SUB_FG = JERFiles::Summer16_23Sep2016_V4_G_L23_AK4PFchs_DATA;
      // JEC_SUB_H = JERFiles::Summer16_23Sep2016_V4_H_L23_AK4PFchs_DATA;
    }
  }

  if(is_mc){
    topjet_corrector_MC.reset(new TopJetCorrector(ctx, JEC_TOP_MC));
    subjet_corrector_MC.reset(new SubJetCorrector(ctx, JEC_SUB_MC));
  }
  else{
    //topjet_corrector_MC.reset(new TopJetCorrector(ctx, JEC_TOP_MC));
    //subjet_corrector_MC.reset(new SubJetCorrector(ctx, JEC_SUB_MC));
    topjet_corrector_B.reset(new TopJetCorrector(ctx, JEC_TOP_B));
    topjet_corrector_C.reset(new TopJetCorrector(ctx, JEC_TOP_C));
    topjet_corrector_D.reset(new TopJetCorrector(ctx, JEC_TOP_D));
    topjet_corrector_E.reset(new TopJetCorrector(ctx, JEC_TOP_E));
    topjet_corrector_F.reset(new TopJetCorrector(ctx, JEC_TOP_F));

    subjet_corrector_B.reset(new SubJetCorrector(ctx, JEC_SUB_B));
    subjet_corrector_C.reset(new SubJetCorrector(ctx, JEC_SUB_C));
    subjet_corrector_D.reset(new SubJetCorrector(ctx, JEC_SUB_D));
    subjet_corrector_E.reset(new SubJetCorrector(ctx, JEC_SUB_E));
    subjet_corrector_F.reset(new SubJetCorrector(ctx, JEC_SUB_F));
  }
}


bool TopJetCorrectionModules::process(uhh2::Event & event) { 

  assert(event.topjets);
 
  if(is_mc){
    topjet_corrector_MC->process(event);
    subjet_corrector_MC->process(event);
  }else{
    //   topjet_corrector_MC->process(event);
    //  subjet_corrector_MC->process(event);
 
    if(event.run <= runnr_B) {
      topjet_corrector_B->process(event);
      subjet_corrector_B->process(event);
    }         
    else if(event.run <= runnr_C) {
      topjet_corrector_C->process(event);
      subjet_corrector_C->process(event);
    } 
    else if(event.run <= runnr_D) {
      topjet_corrector_D->process(event);
      subjet_corrector_D->process(event);
    } 
    else if(event.run <= runnr_E) {
      topjet_corrector_E->process(event);
      subjet_corrector_E->process(event);
    }
    else if(event.run > runnr_E) {
      topjet_corrector_F->process(event);
      subjet_corrector_F->process(event);
    }
   
  }
  return true;
}


bool HOTVRPileupCorrectionModule::process(uhh2::Event & event){

  double rho = event.rho;

  for(auto & topjet : *event.topjets){
    auto subjets = topjet.subjets();
    for(auto & subjet : subjets){
      
      double subjet_pt = subjet.pt(); 
      double subjet_area = subjet.jetArea();

      double pileup_factor = 1. - ((rho*subjet_area)/subjet_pt);

      subjet.set_JEC_factor_raw(1.); //don't set the raw JEC? (L2L3 corrections should be applied on top) 
      if(_area_correction){
	subjet.set_JEC_L1factor_raw(1./pileup_factor);
	subjet.set_v4(subjet.v4()*pileup_factor);
      }
    }
    topjet.set_subjets(move(subjets));
  }

  return true;
}


bool TopJetGroomer::process(uhh2::Event & event) { 

  assert(event.topjets);

  for(auto & topjet : *event.topjets){
    LorentzVector subjet_sum(0,0,0,0);

    for(const auto & subjet : topjet.subjets()){
      if(_corrected) subjet_sum += subjet.v4();
      else subjet_sum += subjet.v4()*subjet.JEC_factor_raw();
    }
    topjet.set_v4(subjet_sum);
    if(!_corrected) topjet.set_JEC_factor_raw(1.);
  }

  return true;
}

bool GetLeadingBjetLepHem(const uhh2::Event &event, Jet &bjet, JetId btag){
  //get the muon
  Muon mu = event.muons->at(0);
  
  //get the bjet
  double max_pt = 0.;
  double pi = 3.14159265359;

  bool b_candidate_found = false;

  for( const auto & ak4jet : *event.jets){
    Jet b_candidate;
    if( btag(ak4jet, event) && (uhh2::deltaPhi(ak4jet,mu) < (2*pi/3)) ){
      b_candidate = ak4jet;
      b_candidate_found = true; 
    }
    if(b_candidate_found){
      if( b_candidate.pt() > max_pt){
	bjet = b_candidate;
	max_pt = b_candidate.pt();
      }
    }
  } 

  if(b_candidate_found) return true;
  else std::cout << "No Bjet candidate found in leptonic hemisphere in event:  " << event.event << std::endl;
  return false;
}

PartonShowerWeight::PartonShowerWeight(uhh2::Context & ctx, std::string sys){

  auto dataset_type = ctx.get("dataset_type");
  bool is_mc = dataset_type == "MC";

  weightIndex = -1;
  if(sys == "ISRup_sqrt2") weightIndex = 2;
  else if(sys == "ISRup_2") weightIndex = 6;
  else if(sys == "ISRup_4") weightIndex = 10;
  
  else if(sys == "ISRdown_sqrt2") weightIndex = 4;
  else if(sys == "ISRdown_2") weightIndex = 8;
  else if(sys == "ISRdown_4") weightIndex = 12;

  else if(sys == "FSRup_sqrt2") weightIndex = 3;
  else if(sys == "FSRup_2") weightIndex = 7;
  else if(sys == "FSRup_4") weightIndex = 11;

  else if(sys == "FSRdown_sqrt2") weightIndex = 5;
  else if(sys == "FSRdown_2") weightIndex = 9;
  else if(sys == "FSRdown_4") weightIndex = 13;
  else
    std::cout << "No parton shower variation sepecified -> module will not have an effect" << std::endl;

  if(!is_mc){
    weightIndex = -1;
    std::cout << "No MC sample. Parton shower weights will have no effec on data" << std::endl;
  }
}

bool PartonShowerWeight::process(uhh2::Event & event){

  if( weightIndex == -1) return true;
  if( event.genInfo->weights().size() == 1) {
    //cout << "no parton shower weights stored. Nothing to be done." << endl;
    return true;
  }
  
  double centralWeight = event.genInfo->weights().at(0);
  double PSweight = event.genInfo->weights().at(weightIndex);
  
  double factor = PSweight/centralWeight;
  event.weight *= factor;
			 
  return true;
}


//AK4
const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_B_L123_AK4PFchs_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017B_V8_DATA/Fall17_17Nov2017B_V8_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017B_V8_DATA/Fall17_17Nov2017B_V8_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017B_V8_DATA/Fall17_17Nov2017B_V8_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017B_V8_DATA/Fall17_17Nov2017B_V8_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_C_L123_AK4PFchs_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017C_V8_DATA/Fall17_17Nov2017C_V8_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017C_V8_DATA/Fall17_17Nov2017C_V8_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017C_V8_DATA/Fall17_17Nov2017C_V8_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017C_V8_DATA/Fall17_17Nov2017C_V8_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_D_L123_AK4PFchs_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017D_V8_DATA/Fall17_17Nov2017D_V8_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017D_V8_DATA/Fall17_17Nov2017D_V8_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017D_V8_DATA/Fall17_17Nov2017D_V8_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017D_V8_DATA/Fall17_17Nov2017D_V8_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_E_L123_AK4PFchs_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017E_V8_DATA/Fall17_17Nov2017E_V8_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017E_V8_DATA/Fall17_17Nov2017E_V8_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017E_V8_DATA/Fall17_17Nov2017E_V8_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017E_V8_DATA/Fall17_17Nov2017E_V8_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_F_L123_AK4PFchs_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017F_V8_DATA/Fall17_17Nov2017F_V8_DATA_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017F_V8_DATA/Fall17_17Nov2017F_V8_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017F_V8_DATA/Fall17_17Nov2017F_V8_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017F_V8_DATA/Fall17_17Nov2017F_V8_DATA_L2L3Residual_AK4PFchs.txt",
};



const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_L123_AK4PFchs_MC = {
  "JECDatabase/textFiles/Fall17_17Nov2017_V8_MC/Fall17_17Nov2017_V8_MC_L1FastJet_AK4PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017_V8_MC/Fall17_17Nov2017_V8_MC_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017_V8_MC/Fall17_17Nov2017_V8_MC_L3Absolute_AK4PFchs.txt",
};


const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_L123_AK4PFPuppi_MC = {
  "JECDatabase/textFiles/Fall17_17Nov2017_V8_MC/Fall17_17Nov2017_V8_MC_L1FastJet_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017_V8_MC/Fall17_17Nov2017_V8_MC_L2Relative_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017_V8_MC/Fall17_17Nov2017_V8_MC_L3Absolute_AK4PFPuppi.txt",
};

const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_L123_AK8PFchs_MC = {
  "JECDatabase/textFiles/Fall17_17Nov2017_V8_MC/Fall17_17Nov2017_V8_MC_L1FastJet_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017_V8_MC/Fall17_17Nov2017_V8_MC_L2Relative_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017_V8_MC/Fall17_17Nov2017_V8_MC_L3Absolute_AK8PFchs.txt",
};

const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_L123_AK8PFPuppi_MC = {
  "JECDatabase/textFiles/Fall17_17Nov2017_V8_MC/Fall17_17Nov2017_V8_MC_L1FastJet_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017_V8_MC/Fall17_17Nov2017_V8_MC_L2Relative_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017_V8_MC/Fall17_17Nov2017_V8_MC_L3Absolute_AK8PFPuppi.txt",
};


//AK4 Puppi
const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_B_L123_AK4PFPuppi_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017B_V8_DATA/Fall17_17Nov2017B_V8_DATA_L1FastJet_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017B_V8_DATA/Fall17_17Nov2017B_V8_DATA_L2Relative_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017B_V8_DATA/Fall17_17Nov2017B_V8_DATA_L3Absolute_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017B_V8_DATA/Fall17_17Nov2017B_V8_DATA_L2L3Residual_AK4PFPuppi.txt",
};

const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_C_L123_AK4PFPuppi_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017C_V8_DATA/Fall17_17Nov2017C_V8_DATA_L1FastJet_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017C_V8_DATA/Fall17_17Nov2017C_V8_DATA_L2Relative_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017C_V8_DATA/Fall17_17Nov2017C_V8_DATA_L3Absolute_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017C_V8_DATA/Fall17_17Nov2017C_V8_DATA_L2L3Residual_AK4PFPuppi.txt",
};

const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_D_L123_AK4PFPuppi_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017D_V8_DATA/Fall17_17Nov2017D_V8_DATA_L1FastJet_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017D_V8_DATA/Fall17_17Nov2017D_V8_DATA_L2Relative_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017D_V8_DATA/Fall17_17Nov2017D_V8_DATA_L3Absolute_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017D_V8_DATA/Fall17_17Nov2017D_V8_DATA_L2L3Residual_AK4PFPuppi.txt",
};

const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_E_L123_AK4PFPuppi_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017E_V8_DATA/Fall17_17Nov2017E_V8_DATA_L1FastJet_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017E_V8_DATA/Fall17_17Nov2017E_V8_DATA_L2Relative_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017E_V8_DATA/Fall17_17Nov2017E_V8_DATA_L3Absolute_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017E_V8_DATA/Fall17_17Nov2017E_V8_DATA_L2L3Residual_AK4PFPuppi.txt",
};

const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_F_L123_AK4PFPuppi_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017F_V8_DATA/Fall17_17Nov2017F_V8_DATA_L1FastJet_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017F_V8_DATA/Fall17_17Nov2017F_V8_DATA_L2Relative_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017F_V8_DATA/Fall17_17Nov2017F_V8_DATA_L3Absolute_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017F_V8_DATA/Fall17_17Nov2017F_V8_DATA_L2L3Residual_AK4PFPuppi.txt",
};


//AK8
const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_B_L123_AK8PFchs_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017B_V8_DATA/Fall17_17Nov2017B_V8_DATA_L1FastJet_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017B_V8_DATA/Fall17_17Nov2017B_V8_DATA_L2Relative_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017B_V8_DATA/Fall17_17Nov2017B_V8_DATA_L3Absolute_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017B_V8_DATA/Fall17_17Nov2017B_V8_DATA_L2L3Residual_AK8PFchs.txt",
};

const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_C_L123_AK8PFchs_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017C_V8_DATA/Fall17_17Nov2017C_V8_DATA_L1FastJet_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017C_V8_DATA/Fall17_17Nov2017C_V8_DATA_L2Relative_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017C_V8_DATA/Fall17_17Nov2017C_V8_DATA_L3Absolute_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017C_V8_DATA/Fall17_17Nov2017C_V8_DATA_L2L3Residual_AK8PFchs.txt",
};

const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_D_L123_AK8PFchs_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017D_V8_DATA/Fall17_17Nov2017D_V8_DATA_L1FastJet_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017D_V8_DATA/Fall17_17Nov2017D_V8_DATA_L2Relative_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017D_V8_DATA/Fall17_17Nov2017D_V8_DATA_L3Absolute_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017D_V8_DATA/Fall17_17Nov2017D_V8_DATA_L2L3Residual_AK8PFchs.txt",
};

const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_E_L123_AK8PFchs_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017E_V8_DATA/Fall17_17Nov2017E_V8_DATA_L1FastJet_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017E_V8_DATA/Fall17_17Nov2017E_V8_DATA_L2Relative_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017E_V8_DATA/Fall17_17Nov2017E_V8_DATA_L3Absolute_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017E_V8_DATA/Fall17_17Nov2017E_V8_DATA_L2L3Residual_AK8PFchs.txt",
};

const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_F_L123_AK8PFchs_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017F_V8_DATA/Fall17_17Nov2017F_V8_DATA_L1FastJet_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017F_V8_DATA/Fall17_17Nov2017F_V8_DATA_L2Relative_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017F_V8_DATA/Fall17_17Nov2017F_V8_DATA_L3Absolute_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017F_V8_DATA/Fall17_17Nov2017F_V8_DATA_L2L3Residual_AK8PFchs.txt",
};



//AK8 Puppi
const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_B_L123_AK8PFPuppi_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017B_V8_DATA/Fall17_17Nov2017B_V8_DATA_L1FastJet_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017B_V8_DATA/Fall17_17Nov2017B_V8_DATA_L2Relative_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017B_V8_DATA/Fall17_17Nov2017B_V8_DATA_L3Absolute_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017B_V8_DATA/Fall17_17Nov2017B_V8_DATA_L2L3Residual_AK8PFPuppi.txt",
};

const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_C_L123_AK8PFPuppi_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017C_V8_DATA/Fall17_17Nov2017C_V8_DATA_L1FastJet_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017C_V8_DATA/Fall17_17Nov2017C_V8_DATA_L2Relative_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017C_V8_DATA/Fall17_17Nov2017C_V8_DATA_L3Absolute_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017C_V8_DATA/Fall17_17Nov2017C_V8_DATA_L2L3Residual_AK8PFPuppi.txt",
};

const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_D_L123_AK8PFPuppi_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017D_V8_DATA/Fall17_17Nov2017D_V8_DATA_L1FastJet_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017D_V8_DATA/Fall17_17Nov2017D_V8_DATA_L2Relative_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017D_V8_DATA/Fall17_17Nov2017D_V8_DATA_L3Absolute_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017D_V8_DATA/Fall17_17Nov2017D_V8_DATA_L2L3Residual_AK8PFPuppi.txt",
};

const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_E_L123_AK8PFPuppi_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017E_V8_DATA/Fall17_17Nov2017E_V8_DATA_L1FastJet_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017E_V8_DATA/Fall17_17Nov2017E_V8_DATA_L2Relative_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017E_V8_DATA/Fall17_17Nov2017E_V8_DATA_L3Absolute_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017E_V8_DATA/Fall17_17Nov2017E_V8_DATA_L2L3Residual_AK8PFPuppi.txt",
};

const std::vector<std::string> JERFiles::Fall17_17Nov2017_V8_F_L123_AK8PFPuppi_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017F_V8_DATA/Fall17_17Nov2017F_V8_DATA_L1FastJet_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017F_V8_DATA/Fall17_17Nov2017F_V8_DATA_L2Relative_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017F_V8_DATA/Fall17_17Nov2017F_V8_DATA_L3Absolute_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017F_V8_DATA/Fall17_17Nov2017F_V8_DATA_L2L3Residual_AK8PFPuppi.txt",
};
