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
  is_mc = ctx.get("dataset_type") == "MC";

  std::vector<std::string> JEC_TOP_MC, JEC_TOP_BCD, JEC_TOP_EF, JEC_TOP_FG, JEC_TOP_H;
  std::vector<std::string> JEC_SUB_MC, JEC_SUB_BCD, JEC_SUB_EF, JEC_SUB_FG, JEC_SUB_H;
  
  if(_type == AK8_PUPPI || _type == CA15_PUPPI){
    if(is_mc){
      JEC_TOP_MC = JERFiles::Summer16_23Sep2016_V4_L123_AK8PFPuppi_MC;
      JEC_SUB_MC = JERFiles::Summer16_23Sep2016_V4_L123_AK4PFPuppi_MC;
    }else{
      JEC_TOP_BCD = JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK8PFPuppi_DATA;
      JEC_TOP_EF = JERFiles::Summer16_23Sep2016_V4_EF_L123_AK8PFPuppi_DATA;
      JEC_TOP_FG = JERFiles::Summer16_23Sep2016_V4_G_L123_AK8PFPuppi_DATA;
      JEC_TOP_H = JERFiles::Summer16_23Sep2016_V4_H_L123_AK8PFPuppi_DATA;
     
      JEC_SUB_BCD = JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK4PFPuppi_DATA;
      JEC_SUB_EF = JERFiles::Summer16_23Sep2016_V4_EF_L123_AK4PFPuppi_DATA;
      JEC_SUB_FG = JERFiles::Summer16_23Sep2016_V4_G_L123_AK4PFPuppi_DATA;
      JEC_SUB_H = JERFiles::Summer16_23Sep2016_V4_H_L123_AK4PFPuppi_DATA;
    }
  }
  
  if(_type == AK8_CHS || _type == CA15_CHS){
    if(is_mc){
      JEC_TOP_MC = JERFiles::Summer16_23Sep2016_V4_L123_AK8PFchs_MC;
      JEC_SUB_MC = JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC;
    }else{
      JEC_TOP_BCD = JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK8PFchs_DATA;
      JEC_TOP_EF = JERFiles::Summer16_23Sep2016_V4_EF_L123_AK8PFchs_DATA;
      JEC_TOP_FG = JERFiles::Summer16_23Sep2016_V4_G_L123_AK8PFchs_DATA;
      JEC_TOP_H = JERFiles::Summer16_23Sep2016_V4_H_L123_AK8PFchs_DATA;
     
      JEC_SUB_BCD = JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK4PFchs_DATA;
      JEC_SUB_EF = JERFiles::Summer16_23Sep2016_V4_EF_L123_AK4PFchs_DATA;
      JEC_SUB_FG = JERFiles::Summer16_23Sep2016_V4_G_L123_AK4PFchs_DATA;
      JEC_SUB_H = JERFiles::Summer16_23Sep2016_V4_H_L123_AK4PFchs_DATA;
    }
  }

  if(_type == HOTVR_CHS){
    if(is_mc){
      JEC_TOP_MC = JERFiles::Summer16_23Sep2016_V4_L23_AK8PFchs_MC;
      JEC_SUB_MC = JERFiles::Summer16_23Sep2016_V4_L23_AK4PFchs_MC;
    }else{
      JEC_TOP_BCD = JERFiles::Summer16_23Sep2016_V4_BCD_L23_AK8PFchs_DATA;
      JEC_TOP_EF = JERFiles::Summer16_23Sep2016_V4_EF_L23_AK8PFchs_DATA;
      JEC_TOP_FG = JERFiles::Summer16_23Sep2016_V4_G_L23_AK8PFchs_DATA;
      JEC_TOP_H = JERFiles::Summer16_23Sep2016_V4_H_L23_AK8PFchs_DATA;
     
      JEC_SUB_BCD = JERFiles::Summer16_23Sep2016_V4_BCD_L23_AK4PFchs_DATA;
      JEC_SUB_EF = JERFiles::Summer16_23Sep2016_V4_EF_L23_AK4PFchs_DATA;
      JEC_SUB_FG = JERFiles::Summer16_23Sep2016_V4_G_L23_AK4PFchs_DATA;
      JEC_SUB_H = JERFiles::Summer16_23Sep2016_V4_H_L23_AK4PFchs_DATA;
    }
  }

  if(is_mc){
    topjet_corrector_MC.reset(new TopJetCorrector(ctx, JEC_TOP_MC));
    subjet_corrector_MC.reset(new SubJetCorrector(ctx, JEC_SUB_MC));
  }
  else{
    topjet_corrector_BCD.reset(new TopJetCorrector(ctx, JEC_TOP_BCD));
    topjet_corrector_EF.reset(new TopJetCorrector(ctx, JEC_TOP_EF));
    topjet_corrector_FG.reset(new TopJetCorrector(ctx, JEC_TOP_FG));
    topjet_corrector_H.reset(new TopJetCorrector(ctx, JEC_TOP_H));

    subjet_corrector_BCD.reset(new SubJetCorrector(ctx, JEC_SUB_BCD));
    subjet_corrector_EF.reset(new SubJetCorrector(ctx, JEC_SUB_EF));
    subjet_corrector_FG.reset(new SubJetCorrector(ctx, JEC_SUB_FG));
    subjet_corrector_H.reset(new SubJetCorrector(ctx, JEC_SUB_H));
  }
}


bool TopJetCorrectionModules::process(uhh2::Event & event) { 

  assert(event.topjets);

  if(is_mc){
    topjet_corrector_MC->process(event);
    subjet_corrector_MC->process(event);
  }else{
    if(event.run <= runnr_BCD) {
      topjet_corrector_BCD->process(event);
      subjet_corrector_BCD->process(event);
    }         
    else if(event.run < runnr_EFearly) {
      topjet_corrector_EF->process(event);
      subjet_corrector_EF->process(event);
    } 
    else if(event.run <= runnr_FlateG) {
      topjet_corrector_FG->process(event);
      subjet_corrector_FG->process(event);
    } 
    else if(event.run > runnr_FlateG) {
      topjet_corrector_H->process(event);
      subjet_corrector_H->process(event);
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

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK4PFPuppi_DATA = {
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFPuppi.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_EF_L123_AK4PFPuppi_DATA = {
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L1FastJet_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2Relative_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L3Absolute_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK4PFPuppi.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_G_L123_AK4PFPuppi_DATA = {
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L1FastJet_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2Relative_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L3Absolute_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2L3Residual_AK4PFPuppi.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_H_L123_AK4PFPuppi_DATA = {
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L1FastJet_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2Relative_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L3Absolute_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PFPuppi.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_L123_AK4PFPuppi_MC = {
  "JECDatabase/textFiles/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L1FastJet_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L2Relative_AK4PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L3Absolute_AK4PFPuppi.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_L123_AK8PFPuppi_MC = {
  "JECDatabase/textFiles/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L1FastJet_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L2Relative_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L3Absolute_AK8PFPuppi.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK8PFPuppi_DATA = {
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK8PFPuppi.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_EF_L123_AK8PFPuppi_DATA = {
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L1FastJet_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2Relative_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L3Absolute_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK8PFPuppi.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_G_L123_AK8PFPuppi_DATA = {
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L1FastJet_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2Relative_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L3Absolute_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2L3Residual_AK8PFPuppi.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_H_L123_AK8PFPuppi_DATA = {
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L1FastJet_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2Relative_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L3Absolute_AK8PFPuppi.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK8PFPuppi.txt",
};




const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_BCD_L23_AK4PFchs_DATA = {
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_EF_L23_AK4PFchs_DATA = {
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_G_L23_AK4PFchs_DATA = {
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_H_L23_AK4PFchs_DATA = {
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L3Absolute_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_L23_AK4PFchs_MC = {
  "JECDatabase/textFiles/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L2Relative_AK4PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L3Absolute_AK4PFchs.txt",
};


const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_BCD_L23_AK8PFchs_DATA = {
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK8PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK8PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016BCDV4_DATA/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK8PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_EF_L23_AK8PFchs_DATA = {
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2Relative_AK8PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L3Absolute_AK8PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016EFV4_DATA/Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK8PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_G_L23_AK8PFchs_DATA = {
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2Relative_AK8PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L3Absolute_AK8PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016GV4_DATA/Summer16_23Sep2016GV4_DATA_L2L3Residual_AK8PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_H_L23_AK8PFchs_DATA = {
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2Relative_AK8PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L3Absolute_AK8PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016HV4_DATA/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK8PFchs.txt",
};

const std::vector<std::string> JERFiles::Summer16_23Sep2016_V4_L23_AK8PFchs_MC = {
  "JECDatabase/textFiles/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L2Relative_AK8PFchs.txt",
  "JECDatabase/textFiles/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC_L3Absolute_AK8PFchs.txt",
};
