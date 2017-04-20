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

  topjet_corrector_MC.reset(new TopJetCorrector(ctx, JEC_TOP_MC));
  subjet_corrector_MC.reset(new SubJetCorrector(ctx, JEC_SUB_MC));

  topjet_corrector_BCD.reset(new TopJetCorrector(ctx, JEC_TOP_BCD));
  topjet_corrector_EF.reset(new TopJetCorrector(ctx, JEC_TOP_EF));
  topjet_corrector_FG.reset(new TopJetCorrector(ctx, JEC_TOP_FG));
  topjet_corrector_H.reset(new TopJetCorrector(ctx, JEC_TOP_H));

  subjet_corrector_BCD.reset(new SubJetCorrector(ctx, JEC_SUB_BCD));
  subjet_corrector_EF.reset(new SubJetCorrector(ctx, JEC_SUB_EF));
  subjet_corrector_FG.reset(new SubJetCorrector(ctx, JEC_SUB_FG));
  subjet_corrector_H.reset(new SubJetCorrector(ctx, JEC_SUB_H));
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

bool TopJetGroomer::process(uhh2::Event & event) { 

  assert(event.topjets);

  for(auto & topjet : *event.topjets){
    LorentzVector subjet_sum;

    for(const auto & subjet : topjet.subjets()){
      if(_corrected) subjet_sum += subjet.v4();
      else subjet_sum += subjet.v4()*subjet.JEC_factor_raw();
    }
    topjet.set_v4(subjet_sum);
    if(!_corrected) topjet.set_JEC_factor_raw(1.);
  }

  return true;
}



