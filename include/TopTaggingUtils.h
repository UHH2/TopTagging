#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/JetCorrections.h"
//#include "UHH2/common/include/Utils.h"
//#include <UHH2/core/include/Utils.h>

class TopJetLeptonDeltaRCleaner : public uhh2::AnalysisModule {
 public:
  explicit TopJetLeptonDeltaRCleaner(float mindr=0.8): minDR_(mindr) {}
  virtual bool process(uhh2::Event&) override;

 private:
  float minDR_;
};


class TopJetCorrectionModules : public uhh2::AnalysisModule {
 public:

  enum jet_type{AK8_PUPPI, AK8_CHS, CA15_PUPPI, CA15_CHS, HOTVR_CHS, HOTVR_PUPPI};

  explicit TopJetCorrectionModules(uhh2::Context & ctx, jet_type type = AK8_CHS);
  virtual bool process(uhh2::Event & event) override;



 private:
 bool is_mc;
 std::unique_ptr<uhh2::AnalysisModule> topjet_corrector_MC, topjet_corrector_BCD, topjet_corrector_EF, topjet_corrector_FG, topjet_corrector_H;
 std::unique_ptr<SubJetCorrector> subjet_corrector_MC, subjet_corrector_BCD, subjet_corrector_EF, subjet_corrector_FG, subjet_corrector_H;

 jet_type _type;

 const int runnr_BCD = 276811;
 const int runnr_EFearly = 278802;
 const int runnr_FlateG = 280385;

};


class Muon_dxydzCut {
public:
 Muon_dxydzCut(float max_dxy_ , float max_dz_): max_dxy(max_dxy_), max_dz(max_dz_){}

    bool operator()(const Muon & muon, const uhh2::Event & ) const{
      return fabs(muon.dxy()) < max_dxy && fabs(muon.dxy()) < max_dz;
    }
    
private:
    float max_dxy, max_dz;
};


class HOTVRPileupCorrectionModule : public uhh2::AnalysisModule {
 public:
  explicit HOTVRPileupCorrectionModule(bool area_correction = true): _area_correction(area_correction) {}
  virtual bool process(uhh2::Event & event) override;

 private:
  bool _area_correction; 
};


class TopJetGroomer : public uhh2::AnalysisModule {
 public:

  explicit TopJetGroomer(bool corrected = true) : _corrected(corrected) {}
  virtual bool process(uhh2::Event & event) override;
 
 private:
  bool _corrected;

};

bool GetLeadingBjetLepHem(const uhh2::Event &event, Jet &bjet, JetId btag);

namespace JERFiles {
 
  extern const std::vector<std::string> Summer16_23Sep2016_V4_BCD_L123_AK4PFPuppi_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_EF_L123_AK4PFPuppi_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_G_L123_AK4PFPuppi_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_H_L123_AK4PFPuppi_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_BCD_L123_AK8PFPuppi_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_EF_L123_AK8PFPuppi_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_G_L123_AK8PFPuppi_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_H_L123_AK8PFPuppi_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_L123_AK4PFPuppi_MC;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_L123_AK8PFPuppi_MC;

  extern const std::vector<std::string> Summer16_23Sep2016_V4_BCD_L23_AK4PFchs_DATA;                                 
  extern const std::vector<std::string> Summer16_23Sep2016_V4_EF_L23_AK4PFchs_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_G_L23_AK4PFchs_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_H_L23_AK4PFchs_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_L23_AK4PFchs_MC;

  extern const std::vector<std::string> Summer16_23Sep2016_V4_BCD_L23_AK8PFchs_DATA;                                 
  extern const std::vector<std::string> Summer16_23Sep2016_V4_EF_L23_AK8PFchs_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_G_L23_AK8PFchs_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_H_L23_AK8PFchs_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_L23_AK8PFchs_MC;


  extern const std::vector<std::string> Summer16_23Sep2016_V4_BCD_L23_AK4PFPuppi_DATA;                                 
  extern const std::vector<std::string> Summer16_23Sep2016_V4_EF_L23_AK4PFPuppi_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_G_L23_AK4PFPuppi_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_H_L23_AK4PFPuppi_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_L23_AK4PFPuppi_MC;

  extern const std::vector<std::string> Summer16_23Sep2016_V4_BCD_L23_AK8PFPuppi_DATA;                                 
  extern const std::vector<std::string> Summer16_23Sep2016_V4_EF_L23_AK8PFPuppi_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_G_L23_AK8PFPuppi_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_H_L23_AK8PFPuppi_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_L23_AK8PFPuppi_MC;


  /*
  extern const std::vector<std::string> Summer16_23Sep2016_V4_BCD_L1RC_AK4PFchs_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_EF_L1RC_AK4PFchs_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_G_L1RC_AK4PFchs_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_H_L1RC_AK4PFchs_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_L1RC_AK4PFchs_MC; 
  */
}

