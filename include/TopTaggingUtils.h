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

  enum jet_type{AK8_PUPPI, AK8_CHS, CA15_PUPPI, CA15_CHS, HOTVR_CHS};

  explicit TopJetCorrectionModules(uhh2::Context & ctx, jet_type type = AK8_CHS);
  virtual bool process(uhh2::Event & event) override;



 private:
 bool is_mc;
 std::unique_ptr<TopJetCorrector> topjet_corrector_MC, topjet_corrector_B, topjet_corrector_C, topjet_corrector_D, topjet_corrector_E, topjet_corrector_F;
 std::unique_ptr<SubJetCorrector> subjet_corrector_MC, subjet_corrector_B, subjet_corrector_C, subjet_corrector_D, subjet_corrector_E, subjet_corrector_F;

 jet_type _type;

 // const int runnr_BCD = 276811;
 // const int runnr_EFearly = 278802;
 // const int runnr_FlateG = 280385;
 const int runnr_B = 299329; 
 const int runnr_C = 302029; 
 const int runnr_D = 303434; 
 const int runnr_E = 304797; 


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

class PartonShowerWeight : public uhh2::AnalysisModule {
 public:

  explicit PartonShowerWeight(uhh2::Context & ctx, std::string sys);
  virtual bool process(uhh2::Event & event) override;

 private:
  int weightIndex;

};

bool GetLeadingBjetLepHem(const uhh2::Event &event, Jet &bjet, JetId btag);

namespace JERFiles {

  /*
  extern const std::vector<std::string> JERFiles::Fall17_17Nov2017_V6_B_L123_AK4PFchs_DATA;
  extern const std::vector<std::string> JERFiles::Fall17_17Nov2017_V6_C_L123_AK4PFchs_DATA;
  extern const std::vector<std::string> JERFiles::Fall17_17Nov2017_V6_D_L123_AK4PFchs_DATA;
  extern const std::vector<std::string> JERFiles::Fall17_17Nov2017_V6_E_L123_AK4PFchs_DATA;
  extern const std::vector<std::string> JERFiles::Fall17_17Nov2017_V6_F_L123_AK4PFchs_DATA;
  */
  //AK4 Puppi
  extern const std::vector<std::string> Fall17_17Nov2017_V6_B_L123_AK4PFPuppi_DATA;
  extern const std::vector<std::string> Fall17_17Nov2017_V6_C_L123_AK4PFPuppi_DATA;
  extern const std::vector<std::string> Fall17_17Nov2017_V6_D_L123_AK4PFPuppi_DATA;
  extern const std::vector<std::string> Fall17_17Nov2017_V6_E_L123_AK4PFPuppi_DATA;
  extern const std::vector<std::string> Fall17_17Nov2017_V6_F_L123_AK4PFPuppi_DATA;
  
  //AK8
  extern const std::vector<std::string> Fall17_17Nov2017_V6_B_L123_AK8PFchs_DATA;
  extern const std::vector<std::string> Fall17_17Nov2017_V6_C_L123_AK8PFchs_DATA;
  extern const std::vector<std::string> Fall17_17Nov2017_V6_D_L123_AK8PFchs_DATA;
  extern const std::vector<std::string> Fall17_17Nov2017_V6_E_L123_AK8PFchs_DATA;
  extern const std::vector<std::string> Fall17_17Nov2017_V6_F_L123_AK8PFchs_DATA;
 
  //AK8 Puppi
  extern const std::vector<std::string> Fall17_17Nov2017_V6_B_L123_AK8PFPuppi_DATA;
  extern const std::vector<std::string> Fall17_17Nov2017_V6_C_L123_AK8PFPuppi_DATA;
  extern const std::vector<std::string> Fall17_17Nov2017_V6_D_L123_AK8PFPuppi_DATA;
  extern const std::vector<std::string> Fall17_17Nov2017_V6_E_L123_AK8PFPuppi_DATA;
  extern const std::vector<std::string> Fall17_17Nov2017_V6_F_L123_AK8PFPuppi_DATA;

  extern const std::vector<std::string> Fall17_17Nov2017_V6_L123_AK4PFPuppi_MC;
  extern const std::vector<std::string> Fall17_17Nov2017_V6_L123_AK8PFchs_MC;
  extern const std::vector<std::string> Fall17_17Nov2017_V6_L123_AK8PFPuppi_MC;

  /*
  extern const std::vector<std::string> Summer16_23Sep2016_V4_BCD_L1RC_AK4PFchs_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_EF_L1RC_AK4PFchs_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_G_L1RC_AK4PFchs_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_H_L1RC_AK4PFchs_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_L1RC_AK4PFchs_MC; 
  */
}



