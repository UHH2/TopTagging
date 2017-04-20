#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/JetCorrections.h"

class TopJetLeptonDeltaRCleaner : public uhh2::AnalysisModule {
 public:
  explicit TopJetLeptonDeltaRCleaner(float mindr=0.8): minDR_(mindr) {}
  virtual bool process(uhh2::Event&) override;

 private:
  float minDR_;
};


class TopJetCorrectionModules : public uhh2::AnalysisModule {
 public:

 enum jet_type{AK8_PUPPI, AK8_CHS, CA15_PUPPI, CA15_CHS};

  explicit TopJetCorrectionModules(uhh2::Context & ctx, jet_type type = AK8_CHS);
  virtual bool process(uhh2::Event & event) override;
  // void init(uhh2::Context & ctx, jet_type type);


 private:
 bool is_mc;
 std::unique_ptr<TopJetCorrector> topjet_corrector_MC, topjet_corrector_BCD, topjet_corrector_EF, topjet_corrector_FG, topjet_corrector_H;
 std::unique_ptr<SubJetCorrector> subjet_corrector_MC, subjet_corrector_BCD, subjet_corrector_EF, subjet_corrector_FG, subjet_corrector_H;

 jet_type _type;

 const int runnr_BCD = 276811;
 const int runnr_EFearly = 278802;
 const int runnr_FlateG = 280385;

};

class TopJetGroomer : public uhh2::AnalysisModule {
 public:

  explicit TopJetGroomer(bool corrected = true) : _corrected(corrected) {}
  virtual bool process(uhh2::Event & event) override;
 
 private:
  bool _corrected;

};


