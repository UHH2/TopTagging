#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include <UHH2/core/include/Utils.h>
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/Utils.h"


#include "UHH2/common/include/ReconstructionHypothesisDiscriminators.h"



namespace uhh2 {
    
/* Select events with at least two jets in which the leading two jets have deltaphi > 2.7 and the third jet pt is
 * below 20% of the average of the leading two jets, where the minimum deltaphi and
 * maximum third jet pt fraction can be changed in the constructor.
 * The jets are assumed to be sorted in pt.
 */
class DijetSelection: public uhh2::Selection {
public:
    DijetSelection(float dphi_min = 2.7f, float third_frac_max = 0.2f);
    virtual bool passes(const uhh2::Event & event) override;
private:
    float dphi_min, third_frac_max;
};


class HTCut: public uhh2::Selection {
 public:
  explicit HTCut(float minHT_, float maxHT_ = uhh2::infinity);
  virtual bool passes(const uhh2::Event&) override;

 private:
  float minHT, maxHT;

};

class PtWSelection: public uhh2::Selection {
 public:
  explicit PtWSelection(double minPt_, double maxPt_ = uhh2::infinity);
  virtual bool passes(const uhh2::Event&) override;

 private:
  double minPt, maxPt;

};

class HTlepCut: public uhh2::Selection {
  public:
    explicit HTlepCut(double minHTLep_, double maxHTLep_, bool useMuons_ = true, bool useElectrons_ = true);
    virtual bool passes(const uhh2::Event&) override; 

  private:
    double minHTLep, maxHTLep; 
    bool useMuons, useElectrons;
};

  class METCut : public uhh2::Selection {
   public:
    explicit METCut(float minMet_, float maxMet_=infinity);
    virtual bool passes(const uhh2::Event & event) override;

   private:
    float minMet, maxMet;
};



  class TwoDCut : public uhh2::Selection {
   public:
    explicit TwoDCut(float min_deltaR_, float min_pTrel_): min_deltaR(min_deltaR_), min_pTrel(min_pTrel_) {}
    virtual bool passes(const uhh2::Event & event) override;

   private:
    float min_deltaR, min_pTrel;
};


class NMuonBTagSelection: public uhh2::Selection {
  public:

    explicit NMuonBTagSelection(int min_nbtag, int max_nbtag=999, JetId btag=CSVBTag(CSVBTag::WP_LOOSE), double ptmin=0., double etamax=infinity );
 
    virtual bool passes(const uhh2::Event &) override;

  private:
    int m_min_nbtag;
    int m_max_nbtag;
    JetId m_btag;
    double m_ptmin;
    double m_etamax;
};


class MergedSelection: public uhh2::Selection {
  public:

    enum mergingOpt{oFullyMerged, oMergedW, oBplusQ, oLight, oBkg, oNotMerged, oSemiMerged};

    explicit MergedSelection( uhh2::Context& ctx,  const std::string ttbarGen_name_, double radius_ , mergingOpt opt_ = oFullyMerged);
    virtual bool passes(const uhh2::Event &) override;
    bool passes_probe(const uhh2::Event &event, const TopJet &probe_jet);

  private:
    uhh2::Event::Handle<TTbarGen> h_ttbarGen;
     const std::string ttbarGen_name;
    double radius;
    mergingOpt opt;
};


class DecayChannelSelection: public uhh2::Selection {
  public:

    explicit DecayChannelSelection( uhh2::Context& ctx,  const std::string ttbarGen_name_, TString channel_ );
    virtual bool passes(const uhh2::Event &) override;

  private:
    uhh2::Event::Handle<TTbarGen> h_ttbarGen;
     const std::string ttbarGen_name;
    TString channel;
};

class MassDiffSelection: public uhh2::Selection{
  public:
    explicit MassDiffSelection( JetId btag=CSVBTag(CSVBTag::WP_MEDIUM));
    virtual bool passes(const uhh2::Event &) override;
    bool passes_probe(const uhh2::Event &event, const TopJet &probeJet);

  private:
    JetId btag_;
};
  

class DPhiMuBSelection: public uhh2::Selection{
  public:
    explicit DPhiMuBSelection( JetId btag=CSVBTag(CSVBTag::WP_MEDIUM), double dPhiMin = 0.);
    virtual bool passes(const uhh2::Event &) override;
    bool passes_probe(const uhh2::Event &event, const TopJet &probeJet);

  private:
    JetId btag_;
    double dPhiMin_;
};


class LeadingAddJetSelection: public uhh2::Selection{
  public:

    explicit LeadingAddJetSelection( JetId btag, double ptMin);
    virtual bool passes(const uhh2::Event &) override;
    bool passes_probe(const uhh2::Event &event, const TopJet &probeJet);

  private:
     JetId btag_;
     double ptMin_;
};

  /*
class HadronicTopSelection: pubic uhh2::Selection{
  public:

    explicit HadronicTopSelection( uhh2::Context& ctx, const std::string ttbarGen_name_);
    virtual bool passes(const uhh2::Event &) override;
    boll passes_jet(const uhh2::Event &event, const TopJet &jet);

  privte:
    uhh2::Event::Handle<TTbarGen> h_ttbarGen;

};
  */

}



