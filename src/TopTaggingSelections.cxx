#include "UHH2/TopTagging/include/TopTaggingSelections.h"
#include "UHH2/TopTagging/include/TopTaggingUtils.h"
#include "UHH2/core/include/Event.h"


#include <stdexcept>

//using namespace uhh2examples;
using namespace std;
using namespace uhh2;


DijetSelection::DijetSelection(float dphi_min_, float third_frac_max_): dphi_min(dphi_min_), third_frac_max(third_frac_max_){}
    
bool DijetSelection::passes(const Event & event){
    assert(event.jets); // if this fails, it probably means jets are not read in
    if(event.jets->size() < 2) return false;
    const auto & jet0 = event.jets->at(0);
    const auto & jet1 = event.jets->at(1);
    auto dphi = deltaPhi(jet0, jet1);
    if(dphi < dphi_min) return false;
    if(event.jets->size() == 2) return true;
    const auto & jet2 = event.jets->at(2);
    auto third_jet_frac = jet2.pt() / (0.5 * (jet0.pt() + jet1.pt()));
    return third_jet_frac < third_frac_max;
}

HTCut::HTCut(float minHT_, float maxHT_):
  minHT(minHT_), maxHT(maxHT_){}

bool HTCut::passes(const Event& event){
  double ht = 0.;
  assert(event.jets);
  for(const auto jet : *event.jets){
    ht += jet.pt();
  }
  if(ht > minHT && ht < maxHT) return true;
  return false;
}


HTlepCut::HTlepCut(double minHTLep_, double maxHTLep_, bool useMuons_, bool useElectrons_):
  minHTLep(minHTLep_), maxHTLep(maxHTLep_), useMuons(useMuons_), useElectrons(useElectrons_){}

bool HTlepCut::passes(const Event& event){
  double ptlep = 0;

  if(useMuons){
    assert(event.muons);
    for(const auto muon : *event.muons){
      if(muon.pt() > ptlep) ptlep = muon.pt();
    }
  }
  if(useElectrons){
    for(const auto ele : *event.electrons){
      if(ele.pt() > ptlep) ptlep = ele.pt();
    }
  }

  double htlep = ptlep + event.met->pt();

  if(htlep > minHTLep && htlep < maxHTLep) return true;
  return false;
}



METCut::METCut(float minMet_, float maxMet_):
  minMet(minMet_), maxMet(maxMet_) {}

bool METCut::passes(const Event & event){

  assert(event.met);

  float MET = event.met->pt();
  return (MET > minMet) && (MET < maxMet);
}


bool TwoDCut::passes(const Event & event){

  assert(event.muons && event.electrons && event.jets);
  /*  if((event.muons->size()+event.electrons->size()) != 1){
    std::cout << "N_elec=" << event.electrons->size() << "N_muon=" << event.muons->size() << "\n @@@ WARNING -- TwoDCut::passes -- unexpected number of muons+electrons in the event (!=1). returning 'false'\n";
    return false;
    }*/

  float drmin, ptrel;  
  if(event.muons->size()) std::tie(drmin, ptrel) = drmin_pTrel(event.muons->at(0), *event.jets);
  else std::tie(drmin, ptrel) = drmin_pTrel(event.electrons->at(0), *event.jets);

  return (drmin > min_deltaR) || (ptrel > min_pTrel);
}



NMuonBTagSelection::NMuonBTagSelection(int min_nbtag, int max_nbtag, JetId btag, double ptmin, double etamax )
{
  m_min_nbtag=min_nbtag;
  m_max_nbtag=max_nbtag;
  m_btag=btag;
  m_ptmin=ptmin;
  m_etamax=etamax;
}

bool NMuonBTagSelection::passes(const Event & event)
{
  int nbtag=0;

  //Assumes to have only one muon                                                                  
  std::vector<Jet>* jets = event.jets;
  std::vector<Muon>* muons = event.muons;
  for(unsigned int i=0; i<event.jets->size(); ++i) {
    int jettagged=0;
    Jet jet=jets->at(i);
    if (m_btag(jet, event)) jettagged=1;
 
    if(muons->size() != 1){
      std::cout << "ATTENTION!!! muon size " << muons->size() << std::endl;
    }

    double deltaphi=deltaPhi(jet,muons->at(0));
    double pi = 3.14159265359;
    if(jettagged&&(deltaphi<(2*pi/3))&&(jet.pt()>m_ptmin)&&(fabs(jet.eta())<m_etamax)){

      nbtag++;

    }
  }

  if(nbtag<m_min_nbtag) return false;
  if(nbtag>m_max_nbtag) return false;
  return true;
}


MergedSelection::MergedSelection( uhh2::Context& ctx, const std::string ttbarGen_name_, double radius_, mergingOpt opt_): ttbarGen_name(ttbarGen_name_), radius(radius_), opt(opt_) {
  h_ttbarGen = ctx.get_handle<TTbarGen>(ttbarGen_name);
 } 

bool MergedSelection::passes(const uhh2::Event &event) {
  return false;
}


bool MergedSelection::passes_probe(const uhh2::Event &event, const TopJet &probe_jet) {

 
  // const TTbarGen& ttbarGen = !ttbarGen_name.empty() ? event.get(h_ttbarGen) : TTbarGen(*event.genparticles,true);
  const auto & ttbarGen = event.get(h_ttbarGen);

  double r = 0.; 
  /* if (radius < 0){
    double Rmin = 0.1;
    double Rmax = 1.5;
    double rho = 600.;
    double reff = probe_jet.pt()/rho;

    if( reff <  Rmin ) r = Rmin;
    else if( reff >  Rmax ) r = Rmax;
    else r = reff;
    //   r = sqrt(probe_jet.jetArea()/3.14159265359); //needed for HOTVR (no fixed radius)
    }*/
  //else 
  r = radius;
  

  if(radius < 0){
    
    //new subjet matching for HOTVR jets

    if(ttbarGen.IsSemiLeptonicDecay()) {
      
      bool q1_matched = false;
      bool q2_matched = false;
      bool b_matched = false;

      GenParticle bHad = ttbarGen.BHad();
      GenParticle q1 = ttbarGen.Q1();
      GenParticle q2 = ttbarGen.Q2();

      for(const auto & subjet: probe_jet.subjets()){
	
	double ri = r = sqrt(subjet.jetArea()/3.14159265359);
	
	if(deltaR(subjet.v4(), q1) < ri) q1_matched = true;
	if(deltaR(subjet.v4(), q2) < ri) q2_matched = true;
	if(deltaR(subjet.v4(), bHad) < ri) b_matched = true;
      }

      int Nmatched = 0;
      if(q1_matched) Nmatched++;
      if(q2_matched) Nmatched++;
      if(b_matched) Nmatched++;

      if(opt == oFullyMerged && Nmatched == 3) return true;
      if(opt == oSemiMerged && Nmatched == 2) return true;
      if(opt == oNotMerged && Nmatched < 2) return true;

      if(opt == oLight && Nmatched < 2) return true;

      if(opt == oMergedW && Nmatched == 2 && b_matched) return true;
      if(opt == oBplusQ && Nmatched == 2 && !b_matched) return true;
      
    }
    else if(opt == oBkg || opt == oNotMerged) {
      return true;
    }    
    
  }else{

    //old matching

    if(ttbarGen.IsSemiLeptonicDecay()) {
      GenParticle bHad = ttbarGen.BHad();
      GenParticle q1 = ttbarGen.Q1();
      GenParticle q2 = ttbarGen.Q2();
      if(opt == oFullyMerged){
	if( deltaR(probe_jet.v4(), bHad.v4()) < r
	    && deltaR(probe_jet.v4(), q1.v4()) < r
	    && deltaR(probe_jet.v4(), q2.v4()) < r) {
	  return true;
	}
      }
      if(opt == oMergedW || opt == oSemiMerged) {
	if( deltaR(probe_jet.v4(), bHad.v4()) > r
	    && deltaR(probe_jet.v4(), q1.v4()) < r
	    && deltaR(probe_jet.v4(), q2.v4()) < r) {
	  return true;
	}
      }
      if(opt == oBplusQ || opt == oSemiMerged) {
	if( deltaR(probe_jet.v4(), bHad.v4()) < r){
	  if (deltaR(probe_jet.v4(), q1.v4()) < r && deltaR(probe_jet.v4(), q2.v4()) > r) return true;
	  if (deltaR(probe_jet.v4(), q1.v4()) > r && deltaR(probe_jet.v4(), q2.v4()) < r) return true;	
	}
      }
      if(opt == oLight|| opt == oNotMerged) {
	unsigned int N = 0;
	if( deltaR(probe_jet.v4(), bHad.v4()) < r ) N++;
	if( deltaR(probe_jet.v4(), q1.v4()) < r) N++;
	if( deltaR(probe_jet.v4(), q2.v4()) < r) N++;
	if( N <= 1) return true;
      }
    }
    else if(opt == oBkg || opt == oNotMerged) {
      return true;
    }
  }
  return false;
}


DecayChannelSelection::DecayChannelSelection( uhh2::Context& ctx,  const std::string ttbarGen_name_, TString channel_): ttbarGen_name(ttbarGen_name_), channel(channel_) {
  h_ttbarGen = ctx.get_handle<TTbarGen>(ttbarGen_name);
} 

bool DecayChannelSelection::passes(const uhh2::Event &event) {
  
  const TTbarGen& ttbarGen = !ttbarGen_name.empty() ? event.get(h_ttbarGen) : TTbarGen(*event.genparticles,false);
  
  if( ttbarGen.DecayChannel() != TTbarGen::e_notfound) {
    if(channel == "dilepton") {
      if (ttbarGen.DecayChannel() == TTbarGen::e_mumu || 
	  ttbarGen.DecayChannel() == TTbarGen::e_ee || 
	  ttbarGen.DecayChannel() == TTbarGen::e_tautau ||
	  ttbarGen.DecayChannel() == TTbarGen::e_emu ||
	  ttbarGen.DecayChannel() == TTbarGen::e_etau ||
	  ttbarGen.DecayChannel() == TTbarGen::e_mutau ) { 
	return true;
      }
    }
    else if(channel == "lepton_jets") {
      if (ttbarGen.IsSemiLeptonicDecay() &&
	  !(ttbarGen.DecayChannel() == TTbarGen::e_tauhad) ) { 
	return true;
      }
    }
    else if(channel == "hadronic") {
      if (ttbarGen.DecayChannel() == TTbarGen::e_had) {
	return true;
      } 
    }
    else if(channel == "tau_jets") {
      if (ttbarGen.DecayChannel() == TTbarGen::e_tauhad) {
	return true;
      } 
    }
    else{
      std::cout << "NO proper decay channel selected: all events will be rejected by the decay channel selection" << std::endl;
    }
  }
  return false;
}

MassDiffSelection::MassDiffSelection(JetId btag):  btag_(btag) {}

bool MassDiffSelection::passes(const uhh2::Event &event){
  std::cout << "passes not used. reject all events!" << std::endl;
  return false;
}

bool MassDiffSelection::passes_probe(const uhh2::Event &event, const TopJet &probeJet){
  //get the muon
  Muon mu = event.muons->at(0);
  
  //get the bjet
  Jet bjet;
 
  bool bjet_found = false; 
  bool b_candidate_found = false;
  
  double max_pt = 0.;
  std::vector<Jet> *ak4jets = event.jets;
  double pi = 3.14159265359;

  for( const auto & ak4jet : *ak4jets){
    Jet b_candidate;
    if( btag_(ak4jet, event) && (deltaPhi(ak4jet,mu) < (2*pi/3)) ){
      b_candidate = ak4jet;
      b_candidate_found = true; 
    }
    if(b_candidate_found){
      if( b_candidate.pt() > max_pt){
	bjet = b_candidate;
	max_pt = b_candidate.pt();
	bjet_found = true; 
      }
    }
  } 

  if(!b_candidate_found) std::cout << "No Bjet candidate found in event:  " << event.event << std::endl;
  if(!bjet_found) std::cout << "No Bjet found in event: " << event.event << std::endl;

  double mProbe_mLep_bjet = probeJet.v4().M()/(bjet.v4()+mu.v4()).M();

  if(bjet_found && mProbe_mLep_bjet > 1.) return true;
  return false;

}


DPhiMuBSelection::DPhiMuBSelection( JetId btag, double dPhiMin):  btag_(btag), dPhiMin_(dPhiMin) {}

bool DPhiMuBSelection::passes(const uhh2::Event &event){
  std::cout << "passes not used. reject all events!" << std::endl;
  return false;
}

bool DPhiMuBSelection::passes_probe(const uhh2::Event &event, const TopJet &probeJet){
  //get the muon
  Muon mu = event.muons->at(0);
  
  //get the bjet
  Jet bjet;
 
  bool bjet_found = false; 
  bool b_candidate_found = false;
  
  double max_pt = 0.;
  std::vector<Jet> *ak4jets = event.jets;
  double pi = 3.14159265359;

  for( const auto & ak4jet : *ak4jets){
    Jet b_candidate;
    if( btag_(ak4jet, event) && (deltaPhi(ak4jet,mu) < (2*pi/3)) ){
      b_candidate = ak4jet;
      b_candidate_found = true; 
    }
    if(b_candidate_found){
      if( b_candidate.pt() > max_pt){
	bjet = b_candidate;
	max_pt = b_candidate.pt();
	bjet_found = true; 
      }
    }
  } 

  if(!b_candidate_found) std::cout << "No Bjet candidate found in event:  " << event.event << std::endl;
  if(!bjet_found) std::cout << "No Bjet found in event: " << event.event << std::endl;

  double dPhi = deltaPhi(probeJet, mu);

  if(bjet_found && dPhi >  dPhiMin_) return true;
  return false;

}

LeadingAddJetSelection::LeadingAddJetSelection( JetId btag, double ptMin): btag_(btag), ptMin_(ptMin) {}

bool LeadingAddJetSelection::passes(const uhh2::Event &event){
  std::cout << "passes not used. reject all events! please use passes_probe instead" << std::endl;
  return false;
}

bool LeadingAddJetSelection::passes_probe(const uhh2::Event &event, const TopJet &probeJet){
  Jet bjet;

  bool bjet_found = GetLeadingBjetLepHem(event, bjet, btag_ );
  if(!bjet_found) return false;

  for( const auto & jet : *event.jets){
    if(deltaR(jet, probeJet) > 0.8 && deltaR(jet, bjet) > 0.1){
      if(jet.pt() > ptMin_) return true;
    }
  }
  return false;
}

/*
HadronicTopSelection::HadronicTopSelection( uhh2::Context& ctx, const std::string ttbarGen_name_): ttbarGen_name(ttbarGen_name_) {
  h_ttbarGen = ctx.get_handle<TTbarGen>(ttbarGen_name);
 } 

bool HadronicTopSelection::passes(const uhh2::Event &event) {
  return false;
}
bool HadronicTopSelection::passes_jet(const uhh2::Event &event, const TopJet &jet) {
 
  // const TTbarGen& ttbarGen = !ttbarGen_name.empty() ? event.get(h_ttbarGen) : TTbarGen(*event.genparticles,true);
  const auto & ttbarGen = event.get(h_ttbarGen);

  vector<GenParticle> tops; 
  double radius = 0.6;

  if(IsTopHadronicDecay()) {
    GenParticle b = ttbarGen.bTop();
    GenParticle q1 = ttbarGen.Wdecay1();
    GenParticle q2 = ttbarGen.Wdecay2();
    if( deltaR(b.v4(),q1.v4()) < radius
	&& deltaR(b.v4(), q2.v4()) < radius
	&& deltaR(q1.v4(), q2.v4()) < radius) {
      tops.push_back(ttbarGen.Top());
    }
  }
  if(IsAntiTopHadronicDecay()) {
    GenParticle b = ttbarGen.bAntitop();
    GenParticle q1 = ttbarGen.WMinusdecay1();
    GenParticle q2 = ttbarGen.WMinusdecay2();
    if( deltaR(b.v4(),q1.v4()) < radius
	&& deltaR(b.v4(), q2.v4()) < radius
	&& deltaR(q1.v4(), q2.v4()) < radius) {
      tops.push_back(ttbarGen.Antitop());
    }
  }

  for( const auto & top : tops){
    if(deltaR(jet.v4(),top.v4()) < radius) return true;
  }

  return false;
}


vector<GenParticles> GetMergedHadronicTops(TTbarGen ttbarGen, 
*/
    
//};
