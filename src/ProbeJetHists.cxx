#include "UHH2/TopTagging/include/ProbeJetHists.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/Utils.h"

#include "TH1F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

ProbeJetHists::ProbeJetHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  fill_PDF = (ctx.get("fill_PDF", "FALSE") == "TRUE");
  if(!(ctx.get("dataset_type") == "MC")) fill_PDF = false;

  mass_scale = "default";
  mass_scale = ctx.get("mass_scale","default");
  
  book<TH1F>("npv", "N_{primary vertices}",  56, -0.5, 55.5);

 //probe jet properties
  book<TH1F>("pt", "Probe jet p_{T} [GeV]", 200, 0., 2000);
  book<TH1F>("eta", "Probe jet #eta", 32, -3.2, 3.2);

  book<TH1F>("pt_subjets", "subjet p_{T} [GeV]", 100,0.,1000.);

  book<TH1F>("mass", "Probe jet mass [GeV]", 50, 0, 500);
  book<TH1F>("mass_sub", "Soft drop mass [GeV]", 100, 0, 500);

  book<TH1F>("mass_SD", "Probe jet soft drop mass raw [GeV]", 50, 0, 500);
  book<TH1F>("mass_SD_Corr", "Probe jet soft drop mass corrected [GeV]", 40, 0, 300);
  book<TH1F>("mass_pruned", "Probe jet pruned mass [GeV]", 50, 0, 500);

  book<TH1F>("tau32", "Probe jet #tau3/#tau2", 40, 0, 1);
  book<TH1F>("tau21", "Probe jet #tau2/#tau1", 40, 0, 1);
  book<TH1F>("tau2", "Probe jet #tau2", 40, 0, 1);
  book<TH1F>("tau3", "Probe jet #tau3", 40, 0, 1);

  book<TH1F>("HTTmass", "Probe jet HTT mass [GeV]", 40, 0, 300);
  book<TH1F>("HTTmass_corr", "Probe jet HTT mass corrrected [GeV]", 40, 0, 300);

  book<TH1F>("fRec", "Probe jet fRec", 40, 0, 1);

  book<TH1F>("m12", "m12", 20, 0., 100.);
  book<TH1F>("m13", "m13", 20, 0., 100.);
  book<TH1F>("m23", "m23", 20, 0., 100.);

  book<TH1F>("Nsub", "Number of subjets", 6, -0.5, 5.5);  

  book<TH1F>("subCSV", "subjets CSV discriminator", 100, 0., 1.);  
  book<TH1F>("subCSV_highest", "highest subjet CSV discriminator", 100, 0., 1.); 
  book<TH1F>("CSV_topjet", "CSV discriminator Top Jet", 100, 0., 1.);

  book<TH1F>("subjetPT", "PT of subjet with highest CSV", 100, 0., 1500.);
  book<TH1F>("subjethadronFlavor", "hadronFlavour of subjet with highest CSV", 100, 0., 1.);

  book<TH1F>("Wmass", "W mass [GeV]", 50, 0., 150.);

  if(fill_PDF){
    for (unsigned i=0; i < 103; ++i) {
      book<TH1F>("mass_sub_PDF_"+std::to_string(i), "Soft drop mass [GeV]", 100, 0, 500);
    }
  }
    
  //additional variables

  //bjet close to muon
  book<TH1F>("mProbe_mLep_bjet1", "m_{probe jet}/m_{b-jet+#mu}", 100, 0., 10.);
  book<TH1F>("mProbe_mLep_bjet_pt", "m_{probe jet}/m_{b-jet+#mu} (leading b-jet)", 100, 0., 10.);
  book<TH1F>("mProbe_mLep_bjet_CSV", "m_{probe jet}/m_{b-jet+#mu} (b-jet with highest CSV)", 100, 0., 10.);
  book<TH1F>("ptProbe_ptLep_bjet1", "p_{T, probe jet}/p_{T, b-jet+#mu}", 100, 0., 10.);
  book<TH1F>("dPhi_mu_b1", "#Delta#Phi(#mu, b-jet)", 50, -1, 4);
  book<TH1F>("pt_leading_add1", "p_{T, leading additional jet} [GeV]" , 30, 0,300);
  book<TH1F>("pt_leading_add_lep1", "p_{T, leading additional jet (lept. hem.)} [GeV]" , 30, 0,300);
  book<TH1F>("pt_leading_addB", "p_{T, leading additional b-jet} [GeV]" , 30, 0,300);
  book<TH1F>("pt_leading_addNoB", "p_{T, leading additional non b-jet} [GeV]" , 30, 0,300);

  //more subjet variables
  book<TH1F>("mSubMax", "m_{subjet} (maximum) [GeV]" , 40, 0, 200);
  book<TH1F>("mSub_ratio", "m_{subjet1}/m_{subjet2} (mass sorted)", 40, 0,40);
  book<TH1F>("rel_mSub_ratio", "#frac{m_{subjet1}}{m_{subjet1}+m_{subjet2}}(mass sorted)", 50, 0, 1); 
  book<TH1F>("ptSub_ratio", "p_{T, subjet1}/p_{T, subjet2}", 40, 0, 10);
  book<TH1F>("rel_ptSub_ratio", "#frac{p_{T, subjet1}}{p_{T, subjet1}+p_{T, subjet2}}", 30, 0.4, 1); 
  
  //scatter plots
  h_ratio_mLB_vs_mB = book<TH2D>("ratio_mLB_vs_mB", ";m_{probe jet}/m_{b-jet+#mu};m_{probe jet}/m_{b-jet}",  100, 0., 10., 1000, 0., 100.);
  h_tau32_vs_pt = book<TH2D>("tau32_vs_pt", ";p_{T} [GeV];#tau_{3}/#tau_{2}",200, 0., 2000,  40, 0, 1);
}

void ProbeJetHists::fill(const Event & event){
}

 void ProbeJetHists::fill_probe(const Event & event, const TopJet & jet){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  
  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  hist("npv")->Fill(event.pvs->size(), weight);
  hist("pt")->Fill(jet.pt(), weight);
  hist("eta")->Fill(jet.eta(), weight);
  
  auto subjets = jet.subjets();
  sort_by_pt(subjets);

  hist("Nsub")->Fill(subjets.size(), weight);

  double m12 = 0; 
  double m13 = 0; 
  double m23 = 0;

  if(subjets.size() == 3) {
    if( (subjets.at(0).v4() + subjets.at(1).v4()).isTimelike() ) {
      m12 = (subjets.at(0).v4() + subjets.at(1).v4()).M();
    }
    if( (subjets.at(0).v4() + subjets.at(2).v4()).isTimelike() ) {
      m13 = (subjets.at(0).v4() + subjets.at(2).v4()).M();
    }
    if( (subjets.at(1).v4() + subjets.at(2).v4()).isTimelike() ) {
      m23 = (subjets.at(1).v4() + subjets.at(2).v4()).M();
    }
  
    hist("m12")->Fill(m12, weight);
    hist("m13")->Fill(m13, weight);
    hist("m23")->Fill(m23, weight);
  }

  LorentzVector subjet_sum(0,0,0,0);
  for (const auto s : subjets) {
    hist("pt_subjets")->Fill(s.pt(), weight);
    subjet_sum += s.v4();
  }

  hist("mass")->Fill(jet.v4().M(), weight);
  //if(subjet_sum.isTimelike()) 

  //hardcoded because of limited time (SAME VALUES IN PostSelectionModule!! ALWAYS CHANGE BOTH!)
  double softdorpmass = subjet_sum.M();
  // if(mass_scale == "up") softdorpmass *= 1.01;
  //if(mass_scale == "down") softdorpmass *= 0.99;

  hist("mass_sub")->Fill(softdorpmass, weight);

  //fill pdf hists
  if(fill_PDF){
    const auto & sys_weights = event.genInfo->systweights();
    float orig_weight = event.genInfo->originalXWGTUP();
    int MY_FIRST_INDEX = 9;
    for (unsigned i=0; i < 103; ++i) {
      TString name = "mass_sub_PDF_"+std::to_string(i);
      hist(name)->Fill(softdorpmass, weight * sys_weights[i+MY_FIRST_INDEX]/orig_weight);
    }
  }

    // else  hist("mass_sub")->Fill(-1, weight);
  hist("mass_SD")->Fill(jet.softdropmass(), weight);
  hist("mass_SD_Corr")->Fill(subjet_sum.M()/jet.JEC_factor_raw(), weight);
  hist("mass_pruned")->Fill(jet.prunedmass(), weight);

  if (jet.has_tag(jet.tagname2tag("fRec"))) hist("fRec")->Fill(jet.get_tag(jet.tagname2tag("fRec")), weight);
  if (jet.has_tag(jet.tagname2tag("mass"))) hist("HTTmass")->Fill(jet.get_tag(jet.tagname2tag("mass")), weight);
  if (jet.has_tag(jet.tagname2tag("mass"))) hist("HTTmass_corr")->Fill(jet.get_tag(jet.tagname2tag("mass"))/jet.JEC_factor_raw(), weight);

  hist("tau32")->Fill(jet.tau3()/jet.tau2(), weight);
  hist("tau21")->Fill(jet.tau2()/jet.tau1(), weight);
  hist("tau3")->Fill(jet.tau3(), weight);
  hist("tau2")->Fill(jet.tau2(), weight);



  hist("CSV_topjet")->Fill(jet.btag_combinedSecondaryVertex(), weight);

  double highestCSV = 0;
  double pt_highestCSV = 0;

  double highestSubjetMass = 0.;
  double lowestSubjetMass = 1000.;

  for (const auto subjet : subjets) {
    hist("subCSV")->Fill(subjet.btag_combinedSecondaryVertex(), weight);
    // hist("subCSV_minus", "CSV discriminator", 100, -12., 12.); 
    if(subjet.btag_combinedSecondaryVertex() > highestCSV){
      highestCSV = subjet.btag_combinedSecondaryVertex();
      pt_highestCSV = subjet.pt();
    }
    if(subjet.v4().M() > highestSubjetMass) highestSubjetMass = subjet.v4().M();
    if(subjet.v4().M() < lowestSubjetMass) lowestSubjetMass = subjet.v4().M();
  }
  hist("subCSV_highest")->Fill(highestCSV, weight);
  hist("subjetPT")->Fill(pt_highestCSV, weight);

  
  hist("mSubMax")->Fill(highestSubjetMass, weight);
  hist("mSub_ratio")->Fill(highestSubjetMass/lowestSubjetMass, weight);
  hist("rel_mSub_ratio")->Fill(highestSubjetMass/subjet_sum.M(), weight);
  if(subjets.size()>1) hist("ptSub_ratio")->Fill(subjets.at(0).pt()/subjets.at(1).pt(), weight);
  if(subjets.size()>1) hist("rel_ptSub_ratio")->Fill(subjets.at(0).pt()/subjet_sum.Pt(), weight);

  JetId checkbtag=CSVBTag(CSVBTag::WP_LOOSE);
  LorentzVector v4_W(0,0,0,0);
  int N_btags = 0;
  if(subjets.size() == 3) {
    for (const auto subjet : subjets) {
      if( checkbtag(subjet, event) ) N_btags++;
      if( !(subjet.btag_combinedSecondaryVertex() == highestCSV) ) v4_W += subjet.v4();
    }
    if(N_btags >= 1) hist("Wmass")->Fill(v4_W.M(), weight);
  }

  h_tau32_vs_pt->Fill(jet.pt(),jet.tau3()/jet.tau2(),weight);

  //======================
  //additoinal hists 
  //======================

  //get the muon
  Muon mu = event.muons->at(0);
  
  //get the bjet
  Jet bjet;
  Jet bjet_max_pt;
  Jet bjet_max_CSV;

  bool bjet_found = false; 
  bool bjet_max_pt_found = false; 
  bool bjet_max_CSV_found = false; 

  bool b_candidate_found = false;
  double min_dphi = 4.;
  double max_pt = 0.;
  double max_CSV = 0.;

  std::vector<Jet> *ak4jets = event.jets;
  //JetId btag = CSVBTag(CSVBTag::WP_MEDIUM);
  JetId btag = DeepCSVBTag(DeepCSVBTag::WP_MEDIUM);
  double pi = 3.14159265359;

  for( const auto & ak4jet : *ak4jets){
    Jet b_candidate;
    if( btag(ak4jet, event) && (deltaPhi(ak4jet,mu) < (2*pi/3)) ){
      b_candidate = ak4jet;
      b_candidate_found = true; 
    }
    if(b_candidate_found){
      if(deltaPhi(b_candidate,mu) < min_dphi ){
	bjet = b_candidate;
	min_dphi = deltaPhi(b_candidate,mu);
	bjet_found = true; 
      }
      if( b_candidate.pt() > max_pt){
	bjet_max_pt = b_candidate;
	max_pt = b_candidate.pt();
	bjet_max_pt_found = true; 
      }
      if( b_candidate.btag_combinedSecondaryVertex() > max_CSV){
	bjet_max_CSV = b_candidate;
	max_CSV = b_candidate.btag_combinedSecondaryVertex();
	bjet_max_CSV_found = true; 
      }
    }
  } 

  bjet = bjet_max_pt; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(!b_candidate_found) cout << "No Bjet candidate found in event:  " << event.event << endl;
  if(!bjet_found) cout << "No Bjet found in event: " << event.event << endl;
  if(!bjet_max_CSV_found) cout << "No leading Bjet found in event: " << event.event << endl;
  if(!bjet_max_pt_found) cout << "No Bjet with highest CSV found in event: " << event.event << endl;

 
  if(bjet_found){
    //get the leding additional jets
    // Jet add_jet;
    // Jet add_jet_lept;
    // bool add_jet_found = false;
    // bool add_jet_lep_found = false;
    double add_jet_maxpt = 0.;
    double add_jet_lept_maxpt = 0.;

    for( const auto & ak4jet : *ak4jets){
      if(deltaR(ak4jet, jet) > 0.8 && deltaR(ak4jet, bjet) > 0.1 && ak4jet.pt() > add_jet_maxpt ) {
	add_jet_maxpt = ak4jet.pt();
      }
      if((deltaPhi(ak4jet, mu) < (2*pi/3)) && deltaR(ak4jet, bjet) > 0.1 && ak4jet.pt() > add_jet_lept_maxpt) 
	add_jet_lept_maxpt = ak4jet.pt();
    }

    //get additional b-jets and non bjets
    std::vector<Jet> add_bjets;
    std::vector<Jet> add_non_bjets;   
    for ( const auto & ak4jet : *ak4jets ) {
      if (deltaR(ak4jet, jet) > 0.8 && deltaR(ak4jet, bjet) > 0.1){
	if (btag(ak4jet, event)) add_bjets.emplace_back(ak4jet);
     	else add_non_bjets.emplace_back(ak4jet);
      }
    }
    sort_by_pt(add_bjets); 
    sort_by_pt(add_non_bjets);

    //fill hists
    double mProbe_mLep_bjet = jet.v4().M()/(bjet.v4()+mu.v4()).M();
    hist("mProbe_mLep_bjet1")->Fill(mProbe_mLep_bjet, weight);

    double mProbe_mLep_bjet_pt = jet.v4().M()/(bjet_max_pt.v4()+mu.v4()).M();
    hist("mProbe_mLep_bjet_pt")->Fill(mProbe_mLep_bjet_pt, weight);

    double mProbe_mLep_bjet_CSV = jet.v4().M()/(bjet_max_CSV.v4()+mu.v4()).M();
    hist("mProbe_mLep_bjet_CSV")->Fill(mProbe_mLep_bjet_CSV, weight);

    double ptProbe_ptLep_bjet = jet.v4().Pt()/(bjet.v4()+mu.v4()).Pt();
    hist("ptProbe_ptLep_bjet1")->Fill(ptProbe_ptLep_bjet, weight);

    hist("dPhi_mu_b1")->Fill(deltaPhi(bjet, mu), weight);
    hist("pt_leading_add1")->Fill(add_jet_maxpt, weight);
    hist("pt_leading_add_lep1")->Fill(add_jet_lept_maxpt, weight);

    if(add_bjets.size()) hist("pt_leading_addB")->Fill(add_bjets.at(0).pt(), weight);
    else hist("pt_leading_addB")->Fill(0., weight);
    if(add_non_bjets.size()) hist("pt_leading_addNoB")->Fill(add_non_bjets.at(0).pt(), weight);
    else hist("pt_leading_addNoB")->Fill(0., weight);

    double mProbe_mbjet = jet.v4().M()/bjet.v4().M();
    h_ratio_mLB_vs_mB->Fill(mProbe_mLep_bjet, mProbe_mbjet, weight);
  }
}


ProbeJetHists::~ProbeJetHists(){}
