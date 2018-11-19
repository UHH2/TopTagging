#include <iostream> 
#include <memory> 

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include <UHH2/common/include/Utils.h>
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/TopTagging/include/TopTaggingUtils.h"

#include <TTree.h>
#include <TFile.h>

using namespace std; 
using namespace uhh2;

class GenObject{
public:
  GenObject(){
    m_size = 0.;
    m_pt = 0.;
  }
  void setGenParticle(GenParticle part) {m_part = part; m_pt = part.pt();}
  GenParticle getGenParticle() {return m_part;}
  double pt() const {return m_pt;}
  void setGenSize(double size) {m_size = size;}
  double gen_size() const { return m_size;}

private:
  GenParticle m_part;
  double m_size;
  double m_pt;
};

class PerformanceModule: public AnalysisModule {

public:
  
  explicit PerformanceModule(Context & ctx);
  virtual bool process(Event & event) override;
  ~PerformanceModule();

private:

  vector<GenObject> GetMergedHadronicTops(TTbarGen ttbarGen, double ptMin, double ptMax, double etaMax, double radius);
  vector<GenObject> GetQuarksAndGluons(vector<GenParticle> genpaticles, double ptMin, double ptMax, double etaMax);
  bool FillTree(const Event & event, GenObject part, int Npart);
  bool SetCMSTTv2Variables(const TopJet jet, bool matched);
  bool SetHOTVRVariables(const TopJet jet, bool matched);
  bool SetHTTVariables(const TopJet jet, bool matched);

  //correctors (JECs)
  std::unique_ptr<TopJetCorrectionModules> topjet_corr;
  std::unique_ptr<uhh2::AnalysisModule> HOTVR_pileupCorrector;
  std::unique_ptr<GenericJetResolutionSmearer> topjetJER_smearer;
  std::unique_ptr<TopJetGroomer> topjet_groomer;
  unique_ptr<uhh2::AnalysisModule> topjet_cleaner;
  
  //ttbar gen 
  unique_ptr<uhh2::AnalysisModule> ttbarGenProducer;
  uhh2::Event::Handle<TTbarGen> h_ttbarGen;
  Event::Handle<vector<TopJet>> h_CHSTopJets;
  Event::Handle<std::vector<Particle>> add_genjet;

  //selections
  unique_ptr<uhh2::Selection> HadronicTopSelection;

  //output tree
  TTree* _flatTree;
 
  //output variables
  int i_evt, i_parton, NPV, type, matched, gen_pdgId;
  double gen_pt, gen_eta, gen_phi, gen_size;

  //CMSTopTaggerV2 vriables
  double clf_softdropmass, clf_tau32, clf_highestSjCSV, clf_softdropmass_uncorrected;

  //HOTVR variables
  int clf_Nsubjets;
  double clf_jetmass, clf_fpt, clf_mpair, clf_tau32_groomed;
  double dr_jet;

  //HTT variables
  double clf_mass, clf_mass_raw, clf_mass_tag, clf_fRec, clf_fRec_raw, clf_fRec_tag;

  TString version;
  TString fileName;
  bool usePUPPI, useHOTVR, useHTT;
  bool useJetFlavor;
 
};

PerformanceModule::PerformanceModule(Context & ctx){

  //get configuration values
  version = (TString)ctx.get("dataset_version", "");
  TString path = (TString)ctx.get("OutputPath", "");
  TString postFix = (TString)ctx.get("PFix", "");
  cout << "postfix:  " << postFix << endl;
  fileName = path+"PerformanceTree_"+version+postFix+".root";

  usePUPPI = (ctx.get("usePUPPI", "<not set>") == "TRUE");
  useHOTVR = (ctx.get("useHOTVR", "<not set>") == "TRUE");
  useHTT = (ctx.get("useHTT", "<not set>") == "TRUE");

  useJetFlavor = (ctx.get("useJetFlavor", "<not set>") == "TRUE");

  topjet_cleaner.reset(new TopJetCleaner(ctx, TopJetId(PtEtaCut(0., 5.0))));

  //JEC and JER settings 
 
  //if(useHOTVR && !usePUPPI) topjet_corr.reset(new TopJetCorrectionModules(ctx, TopJetCorrectionModules::HOTVR_CHS));
  //else if(useHOTVR && usePUPPI) topjet_corr.reset(new TopJetCorrectionModules(ctx, TopJetCorrectionModules::HOTVR_PUPPI));
  //else if(usePUPPI) topjet_corr.reset(new TopJetCorrectionModules(ctx, TopJetCorrectionModules::AK8_PUPPI));
  if(usePUPPI) topjet_corr.reset(new TopJetCorrectionModules(ctx, TopJetCorrectionModules::AK8_PUPPI));
  else topjet_corr.reset(new TopJetCorrectionModules(ctx, TopJetCorrectionModules::AK8_CHS));

  if(useHOTVR) topjetJER_smearer.reset(new GenericJetResolutionSmearer(ctx, "topjets", "gentopjets", false));
  else{
    add_genjet = ctx.get_handle<std::vector<Particle>>("slimmedGenJetsAK8");
    topjetJER_smearer.reset(new GenericJetResolutionSmearer(ctx, "topjets", "slimmedGenJetsAK8", false));
  }
  
  topjet_groomer.reset(new TopJetGroomer()); // just take the subjet sum (make sure that the subjets are corrected properly)

  //ttbar gen settings
  ttbarGenProducer.reset(new TTbarGenProducer(ctx, "ttbargen", true));
  h_ttbarGen = ctx.get_handle<TTbarGen>("ttbargen");

  //set tree branches
  type = 0; 
  if( version.Contains("Zprime") ) type = 2;

  _flatTree = new TTree("tree", "tree");

  //ctx.do_undeclare_all_event_output();

  _flatTree->Branch("i_evt", &i_evt, "i_evt/I");
  _flatTree->Branch("i_parton", &i_parton, "i_parton/I");
  _flatTree->Branch("npv", &NPV, "npv/I");
  _flatTree->Branch("type", &type, "type/I");

  _flatTree->Branch("gen_pt", &gen_pt, "gen_pt/D");
  _flatTree->Branch("gen_eta", &gen_eta, "gen_eta/D");
  _flatTree->Branch("gen_phi", &gen_phi, "gen_phi/D");
  _flatTree->Branch("gen_size", &gen_size, "gen_size/D");

  _flatTree->Branch("gen_pdgId", &gen_pdgId, "gen_pdgId/I");

  if(useHOTVR){
    // _flatTree->Branch("dr_jet" , &dr_jet, "dr_jet/D");
    _flatTree->Branch("clf_Nsubjets" , &clf_Nsubjets, "clf_Nsubjets/I");
    _flatTree->Branch("clf_jetmass" , &clf_jetmass, "clf_jetmass/D");
    _flatTree->Branch("clf_fpt" , &clf_fpt, "clf_fpt/D");
    _flatTree->Branch("clf_mpair" , &clf_mpair , "clf_mpair/D");
    _flatTree->Branch("clf_tau32" , &clf_tau32_groomed , "clf_tau32/D");
  }
  else if(useHTT){
    //_flatTree->Branch("clf_Nsubjets", &clf_Nsubjets,       "clf_Nsubjets/I");
    _flatTree->Branch("mass",         &clf_mass,           "clf_mass/D");
    //  _flatTree->Branch("mass_raw",     &clf_mass_raw,       "clf_mass_raw/D");
    //  _flatTree->Branch("mass_tag",     &clf_mass_tag,       "clf_mass_tag/D");
    _flatTree->Branch("fRec",         &clf_fRec,           "clf_fRec/D");
    //  _flatTree->Branch("fRec_raw",     &clf_fRec_raw,       "clf_fRec_raw/D");
    //  _flatTree->Branch("fRec_tag",     &clf_fRec_tag,       "clf_fRec_tag/D");
    _flatTree->Branch("clf_tau32" ,   &clf_tau32_groomed , "clf_tau32/D");
  }else{
    _flatTree->Branch("clf_softdropmass", &clf_softdropmass, "clf_softdropmass/D");
    // _flatTree->Branch("clf_softdropmass_uncorrected", &clf_softdropmass_uncorrected, "clf_softdropmass_uncorrected/D");
    _flatTree->Branch("clf_tau32", &clf_tau32, "clf_tau32/D");
    _flatTree->Branch("clf_highestSubjetCSV", &clf_highestSjCSV, "clf_highestSubjetCSV/D");
  }  

  //_flatTree->Branch("matched", &matched, "matched/I");

}



bool PerformanceModule::process(Event & event){


  for(const auto & topjet: *event.topjets){
    cout << "event: "<< event.event << "   topjet pt: " << topjet.pt() << endl;

    LorentzVector subjet_sum(0,0,0,0);
    for (const auto s : topjet.subjets()) {
      subjet_sum += s.v4();
    }
    cout << "   subjet sum pt: " << subjet_sum.Pt() << endl;
    
    for(const auto & subjet: topjet.subjets()){
      cout << "            subjet pt: " << subjet.pt() << endl;
    }
  }
 
  assert(event.topjets);

  //Jet corrections
  if(useHOTVR){
    topjet_corr->process(event); 
    topjet_groomer->process(event);
    //topjet_cleaner->process(event);
    //topjetJER_smearer->process(event);
  }
  else if(useHTT) {
    topjet_corr->process(event); //apply topjet and subjet corrections
    topjet_groomer->process(event); //now get the subjet sum from corrected subjets
  }else{
    topjet_corr->process(event); //apply AK8 corrections on the full jet and AK4 corrections on the subjets
    topjetJER_smearer->process(event);
  }

  sort_by_pt(*event.topjets);

  // particles
  vector<GenObject> parts;

  if(version.Contains("Zprime")){
    ttbarGenProducer->process(event);
    const auto & ttbarGen = event.get(h_ttbarGen);
    //if(version.Contains("highPt")) parts = GetMergedHadronicTops(ttbarGen, 1000., 1400., 1.5, 0.6);
    //else parts = GetMergedHadronicTops(ttbarGen, 300., 470., 2.4, 1.2);
    parts = GetMergedHadronicTops(ttbarGen, 200., uhh2::infinity, 2.4, 99.);
  }
  if(version.Contains("QCD")){
    //if(version.Contains("highPt")) parts = GetQuarksAndGluons(*event.genparticles, 1000., 1400., 1.5);
    //else parts = GetQuarksAndGluons(*event.genparticles, 300., 470., 2.4);
    parts = GetQuarksAndGluons(*event.genparticles, 200., uhh2::infinity, 2.4);
  }

  //fill output tree
  int Npart = 1; 
  for( const auto & part : parts) {
    FillTree(event, part, Npart);
    Npart++; 
  }
 
  return true;
}

PerformanceModule::~PerformanceModule(){
  TFile* _file = new TFile(fileName, "Recreate");
   _flatTree->Write("PerformanceTree");
   _file->Close();
}



bool PerformanceModule::FillTree(const Event & event, GenObject g_part, int Npart){

  //-------------------------------------------
  //get TopJets and match them to the particle
  //-------------------------------------------
  vector<TopJet> topjets = *event.topjets;
  GenParticle part = g_part.getGenParticle();

  bool matchedJet = false;
  double radius = 0.6;
  TopJet jet; 

  double mindR = 5.;
  // double maxPt = 1.;
  for( const auto & topjet : topjets){
    double dR = deltaR( topjet.v4(), part.v4());
    /* if(useHOTVR){
      radius = 600./topjet.pt();
      if(radius < 0.1) radius = 0.1;
      if(radius > 1.5) radius = 1.5;
      radius *= 0.75;
    }
    if(useHOTVR){
      if(dR < mindR){
	jet = topjet;
	mindR = dR;
	//maxPt = topjet.pt();
	matchedJet = true;
      }
      }else{*/
    if(dR < radius){
      if(dR < mindR){
	//if(topjet.pt() > maxPt){
	if(matchedJet) cout << "more than one jet to match" << endl;
	jet = topjet;
	mindR = dR;
	//maxPt = topjet.pt();
	matchedJet = true;
      }
    }
    // } 
  }
 
  // if(!matchedJet){
  //  cout << "NO JET MATCHED to the particle" << endl;
  //   return false;
  //}

  //----------------------------------------------
  //set the output variables
  //----------------------------------------------
  i_evt = event.event;
  i_parton = Npart;
  
  gen_pt = part.pt();
  gen_eta = part.eta();
  gen_phi = part.phi();
  gen_size = g_part.gen_size();

  dr_jet = mindR;

  gen_pdgId = part.pdgId();

  NPV = event.pvs->size();
 
  if(useHOTVR && !useHTT) SetHOTVRVariables(jet, matchedJet);
  else if (useHTT && !useHOTVR) SetHTTVariables(jet, matchedJet);
  else if (!useHTT && !useHOTVR) SetCMSTTv2Variables(jet, matchedJet);
  else cout << "no tagger output specified" << endl;

  if(matchedJet) matched = 1;
  else matched = 0; 

  _flatTree->Fill();

  return true;
}


vector<GenObject> PerformanceModule::GetQuarksAndGluons(vector<GenParticle> genpaticles, double ptMin, double ptMax, double etaMax){

  vector<GenObject> parts;

  for(const auto & part : genpaticles){
    int absId = fabs(part.pdgId());

    if( (absId < 6 || absId == 21) 
	&& part.status() == 23
	&& part.pt() > ptMin 
	&& part.pt() < ptMax 
	&& fabs(part.eta()) < etaMax){
      GenObject g_part;
      g_part.setGenParticle(part);
      parts.push_back(g_part);
    }
  }

  sort_by_pt(parts);
  return parts;

}


vector<GenObject> PerformanceModule::GetMergedHadronicTops(TTbarGen ttbarGen, double ptMin, double ptMax, double etaMax, double radius){

  vector<GenObject> tops; 

  if(ttbarGen.IsTopHadronicDecay()) {
    GenParticle top = ttbarGen.Top_beforeDecay();
    vector<GenParticle> quarks;
    quarks.push_back(ttbarGen.bTop());
    quarks.push_back(ttbarGen.Wdecay1());
    quarks.push_back(ttbarGen.Wdecay2());
    bool merged = true;
    double gen_size = 0.;
    for( const auto & q : quarks){
      if( deltaR(top.v4(), q.v4()) > radius ) merged = false;
      if( deltaR(top.v4(), q.v4()) > gen_size) gen_size = deltaR(top.v4(), q.v4());
    }
    if( merged 
	&& top.pt() > ptMin 
	&& top.pt() < ptMax 
	&& fabs(top.eta()) < etaMax){
      GenObject g_top;
      g_top.setGenParticle(top);
      g_top.setGenSize(gen_size);
      tops.push_back(g_top);
    }
  }
  if(ttbarGen.IsAntiTopHadronicDecay()) {
    GenParticle top = ttbarGen.Antitop_beforeDecay();
    vector<GenParticle> quarks;
    quarks.push_back(ttbarGen.bAntitop());
    quarks.push_back(ttbarGen.WMinusdecay1());
    quarks.push_back(ttbarGen.WMinusdecay2());
    bool merged = true;
    double gen_size = 0.;
    for( const auto & q : quarks){
      if( deltaR(top.v4(), q.v4()) > radius ) merged = false;
      if( deltaR(top.v4(), q.v4()) > gen_size) gen_size = deltaR(top.v4(), q.v4()); 
    }
    if( merged 
	&& top.pt() > ptMin 
	&& top.pt() < ptMax 
	&& fabs(top.eta()) < etaMax){ 
      GenObject g_top;
      g_top.setGenParticle(top);
      g_top.setGenSize(gen_size);
      tops.push_back(g_top);
    }
  }
  
  sort_by_pt(tops);

  return tops;

}


bool PerformanceModule::SetHTTVariables(const TopJet jet, bool matched){

  auto subjets = jet.subjets();
  sort_by_pt(subjets);
  
  clf_Nsubjets = subjets.size();

  LorentzVector subjet_sum(0,0,0,0);
  for (const auto s : subjets) {
    subjet_sum += s.v4();
  }

  LorentzVector subjet_sumUNCORR(0,0,0,0);
  for (const auto s : subjets) {
    subjet_sumUNCORR += s.v4()*s.JEC_factor_raw();
  }
  
  double mass = subjet_sum.M();
  double mass_raw = subjet_sumUNCORR.M();

  double m12 = 0., m13 = 0., m23 = 0.;
  if(subjets.size() > 2){
    m12 = (subjets.at(0).v4() + subjets.at(1).v4()).M();
    m13 = (subjets.at(0).v4() + subjets.at(2).v4()).M();
    m23 = (subjets.at(1).v4() + subjets.at(2).v4()).M();
  }

  double m12_raw = 0., m13_raw = 0., m23_raw = 0.;
  if(subjets.size() > 2){
    m12_raw = (subjets.at(0).v4()*subjets.at(0).JEC_factor_raw() + subjets.at(1).v4()*subjets.at(1).JEC_factor_raw()).M();
    m13_raw = (subjets.at(0).v4()*subjets.at(0).JEC_factor_raw() + subjets.at(2).v4()*subjets.at(2).JEC_factor_raw()).M();
    m23_raw = (subjets.at(1).v4()*subjets.at(1).JEC_factor_raw() + subjets.at(2).v4()*subjets.at(2).JEC_factor_raw()).M();
  }

  double mW_mt = 80.4/172.3;

  double fRec = 0.;
  if(subjets.size() > 2) fRec = min(fabs(((m12/mass)/mW_mt)-1), min( fabs(((m23/mass)/mW_mt)-1), fabs(((m13/mass)/mW_mt)-1)));

  double fRec_raw = 0.;
  if(subjets.size() > 2)fRec_raw = min(fabs(((m12_raw/mass_raw)/mW_mt)-1), min( fabs(((m23_raw/mass_raw)/mW_mt)-1), fabs(((m13_raw/mass_raw)/mW_mt)-1)));

  double fRec_tag = 0.;
  if (jet.has_tag(jet.tagname2tag("fRec"))) fRec_tag = jet.get_tag(jet.tagname2tag("fRec"));

  double mass_tag = 0.;
  if (jet.has_tag(jet.tagname2tag("mass"))) mass_tag = jet.get_tag(jet.tagname2tag("mass"));			
  
  clf_Nsubjets = subjets.size();
  clf_mass = mass;
  clf_mass_raw = mass_raw;
  clf_mass_tag = mass_tag;
  clf_fRec = fRec;
  clf_fRec_raw = fRec_raw;
  clf_fRec_tag = fRec_tag;
  clf_tau32_groomed = jet.tau3_groomed()/jet.tau2_groomed();
  
  if(!matched) clf_tau32_groomed = 99.;
  
  return true;
}

bool PerformanceModule::SetCMSTTv2Variables(const TopJet jet, bool matched){

  auto subjets = jet.subjets();

  LorentzVector subjet_sum(0,0,0,0);
  for (const auto s : subjets) {
    subjet_sum += s.v4();
  }

  LorentzVector subjet_sumUNCORR(0,0,0,0);
  for (const auto s : subjets) {
    subjet_sumUNCORR += s.v4()*s.JEC_factor_raw();
  }

  clf_softdropmass = subjet_sum.M();
  clf_softdropmass_uncorrected = subjet_sumUNCORR.M();
  clf_tau32 = jet.tau3()/jet.tau2();

  if(!matched) clf_tau32 = 99.;

  double highestCSV = 0.;
  for (const auto s : subjets) {
    double CSV = s.btag_combinedSecondaryVertex();
    if(CSV > highestCSV) highestCSV = CSV;
  }

  clf_highestSjCSV = highestCSV;

  return true;
} 

bool PerformanceModule::SetHOTVRVariables(const TopJet jet, bool matched){

  double m12 = 0., m13 = 0., m23 = 0.;

  if(jet.subjets().size()){
    auto subjets = jet.subjets();
    sort_by_pt(subjets);
  
    clf_Nsubjets = subjets.size();
    
    clf_fpt = subjets.at(0).pt()/jet.pt();

    if(clf_Nsubjets > 2){
      m12 = (subjets.at(0).v4() + subjets.at(1).v4()).M();
      m13 = (subjets.at(0).v4() + subjets.at(2).v4()).M();
      m23 = (subjets.at(1).v4() + subjets.at(2).v4()).M();
    }
  }else{
    clf_Nsubjets = 0;
    clf_fpt = 0;
  }  
  
  clf_mpair = min(min(m12, m13), m23);
    
  clf_jetmass = jet.v4().M();

  clf_tau32_groomed = jet.tau3_groomed()/jet.tau2_groomed();
  
  if(!matched) clf_tau32_groomed = 99.;

  return true;
} 

UHH2_REGISTER_ANALYSIS_MODULE(PerformanceModule)

  
