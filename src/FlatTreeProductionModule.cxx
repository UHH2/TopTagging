#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/TTbarGen.h"

#include "UHH2/TopTagging/include/TopTaggingSelections.h"
#include "UHH2/TopTagging/include/ProbeJetHists.h"


using namespace std;
using namespace uhh2;
using namespace uhh2examples;


class FlatTreeProductionModule: public AnalysisModule {
public:
    
  explicit FlatTreeProductionModule(Context & ctx);
  virtual bool process(Event & event) override;


private:

  //correctors
  std::unique_ptr<CommonModules> common;

  //reweighting and scale factors
  std::vector<std::unique_ptr<AnalysisModule>> reweighting_modules;
  std::unique_ptr<uhh2::AnalysisModule> muo_tight_noniso_SF, muo_trigger_SF;
  std::unique_ptr<uhh2::AnalysisModule> btagwAK8, subjet_btagwAK8;
  std::unique_ptr<uhh2::AnalysisModule> scale_variation;
  std::unique_ptr<GenericJetResolutionSmearer> topjetJER_smearer;

  std::unique_ptr<uhh2::AnalysisModule> ttGenProducer;

  std::unique_ptr<TopJetCleaner> topjet_cleaner;

  //selections
  std::unique_ptr<MergedSelection> mergedEvent_selection; 

  //variables and functions
  bool fill_PDF;
  bool isMC;
  bool useHTT, usePUPPI;
  bool merged_category;
  string version;

  uhh2::Event::Handle<TTbarGen> h_ttbarGen;

  //output branches
  Event::Handle<double> h_weight;

  Event::Handle<double> h_jet_mass;
  Event::Handle<double> h_jet_pt;
  Event::Handle<double> h_jet_eta;
  Event::Handle<double> h_m_tt;

  //hists
  std::unique_ptr<ProbeJetHists> hists, hists200, hists300;

};


FlatTreeProductionModule::FlatTreeProductionModule(Context & ctx){

  //=============================
  // access configuaration values
  //=============================

  isMC = (ctx.get("dataset_type") == "MC");
  
  useHTT = (ctx.get("useHTT", "<not set>") == "TRUE");
  usePUPPI = (ctx.get("usePUPPI", "<not set>") == "TRUE");

  fill_PDF = (ctx.get("fill_PDF", "FALSE") == "TRUE");

  if(usePUPPI) cout << "use PUPPI topjets" << endl;
  else cout << "use CHS topjets" << endl;

  if(useHTT) cout << "run the HTT" << endl;
  else cout << "run CMS tagger" << endl;
  
  version = ctx.get("dataset_version", "");

  string merged = ctx.get("MergedSelection", "All");
  TString vers = (TString)version;
  if(vers.Contains("mergedTop")) merged = "mergedTop";
  if(vers.Contains("mergedW")){
    merged = "mergedW";
    cout << "merged = mergedW" << endl;
  }
  if(vers.Contains("mergedQB")) merged = "mergedQB";
  if(vers.Contains("light")) merged = "light";
  if(vers.Contains("bkg")) merged = "bkg";
  if(vers.Contains("notmerged")) merged = "notmerged";

  JetId jetid = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0, 2.4));
  TopJetId topjetid = AndId<TopJet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(150., 2.4));
  topjet_cleaner.reset(new TopJetCleaner(ctx,topjetid));


  string PU_variation = "central";
  PU_variation = ctx.get("PU_variation","central");

  string BTag_variation = "central";
  BTag_variation = ctx.get("BTag_variation","central");

  string SubjetBTag_variation = "central";
  SubjetBTag_variation = ctx.get("SubjetBTag_variation","central");

  string MuonID_variation = "none";
  MuonID_variation = ctx.get("MuonID_variation","none");

  string MuonTrigger_variation = "none";
  MuonTrigger_variation = ctx.get("MuonTrigger_variation","none");


  //===========================
  //setup corrections
  //===========================

  common.reset(new CommonModules());
  common->set_HTjetid(jetid);
  common->disable_jec();
  common->disable_jersmear();
  common->disable_lumisel();
  common->disable_metfilters();
  common->disable_pvfilter();
  common->disable_jetpfidfilter();
  common->init(ctx, PU_variation);
  cout << "common init" <<endl;


  //==============================
  //reweighting and scale factors
  //==============================

  btagwAK8.reset(new MCBTagScaleFactor(ctx, CSVBTag::WP_MEDIUM, "jets", BTag_variation,"mujets","incl","MCBtagEfficiencies"));

  subjet_btagwAK8.reset(new MCBTagScaleFactor(ctx, CSVBTag::WP_LOOSE, "topjets", SubjetBTag_variation,"lt","incl","MCSubjetBtagEfficiencies","", "SubjetBTagCalibration"));//,"SubjetBTagCalibration"));

  muo_tight_noniso_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/dreyert/CMSSW_8_0_24_patch1/src/UHH2/common/data//MuonID_EfficienciesAndSF_average_RunBtoH.root","MC_NUM_TightID_DEN_genTracks_PAR_pt_eta",1, "tightID", true, MuonID_variation));
  muo_trigger_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/dreyert/CMSSW_8_0_24_patch1/src/UHH2/common/data//MuonTrigger_EfficienciesAndSF_average_RunBtoH.root","IsoMu50_OR_IsoTkMu50_PtEtaBins",1, "muonTrigger", true, MuonTrigger_variation));
 
  scale_variation.reset(new MCScaleVariation(ctx));

  ttGenProducer.reset(new TTbarGenProducer(ctx, "ttbargen", true));
  h_ttbarGen = ctx.get_handle<TTbarGen>("ttbargen");
 
  //==========
  //matching
  //=========
  double jet_radius =  0.8;
  //if (useHTT) jet_radius = 1.5; 
  //else jet_radius = 0.8; 

  merged_category = true;
  if (merged == "mergedTop") mergedEvent_selection.reset( new MergedSelection(ctx, "ttbargen", jet_radius));
  else if (merged == "mergedW") mergedEvent_selection.reset( new MergedSelection(ctx, "ttbargen", jet_radius, MergedSelection::oMergedW));
  else if (merged == "mergedQB") mergedEvent_selection.reset( new MergedSelection(ctx, "ttbargen", jet_radius, MergedSelection::oBplusQ));
  else if (merged == "light") mergedEvent_selection.reset( new MergedSelection(ctx, "ttbargen", jet_radius, MergedSelection::oLight));
  else if (merged == "bkg") mergedEvent_selection.reset( new MergedSelection(ctx, "ttbargen", jet_radius, MergedSelection::oBkg));
  else if (merged == "notmerged") mergedEvent_selection.reset( new MergedSelection(ctx, "ttbargen", jet_radius, MergedSelection::oNotMerged));
  else merged_category = false;

  //================
  //declare new output 
  //================
  ctx.undeclare_all_event_output();

  h_weight = ctx.declare_event_output<double>("weight");
  h_jet_mass = ctx.declare_event_output<double>("jet_mass");
  h_jet_pt = ctx.declare_event_output<double>("jet_pt");
  h_jet_eta = ctx.declare_event_output<double>("jet_eta");
  h_m_tt = ctx.declare_event_output<double>("mttbar");

  hists.reset(new ProbeJetHists(ctx, "histograms"));
  hists200.reset(new ProbeJetHists(ctx, "histograms200"));
  hists300.reset(new ProbeJetHists(ctx, "histograms300"));
}

bool FlatTreeProductionModule::process(Event &event) {

  //========================================
  //corrections, reweighting, scale factors
  //========================================

  bool ok = common->process(event);
  if(!ok) return false;

  //  topjet_cleaner->process(event);
  
  //apply muon b tagging scale factors 
  muo_tight_noniso_SF->process(event);
  muo_trigger_SF->process(event);
  btagwAK8->process(event);

  scale_variation->process(event);

  if( ((TString)version).Contains("TTbar")) ttGenProducer->process(event);

  //=====================
  //get the probe jet
  //=====================

  std::vector<TopJet>* topjets = event.topjets;
  std::vector<Muon>* muons = event.muons;
  TopJet probe_jet; 
  bool probejet_found = false;

  for(auto & topjet : *topjets){
    double pi = 3.14159265359;
    double delPhi = deltaPhi(topjet, muons->at(0));
    if(delPhi > ((2./3.)*pi)){
      probe_jet = topjet;
      probejet_found = true;
      break;
    }
  }

  if(!probejet_found){
    //   cout << "WARNING: no probe jet found"<< endl;
    return false;
  }
  
  if( ((TString)version).Contains("TTbar") && merged_category && !mergedEvent_selection->passes_probe(event, probe_jet)) return false;

  //====================
  //Merged W Selection
  //====================
  double probejet_tau32 = probe_jet.tau3()/probe_jet.tau2();
  double probejet_tau21 = probe_jet.tau2()/probe_jet.tau1();

  //get a subjet b tag on the probe jet and the probe jet mass from the subjets
  bool subjet_btag = false;

  for(const auto & subjet : probe_jet.subjets()){
    JetId btag = CSVBTag(CSVBTag::WP_LOOSE);
    if( btag(subjet, event) ) subjet_btag = true;
  }

  LorentzVector subjet_sum(0,0,0,0);
  for(const auto & subjet : probe_jet.subjets()){
    subjet_sum += subjet.v4();
  }

  double probejet_mass = subjet_sum.M();

  //selection 
  if(subjet_btag) return false; //subjet b-tag veto
  if(probejet_tau32 < 0.8 ) return false; // tau32 > 0.8 (loose top tagging wp)
  if(probejet_tau21 > 0.4 ) return false; // tau21 < 0.4 (medium W tag)

  double mtt = 0.;
  if( ((TString)version).Contains("TTbar")){
    const TTbarGen& ttbarGen = event.get(h_ttbarGen);
    mtt = (ttbarGen.Top().v4()+ttbarGen.Antitop().v4()).M();
  }

  //================
  //set output 
  //================
  event.set(h_weight, event.weight);
  event.set(h_jet_mass, probejet_mass);
  event.set(h_jet_pt, probe_jet.pt());
  event.set(h_jet_eta, probe_jet.eta());
  event.set(h_m_tt, mtt);

  hists->fill_probe(event, probe_jet);
  if(probe_jet.pt() > 200) hists200->fill_probe(event, probe_jet);
  if(probe_jet.pt() > 300) hists300->fill_probe(event, probe_jet);

  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(FlatTreeProductionModule)
