#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/AdditionalSelections.h"
#include <UHH2/common/include/JetCorrections.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/ElectronIds.h>
#include <UHH2/common/include/JetIds.h>
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/TopTagging/include/TopTaggingSelections.h"
#include "UHH2/TopTagging/include/ProbeJetHists.h"
#include <UHH2/common/include/Utils.h>
#include "UHH2/TopTagging/include/TopTaggingUtils.h"
#include "UHH2/HOTVR/include/HOTVRJetCorrector.h"


using namespace std;
using namespace uhh2;
using namespace uhh2examples;

class MistagSelectionModule: public AnalysisModule {
public:

  explicit MistagSelectionModule(Context & ctx);
  virtual bool process(Event &  event) override;

private:

  //correctors
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<uhh2::AnalysisModule> topjet_corr;
  std::unique_ptr<GenericJetResolutionSmearer> topjetJER_smearer;
  std::unique_ptr<TopJetGroomer> topjet_groomer;
  
  //cleaners
  std::unique_ptr<JetCleaner> jet_cleaner;
  std::unique_ptr<TopJetCleaner> topjet_cleaner;

  //selections
  std::unique_ptr<uhh2::Selection> trigger_selection;
  std::unique_ptr<AndSelection> first_selection;
  std::unique_ptr<uhh2::Selection> ht_selection;

  //histograms
  std::unique_ptr<ProbeJetHists> hists_probe, hists_probe_200, hists_probe_300to400, hists_probe_600to1200;
  std::vector<std::unique_ptr<uhh2::Hists>> hists_after_sel;

  bool useHTT, usePUPPI, useHOTVR;
  string version;

  bool isMC;

  Event::Handle<std::vector<Particle>> add_genjet;

};

MistagSelectionModule::MistagSelectionModule(Context & ctx){

  //=============================
  // access configuaration values
  //=============================

  isMC = (ctx.get("dataset_type") == "MC");
  
  usePUPPI = (ctx.get("usePUPPI", "<not set>") == "TRUE");

  useHOTVR = (ctx.get("useHOTVR", "<not set>") == "TRUE");

  if(usePUPPI) cout << "use PUPPI topjets" << endl;
  else cout << "use CHS topjets" << endl;

  if(useHOTVR) cout << "run HOTVR" << endl;
  else cout << "run CMS tagger" << endl;

  version = ctx.get("dataset_version", "");

  //=============================
  // IDs and corrections
  //=============================

  MuonId muid = AndId<Muon>(MuonIDTight(), Muon_MINIIso(0.2,"delta-beta"), Muon_dxydzCut(0.2, 0.5), PtEtaCut(5., 2.4));
  ElectronId eleid = AndId<Electron>(ElectronID_Spring16_veto, Electron_MINIIso(0.1, "delta-beta"), PtEtaCut(5., 2.4));

  JetId jetid = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(20.0, 2.4));

  string PU_variation ="central";
  PU_variation = ctx.get("PU_variation","central");

  common.reset(new CommonModules());
  common->set_muon_id(muid);
  common->set_electron_id(eleid);
  common->set_jet_id(jetid);
  common->set_HTjetid(jetid);
  common->switch_jetlepcleaner(true);
  common->switch_jetPtSorter();
  common->init(ctx, PU_variation);

  if(useHOTVR){
    if(usePUPPI) topjet_corr.reset(new TopJetCorrectionModules(ctx, TopJetCorrectionModules::AK8_PUPPI));
    else topjet_corr.reset(new TopJetCorrectionModules(ctx, TopJetCorrectionModules::AK8_CHS));
  }else{
    if(usePUPPI) topjet_corr.reset(new TopJetCorrectionModules(ctx, TopJetCorrectionModules::AK8_PUPPI));
    else topjet_corr.reset(new TopJetCorrectionModules(ctx, TopJetCorrectionModules::AK8_CHS));
  }

  if(isMC) {
    if(useHOTVR) topjetJER_smearer.reset(new GenericJetResolutionSmearer(ctx, "topjets", "gentopjets", false));
    else{
      add_genjet = ctx.get_handle<std::vector<Particle>>("slimmedGenJetsAK8");
      topjetJER_smearer.reset(new GenericJetResolutionSmearer(ctx, "topjets", "slimmedGenJetsAK8", false));
    }
  }

  topjet_groomer.reset(new TopJetGroomer()); //take the subjet sum (make sure that the subjets are corrected properly)
  topjet_cleaner.reset(new TopJetCleaner(ctx, TopJetId(PtEtaCut(170., 2.4))));
  
  first_selection.reset(new AndSelection(ctx,"first selection"));  
  first_selection->add<NMuonSelection>("veto on muons",0,0);
  first_selection->add<NElectronSelection>("veto on electrons",0,0);
  first_selection->add<NJetSelection>("Number of jets >= 2",2);

  ht_selection.reset(new HTCut(ctx,1000));

  trigger_selection.reset(new TriggerSelection("HLT_PFHT900_v*"));

  hists_probe.reset(new ProbeJetHists(ctx, "probe_jet_hists_post_sel"));
  hists_probe_200.reset(new ProbeJetHists(ctx, "probe_jet_hists_post_sel_pt200"));
  hists_probe_300to400.reset(new ProbeJetHists(ctx, "probe_jet_hists_post_sel_pt300to400"));
  hists_probe_600to1200.reset(new ProbeJetHists(ctx, "probe_jet_hists_post_sel_pt600to1200"));

  hists_after_sel.emplace_back(new EventHists(ctx, "Event_sel"));
  hists_after_sel.emplace_back(new MuonHists(ctx, "Muon_sel"));
  hists_after_sel.emplace_back(new ElectronHists(ctx, "Electron_sel"));
  hists_after_sel.emplace_back(new JetHists(ctx, "Jet_sel"));
  hists_after_sel.emplace_back(new TopJetHists(ctx, "TopJet_sel"));
}

bool MistagSelectionModule::process(Event & event) {

  //=============================         
  //run corrections
  //=============================    
  bool ok = common->process(event);
  if(!ok) return false;

  topjet_corr->process(event); 
  if(useHOTVR) topjet_groomer->process(event);
  if(isMC && !useHOTVR) topjetJER_smearer->process(event);

  //=============================         
  //run cleaners
  //=============================       
  topjet_cleaner->process(event);

  //======================================        
  //apply selections
  //======================================  
  if(!trigger_selection->passes(event)) return false;
  if(!first_selection->passes(event)) return false;
  if(!ht_selection->passes(event)) return false; 

  //======================================        
  //fill histograms after selection
  //====================================== 
  for(auto & h : hists_after_sel){
    h->fill(event);
  }
  if(event.topjets->size() > 0){
    sort_by_pt(*event.topjets);
    hists_probe->fill_probe(event, event.topjets->at(0));
    if(event.topjets->at(0).pt() > 200) hists_probe_200->fill_probe(event, event.topjets->at(0));
    if(event.topjets->at(0).pt() > 300 && event.topjets->at(0).pt() < 400) hists_probe_300to400->fill_probe(event, event.topjets->at(0));
    if(event.topjets->at(0).pt() > 600 && event.topjets->at(0).pt() < 1200) hists_probe_600to1200->fill_probe(event, event.topjets->at(0));
  }

  return true;

}

UHH2_REGISTER_ANALYSIS_MODULE(MistagSelectionModule)
