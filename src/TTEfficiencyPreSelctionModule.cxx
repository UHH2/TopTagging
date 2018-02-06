#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/LuminosityHists.h"
//#include <UHH2/common/include/JetCorrections.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/ElectronIds.h>
#include <UHH2/common/include/JetIds.h>
#include "UHH2/common/include/NSelections.h"
#include "UHH2/TopTagging/include/TopTaggingSelections.h"
#include <UHH2/common/include/Utils.h>


using namespace std;
using namespace uhh2;
//using namespace uhh2examples;

class TTEfficiencyPreSelectionModule: public AnalysisModule {
public:
    
  explicit TTEfficiencyPreSelectionModule(Context & ctx);
  virtual bool process(Event & event) override;

private:
    
  //correctors
  std::unique_ptr<CommonModules> common;
  
  //cleaners
  std::unique_ptr<JetCleaner> jet_cleaner;

  //selections
  // std::unique_ptr<uhh2::Selection> muon_sel, pv_sel, jet_sel;
  std::unique_ptr<uhh2::Selection> htlep_sel, met_sel;
  std::unique_ptr<AndSelection> pre_selection;

  //histograms
  std::vector<std::unique_ptr<uhh2::Hists>> hists_before_presel;
  std::vector<std::unique_ptr<uhh2::Hists>> hists_after_presel; 
};


TTEfficiencyPreSelectionModule::TTEfficiencyPreSelectionModule(Context & ctx){

  const bool isMC = (ctx.get("dataset_type") == "MC");

  MuonId muid = AndId<Muon>(MuonIDTight(), PtEtaCut(55., 2.4));
  ElectronId eleid = AndId<Electron>(ElectronID_Spring16_medium_noIso, PtEtaCut(55., 2.4));

  common.reset(new CommonModules());

  common->set_muon_id(muid);
  common->set_electron_id(eleid);

  common->switch_jetlepcleaner(true);
  common->switch_jetPtSorter();
  common->switch_metcorrection();

  common->disable_mcpileupreweight();
  //common->disable_jersmear();

  common->init(ctx);
  cout << "common init" <<endl;

  jet_cleaner.reset(new JetCleaner(ctx, 15., 2.4));

  //selections
  met_sel.reset(new METCut(20., std::numeric_limits<double>::infinity()));
  htlep_sel.reset(new HTlepCut(70., std::numeric_limits<double>::infinity()));
 
  pre_selection.reset(new AndSelection(ctx,"first selection"));  
  pre_selection->add<NMuonSelection>("Number of muons == 1",1,1);
  pre_selection->add<NJetSelection>("Number of jets >= 2",2);
 
  //histograms
  hists_before_presel.emplace_back(new EventHists(ctx, "Event_before_presel"));
  hists_before_presel.emplace_back(new MuonHists(ctx, "Muon_before_presel"));
  hists_before_presel.emplace_back(new ElectronHists(ctx, "Electron_before_presel"));
  hists_before_presel.emplace_back(new JetHists(ctx, "Jet_before_presel"));
  hists_before_presel.emplace_back(new TopJetHists(ctx, "TopJet_before_presel"));
  hists_before_presel.emplace_back(new LuminosityHists(ctx, "Lumi_presel"));

  hists_after_presel.emplace_back(new EventHists(ctx, "Event_after_presel"));
  hists_after_presel.emplace_back(new MuonHists(ctx, "Muon_after_presel"));
  hists_after_presel.emplace_back(new ElectronHists(ctx, "Electron_after_presel"));
  hists_after_presel.emplace_back(new JetHists(ctx, "Jet_after_presel"));
  hists_after_presel.emplace_back(new TopJetHists(ctx, "TopJet_after_presel"));
  hists_after_presel.emplace_back(new LuminosityHists(ctx, "Lumi_after_presel"));

}



bool TTEfficiencyPreSelectionModule::process(Event & event) {
   
  //cout << "TTEfficiencyPreSelectionModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;

  //fill histograms before preselection
  for(auto & h : hists_before_presel){
    h->fill(event);
  }

  //keep uncleaned jets and MET
  std::unique_ptr< std::vector<Jet> >    uncleaned_jets   (new std::vector<Jet>   (*event.jets));
  std::unique_ptr<MET> uncleaned_met(new MET(*event.met));

  //run corrections
  bool ok = common->process(event);
  if(!ok) return false;

  //run jet cleaners 
  jet_cleaner->process(event);
    
  //apply selections
  if(!pre_selection->passes(event)) return false;
  if(!met_sel->passes(event)) return false;
  if(!htlep_sel->passes(event)) return false;
 
  // store uncleaned jets and MET
  event.jets->clear();
  event.jets->reserve(uncleaned_jets->size());
  for(const auto & j : *uncleaned_jets) event.jets->push_back(j); 
  sort_by_pt<Jet>(*event.jets);

  event.met->set_pt(uncleaned_met->pt());
  event.met->set_phi(uncleaned_met->phi());

  // fill histograms after preselection
  for(auto & h : hists_after_presel){
    h->fill(event);
  }
  
  return true;

}

UHH2_REGISTER_ANALYSIS_MODULE(TTEfficiencyPreSelectionModule)


