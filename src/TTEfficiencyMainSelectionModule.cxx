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
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/common/include/AdditionalSelections.h"
#include <UHH2/common/include/JetCorrections.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/ElectronIds.h>
#include <UHH2/common/include/JetIds.h>
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/TopPtReweight.h"
#include "UHH2/TopTagging/include/TopTaggingSelections.h"
#include "UHH2/TopTagging/include/ProbeJetHists.h"
#include <UHH2/common/include/Utils.h>
#include "UHH2/TopTagging/include/TopTaggingUtils.h"
#include "UHH2/HOTVR/include/HOTVRJetCorrector.h"

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

class TTEfficiencyMainSelectionModule: public AnalysisModule {
public:
    
  explicit TTEfficiencyMainSelectionModule(Context & ctx);
  virtual bool process(Event & event) override;

private:
    
  //correctors
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<uhh2::AnalysisModule> topjet_corr;
  std::unique_ptr<GenericJetResolutionSmearer> topjetJER_smearer;
  std::unique_ptr<TopJetGroomer> topjet_groomer;
  
  //cleaners
  std::unique_ptr<JetCleaner> jet_cleaner;
  std::unique_ptr<TopJetCleaner> topjet_cleaner;

  //reweighting and scale factors
  std::vector<std::unique_ptr<AnalysisModule>> reweighting_modules;
  std::unique_ptr<uhh2::AnalysisModule> muo_tight_noniso_SF, muo_trigger_SF;
  std::unique_ptr<AnalysisModule> btagwAK8;

  //selections
  std::unique_ptr<uhh2::Selection>  met_sel, ptW_sel, bjetCloseToLepton_sel, mttgen_sel, twoDcut;
  std::unique_ptr<AndSelection> first_selection;
  std::unique_ptr<uhh2::Selection> hadronic_selection, lepton_jets_seletion, dilepton_selection, tau_jets_selection; 

  //histograms
  std::unique_ptr<BTagMCEfficiencyHists>  hists_btag_medium_eff;
  std::vector<std::unique_ptr<uhh2::Hists>> hists_before_sel;
  std::vector<std::unique_ptr<uhh2::Hists>> hists_after_sel; 

  bool useHTT, usePUPPI, useHOTVR;
  string version;

  bool isMC;

  Event::Handle<std::vector<Particle>> add_genjet;
  
};


TTEfficiencyMainSelectionModule::TTEfficiencyMainSelectionModule(Context & ctx){

  //=============================
  // access configuaration values
  //=============================

  isMC = (ctx.get("dataset_type") == "MC");
  
  useHTT = (ctx.get("useHTT", "<not set>") == "TRUE");
  usePUPPI = (ctx.get("usePUPPI", "<not set>") == "TRUE");

  useHOTVR = (ctx.get("useHOTVR", "<not set>") == "TRUE");

  if(usePUPPI) cout << "use PUPPI topjets" << endl;
  else cout << "use CHS topjets" << endl;

  if(useHTT) cout << "run the HTT" << endl;
  else if(useHOTVR) cout << "run HOTVR" << endl;
  else cout << "run CMS tagger" << endl;

  version = ctx.get("dataset_version", "");


  bool TopPtReweighting = false;
  TopPtReweighting = (ctx.get("TopPtReweight","FALSE")== "TRUE");


  //=============================
  // IDs and corrections
  //=============================

  // MuonId muid = AndId<Muon>(MuonIDTight(), MuonIso(0.15), Muon_dxydzCut(0.2, 0.5), PtEtaCut(45., 2.4));
  MuonId muid = AndId<Muon>(MuonIDTight(), PtEtaCut(55., 2.4));
  ElectronId eleid = AndId<Electron>(ElectronID_Spring16_medium_noIso, PtEtaCut(55., 2.4));

  JetId jetid = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0, 2.4));

  string PU_variation ="central";
  PU_variation = ctx.get("PU_variation","central");

  //common modules
  common.reset(new CommonModules());
  common->set_jet_id(jetid);
  common->set_muon_id(muid);
  common->set_electron_id(eleid);
  common->switch_jetlepcleaner(true);
  common->switch_jetPtSorter();
  common->switch_metcorrection();
  common->set_HTjetid(jetid);
  common->init(ctx, PU_variation);
  cout << "common init" <<endl;
	      
  if(useHOTVR){
    // if(usePUPPI) topjet_corr.reset(new TopJetCorrectionModules(ctx, TopJetCorrectionModules::HOTVR_PUPPI));
    // else topjet_corr.reset(new TopJetCorrectionModules(ctx, TopJetCorrectionModules::HOTVR_CHS));
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

  //=============================
  // cleaners
  //=============================

  jet_cleaner.reset(new JetCleaner(ctx, 30., 2.4));

  if(useHTT) topjet_cleaner.reset(new TopJetCleaner(ctx, TopJetId(PtEtaCut(150., 2.4))));
  else topjet_cleaner.reset(new TopJetCleaner(ctx, TopJetId(PtEtaCut(170., 2.4))));


  //=============================
  //reweighting and scale factors
  //=============================
  if (version == "TTbar_Incl" || version ==  "TTbar_700to1000" || version ==  "TTbar_1000toInf") {
    reweighting_modules.emplace_back(new TTbarGenProducer(ctx, "ttbargen", true));
    // reweighting_modules.emplace_back(new TopPtReweight(ctx, 0.156, -0.00137, "ttbargen", "weight_ttbar", true)); //8TeV
    if(TopPtReweighting) reweighting_modules.emplace_back(new TopPtReweight(ctx, 0.0615, -0.0005, "ttbargen", "weight_ttbar", true)); //13TeV
  }

  if (version == "TTbar_Incl" ) mttgen_sel.reset(new MttbarGenSelection(0., 700.));

  btagwAK8.reset(new MCBTagScaleFactor(ctx, CSVBTag::WP_MEDIUM, "jets","central","mujets","incl","MCBtagEfficiencies"));
  muo_tight_noniso_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/dreyert/CMSSW_8_0_24_patch1/src/UHH2/common/data//MuonID_EfficienciesAndSF_average_RunBtoH.root","MC_NUM_TightID_DEN_genTracks_PAR_pt_eta",1, "tightID"));
  muo_trigger_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/dreyert/CMSSW_8_0_24_patch1/src/UHH2/common/data//MuonTrigger_EfficienciesAndSF_average_RunBtoH.root","IsoMu50_OR_IsoTkMu50_PtEtaBins",1, "muonTrigger"));
 
   
  //=============================
  //selections 
  //=============================
 
  first_selection.reset(new AndSelection(ctx,"first selection"));  
  first_selection->add<NMuonSelection>("Number of muons == 1",1,1);
  first_selection->add<NJetSelection>("Number of jets >= 2",2);
  
  twoDcut.reset(new TwoDCut(.4, 25.));
  met_sel.reset(new METCut(50., std::numeric_limits<double>::infinity()));
  // ptW_sel.reset(new PtWSelection(250.));
  ptW_sel.reset(new PtWSelection(150.));
  bjetCloseToLepton_sel.reset(new NMuonBTagSelection(1,999,CSVBTag(CSVBTag::WP_MEDIUM)));
 
  //partonlevel decay channel slections
  hadronic_selection.reset( new DecayChannelSelection(ctx, "ttbargen", "hadronic"));
  lepton_jets_seletion.reset( new DecayChannelSelection(ctx, "ttbargen", "lepton_jets"));
  dilepton_selection.reset( new DecayChannelSelection(ctx, "ttbargen", "dilepton"));
  tau_jets_selection.reset( new DecayChannelSelection(ctx, "ttbargen", "tau_jets"));

  //=============================         
  //histograms
  //=============================
  hists_btag_medium_eff.reset(new BTagMCEfficiencyHists(ctx,"BTagMedium",CSVBTag::WP_MEDIUM));

  hists_before_sel.emplace_back(new EventHists(ctx, "Event_presel"));
  hists_before_sel.emplace_back(new MuonHists(ctx, "Muon_presel"));
  hists_before_sel.emplace_back(new ElectronHists(ctx, "Electron_presel"));
  hists_before_sel.emplace_back(new JetHists(ctx, "Jet_presel"));
  hists_before_sel.emplace_back(new TopJetHists(ctx, "TopJet_presel"));

  hists_after_sel.emplace_back(new EventHists(ctx, "Event_sel"));
  hists_after_sel.emplace_back(new MuonHists(ctx, "Muon_sel"));
  hists_after_sel.emplace_back(new ElectronHists(ctx, "Electron_sel"));
  hists_after_sel.emplace_back(new JetHists(ctx, "Jet_sel"));
  hists_after_sel.emplace_back(new TopJetHists(ctx, "TopJet_sel"));

}



bool TTEfficiencyMainSelectionModule::process(Event & event) {
 
  if (version == "TTbar_Incl" ) {
    if(!mttgen_sel->passes(event)) return false;
  }

 //============================================         
  //apply top pt reweighting and scale factors
 //============================================        
  for (auto & rew : reweighting_modules) {
    rew->process(event);
  }
  muo_tight_noniso_SF->process(event);
  muo_trigger_SF->process(event);

  //=============================         
  //run corrections
  //=============================         
  bool ok = common->process(event);
  if(!ok) return false;
  
  topjet_corr->process(event); //apply AK8 corrections on the full jet and AK4 corrections on the subjets OR apply HOTVR corrections
  if(useHOTVR) topjet_groomer->process(event);
  if(isMC && !useHOTVR) topjetJER_smearer->process(event);

  //=============================         
  //run cleaners
  //=============================       
  jet_cleaner->process(event);
  topjet_cleaner->process(event);
 
  //======================================        
  //fill histograms before selection
  //====================================== 
  for(auto & h : hists_before_sel){
    h->fill(event);
  }

  //======================================        
  //apply selections
  //======================================        
  if(!first_selection->passes(event)) return false; 
  if(!twoDcut->passes(event)) return false; 
  if(!met_sel->passes(event)) return false; 
  if(!ptW_sel->passes(event)) return false; 
  hists_btag_medium_eff->fill(event);
  if(!bjetCloseToLepton_sel->passes(event)) return false; 
  
  //====================================================         
  //apply b tagging scale factors after the selection
  //==================================================== 
  btagwAK8->process(event);
  
  //======================================        
  //fill histograms after selection
  //====================================== 
  for(auto & h : hists_after_sel){
    h->fill(event);
  }
   
  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(TTEfficiencyMainSelectionModule)


