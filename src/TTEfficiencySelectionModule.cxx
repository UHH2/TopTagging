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


using namespace std;
using namespace uhh2;
using namespace uhh2examples;

class TTEfficiencySelectionModule: public AnalysisModule {
public:
    
  explicit TTEfficiencySelectionModule(Context & ctx);
  virtual bool process(Event & event) override;

private:
    
  //correctors
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<TopJetCorrectionModules> topjet_corr;
  std::unique_ptr<GenericJetResolutionSmearer> topjetJER_smearer;
  std::unique_ptr<TopJetGroomer> topjet_groomer;
  std::unique_ptr<uhh2::AnalysisModule> pileupRW;
  
  //cleaners
  std::unique_ptr<JetLeptonCleaner> jet_lepton_cleaner;
  std::unique_ptr<TopJetLeptonDeltaRCleaner> topjet_cleaner_dRlep;
  std::unique_ptr<JetCleaner> jet_cleaner;
  std::unique_ptr<TopJetCleaner> topjet_cleaner, topjet_cleaner_400;

  //reweighting and scale factors
  std::vector<std::unique_ptr<AnalysisModule>> reweighting_modules;
  std::unique_ptr<uhh2::AnalysisModule> muo_tight_noniso_SF, muo_trigger_SF;
  //std::unique_ptr<uhh2::AnalysisModule> muo_medium_noniso_SF;
  std::unique_ptr<AnalysisModule> btagwAK8, subjet_btagwAK8;

  //selections
  std::unique_ptr<uhh2::Selection> topjet_sel, htlep_sel, trigger_sel,trigger_sel2,  met_sel, bjetCloseToLepton_sel, twoDcut, mttgen_sel;

  std::unique_ptr<AndSelection> first_selection;

  std::unique_ptr<uhh2::Selection> hadronic_selection, lepton_jets_seletion, dilepton_selection, tau_jets_selection; 
  std::unique_ptr<MergedSelection> merged_selection, mergedW_selection, mergedQB_selection; 
  std::unique_ptr<MassDiffSelection> massDiff_selection;

  //histograms
  std::unique_ptr<BTagMCEfficiencyHists> hists_btag_eff, hists_btag_medium_eff, hists_subjet_btag_eff;
  std::vector<std::unique_ptr<uhh2::Hists>> hists_before_sel;
  std::vector<std::unique_ptr<uhh2::Hists>> hists_after_sel; 
  std::unique_ptr<uhh2::Hists> hists_notrigger, hists_trigger;
  std::unique_ptr<ProbeJetHists> hists_all, hists_all_400, hists_all_400to550, hists_all_550 ,hists_tagged;
  std::unique_ptr<ProbeJetHists> hists_noTag, hists_tau32_wp1, hists_tau32_wp2, hists_tagged_wp1, hists_tagged_wp2;
  std::unique_ptr<ProbeJetHists> hists_merged, hists_mergedW, hists_mergedQB, hists_lightjets, hists_unmerged, hists_dilepton, hists_lepton_jets, hists_taujets, hists_hadronic; 
  std::unique_ptr<ProbeJetHists> hists_merged_400, hists_mergedW_400, hists_mergedQB_400, hists_lightjets_400, hists_unmerged_400, hists_dilepton_400, hists_lepton_jets_400, hists_taujets_400, hists_hadronic_400; 

  std::vector<std::unique_ptr<ProbeJetHists>> h_tau32, h_tau32_mass, h_tau32_btag, h_tau32_mass_btag, h_tau32_mass_btag_pt400to550, h_tau32_mass_btag_pt550, h_tau32_mass_btag_pt400;

  bool useHTT, usePUPPI;
  bool subjet_correction;
  string version;

  bool isMC;

};


TTEfficiencySelectionModule::TTEfficiencySelectionModule(Context & ctx){

    
  cout << "Hello World from TopTaggingStudiesPreSelectionModule_Eff!" << endl;
    
  // If needed, access the configuration of the module here, e.g.:
  string testvalue = ctx.get("TestKey", "<not set>");
  cout << "TestKey in the configuration was: " << testvalue << endl;

  string triggerName = ctx.get("Trigger", "<not set>");

  isMC = (ctx.get("dataset_type") == "MC");
  
  useHTT = (ctx.get("useHTT", "<not set>") == "TRUE");
  usePUPPI = (ctx.get("usePUPPI", "<not set>") == "TRUE");

  if(usePUPPI) cout << "use PUPPI topjets" << endl;
  else cout << "use CHS topjets" << endl;

  if(useHTT) cout << "run the HTT" << endl;
  else cout << "run CMS tagger" << endl;

  version = ctx.get("dataset_version", "");

  subjet_correction = false; 
  string correction_mode = ctx.get("TopJetCorrectionMode", "<not set>");
  if (correction_mode == "SUB") {
    subjet_correction = true; 
    cout << "use AK4 correction on subjets to get the groomed topjet." << endl;
  }



  //MuonId muid = AndId<Muon>(MuonIDMedium_ICHEP(), PtEtaCut(50., 2.1));
  MuonId muid = AndId<Muon>(MuonIDTight(), PtEtaCut(55., 2.4));

  ElectronId eleid = AndId<Electron>(ElectronID_Spring16_medium_noIso, PtEtaCut(55., 2.4));
  //ElectronId eleid = AndId<Electron>(ElectronID_MVAGeneralPurpose_Spring16_tight, PtEtaCut(50., 2.1));

  JetId jetid = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0, 2.4));

  string PU_variation ="central";
  PU_variation = ctx.get("PU_variation","central");


  //common modules
  common.reset(new CommonModules());
  common->set_jet_id(jetid);
  common->set_muon_id(muid);
  common->set_electron_id(eleid);
  common->set_HTjetid(jetid);
  common->switch_jetlepcleaner(true);
  common->switch_jetPtSorter();
  common->switch_metcorrection();
  //common->disable_mcpileupreweight();
  common->init(ctx, PU_variation);
  cout << "common init" <<endl;
 
  if(usePUPPI) topjet_corr.reset(new TopJetCorrectionModules(ctx, TopJetCorrectionModules::AK8_PUPPI));
  else topjet_corr.reset(new TopJetCorrectionModules(ctx, TopJetCorrectionModules::AK8_CHS));

  topjetJER_smearer.reset(new GenericJetResolutionSmearer(ctx, "topjets", "gentopjets", false));

  if(subjet_correction) topjet_groomer.reset(new TopJetGroomer()); // just take the subjet sum (make sure that the subjets are corrected properly)
  else topjet_groomer.reset(new TopJetGroomer(false)); //undo the subjet JEC corrections in case you want to correct the resulting subjet sum


//cleaners
//  jet_lepton_cleaner.reset(new JetLeptonCleaner(ctx, JEC_AK4));
  //jet_lepton_cleaner->set_drmax(.4);
  //jet_cleaner1.reset(new JetCleaner(ctx, 15., 3.0));
  jet_cleaner.reset(new JetCleaner(ctx, 30., 2.4));
  if (useHTT) topjet_cleaner_dRlep.reset(new TopJetLeptonDeltaRCleaner(1.5));
  else topjet_cleaner_dRlep.reset(new TopJetLeptonDeltaRCleaner(0.8));

  if(useHTT) topjet_cleaner.reset(new TopJetCleaner(ctx, TopJetId(PtEtaCut(150., 2.4))));
  else topjet_cleaner.reset(new TopJetCleaner(ctx, TopJetId(PtEtaCut(170., 2.4))));
  topjet_cleaner_400.reset(new TopJetCleaner(ctx, TopJetId(PtEtaCut(400., 2.4))));


  //reweighting and scale factors
  if (version == "TTbar_Incl" || version ==  "TTbar_700to1000" || version ==  "TTbar_1000toInf") {
    reweighting_modules.emplace_back(new TTbarGenProducer(ctx, "ttbargen", true));
    // reweighting_modules.emplace_back(new TopPtReweight(ctx, 0.156, -0.00137, "ttbargen", "weight_ttbar", true)); //8TeV
    // reweighting_modules.emplace_back(new TopPtReweight(ctx, 0.0615, -0.0005, "ttbargen", "weight_ttbar", true)); //13TeV
  }

  if (version == "TTbar_Incl" ) mttgen_sel.reset(new MttbarGenSelection(0., 700.));

  //  btagwAK8.reset(new MCBTagScaleFactor(ctx, CSVBTag::WP_LOOSE, "jets","central","mujets","incl","MCBtagEfficiencies"));
  btagwAK8.reset(new MCBTagScaleFactor(ctx, CSVBTag::WP_MEDIUM, "jets","central","mujets","incl","MCBtagEfficiencies"));

  subjet_btagwAK8.reset(new MCBTagScaleFactor(ctx, CSVBTag::WP_LOOSE, "topjets","central","lt","incl","MCSubjetBtagEfficiencies","", "SubjetBTagCalibration"));//,"SubjetBTagCalibration"));

  muo_tight_noniso_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/dreyert/CMSSW_8_0_24_patch1/src/UHH2/common/data//MuonID_EfficienciesAndSF_average_RunBtoH.root","MC_NUM_TightID_DEN_genTracks_PAR_pt_eta",1, "tightID"));
  muo_trigger_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/dreyert/CMSSW_8_0_24_patch1/src/UHH2/common/data//MuonTrigger_EfficienciesAndSF_average_RunBtoH.root","IsoMu50_OR_IsoTkMu50_PtEtaBins",1, "muonTrigger"));
  //muo_medium_noniso_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/dreyert/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonID_EfficienciesAndSF_average_RunBtoH.root","MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta",1, "mediumID"));
  
 

  //selections 
  topjet_sel.reset(new NTopJetSelection(2));

  trigger_sel = uhh2::make_unique<TriggerSelection>("HLT_Mu50_v*");
  trigger_sel2 = uhh2::make_unique<TriggerSelection>("HLT_TkMu50_v*");

  twoDcut.reset(new TwoDCut(.4, 25.));
  met_sel.reset(new METCut(50., std::numeric_limits<double>::infinity()));
  htlep_sel.reset(new HTlepCut(150., std::numeric_limits<double>::infinity()));
  bjetCloseToLepton_sel.reset(new NMuonBTagSelection(1,999,CSVBTag(CSVBTag::WP_MEDIUM)));
 
  first_selection.reset(new AndSelection(ctx,"first selection"));  
  first_selection->add<NMuonSelection>("Number of muons == 1",1,1);
  first_selection->add<NJetSelection>("Number of jets >= 2",2);
  first_selection->add<NTopJetSelection>("Number of topjets >= 1",1);

  //partonlevel decay channel slections for 
  hadronic_selection.reset( new DecayChannelSelection(ctx, "ttbargen", "hadronic"));
  lepton_jets_seletion.reset( new DecayChannelSelection(ctx, "ttbargen", "lepton_jets"));
  dilepton_selection.reset( new DecayChannelSelection(ctx, "ttbargen", "dilepton"));
  tau_jets_selection.reset( new DecayChannelSelection(ctx, "ttbargen", "tau_jets"));


  double jet_radius = 0.;
  if (useHTT) jet_radius = 1.5; 
  else jet_radius = 0.8; 

  //merged selections
  merged_selection.reset( new MergedSelection(ctx, "ttbargen", jet_radius));
  mergedW_selection.reset( new MergedSelection(ctx, "ttbargen", jet_radius, MergedSelection::oMergedW));
  mergedQB_selection.reset( new MergedSelection(ctx, "ttbargen", jet_radius, MergedSelection::oBplusQ));

  //mass diff selection
   massDiff_selection.reset( new MassDiffSelection(ctx) );
                  
  //histograms
  hists_btag_eff.reset(new BTagMCEfficiencyHists(ctx,"BTagLoose",CSVBTag::WP_LOOSE));
  hists_btag_medium_eff.reset(new BTagMCEfficiencyHists(ctx,"BTagMedium",CSVBTag::WP_MEDIUM));
  hists_subjet_btag_eff.reset(new BTagMCEfficiencyHists(ctx,"SubjetBTag",CSVBTag::WP_LOOSE, "topjets") );


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

  TString name = "ProbeJet";
  for(int i = 1; i <= 9; i++){
     h_tau32.emplace_back(new ProbeJetHists(ctx, (name+"_tau32_wp"+TString::Format("%1.0i",i)).Data() ) );
     h_tau32_mass.emplace_back(new ProbeJetHists(ctx, (name+"_tau32_mass_wp"+TString::Format("%1.0i",i)).Data() ) );
     h_tau32_btag.emplace_back(new ProbeJetHists(ctx, (name+"_tau32_btag_wp"+TString::Format("%1.0i",i)).Data() ) );
     h_tau32_mass_btag.emplace_back(new ProbeJetHists(ctx, (name+"_tau32_mass_btag_wp"+TString::Format("%1.0i",i)).Data() ) );
     h_tau32_mass_btag_pt400.emplace_back(new ProbeJetHists(ctx, (name+"_tau32_mass_btag_pt400_wp"+TString::Format("%1.0i",i)).Data() ) );
     h_tau32_mass_btag_pt400to550.emplace_back(new ProbeJetHists(ctx, (name+"_tau32_mass_btag_pt400to550_wp"+TString::Format("%1.0i",i)).Data() ) );
     h_tau32_mass_btag_pt550.emplace_back(new ProbeJetHists(ctx, (name+"_tau32_mass_btag_pt550_wp"+TString::Format("%1.0i",i)).Data() ) );
  }

  hists_notrigger.reset(new EventHists(ctx, "Event_noTrigger"));
  hists_trigger.reset(new EventHists(ctx, "Event_Trigger"));

  hists_all.reset(new ProbeJetHists(ctx, "ProbeJet_All"));
  hists_all_400.reset(new ProbeJetHists(ctx, "ProbeJet_All_Pt400"));
  hists_all_400to550.reset(new ProbeJetHists(ctx, "ProbeJet_All_Pt400to550"));
  hists_all_550.reset(new ProbeJetHists(ctx, "ProbeJet_All_Pt550"));

  hists_tagged.reset(new ProbeJetHists(ctx, "ProbeJet_Tagged"));

  hists_noTag.reset(new ProbeJetHists(ctx, "ProbeJet_noTag"));

  hists_tagged_wp1.reset(new ProbeJetHists(ctx, "ProbeJet_tau32_wp1_mass"));
  hists_tagged_wp2.reset(new ProbeJetHists(ctx, "ProbeJet_tau32_wp2_mass"));

  hists_merged.reset(new ProbeJetHists(ctx, "ProbeJet_all_merged"));
  hists_mergedW.reset(new ProbeJetHists(ctx, "ProbeJet_all_mergedW"));
  hists_mergedQB.reset(new ProbeJetHists(ctx, "ProbeJet_all_mergedQB"));
  hists_lightjets.reset(new ProbeJetHists(ctx, "ProbeJet_all_lightjets"));
  hists_unmerged.reset(new ProbeJetHists(ctx, "ProbeJet_all_unmerged"));
  hists_dilepton.reset(new ProbeJetHists(ctx, "ProbeJet_all_dilepton"));
  hists_lepton_jets.reset(new ProbeJetHists(ctx, "ProbeJet_all_lepton_jets"));
  hists_taujets.reset(new ProbeJetHists(ctx, "ProbeJet_all_tau_jets"));
  hists_hadronic.reset(new ProbeJetHists(ctx, "ProbeJet_all_hadronic"));

  hists_merged_400.reset(new ProbeJetHists(ctx, "ProbeJet_all_merged_Pt400"));
  hists_mergedW_400.reset(new ProbeJetHists(ctx, "ProbeJet_all_mergedW_Pt400"));
  hists_mergedQB_400.reset(new ProbeJetHists(ctx, "ProbeJet_all_mergedQB_Pt400"));
  hists_lightjets_400.reset(new ProbeJetHists(ctx, "ProbeJet_all_lightjets_Pt400"));
  hists_unmerged_400.reset(new ProbeJetHists(ctx, "ProbeJet_all_unmerged_Pt400"));
  hists_dilepton_400.reset(new ProbeJetHists(ctx, "ProbeJet_all_dilepton_Pt400"));
  hists_lepton_jets_400.reset(new ProbeJetHists(ctx, "ProbeJet_all_lepton_jets_Pt400"));
  hists_taujets_400.reset(new ProbeJetHists(ctx, "ProbeJet_all_tau_jets_Pt400"));
  hists_hadronic_400.reset(new ProbeJetHists(ctx, "ProbeJet_all_hadronic_Pt400"));

}



bool TTEfficiencySelectionModule::process(Event & event) {
   
  //  cout << "TTEfficiencySelectionModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;

  if(version == "TTbar_Incl" && event.run == 1 && event.event == 67315776) return false;
  if(version == "TTbar_1000toInf" && event.run == 1 && event.event == 185748260) return false;
  /*
  cout <<"=============================================" << endl;
  std::vector<GenParticle>* gens = event.genparticles;
  for( auto & gen : *gens){
    cout << "index: " << gen.index() << "  id: " << gen.pdgId() << "  status: "<< gen.status()  << "  mother1: " << gen.mother1() << "  daughter1: " << gen.daughter1() << "  daughter2: " << gen.daughter2() << endl;
  }
  */
  
  

 //apply top pt reweighting
  for (auto & rew : reweighting_modules) {
    rew->process(event);
  }
  if (version == "TTbar_Incl" ) {
    if(!mttgen_sel->passes(event)) return false;
  }

   //run corrections
  bool ok = common->process(event);
  if(!ok) return false;
  
  if(subjet_correction) {
    topjet_corr->process(event); //apply topjet and subjet corrections first 
    topjet_groomer->process(event); //now get the subjet sum from corrected subjets
  }else{
    topjet_groomer->process(event); //get the subjet sum from uncorrected subjets (sets the JEC_Raw for the topjet to 1)
    topjet_corr->process(event); //now apply corrections to the resulting groomed jet
  }

  // topjetJER_smearer->process(event);
 
  //run jet cleaners 
  // jet_lepton_cleaner->process(event);
 
  jet_cleaner->process(event);
  topjet_cleaner->process(event);
  topjet_cleaner_dRlep->process(event); //remove topjets that overlap with a lepton (e/mu)
  

   // muo_medium_noniso_SF->process(event);
  muo_tight_noniso_SF->process(event);
  
  //fill histograms before selection
  for(auto & h : hists_before_sel){
    h->fill(event);
  }
  
  //apply selections
  if( !isMC && event.run < 274954) {
    if(!trigger_sel->passes(event)) return false;
  }else{
    if( !(trigger_sel->passes(event) || trigger_sel2->passes(event)) ) return false;
  }
  if(!first_selection->passes(event)) return false; 
  if(!twoDcut->passes(event)) return false; 
  if(!met_sel->passes(event)) return false; 
  if(!htlep_sel->passes(event)) return false; 
  hists_btag_eff->fill(event);
  hists_btag_medium_eff->fill(event);
  if(!bjetCloseToLepton_sel->passes(event)) return false; 
  
  
  //apply b tagging scale factors after the selection
  muo_trigger_SF->process(event);
  btagwAK8->process(event);
  
  for(auto & h : hists_after_sel){
    h->fill(event);
  }
  
  //get the probe jet
  std::vector<TopJet>* topjets = event.topjets;
  std::vector<Muon>* muons = event.muons;
  TopJet probe_jet; 
  bool probejet_found = false;

  for(auto & topjet : *topjets){
    double pi = 3.14159265359;
    double delPhi = deltaPhi(topjet, muons->at(0));
    if(delPhi > ((2/3)*pi)){
      probe_jet = topjet;
      probejet_found = true;
      break;
    }
  }

  if(!probejet_found){
    //   cout << "WARNING: no probe jet found"<< endl;
    return false;
  }

 
  //topjet_corrector->process(event);
  //=================================
  //define working points and tags
  //=================================
  double probejet_tau32 = probe_jet.tau3()/probe_jet.tau2();

  //get a subjet b tag on the probe jet and the probe jet mass from the subjets
  bool subjet_btag = false;
  LorentzVector subjet_sum;
  LorentzVector subjet_sum_uncorr;

  for(const auto & subjet : probe_jet.subjets()){
    JetId btag = CSVBTag(CSVBTag::WP_LOOSE);
    if( btag(subjet, event) ) subjet_btag = true;
  }

  /*
    for(const auto & subjet : probe_jet.subjets()){
    //  subjet_sum += subjet.v4();

    if(subjet.btag_combinedSecondaryVertex() > 0.5426) 
    subjet_btag = true;
    //  cout << "subjet JEC_23: " << 1/subjet.JEC_factor_raw() << endl;
    }
  */

  double probejet_mass = probe_jet.v4().M();
  double probejet_pt = probe_jet.pt();

  bool mass_cut = false;
  if(usePUPPI && probejet_mass > 105 && probejet_mass < 210) mass_cut = true;
  else if(!usePUPPI && probejet_mass > 105 && probejet_mass < 220) mass_cut = true;

  std::vector<double> tau32_wps_CHS           {0.45, 0.51, 0.60, 0.69, 0.84, 0.50, 0.57, 0.67, 0.81};
  std::vector<double> tau32_wps_PUPPI         {0.41, 0.48, 0.57, 0.67, 0.84, 0.46, 0.54, 0.65, 0.80};

  std::vector<double> tau32_wps_HTT_CHS       {0.49, 0.54, 0.65, 0.97, 0.98, 0.55, 0.62, 0.93, 0.97};
  std::vector<double> fRec_wps_CHS            {0.14, 0.20, 0.25, 0.22, 0.50, 0.17, 0.27, 0.20, 0.47};
  std::vector<double> tau32_wps_HTT_PUPPI     {0.39, 0.48, 0.60, 0.96, 0.98, 0.50, 0.60, 0.99, 0.97};
  std::vector<double> fRec_wps_PUPPI          {0.36, 0.43, 0.22, 0.25, 0.48, 0.16, 0.25, 0.22, 0.46};

  double fRec = 0;
  if(useHTT) fRec = probe_jet.get_tag( probe_jet.tagname2tag("fRec"));

  std::vector<double> tau32_wps; 
  std::vector<double> fRec_wps; 
  if(!useHTT && !usePUPPI) tau32_wps = tau32_wps_CHS;
  else if(!useHTT && usePUPPI) tau32_wps = tau32_wps_PUPPI;
  else if(useHTT && !usePUPPI) {
    tau32_wps = tau32_wps_HTT_CHS; 
    fRec_wps = fRec_wps_CHS;
  }
  else if(useHTT && usePUPPI) {
    tau32_wps = tau32_wps_HTT_PUPPI; 
    fRec_wps = fRec_wps_PUPPI;
  }

  if(probe_jet.v4().M()<10) return false; //mild mass cut on the denominator 
  //if(!massDiff_selection->passes_probe(event, probe_jet)) return false;
  //===================
  //fill the histograms
  //===================
  hists_all->fill_probe(event, probe_jet);
  hists_subjet_btag_eff->fill(event);

  if(probejet_pt > 400) hists_all_400->fill_probe(event, probe_jet);
  if(probejet_pt > 400 && probejet_pt < 550) hists_all_400to550->fill_probe(event, probe_jet);
  if(probejet_pt > 550) hists_all_550->fill_probe(event, probe_jet);


  if (version == "TTbar_Incl" || version ==  "TTbar_700to1000" || version ==  "TTbar_1000toInf") {
    bool fullymerged = merged_selection->passes_probe(event, probe_jet);
    bool mergedW = mergedW_selection->passes_probe(event, probe_jet);;
    bool mergedQB = mergedQB_selection->passes_probe(event, probe_jet);; 
    if( fullymerged ){
      hists_merged->fill_probe(event, probe_jet);
      if(probejet_pt > 400) hists_merged_400->fill_probe(event, probe_jet);
    }else{
      hists_unmerged->fill_probe(event, probe_jet); 
      if(probejet_pt > 400) hists_unmerged_400->fill_probe(event, probe_jet);
    }
    if (mergedW) {
      hists_mergedW->fill_probe(event, probe_jet);
      if(probejet_pt > 400)  hists_mergedW_400->fill_probe(event, probe_jet);
    }
    if (mergedQB) {
      hists_mergedQB->fill_probe(event, probe_jet);
      if(probejet_pt > 400) hists_mergedQB_400->fill_probe(event, probe_jet);
    }
    if (!fullymerged && !mergedW && !mergedQB){
      hists_lightjets->fill_probe(event, probe_jet);
      if(probejet_pt > 400) hists_lightjets_400->fill_probe(event, probe_jet);
    }
    if( tau_jets_selection->passes(event) ) {
      hists_taujets->fill_probe(event, probe_jet); 
      if(probejet_pt > 400) hists_taujets_400->fill_probe(event, probe_jet);
    }
    if( lepton_jets_seletion->passes(event) ) {
      hists_lepton_jets->fill_probe(event, probe_jet);
      if(probejet_pt > 400) hists_lepton_jets_400->fill_probe(event, probe_jet);
    } 
    if( dilepton_selection->passes(event) ) {
      hists_dilepton->fill_probe(event, probe_jet);
      if(probejet_pt > 400) hists_dilepton_400->fill_probe(event, probe_jet);
    }
    if( hadronic_selection->passes(event) ) {
      hists_hadronic->fill_probe(event, probe_jet); 
      if(probejet_pt > 400)  hists_hadronic_400->fill_probe(event, probe_jet);
    }
  }

  for(unsigned int i = 0; i < tau32_wps.size(); i++){
    if(probejet_tau32 < tau32_wps.at(i) && (!useHTT || fRec < fRec_wps.at(i))){
      if(i < h_tau32.size()) h_tau32.at(i)->fill_probe(event, probe_jet);
      else cout << "More working points than histogram collections (" << tau32_wps.size() << " wps and " << h_tau32.size() << " tau32 collections)" << endl;
      if(mass_cut){
	if(i < h_tau32_mass.size()) h_tau32_mass.at(i)->fill_probe(event, probe_jet);
	else cout << "More working points than histogram collections (" << tau32_wps.size() << " wps and " << h_tau32_mass.size() << " tau32_mass collections)" << endl;
      }
    }
  }
  
  subjet_btagwAK8->process(event);
  if(subjet_btag){

    for(unsigned int i = 0; i < tau32_wps.size(); i++){
      if(probejet_tau32 < tau32_wps.at(i) && (!useHTT || fRec < fRec_wps.at(i))){
    
	if(i < h_tau32_btag.size()) h_tau32_btag.at(i)->fill_probe(event, probe_jet);
	else cout << "More working points than histogram collections (" << tau32_wps.size() << " wps and " << h_tau32_btag.size() << " tau32btag collections)" << endl;
  
	if(mass_cut){
	  if(i < h_tau32_mass_btag.size()) h_tau32_mass_btag.at(i)->fill_probe(event, probe_jet);
	  else cout << "More working points than histogram collections (" << tau32_wps.size() << " wps and " << h_tau32_mass_btag.size() << " tau32_mass_btag collections)" << endl;
	  if(probejet_pt > 400) {
	    if(i < h_tau32_mass_btag_pt400.size()) h_tau32_mass_btag_pt400.at(i)->fill_probe(event, probe_jet);
	    else cout << "More working points than histogram collections (" << tau32_wps.size() << " wps and " << h_tau32_mass_btag_pt400.size() << " tau32_mass_btag_pt400 collections)" << endl;
	  }
	  if(probejet_pt > 400 && probejet_pt < 550) {
	    if(i < h_tau32_mass_btag_pt400to550.size()) h_tau32_mass_btag_pt400to550.at(i)->fill_probe(event, probe_jet);
	    else cout << "More working points than histogram collections (" << tau32_wps.size() << " wps and " << h_tau32_mass_btag_pt400to550.size() << " tau32_mass_btag_pt400to550 collections)" << endl;
	  }
	  if(probejet_pt > 550) {
	    if(i < h_tau32_mass_btag_pt550.size()) h_tau32_mass_btag_pt550.at(i)->fill_probe(event, probe_jet);
	    else cout << "More working points than histogram collections (" << tau32_wps.size() << " wps and " << h_tau32_mass_btag_pt550.size() << " tau32_mass_btag_550 collections)" << endl;
	  }
	}
	
      }

    }
  }
     
   
  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(TTEfficiencySelectionModule)


