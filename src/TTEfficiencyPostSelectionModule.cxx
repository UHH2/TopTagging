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

class TTEfficiencyPostSelectionModule: public AnalysisModule {
public:
    
  explicit TTEfficiencyPostSelectionModule(Context & ctx);
  virtual bool process(Event & event) override;

private:
    
  //correctors
  std::unique_ptr<CommonModules> common;
  std::unique_ptr<uhh2::AnalysisModule> pileupRW;

  //reweighting and scale factors
  std::vector<std::unique_ptr<AnalysisModule>> reweighting_modules;
  std::unique_ptr<uhh2::AnalysisModule> muo_tight_noniso_SF, muo_trigger_SF;
  std::unique_ptr<AnalysisModule> btagwAK8, subjet_btagwAK8;
  std::unique_ptr<uhh2::AnalysisModule> scale_variation;
  std::unique_ptr<GenericJetResolutionSmearer> topjetJER_smearer;

  //selections
  std::unique_ptr<uhh2::Selection> hadronic_selection, lepton_jets_seletion, dilepton_selection, tau_jets_selection; 
  std::unique_ptr<MergedSelection> merged_selection, mergedW_selection, mergedQB_selection, mergedEvent_selection; 
  std::unique_ptr<MassDiffSelection> massDiff_selection;
  std::unique_ptr<DPhiMuBSelection> dphi_selection;
  std::unique_ptr<LeadingAddJetSelection> addJet_selection;

  //histograms
  std::unique_ptr<BTagMCEfficiencyHists> hists_btag_eff, hists_btag_medium_eff, hists_subjet_btag_eff;
  std::vector<std::unique_ptr<uhh2::Hists>> hists_before_sel;
  std::vector<std::unique_ptr<uhh2::Hists>> hists_after_sel; 
  std::unique_ptr<uhh2::Hists> hists_notrigger, hists_trigger;
  std::unique_ptr<ProbeJetHists> hists_all, hists_all_400, hists_all_400to550, hists_all_550 ,hists_tagged;
  std::unique_ptr<ProbeJetHists>  hists_dilepton, hists_lepton_jets, hists_taujets, hists_hadronic, hists_dilepton_400, hists_lepton_jets_400, hists_taujets_400, hists_hadronic_400; 

  std::vector<std::unique_ptr<ProbeJetHists>> h_tau32, h_tau32_mass, h_tau32_btag, h_tau32_mass_btag, h_tau32_mass_btag_pt400to550, h_tau32_mass_btag_pt550, h_tau32_mass_btag_pt400, h_tau32_btag_pt400;

  std::vector<std::vector<std::unique_ptr<ProbeJetHists>>> h_probe_all_pass, h_probe_mass_pass, h_probe_btag_pass, h_probe_mass_btag_pass;
  std::vector<std::vector<std::unique_ptr<ProbeJetHists>>> h_probe_all_fail, h_probe_mass_fail, h_probe_btag_fail, h_probe_mass_btag_fail;

  //std::vector<std::vector<std::unique_ptr<ProbeJetHists>>>  h_probe_btag_NoSF_pass, h_probe_btag_NoSF_fail;

  //std::vector<std::vector<std::unique_ptr<ProbeJetHists>>> h_probeNPV_all_pass, h_probeNPV_mass_pass, h_probeNPV_btag_pass, h_probeNPV_mass_btag_pass;
  //std::vector<std::vector<std::unique_ptr<ProbeJetHists>>> h_probeNPV_all_fail, h_probeNPV_mass_fail, h_probeNPV_btag_fail, h_probeNPV_mass_btag_fail;

 
  //variables and functions
  string mass_scale;
  bool useHTT, usePUPPI;
  string version;
  bool fill_PDF;
  bool isMC;
  bool invert_merged_selection;
  bool merged_category;

  const std::vector<double> pt_bins{-1,-400 , 300, 400, 480, 600};
  const std::vector<int> npv_bins{0,10,15,20,25};

  bool get_pt_cut(unsigned int bin, const double pt);
  bool get_npv_cut(unsigned int bin, const int npv);
  bool get_tau32_cut(unsigned int bin, const double tau32, const std::vector<double> wps); 
};


TTEfficiencyPostSelectionModule::TTEfficiencyPostSelectionModule(Context & ctx){

  //=============================
  // access configuaration values
  //=============================

  string testvalue = ctx.get("TestKey", "<not set>");
  cout << "TestKey in the configuration was: " << testvalue << endl;

  string triggerName = ctx.get("Trigger", "<not set>");

  isMC = (ctx.get("dataset_type") == "MC");
  
  useHTT = (ctx.get("useHTT", "<not set>") == "TRUE");
  usePUPPI = (ctx.get("usePUPPI", "<not set>") == "TRUE");

  fill_PDF = (ctx.get("fill_PDF", "FALSE") == "TRUE");

  if(usePUPPI) cout << "use PUPPI topjets" << endl;
  else cout << "use CHS topjets" << endl;

  if(useHTT) cout << "run the HTT" << endl;
  else cout << "run CMS tagger" << endl;

  version = ctx.get("dataset_version", "");

  string merged = ctx.get("MergedSelection", "<not set>");
  TString vers = (TString)version;
  if(vers.Contains("mergedTop")) merged = "mergedTop";
  if(vers.Contains("mergedW")) merged = "mergedW";
  if(vers.Contains("mergedQB")) merged = "mergedQB";
  if(vers.Contains("light")) merged = "light";
  if(vers.Contains("bkg")) merged = "bkg";
  if(vers.Contains("notmerged")) merged = "notmerged";
  if(vers.Contains("semimerged")) merged = "semimerged";

  JetId jetid = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0, 2.4));

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

  bool TopPtReweighting = false;
  TopPtReweighting = (ctx.get("TopPtReweight","FALSE")== "TRUE");

  mass_scale = "default";
  mass_scale = ctx.get("mass_scale","default");


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

  // if (version == "TTbar_Incl" || version ==  "TTbar_700to1000" || version ==  "TTbar_1000toInf") {
  if ( vers.Contains("TTbar") ) {
    reweighting_modules.emplace_back(new TTbarGenProducer(ctx, "ttbargen", true));
    // reweighting_modules.emplace_back(new TopPtReweight(ctx, 0.156, -0.00137, "ttbargen", "weight_ttbar", true)); //8TeV
    if(TopPtReweighting) reweighting_modules.emplace_back(new TopPtReweight(ctx, 0.0615, -0.0005, "ttbargen", "weight_ttbar", true)); //13TeV
  }

  //  btagwAK8.reset(new MCBTagScaleFactor(ctx, CSVBTag::WP_LOOSE, "jets","central","mujets","incl","MCBtagEfficiencies"));
  btagwAK8.reset(new MCBTagScaleFactor(ctx, CSVBTag::WP_MEDIUM, "jets", BTag_variation,"mujets","incl","MCBtagEfficiencies"));

  subjet_btagwAK8.reset(new MCBTagScaleFactor(ctx, CSVBTag::WP_LOOSE, "topjets", SubjetBTag_variation,"lt","incl","MCSubjetBtagEfficiencies","", "SubjetBTagCalibration"));//,"SubjetBTagCalibration"));

  muo_tight_noniso_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/dreyert/CMSSW_8_0_24_patch1/src/UHH2/common/data//MuonID_EfficienciesAndSF_average_RunBtoH.root","MC_NUM_TightID_DEN_genTracks_PAR_pt_eta",1, "tightID", true, MuonID_variation));
  muo_trigger_SF.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/dreyert/CMSSW_8_0_24_patch1/src/UHH2/common/data//MuonTrigger_EfficienciesAndSF_average_RunBtoH.root","IsoMu50_OR_IsoTkMu50_PtEtaBins",1, "muonTrigger", true, MuonTrigger_variation));
 
  scale_variation.reset(new MCScaleVariation(ctx));


  //==========================================
  //partonlevel decay channel slections 
  //==========================================
  hadronic_selection.reset( new DecayChannelSelection(ctx, "ttbargen", "hadronic"));
  lepton_jets_seletion.reset( new DecayChannelSelection(ctx, "ttbargen", "lepton_jets"));
  dilepton_selection.reset( new DecayChannelSelection(ctx, "ttbargen", "dilepton"));
  tau_jets_selection.reset( new DecayChannelSelection(ctx, "ttbargen", "tau_jets"));

  //---------
  //matching
  //---------
  double jet_radius = 0.;
  if (useHTT) jet_radius = 1.5; 
  else jet_radius = 0.8; 

  merged_selection.reset( new MergedSelection(ctx, "ttbargen", jet_radius));
  mergedW_selection.reset( new MergedSelection(ctx, "ttbargen", jet_radius, MergedSelection::oMergedW));
  mergedQB_selection.reset( new MergedSelection(ctx, "ttbargen", jet_radius, MergedSelection::oBplusQ));
 
  merged_category = true;
  if (merged == "mergedTop" || merged == "notmergedTop") mergedEvent_selection.reset( new MergedSelection(ctx, "ttbargen", jet_radius));
  else if (merged == "mergedW") mergedEvent_selection.reset( new MergedSelection(ctx, "ttbargen", jet_radius, MergedSelection::oMergedW));
  else if (merged == "mergedQB") mergedEvent_selection.reset( new MergedSelection(ctx, "ttbargen", jet_radius, MergedSelection::oBplusQ));
  else if (merged == "light") mergedEvent_selection.reset( new MergedSelection(ctx, "ttbargen", jet_radius, MergedSelection::oLight));
  else if (merged == "bkg") mergedEvent_selection.reset( new MergedSelection(ctx, "ttbargen", jet_radius, MergedSelection::oBkg));
  else if (merged == "notmerged") mergedEvent_selection.reset( new MergedSelection(ctx, "ttbargen", jet_radius, MergedSelection::oNotMerged));
  else if (merged == "semimerged") mergedEvent_selection.reset( new MergedSelection(ctx, "ttbargen", jet_radius, MergedSelection::oSemiMerged));
  else merged_category = false;

  invert_merged_selection = false;
  if (merged == "notmergedTop") {invert_merged_selection = true; merged_category = false;}

  if (!(vers.Contains("TTbar"))){
    merged_category = false;
    invert_merged_selection = false;
  }


  //===========================
  //additional selections
  //===========================

  massDiff_selection.reset( new MassDiffSelection(ctx) );
  dphi_selection.reset(new DPhiMuBSelection(ctx, CSVBTag(CSVBTag::WP_MEDIUM), 1.2));
  addJet_selection.reset(new LeadingAddJetSelection(ctx, CSVBTag(CSVBTag::WP_MEDIUM), 50));


  //===========================
  //histograms
  //===========================

  hists_subjet_btag_eff.reset(new BTagMCEfficiencyHists(ctx,"SubjetBTag",CSVBTag::WP_LOOSE, "topjets") );

  hists_after_sel.emplace_back(new EventHists(ctx, "Event_sel"));
  hists_after_sel.emplace_back(new MuonHists(ctx, "Muon_sel"));
  hists_after_sel.emplace_back(new ElectronHists(ctx, "Electron_sel"));
  hists_after_sel.emplace_back(new JetHists(ctx, "Jet_sel"));
  hists_after_sel.emplace_back(new TopJetHists(ctx, "TopJet_sel"));

  std::vector<TString> wps;
  std::vector<TString> wps_PUPPI{"", "_wp1", "_wp2", "_wp3", "_wp4","_wp5"};
  std::vector<TString> wps_CHS{"", "_wp1", "_wp2", "_wp3", "_wp4","_wp5"};

  if(usePUPPI) wps = wps_PUPPI;
  else wps = wps_CHS;

  TString name = "ProbeJet";
  for(int i = 1; i <= 4; i++){
    h_tau32.emplace_back(new ProbeJetHists(ctx, (name+"_tau32_wp"+TString::Format("%1.0i",i)).Data() ) );
    h_tau32_mass.emplace_back(new ProbeJetHists(ctx, (name+"_tau32_mass_wp"+TString::Format("%1.0i",i)).Data() ) );
    h_tau32_btag.emplace_back(new ProbeJetHists(ctx, (name+"_tau32_btag_wp"+TString::Format("%1.0i",i)).Data() ) );
    h_tau32_mass_btag.emplace_back(new ProbeJetHists(ctx, (name+"_tau32_mass_btag_wp"+TString::Format("%1.0i",i)).Data() ) );
    h_tau32_btag_pt400.emplace_back(new ProbeJetHists(ctx, (name+"_tau32_btag_pt400_wp"+TString::Format("%1.0i",i)).Data() ) );
    h_tau32_mass_btag_pt400.emplace_back(new ProbeJetHists(ctx, (name+"_tau32_mass_btag_pt400_wp"+TString::Format("%1.0i",i)).Data() ) );
    h_tau32_mass_btag_pt400to550.emplace_back(new ProbeJetHists(ctx, (name+"_tau32_mass_btag_pt400to550_wp"+TString::Format("%1.0i",i)).Data() ) );
    h_tau32_mass_btag_pt550.emplace_back(new ProbeJetHists(ctx, (name+"_tau32_mass_btag_pt550_wp"+TString::Format("%1.0i",i)).Data() ) );
  }

  h_probe_all_pass.resize(pt_bins.size());
  h_probe_mass_pass.resize(pt_bins.size());
  h_probe_btag_pass.resize(pt_bins.size());
  h_probe_mass_btag_pass.resize(pt_bins.size()); 
  // h_probe_btag_NoSF_pass.resize(pt_bins.size()); 

  h_probe_all_fail.resize(pt_bins.size());
  h_probe_mass_fail.resize(pt_bins.size());
  h_probe_btag_fail.resize(pt_bins.size());
  h_probe_mass_btag_fail.resize(pt_bins.size());
  // h_probe_btag_NoSF_fail.resize(pt_bins.size()); 
  
  for(unsigned int bin = 0; bin < pt_bins.size(); ++bin){   
    TString ptString = "";
    if(pt_bins.at(bin) == -1) ptString = "ptInclusive";
    else if(pt_bins.at(bin) >= 0 && bin+1 < pt_bins.size() && pt_bins.at(bin+1) >= 0) ptString = "pt"+TString::Format("%3.0f",pt_bins.at(bin))+"to"+TString::Format("%3.0f",pt_bins.at(bin+1));
    else ptString = "pt"+TString::Format("%3.0f",fabs(pt_bins.at(bin)));

    for(unsigned int b = 0; b < wps.size(); ++b){
      h_probe_all_pass.at(bin).emplace_back(new ProbeJetHists(ctx, (name+"_"+ptString+wps.at(b)+"_all_pass").Data()));
      h_probe_mass_pass.at(bin).emplace_back(new ProbeJetHists(ctx, (name+"_"+ptString+wps.at(b)+"_mass_pass").Data()));
      h_probe_btag_pass.at(bin).emplace_back(new ProbeJetHists(ctx,  (name+"_"+ptString+wps.at(b)+"_btag_pass").Data()));
      h_probe_mass_btag_pass.at(bin).emplace_back(new ProbeJetHists(ctx,  (name+"_"+ptString+wps.at(b)+"_mass_btag_pass").Data()));
      //    h_probe_btag_NoSF_pass.at(bin).emplace_back(new ProbeJetHists(ctx,  (name+"_"+ptString+wps.at(b)+"_mass_btag_NoSF_pass").Data()));

      h_probe_all_fail.at(bin).emplace_back(new ProbeJetHists(ctx, (name+"_"+ptString+wps.at(b)+"_all_fail").Data()));
      h_probe_mass_fail.at(bin).emplace_back(new ProbeJetHists(ctx, (name+"_"+ptString+wps.at(b)+"_mass_fail").Data()));
      h_probe_btag_fail.at(bin).emplace_back(new ProbeJetHists(ctx,  (name+"_"+ptString+wps.at(b)+"_btag_fail").Data()));
      h_probe_mass_btag_fail.at(bin).emplace_back(new ProbeJetHists(ctx,  (name+"_"+ptString+wps.at(b)+"_mass_btag_fail").Data()));
      //  h_probe_btag_NoSF_fail.at(bin).emplace_back(new ProbeJetHists(ctx,  (name+"_"+ptString+wps.at(b)+"_mass_btag_NoSF_fail").Data()));
    }  
  }
  /*
  h_probeNPV_all_pass.resize(npv_bins.size());
  h_probeNPV_mass_pass.resize(npv_bins.size());
  h_probeNPV_btag_pass.resize(npv_bins.size());
  h_probeNPV_mass_btag_pass.resize(npv_bins.size()); 

  h_probeNPV_all_fail.resize(npv_bins.size());
  h_probeNPV_mass_fail.resize(npv_bins.size());
  h_probeNPV_btag_fail.resize(npv_bins.size());
  h_probeNPV_mass_btag_fail.resize(npv_bins.size()); 

  for(unsigned int bin = 0; bin < npv_bins.size(); ++bin){   
    TString npvString = "";
    if(npv_bins.at(bin) >= 0 && bin+1 < npv_bins.size() && npv_bins.at(bin+1) >= 0) npvString = "npv"+TString::Format("%2.0i",npv_bins.at(bin))+"to"+TString::Format("%2.0i",npv_bins.at(bin+1));
    else npvString = "npv"+TString::Format("%2.0i",abs(npv_bins.at(bin)));

    for(unsigned int b = 0; b < wps.size(); ++b){
      h_probeNPV_all_pass.at(bin).emplace_back(new ProbeJetHists(ctx, (name+"_"+npvString+wps.at(b)+"_all_pass").Data()));
      h_probeNPV_mass_pass.at(bin).emplace_back(new ProbeJetHists(ctx, (name+"_"+npvString+wps.at(b)+"_mass_pass").Data()));
      h_probeNPV_btag_pass.at(bin).emplace_back(new ProbeJetHists(ctx,  (name+"_"+npvString+wps.at(b)+"_btag_pass").Data()));
      h_probeNPV_mass_btag_pass.at(bin).emplace_back(new ProbeJetHists(ctx,  (name+"_"+npvString+wps.at(b)+"_mass_btag_pass").Data()));

      h_probeNPV_all_fail.at(bin).emplace_back(new ProbeJetHists(ctx, (name+"_"+npvString+wps.at(b)+"_all_fail").Data()));
      h_probeNPV_mass_fail.at(bin).emplace_back(new ProbeJetHists(ctx, (name+"_"+npvString+wps.at(b)+"_mass_fail").Data()));
      h_probeNPV_btag_fail.at(bin).emplace_back(new ProbeJetHists(ctx,  (name+"_"+npvString+wps.at(b)+"_btag_fail").Data()));
      h_probeNPV_mass_btag_fail.at(bin).emplace_back(new ProbeJetHists(ctx,  (name+"_"+npvString+wps.at(b)+"_mass_btag_fail").Data()));
    }  
  }
  */
  hists_all.reset(new ProbeJetHists(ctx, "ProbeJet_All"));
  hists_all_400.reset(new ProbeJetHists(ctx, "ProbeJet_All_Pt400"));
  hists_all_400to550.reset(new ProbeJetHists(ctx, "ProbeJet_All_Pt400to550"));
  hists_all_550.reset(new ProbeJetHists(ctx, "ProbeJet_All_Pt550"));

  hists_dilepton.reset(new ProbeJetHists(ctx, "ProbeJet_all_dilepton"));
  hists_lepton_jets.reset(new ProbeJetHists(ctx, "ProbeJet_all_lepton_jets"));
  hists_taujets.reset(new ProbeJetHists(ctx, "ProbeJet_all_tau_jets"));
  hists_hadronic.reset(new ProbeJetHists(ctx, "ProbeJet_all_hadronic"));

  hists_dilepton_400.reset(new ProbeJetHists(ctx, "ProbeJet_all_dilepton_Pt400"));
  hists_lepton_jets_400.reset(new ProbeJetHists(ctx, "ProbeJet_all_lepton_jets_Pt400"));
  hists_taujets_400.reset(new ProbeJetHists(ctx, "ProbeJet_all_tau_jets_Pt400"));
  hists_hadronic_400.reset(new ProbeJetHists(ctx, "ProbeJet_all_hadronic_Pt400"));
}



bool TTEfficiencyPostSelectionModule::process(Event & event) {
  
  //========================================
  //corrections, reweighting, scale factors
  //========================================

  //run corrections
  bool ok = common->process(event);
  if(!ok) return false;

  //apply top pt reweighting
  for (auto & rew : reweighting_modules) {
    rew->process(event);
  }
  
  //apply muon b tagging scale factors 
  muo_tight_noniso_SF->process(event);
  muo_trigger_SF->process(event);
  btagwAK8->process(event);

  scale_variation->process(event);

  for(auto & h : hists_after_sel){
    h->fill(event);
  }
  

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
  
  if(invert_merged_selection){ 
    if(mergedEvent_selection->passes_probe(event, probe_jet)) return false;
  }
  else if (merged_category){
    if(!mergedEvent_selection->passes_probe(event, probe_jet)) return false;
  }
  
 
  //=================================
  //define working points and tags
  //=================================

  double probejet_tau32 = probe_jet.tau3()/probe_jet.tau2();

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
  double probejet_pt = probe_jet.pt();

  bool mass_cut = false;
  if(usePUPPI && probejet_mass > 105 && probejet_mass < 210) mass_cut = true;
  else if(!usePUPPI && probejet_mass > 105 && probejet_mass < 220) mass_cut = true;

  std::vector<double> tau32_wps_CHS           {0.50, 0.57, 0.67, 0.81};
  std::vector<double> tau32_wps_PUPPI         {0.46, 0.54, 0.65, 0.80};

  std::vector<double> tau32_wps_HTT_CHS       {0.49, 0.54, 0.65, 0.97, 0.98, 0.55, 0.62, 0.93, 0.97};
  std::vector<double> fRec_wps_CHS            {0.14, 0.20, 0.25, 0.22, 0.50, 0.17, 0.27, 0.20, 0.47};
  std::vector<double> tau32_wps_HTT_PUPPI     {0.39, 0.48, 0.60, 0.96, 0.98, 0.50, 0.60, 0.99, 0.97};
  std::vector<double> fRec_wps_PUPPI          {0.36, 0.43, 0.22, 0.25, 0.48, 0.16, 0.25, 0.22, 0.46};

  //  double fRec = 0;
  //if(useHTT) fRec = probe_jet.get_tag( probe_jet.tagname2tag("fRec"));

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

  if(probejet_mass < 10) return false; //mild mass cut on the denominator 
 

  //===================
  //fill histograms
  //===================

  hists_all->fill_probe(event, probe_jet);
  hists_subjet_btag_eff->fill(event);

  if(probejet_pt > 400) hists_all_400->fill_probe(event, probe_jet);
  if(probejet_pt > 400 && probejet_pt < 550) hists_all_400to550->fill_probe(event, probe_jet);
  if(probejet_pt > 550) hists_all_550->fill_probe(event, probe_jet);


  // if (version == "TTbar_Incl" || version ==  "TTbar_700to1000" || version ==  "TTbar_1000toInf") {
  if ( ((TString)version).Contains("TTbar") ){
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
 
 
  //===================
  //new working points
  //===================

  std::vector<bool> toptag, toptag_mass, toptag_btag, toptag_mass_btag;

  std::vector<double> tau32_wps_CHS_new           {-1, 0.40, 0.50, 0.57, 0.67, 0.81};
  std::vector<double> tau32_wps_PUPPI_new         {-1, 0.40, 0.46, 0.54, 0.65, 0.80};

  std::vector<double> tau32_wps_new; 
  std::vector<double> fRec_wps_new; 
  if(!useHTT && !usePUPPI) tau32_wps_new = tau32_wps_CHS_new;
  else if(!useHTT && usePUPPI) tau32_wps_new = tau32_wps_PUPPI_new;
  else if(useHTT && !usePUPPI) {
    tau32_wps_new = tau32_wps_HTT_CHS; 
    fRec_wps = fRec_wps_CHS;
  }
  else if(useHTT && usePUPPI) {
    tau32_wps = tau32_wps_HTT_PUPPI; 
    fRec_wps = fRec_wps_PUPPI;
  }
    
  for(unsigned int wp = 0; wp < tau32_wps_new.size(); ++wp){
  
    bool tau_cut = get_tau32_cut(wp, probejet_tau32, tau32_wps_new);
    
    if(tau_cut) toptag.push_back(true); 
    else toptag.push_back(false);
      
    if(tau_cut && mass_cut) toptag_mass.push_back(true); 
    else toptag_mass.push_back(false);
    
    if(tau_cut && subjet_btag) toptag_btag.push_back(true); 
    else toptag_btag.push_back(false);
      
    if(tau_cut && mass_cut && subjet_btag) toptag_mass_btag.push_back(true); 
    else toptag_mass_btag.push_back(false);
  }
  
  
  //=======================
  //fill the new histograms
  //=======================

  for(unsigned int bin = 0; bin < pt_bins.size(); ++bin){
    for(unsigned int wp = 0; wp < tau32_wps_new.size(); ++wp){

      bool pt_cut = get_pt_cut(bin, probejet_pt);
      if(pt_cut){
	if(toptag.at(wp)) h_probe_all_pass.at(bin).at(wp)->fill_probe(event, probe_jet);
        else h_probe_all_fail.at(bin).at(wp)->fill_probe(event, probe_jet);

	if(toptag_mass.at(wp)) h_probe_mass_pass.at(bin).at(wp)->fill_probe(event, probe_jet);
	else h_probe_mass_fail.at(bin).at(wp)->fill_probe(event, probe_jet);
      }

    }
  }

  /*//fill NPV hists
  if(probejet_pt > 400){
    for(unsigned int bin = 0; bin < npv_bins.size(); ++bin){
      for(unsigned int wp = 0; wp < tau32_wps_new.size(); ++wp){
	
	bool npv_cut = get_npv_cut(bin,  event.pvs->size());
	if(npv_cut){
	  if(toptag.at(wp)) h_probeNPV_all_pass.at(bin).at(wp)->fill_probe(event, probe_jet);
	  else h_probeNPV_all_fail.at(bin).at(wp)->fill_probe(event, probe_jet);
	  
	  if(toptag_mass.at(wp)) h_probeNPV_mass_pass.at(bin).at(wp)->fill_probe(event, probe_jet);
	  else h_probeNPV_mass_fail.at(bin).at(wp)->fill_probe(event, probe_jet);
	}
	
      }
    }
    
  }
  */
  /*
 for(unsigned int bin = 0; bin < pt_bins.size(); ++bin){
    for(unsigned int wp = 0; wp < tau32_wps_new.size(); ++wp){

      bool pt_cut = get_pt_cut(bin, probejet_pt);
      if(pt_cut){
	if(toptag_btag.at(wp)) h_probe_btag_NoSF_pass.at(bin).at(wp)->fill_probe(event, probe_jet);
	else h_probe_btag_NoSF_fail.at(bin).at(wp)->fill_probe(event, probe_jet);
      }

    }
  }
  */
 subjet_btagwAK8->process(event);
 
  for(unsigned int bin = 0; bin < pt_bins.size(); ++bin){
    for(unsigned int wp = 0; wp < tau32_wps_new.size(); ++wp){

      bool pt_cut = get_pt_cut(bin, probejet_pt);
      if(pt_cut){
	if(toptag_btag.at(wp)) h_probe_btag_pass.at(bin).at(wp)->fill_probe(event, probe_jet);
	else h_probe_btag_fail.at(bin).at(wp)->fill_probe(event, probe_jet);

	if(toptag_mass_btag.at(wp)) h_probe_mass_btag_pass.at(bin).at(wp)->fill_probe(event, probe_jet);
	else h_probe_mass_btag_fail.at(bin).at(wp)->fill_probe(event, probe_jet);
      }

    }
  }
 
  /*
  //Fill NPV hists
  if(probejet_pt > 400){

    for(unsigned int bin = 0; bin < npv_bins.size(); ++bin){
      for(unsigned int wp = 0; wp < tau32_wps_new.size(); ++wp){
	
	bool npv_cut = get_npv_cut(bin, event.pvs->size());
	if(npv_cut){
	  if(toptag_btag.at(wp)) h_probeNPV_btag_pass.at(bin).at(wp)->fill_probe(event, probe_jet);
	  else h_probeNPV_btag_fail.at(bin).at(wp)->fill_probe(event, probe_jet);
	  
	  if(toptag_mass_btag.at(wp)) h_probeNPV_mass_btag_pass.at(bin).at(wp)->fill_probe(event, probe_jet);
	  else h_probeNPV_mass_btag_fail.at(bin).at(wp)->fill_probe(event, probe_jet);
	}
	
      }
    }
  }
  */
  return true;
}

bool TTEfficiencyPostSelectionModule::get_pt_cut(unsigned int bin, const double pt){
  if(bin > pt_bins.size()) return false;

  bool m_pt_cut = false; 
  if(pt_bins.at(bin) == -1) m_pt_cut = true; 
  else if(pt_bins.at(bin) >= 0 && bin+1 < pt_bins.size() && pt_bins.at(bin+1) >= 0){
    if( pt > pt_bins.at(bin) && pt < pt_bins.at(bin+1) ) m_pt_cut = true; 
  }
  else{
    if( pt > fabs(pt_bins.at(bin)) ) m_pt_cut = true; 
  }
  return m_pt_cut; 
}


bool TTEfficiencyPostSelectionModule::get_npv_cut(unsigned int bin, const int npv){
  if(bin > npv_bins.size()) return false;

  bool m_npv_cut = false; 
  if(npv_bins.at(bin) == -1) m_npv_cut = true; 
  else if(npv_bins.at(bin) >= 0 && bin+1 < npv_bins.size() && npv_bins.at(bin+1) >= 0){
    if( npv > npv_bins.at(bin) && npv <= npv_bins.at(bin+1) ) m_npv_cut = true; 
  }
  else{
    if( npv > fabs(npv_bins.at(bin)) ) m_npv_cut = true; 
  }
  return m_npv_cut; 
}

bool TTEfficiencyPostSelectionModule::get_tau32_cut(unsigned int bin, const double tau32, const std::vector<double> wps){
  if(bin > wps.size()) return false;

  bool m_tau_cut = false; 
  if(wps.at(bin) < 0) m_tau_cut = true; 
  else{
    if( tau32 < wps.at(bin)) m_tau_cut = true; 
  }
  return m_tau_cut;
}

UHH2_REGISTER_ANALYSIS_MODULE(TTEfficiencyPostSelectionModule)


