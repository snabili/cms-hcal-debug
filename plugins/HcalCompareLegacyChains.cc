// -*- C++ -*-
//
// Package:    HcalCompareLegacyChains
// Class:      HcalCompareLegacyChains
// 
/**\class HcalCompareLegacyChains HcalCompareLegacyChains.cc HcalDebug/CompareChans/src/HcalCompareLegacyChains.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  matthias wolf
//         Created:  Fri Aug 26 11:37:21 CDT 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CondFormats/DataRecord/interface/HcalChannelQualityRcd.h"

#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"

//HCAL
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"

//ECAL
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "Geometry/CaloTopology/interface/EcalEndcapTopology.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputer.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputerRcd.h"

//PF Algo
#include "RecoParticleFlow/Configuration/test/PFChargedHadronAnalyzer.h"
#include "RecoParticleFlow/PFTracking/interface/PFTrackAlgoTools.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFLayer.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFClusterWidthAlgo.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

//Calo MET
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"

//L1T Offline
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1Upgrade.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoMetDataFormat.h"

//Primary Vertex 
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFinding.h"
#include "RecoVertex/ConfigurableVertexReco/interface/ConfigurableVertexReconstructor.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"


#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TTree.h"
//
// class declaration
//

using namespace std;
using namespace edm;
using namespace reco;

class HcalCompareLegacyChains : public edm::EDAnalyzer {
   public:
      explicit HcalCompareLegacyChains(const edm::ParameterSet&);
      ~HcalCompareLegacyChains();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void beginLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override;

      double get_cosh(const HcalDetId&);

      // ----------member data ---------------------------
      bool first_;


      std::vector<edm::InputTag> frames_;
      edm::InputTag digis_;
      std::vector<edm::InputTag> rechits_;

      edm::ESHandle<CaloGeometry> gen_geo_;
      edm::ESHandle<CaloGeometry> EcalBarrelGeometry_;
      edm::ESHandle<CaloGeometry> EcalEndcapGeometry_;

      edm::EDGetTokenT<l1t::EtSumBxCollection> sumToken_;

      edm::EDGetTokenT<reco::PFMETCollection>      metToken_;
      edm::EDGetTokenT<reco::VertexCollection>    pvToken_;

      bool swap_iphi_;

      TH2D *df_multiplicity_;
      TH2D *tp_multiplicity_;

      //====== 1D Histo Declaration =========
      TH1D *h_PFRecHit;
      TH1D *h_CHPVMPF_noPU;
      TH1D *h_CHPVMPF_PU;
      TH1D *h_CHPVMPF;
      TH1D *h_RecMET_noPU;
      TH1D *h_RecMET_PU;
      TH1D *h_RecMET;
      TH1D *h_RecHEdepth1_noPU;
      TH1D *h_RecHEdepth2_noPU;
      TH1D *h_RecHEdepth3_noPU;
      TH1D *h_RecHEdepth4_noPU;
      TH1D *h_RecHEdepth5_noPU;
      TH1D *h_RecHEdepth6_noPU;
      TH1D *h_RecHEdepth7_noPU;
      TH1D *h_RecHEdepth1_PU;
      TH1D *h_RecHEdepth2_PU;
      TH1D *h_RecHEdepth3_PU;
      TH1D *h_RecHEdepth4_PU;
      TH1D *h_RecHEdepth5_PU;
      TH1D *h_RecHEdepth6_PU;
      TH1D *h_RecHEdepth7_PU;
      TH1D *h_PFRecHit_Et;

      //====== 2D Histo Declaration =========
      TH2D *h_PFMET_PFMETCandidate;
      TH2D *h_PFMET_RecMET_noPU;
      TH2D *h_PFMET_RecMET_depthHE_noPU;
      TH2D *h_PFMET_RecMET_PU;
      TH2D *h_PFMET_RecMET_depthHE_PU;
      TH2D *h_PFMET_RecMET_nodepthHE;
      TH2D *h_depth_RecHE;
      TH2D *h_depth_RecHE_noPU;
      TH2D *h_depth_RecHE_PU;
      TH2D *h_depth_fracRecHE_noPU;
      TH2D *h_depth_fracRecHE_PU;
      TH2D *h_PFRecMET_PFMET_noPU;
      TH2D *h_PFRecMET_PFMET_PU;
      TH2D *h_PFRecMET_PFMET;

      TTree *tps_;
      TTree *tpsplit_;
      TTree *events_;
      TTree *matches_;

      double tp_energy_;
      int tp_ieta_;
      int tp_iphi_;
      int tp_soi_;

      double tpsplit_energy_;
      double tpsplit_oot_;
      int tpsplit_ieta_;
      int tpsplit_iphi_;
      int tpsplit_depth_;
      double tpsplit_ettot_;
      double tpsplit_rise_avg_;
      double tpsplit_rise_rms_;
      double tpsplit_fall_avg_;
      double tpsplit_fall_rms_;

      double ev_rh_energy0_;
      double ev_rh_energy2_;
      double ev_rh_energy3_;
      double ev_tp_energy_;
      int ev_rh_unmatched_;
      int ev_tp_unmatched_;

      double mt_rh_energy0_;
      double mt_rh_energy2_;
      double mt_rh_energy3_;
      double mt_tp_energy_;

      int mt_ieta_;
      int mt_iphi_;
      int mt_version_;
      int mt_tp_soi_;

      int max_severity_;
      L1Analysis::L1AnalysisRecoMetDataFormat* met_data;

      const HcalChannelQuality* status_;
      const HcalSeverityLevelComputer* comp_;

      edm::EDGetTokenT<EcalRecHitCollection>   tok_EB_;
      edm::EDGetTokenT<EcalRecHitCollection>   tok_EE_;
      edm::EDGetTokenT<HBHERecHitCollection>   hRhToken;
      edm::EDGetTokenT<PFCandidateCollection>  inputTagPFCandidates_;

      const reco::Vertex* primary_vertex_;
};

HcalCompareLegacyChains::HcalCompareLegacyChains(const edm::ParameterSet& config) :
   edm::EDAnalyzer(),
   first_(true),
   frames_(config.getParameter<std::vector<edm::InputTag>>("dataFrames")),
   digis_(config.getParameter<edm::InputTag>("triggerPrimitives")),
   rechits_(config.getParameter<std::vector<edm::InputTag>>("recHits")),
   swap_iphi_(config.getParameter<bool>("swapIphi")),
   max_severity_(config.getParameter<int>("maxSeverity"))
{
   consumes<HcalTrigPrimDigiCollection>(digis_);
   consumes<HBHEDigiCollection>(frames_[0]);
   consumes<HFDigiCollection>(frames_[1]);
   consumes<edm::SortedCollection<HBHERecHit>>(rechits_[0]);
   consumes<edm::SortedCollection<HFRecHit>>(rechits_[1]);

   tok_EB_     = consumes<EcalRecHitCollection>(edm::InputTag("ecalRecHit","EcalRecHitsEB"));
   tok_EE_     = consumes<EcalRecHitCollection>(edm::InputTag("ecalRecHit","EcalRecHitsEE"));

   inputTagPFCandidates_ = consumes<PFCandidateCollection>(edm::InputTag("particleFlow"));

   consumes<reco::PFMETCollection> (edm::InputTag("pfMet"));
   sumToken_ = consumes<l1t::EtSumBxCollection>(config.getUntrackedParameter<edm::InputTag>("sumToken"));
   metToken_ = consumes<reco::PFMETCollection>(config.getUntrackedParameter<edm::InputTag>("metToken"));

   hRhToken = consumes<HBHERecHitCollection >(config.getUntrackedParameter<string>("HBHERecHits","hbheprereco"));

   pvToken_ = consumes<reco::VertexCollection>(config.getParameter<edm::InputTag>("pvToken"));

   edm::Service<TFileService> fs;
   met_data = new L1Analysis::L1AnalysisRecoMetDataFormat();

   df_multiplicity_ = fs->make<TH2D>("df_multiplicity", "DataFrame multiplicity;ieta;iphi", 65, -32.5, 32.5, 72, 0.5, 72.5);
   tp_multiplicity_ = fs->make<TH2D>("tp_multiplicity", "TrigPrim multiplicity;ieta;iphi", 65, -32.5, 32.5, 72, 0.5, 72.5);

   tps_ = fs->make<TTree>("tps", "Trigger primitives");
   tps_->Branch("et", &tp_energy_);
   tps_->Branch("ieta", &tp_ieta_);
   tps_->Branch("iphi", &tp_iphi_);
   tps_->Branch("soi", &tp_soi_);

   tpsplit_ = fs->make<TTree>("tpsplit", "Trigger primitives");
   tpsplit_->Branch("et", &tpsplit_energy_);
   tpsplit_->Branch("oot", &tpsplit_oot_);
   tpsplit_->Branch("ieta", &tpsplit_ieta_);
   tpsplit_->Branch("iphi", &tpsplit_iphi_);
   tpsplit_->Branch("depth", &tpsplit_depth_);
   tpsplit_->Branch("etsum", &tpsplit_ettot_);
   tpsplit_->Branch("rise_avg", &tpsplit_rise_avg_);
   tpsplit_->Branch("rise_rms", &tpsplit_rise_rms_);
   tpsplit_->Branch("fall_avg", &tpsplit_fall_avg_);
   tpsplit_->Branch("fall_rms", &tpsplit_fall_rms_);

   events_ = fs->make<TTree>("events", "Event quantities");
   events_->Branch("RH_energyM0", &ev_rh_energy0_);
   events_->Branch("RH_energyM2", &ev_rh_energy2_);
   events_->Branch("RH_energyM3", &ev_rh_energy3_);
   events_->Branch("TP_energy", &ev_tp_energy_);
   events_->Branch("RH_unmatched", &ev_rh_unmatched_);
   events_->Branch("TP_unmatched", &ev_tp_unmatched_);

   matches_ = fs->make<TTree>("matches", "Matched RH and TP");
   matches_->Branch("RH_energyM0", &mt_rh_energy0_);
   matches_->Branch("RH_energyM2", &mt_rh_energy2_);
   matches_->Branch("RH_energyM3", &mt_rh_energy3_);
   matches_->Branch("TP_energy", &mt_tp_energy_);
   matches_->Branch("ieta", &mt_ieta_);
   matches_->Branch("iphi", &mt_iphi_);
   matches_->Branch("tp_version", &mt_version_);
   matches_->Branch("tp_soi", &mt_tp_soi_);

   //====MET 1D Histograms=================================================================================
   h_PFRecHit = fs->make<TH1D>("h_PFRecHit","Computed PFRecHit from PFCandidate for Charged Hadrons",100,0,100);
   // Charged Hadron Vertex 
   h_CHPVMPF_noPU = fs->make<TH1D>("h_CHPVMPF_noPU","Charged Hadrons Primary Vertex Z Position", 100,-15,15);
   h_CHPVMPF_PU = fs->make<TH1D>("h_CHPVMPF_PU" ,"Charged Hadrons PFCandidate Vertex Z Position", 100,-15,15);
   h_CHPVMPF = fs->make<TH1D>("h_CHPVMPF","Charged Hadrons Vertex Minus PFCandidate Z Position",500,-15,15);

   //====MET 2D Histograms=================================================================================
   h_PFMET_PFMETCandidate = fs->make<TH2D>("h_PFMET_PFMETCandidate","PFMET vs PFMET Computed from PFCandidate",100,0,120,100,0,120);

   //====Depth Study=================================================================================
   h_PFMET_RecMET_noPU = fs->make<TH2D>("h_PFMET_RecMET_noPU","no Pileup PF MET vs Rechit MET",100,0,150,100,0,150);
   h_PFMET_RecMET_depthHE_noPU = fs->make<TH2D>("h_PFMET_RecMET_depthHE_noPU","no Pileup PF MET vs Rechit MET with Depth Threshold",100,0,150,100,0,150);
   h_PFMET_RecMET_PU = fs->make<TH2D>("h_PFMET_RecMET_PU","Pileup PF MET vs Rechit MET",100,0,150,100,0,150);
   h_PFMET_RecMET_depthHE_PU = fs->make<TH2D>("h_PFMET_RecMET_depthHE_PU","Pileup PF MET vs Rechit MET with Depth Threshold",100,0,150,100,0,150);
   h_PFMET_RecMET_nodepthHE = fs->make<TH2D>("h_PFMET_RecMET_nodepthHE","PF MET vs Rechit MET without Depth Threshold",100,0,150,100,0,150);
   h_RecMET_noPU = fs->make<TH1D>("h_RecMET_noPU","RecHit MET from no Pileup",100,-10,300);
   h_RecMET_PU = fs->make<TH1D>("h_RecMET_PU","RecHit MET from Pileup",100,-10,300);
   h_RecMET = fs->make<TH1D>("h_RecMET","RecHit MET",100,-10,300);
   h_depth_RecHE = fs->make<TH2D>("h_depth_RecHE","RechitHE vs Depth",1000,0,10,1000,0,7.5);
   h_depth_RecHE_noPU = fs->make<TH2D>("h_depth_RecHE_noPU","RechitHE vs Depth for Dz<0.1",1000,0,10,1000,0,7.5);
   h_depth_RecHE_PU = fs->make<TH2D>("h_depth_RecHE_PU","RechitHE vs Depth for Dz>0.1",1000,0,10,1000,0,7.5);

   h_RecHEdepth1_noPU = fs->make<TH1D>("h_RecHEdepth1_noPU","Rechit depth == 1 no Pileup",500,-1,10);
   h_RecHEdepth2_noPU = fs->make<TH1D>("h_RecHEdepth2_noPU","Rechit depth == 2 no Pileup",500,-1,10);
   h_RecHEdepth3_noPU = fs->make<TH1D>("h_RecHEdepth3_noPU","Rechit depth == 3 no Pileup",500,-1,10);
   h_RecHEdepth4_noPU = fs->make<TH1D>("h_RecHEdepth4_noPU","Rechit depth == 4 no Pileup",500,-1,10);
   h_RecHEdepth5_noPU = fs->make<TH1D>("h_RecHEdepth5_noPU","Rechit depth == 5 no Pileup",500,-1,10);
   h_RecHEdepth6_noPU = fs->make<TH1D>("h_RecHEdepth6_noPU","Rechit depth == 6 no Pileup",500,-1,10);
   h_RecHEdepth7_noPU = fs->make<TH1D>("h_RecHEdepth7_noPU","Rechit depth == 7 no Pileup",500,-1,10);

   h_RecHEdepth1_PU = fs->make<TH1D>("h_RecHEdepth1_PU","Rechit depth == 1 Pileup",500,-1,10);
   h_RecHEdepth2_PU = fs->make<TH1D>("h_RecHEdepth2_PU","Rechit depth == 2 Pileup",500,-1,10);
   h_RecHEdepth3_PU = fs->make<TH1D>("h_RecHEdepth3_PU","Rechit depth == 3 Pileup",500,-1,10);
   h_RecHEdepth4_PU = fs->make<TH1D>("h_RecHEdepth4_PU","Rechit depth == 4 Pileup",500,-1,10);
   h_RecHEdepth5_PU = fs->make<TH1D>("h_RecHEdepth5_PU","Rechit depth == 5 Pileup",500,-1,10);
   h_RecHEdepth6_PU = fs->make<TH1D>("h_RecHEdepth6_PU","Rechit depth == 6 Pileup",500,-1,10);
   h_RecHEdepth7_PU = fs->make<TH1D>("h_RecHEdepth7_PU","Rechit depth == 7 Pileup",500,-1,10);

   h_depth_fracRecHE_noPU = fs->make<TH2D>("h_depth_fracRecHE_noPU","No Pileup Fractio RecHit energy vs Depth",1000,0,10,100,0,2);
   h_depth_fracRecHE_PU = fs->make<TH2D>("h_depth_fracRecHE_PU","Pileup Fractio RecHit energy vs Depth",1000,0,10,100,0,2);
   h_PFRecHit_Et = fs->make<TH1D>("h_PFRecHit_Et","Sum Et of PF hits in ECAL EB, HCAL HB",100,0,250);
   h_PFRecMET_PFMET_noPU = fs->make<TH2D>("h_PFRecMET_PFMET_noPU","noPU PFRecMET Vs PFMET",100,0,200,100,0,200);
   h_PFRecMET_PFMET_PU = fs->make<TH2D>("h_PFRecMET_PFMET_PU","PU PFRecMET Vs PFMET",100,0,200,100,0,200);
   h_PFRecMET_PFMET = fs->make<TH2D>("h_PFRecMET_PFMET","PFRecMET Vs PFMET",100,0,200,100,0,200);

}

HcalCompareLegacyChains::~HcalCompareLegacyChains() {}

double
HcalCompareLegacyChains::get_cosh(const HcalDetId& id)
{
   const auto *sub_geo = dynamic_cast<const HcalGeometry*>(gen_geo_->getSubdetectorGeometry(id));
   auto eta = sub_geo->getPosition(id).eta();
   return cosh(eta);
}

void
HcalCompareLegacyChains::beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& setup)
{
   edm::ESHandle<HcalChannelQuality> status;
   setup.get<HcalChannelQualityRcd>().get("withTopo", status);
   status_ = status.product();
   edm::ESHandle<HcalSeverityLevelComputer> comp;
   setup.get<HcalSeverityLevelComputerRcd>().get(comp);
   comp_ = comp.product();
}
unsigned int ievt = 0;
void
HcalCompareLegacyChains::analyze(const edm::Event& event, const edm::EventSetup& setup)
{

   ievt++;
   std::cout<<"Event number: "<<ievt<<" started"<<std::endl;

   //============Offline PFMET ===================================================================== 
   edm::Handle<reco::PFMETCollection> metLabel_;
   event.getByToken(metToken_, metLabel_);
   const reco::PFMETCollection *metCol = metLabel_.product();
   const reco::PFMET theMet = metCol->front();
   met_data->met     = theMet.et();
   met_data->metPhi  = theMet.phi();
   met_data->sumEt   = theMet.sumEt();

   //================ Primary Vertex Position =============================================
   edm::Handle<reco::VertexCollection> PVs;
   event.getByToken(pvToken_, PVs);
   math::XYZPoint leadPV(0,0,0);

   //================ PFRecHits from PFCandidates->PFBlocks->PFClusters->PFRecHitFractions =============================================
   edm::Handle<reco::PFCandidateCollection> pfCandidate;
   event.getByToken(inputTagPFCandidates_, pfCandidate);

   double PFRecHit[7] = {0.0};
   double PFRecHit_tot[7] = {0.0};
   double PFRecHit_Et = 0;
   double HEx_D1, HEx_D2, HEx_D3, HEx_D4, HEx_D5, HEx_D6, HEx_D7 = 0;
   double HEy_D1, HEy_D2, HEy_D3, HEy_D4, HEy_D5, HEy_D6, HEy_D7 = 0;
   double HEx_d1, HEx_d2, HEx_d3, HEx_d4, HEx_d5, HEx_d6, HEx_d7 = 0;
   double HEy_d1, HEy_d2, HEy_d3, HEy_d4, HEy_d5, HEy_d6, HEy_d7 = 0;
   double HB_Ex, HB_Ey, Ecal_Ex, Ecal_Ey = 0;
   double HE_Ex, HE_Ey, Ex, Ey = 0;
   double HEx_PU, HEy_PU, Ex_PU, Ey_PU = 0;
   double HEx_noPU, HEy_noPU, Ex_noPU, Ey_noPU = 0;
   double PFRecMET, PFRecMET_PU, PFRecMET_noPU = 0;
   double posX,posY,posZ = 0;
   for( reco::PFCandidateCollection::const_iterator ci  = pfCandidate->begin(); ci!=pfCandidate->end(); ++ci)  { //PFCandidate Loop
     const reco::PFCandidate& pfc = *ci;
     if ( pfc.particleId() != 1 ) continue;// Charged Hadron Particles
     math::XYZPoint vtx = pfc.vertex();
     if (!PVs->empty() && !((*PVs)[0].isFake())) {
       leadPV = math::XYZPoint((*PVs)[0].x(),(*PVs)[0].y(),(*PVs)[0].z());
     }
     double dz = vtx.z()-leadPV.z();
     h_CHPVMPF->Fill(dz);
     const reco::PFCandidate::ElementsInBlocks& pfBlocks = ci->elementsInBlocks();
     for ( reco::PFCandidate::ElementsInBlocks::const_iterator pfBlock = pfBlocks.begin(); pfBlock != pfBlocks.end(); ++pfBlock ) {//PFBlock Loop
       const edm::OwnVector<reco::PFBlockElement>& pfBlockElements = pfBlock->first->elements();
       for ( edm::OwnVector<reco::PFBlockElement>::const_iterator pfBlockElement = pfBlockElements.begin(); pfBlockElement != pfBlockElements.end(); ++pfBlockElement ) {//PFBlock Elements Loop
	 HE_Ex = 0;
         HE_Ey = 0;
	 HEx_D1 = 0;
         HEx_D2 = 0;
         HEx_D3 = 0;
         HEx_D4 = 0;
         HEx_D5 = 0;
         HEx_D6 = 0;
         HEx_D7 = 0;

         HEy_D1 = 0;
         HEy_D2 = 0;
         HEy_D3 = 0;
         HEy_D4 = 0;
         HEy_D5 = 0;
         HEy_D6 = 0;
         HEy_D7 = 0;

         HEx_d1 = 0;
         HEx_d2 = 0;
         HEx_d3 = 0;
         HEx_d4 = 0;
         HEx_d5 = 0;
         HEx_d6 = 0;
         HEx_d7 = 0;

         HEy_d1 = 0;
         HEy_d2 = 0;
         HEy_d3 = 0;
         HEy_d4 = 0;
         HEy_d5 = 0;
         HEy_d6 = 0;
         HEy_d7 = 0;
	 Ecal_Ex = 0;
         Ecal_Ey = 0;
	 HB_Ex = 0;
	 HB_Ey = 0;

         if ( pfBlockElement->clusterRef().isNonnull() ) {
           reco::PFClusterRef pfCluster = pfBlockElement->clusterRef();
           const std::vector<reco::PFRecHitFraction>& pfRecHitFractions = pfCluster->recHitFractions();
           for ( std::vector<reco::PFRecHitFraction>::const_iterator it = pfRecHitFractions.begin(); it != pfRecHitFractions.end(); ++it ) {//PFCluster Loop
             const reco::PFRecHitRef& pfRecHits = it->recHitRef();
             posX = pfRecHits->position().x();
             posY = pfRecHits->position().y();
             posZ = pfRecHits->position().z();
             math::XYZPoint pflowPos(posX,posY,posZ);
             double Eta    = pflowPos.eta();
             double Phi    = pflowPos.phi(); 
             unsigned int depth = pfRecHits->depth(); 
             PFRecHit[depth] = pfRecHits->energy()/cosh(Eta); 
             PFRecHit_tot[depth] += pfRecHits->energy()/cosh(Eta);
	     if (pfCluster->layer() == PFLayer::HCAL_ENDCAP && 1.3<Eta && Eta<=3) {
               h_depth_RecHE->Fill(depth,PFRecHit[depth]);
               h_depth_RecHE->SetMarkerStyle(3);
               h_depth_RecHE->SetMarkerSize(1.5);
               h_depth_RecHE->GetXaxis()->SetTitle("depth");
               h_depth_RecHE->GetYaxis()->SetTitle("PFRecHit[depth] (GeV)");
	       HE_Ex += cos(Phi)*pfRecHits->energy()/cosh(Eta);
               HE_Ey += sin(Phi)*pfRecHits->energy()/cosh(Eta);

       	       if(fabs(dz)<0.1){//no Pileup
	         h_CHPVMPF_noPU->Fill(dz);
                 h_depth_RecHE_noPU->Fill(depth,PFRecHit[depth]);
                 h_depth_RecHE_noPU->SetMarkerStyle(3);
                 h_depth_RecHE_noPU->SetMarkerColor(2);
                 h_depth_RecHE_noPU->SetMarkerSize(1.5);
                 h_depth_RecHE_noPU->GetXaxis()->SetTitle("depth");
                 h_depth_RecHE_noPU->GetYaxis()->SetTitle("PFRecHit (GeV)");

                 h_depth_fracRecHE_noPU->Fill(depth,PFRecHit[depth]/PFRecHit_tot[depth]);
                 h_depth_fracRecHE_noPU->SetMarkerStyle(3);
                 h_depth_fracRecHE_noPU->SetMarkerColor(6);
                 h_depth_fracRecHE_noPU->SetMarkerSize(1.5);
                 h_depth_fracRecHE_noPU->GetXaxis()->SetTitle("depth");
                 h_depth_fracRecHE_noPU->GetYaxis()->SetTitle("RechitHE (GeV)");
                 if(depth == 1 && (pfRecHits->energy()/cosh(Eta))>0.38){
                   h_RecHEdepth1_noPU->Fill(pfRecHits->energy()/cosh(Eta));
                    HEx_D1 += cos(Phi)*pfRecHits->energy()/cosh(Eta);
                    HEy_D1 += sin(Phi)*pfRecHits->energy()/cosh(Eta);
                 }
                 if(depth == 2 && (pfRecHits->energy()/cosh(Eta))>0.54){
                    h_RecHEdepth2_noPU->Fill(pfRecHits->energy()/cosh(Eta));
	            HEx_D2 += cos(Phi)*pfRecHits->energy()/cosh(Eta);
                    HEy_D2 += sin(Phi)*pfRecHits->energy()/cosh(Eta);
                 }
                 if(depth == 3 && (pfRecHits->energy()/cosh(Eta))>0.3){
                    h_RecHEdepth3_noPU->Fill(pfRecHits->energy()/cosh(Eta));
                    HEx_D3 += cos(Phi)*pfRecHits->energy()/cosh(Eta);
                    HEy_D3 += sin(Phi)*pfRecHits->energy()/cosh(Eta);
                 }
                 if(depth == 4 && (pfRecHits->energy()/cosh(Eta))>0.2){
                    h_RecHEdepth4_noPU->Fill(pfRecHits->energy()/cosh(Eta));
                    HEx_D4 += cos(Phi)*pfRecHits->energy()/cosh(Eta);
                    HEy_D4 += sin(Phi)*pfRecHits->energy()/cosh(Eta);
                 }
                 if(depth == 5 && (pfRecHits->energy()/cosh(Eta))>0.17){
                    h_RecHEdepth5_noPU->Fill(pfRecHits->energy()/cosh(Eta));
                    HEx_D5 += cos(Phi)*pfRecHits->energy()/cosh(Eta);
                    HEy_D5 += sin(Phi)*pfRecHits->energy()/cosh(Eta);
                 }
                 if(depth == 6 && (pfRecHits->energy()/cosh(Eta))>0.167){
                    h_RecHEdepth6_noPU->Fill(pfRecHits->energy()/cosh(Eta));
                    HEx_D6 += cos(Phi)*pfRecHits->energy()/cosh(Eta);
                    HEy_D6 += sin(Phi)*pfRecHits->energy()/cosh(Eta);
                 }
                 if(depth == 7 && (pfRecHits->energy()/cosh(Eta))>0.125){
                    h_RecHEdepth7_noPU->Fill(pfRecHits->energy()/cosh(Eta));
                    HEx_D7 += cos(Phi)*pfRecHits->energy()/cosh(Eta);
                    HEy_D7 += sin(Phi)*pfRecHits->energy()/cosh(Eta);
                 }
	       }
               if(fabs(dz)>0.1){//Pileup
                 h_CHPVMPF_PU->Fill(dz);
                 h_depth_RecHE_PU->Fill(depth,PFRecHit[depth]);
                 h_depth_RecHE_PU->SetMarkerStyle(3);
                 h_depth_RecHE_PU->SetMarkerColor(4);
                 h_depth_RecHE_PU->SetMarkerSize(1.5);
                 h_depth_RecHE_PU->GetXaxis()->SetTitle("depth");
                 h_depth_RecHE_PU->GetYaxis()->SetTitle("PFRecHit (GeV)");

                 h_depth_fracRecHE_PU->Fill(depth,PFRecHit[depth]/PFRecHit_tot[depth]);
                 h_depth_fracRecHE_PU->SetMarkerStyle(3);
                 h_depth_fracRecHE_PU->SetMarkerColor(9);
                 h_depth_fracRecHE_PU->SetMarkerSize(1.5);
                 h_depth_fracRecHE_PU->GetXaxis()->SetTitle("depth");
                 h_depth_fracRecHE_PU->GetYaxis()->SetTitle("RechitHE (GeV)");
                 if(depth == 1 && (pfRecHits->energy()/cosh(Eta))>0.38){
                   h_RecHEdepth1_PU->Fill(pfRecHits->energy()/cosh(Eta));
                   HEx_d1 += cos(Phi)*pfRecHits->energy()/cosh(Eta);
                   HEy_d1 += sin(Phi)*pfRecHits->energy()/cosh(Eta);
                 }
                 if(depth == 2 && (pfRecHits->energy()/cosh(Eta))>0.54){
                    h_RecHEdepth2_PU->Fill(pfRecHits->energy()/cosh(Eta));
                    HEx_d2 += cos(Phi)*pfRecHits->energy()/cosh(Eta);
                    HEy_d2 += sin(Phi)*pfRecHits->energy()/cosh(Eta);
                 }
                 if(depth == 3 && (pfRecHits->energy()/cosh(Eta))>0.3){
                    h_RecHEdepth3_PU->Fill(pfRecHits->energy()/cosh(Eta));
                    HEx_d3 += cos(Phi)*pfRecHits->energy()/cosh(Eta);
                    HEy_d3 += sin(Phi)*pfRecHits->energy()/cosh(Eta);
                 }
                 if(depth == 4 && (pfRecHits->energy()/cosh(Eta))>0.2){
                    h_RecHEdepth4_PU->Fill(pfRecHits->energy()/cosh(Eta));
                    HEx_d4 += cos(Phi)*pfRecHits->energy()/cosh(Eta);
                    HEy_d4 += sin(Phi)*pfRecHits->energy()/cosh(Eta);
                 }
                 if(depth == 5 && (pfRecHits->energy()/cosh(Eta))>0.17){
                    h_RecHEdepth5_PU->Fill(pfRecHits->energy()/cosh(Eta));
                    HEx_d5 += cos(Phi)*pfRecHits->energy()/cosh(Eta);
                    HEy_d5 += sin(Phi)*pfRecHits->energy()/cosh(Eta);
                 }
                 if(depth == 6 && (pfRecHits->energy()/cosh(Eta))>0.167){
                    h_RecHEdepth6_PU->Fill(pfRecHits->energy()/cosh(Eta));
                    HEx_d6 += cos(Phi)*pfRecHits->energy()/cosh(Eta);
                    HEy_d6 += sin(Phi)*pfRecHits->energy()/cosh(Eta);
                 }
                 if(depth == 7 && (pfRecHits->energy()/cosh(Eta))>0.125){
                    h_RecHEdepth7_PU->Fill(pfRecHits->energy()/cosh(Eta));
                    HEx_d7 += cos(Phi)*pfRecHits->energy()/cosh(Eta);
                    HEy_d7 += sin(Phi)*pfRecHits->energy()/cosh(Eta);
                 }
               }
	     }	 
	     if(pfCluster->layer() == PFLayer::HCAL_BARREL2 && Eta<1.3){
               HB_Ex += cos(Phi)*pfRecHits->energy()/cosh(Eta);
               HB_Ey += sin(Phi)*pfRecHits->energy()/cosh(Eta);
               PFRecHit_Et = pfRecHits->energy()/cosh(Eta);
               h_PFRecHit_Et->Fill(PFRecHit_Et);
             }
	     if (pfCluster->layer() == PFLayer::ECAL_ENDCAP || pfCluster->layer() == PFLayer::ECAL_BARREL) {
	       Ecal_Ex += cos(Phi)*pfRecHits->energy()/cosh(Eta);
               Ecal_Ey += sin(Phi)*pfRecHits->energy()/cosh(Eta);
	     }
             HEx_noPU = HEx_d1 + HEx_d2 + HEx_d3 + HEx_d4 + HEx_d5 + HEx_d6 + HEx_d7;
             HEy_noPU = HEy_d1 + HEy_d2 + HEy_d3 + HEy_d4 + HEy_d5 + HEy_d6 + HEy_d7;
             Ex_noPU = HEx_noPU + HB_Ex + Ecal_Ex;
             Ey_noPU = HEy_noPU + HB_Ey + Ecal_Ey;
             PFRecMET_noPU = sqrt(Ex_noPU*Ex_noPU+Ey_noPU*Ey_noPU);

             HEx_PU = HEx_D1 + HEx_D2 + HEx_D3 + HEx_D4 + HEx_D5 + HEx_D6 + HEx_D7;
             HEy_PU = HEy_D1 + HEy_D2 + HEy_D3 + HEy_D4 + HEy_D5 + HEy_D6 + HEy_D7;
             Ex_PU = HEx_PU + HB_Ex + Ecal_Ex;
             Ey_PU = HEy_PU + HB_Ey + Ecal_Ey;
	     PFRecMET_PU = sqrt(Ex_PU*Ex_PU+Ey_PU*Ey_PU);

             Ex = HE_Ex + HB_Ex + Ecal_Ex;
             Ey = HE_Ey + HB_Ey + Ecal_Ey;
	     PFRecMET = sqrt(Ex*Ex+Ey*Ey);
           }   //End loop over PFCluster
         }
       }  //Loop over PFBlockElements
     }  //Loop over PFBlock
   }  //Loop over PFCandidate
                h_PFRecMET_PFMET->Fill(PFRecMET,theMet.et());
   h_PFRecMET_PFMET->SetMarkerStyle(3);
   h_PFRecMET_PFMET->SetMarkerColor(4);
   h_PFRecMET_PFMET->SetMarkerSize(1.5);
   h_PFRecMET_PFMET->GetXaxis()->SetTitle("PFRecMET (GeV)");
   h_PFRecMET_PFMET->GetYaxis()->SetTitle("PFMET (GeV)");
   h_PFRecMET_PFMET_PU->Fill(PFRecMET_PU,theMet.et());
   h_PFRecMET_PFMET_PU->SetMarkerStyle(3);
   h_PFRecMET_PFMET_PU->SetMarkerColor(4);
   h_PFRecMET_PFMET_PU->SetMarkerSize(1.5);
   h_PFRecMET_PFMET_PU->GetXaxis()->SetTitle("PFRecMET_PU (GeV)");
   h_PFRecMET_PFMET_PU->GetYaxis()->SetTitle("PFMET (GeV)");
   h_PFRecMET_PFMET_noPU->Fill(PFRecMET_noPU,theMet.et());
   h_PFRecMET_PFMET_noPU->SetMarkerStyle(3);
   h_PFRecMET_PFMET_noPU->SetMarkerColor(4);
   h_PFRecMET_PFMET_noPU->SetMarkerSize(1.5);
   h_PFRecMET_PFMET_noPU->GetXaxis()->SetTitle("PFRecMET_noPU (GeV)");
   h_PFRecMET_PFMET_noPU->GetYaxis()->SetTitle("PFMET (GeV)");


   edm::ESHandle<HcalTrigTowerGeometry> tpd_geo_h;
   setup.get<CaloGeometryRecord>().get(tpd_geo_h);
   edm::ESHandle<HcalDbService> conditions;
   setup.get<HcalDbRecord>().get(conditions);
   const HcalTrigTowerGeometry& tpd_geo = *tpd_geo_h;

   // ==========
   // Dataframes
   // ==========

   Handle<HBHEDigiCollection> frames;
   Handle<HFDigiCollection> hfframes;
   if (frames_.size() == 2) {
      if (first_ && event.getByLabel(frames_[0], frames) && event.getByLabel(frames_[1], hfframes)) {
         std::set<HcalTrigTowerDetId> ids;

         for (const auto& frame: *(frames.product())) {
            auto mapped = tpd_geo_h->towerIds(frame.id());

            for (const auto& id: mapped) {
               df_multiplicity_->Fill(id.ieta(), id.iphi());
               ids.insert(id);
            }
         }

         for (const auto& frame: *(hfframes.product())) {
            auto mapped = tpd_geo_h->towerIds(frame.id());

            for (const auto& id: mapped) {
               df_multiplicity_->Fill(id.ieta(), id.iphi());
               ids.insert(id);
            }
         }

         for (const auto& id: ids) {
            tp_multiplicity_->Fill(id.ieta(), id.iphi());
         }

         first_ = false;
      }
   }

   // ==============
   // Matching stuff
   // ==============

   ev_rh_energy0_ = 0.;
   ev_rh_energy2_ = 0.;
   ev_rh_energy3_ = 0.;
   ev_rh_unmatched_ = 0.;
   ev_tp_energy_ = 0.;
   ev_tp_unmatched_ = 0.;

   std::map<HcalTrigTowerDetId, std::vector<HBHERecHit>> rhits;
   std::map<HcalTrigTowerDetId, std::vector<HFRecHit>> fhits;
   std::map<HcalTrigTowerDetId, std::vector<HcalTriggerPrimitiveDigi>> tpdigis;

   Handle<HcalTrigPrimDigiCollection> digis;
   if (!event.getByLabel(digis_, digis)) {
      LogError("HcalTrigPrimDigiCleaner") <<
         "Can't find hcal trigger primitive digi collection with tag '" <<
         digis_ << "'" << std::endl;
      return;
   }

   edm::Handle< edm::SortedCollection<HBHERecHit> > hits;
   if (!event.getByLabel(rechits_[0], hits)) {
      edm::LogError("HcalCompareLegacyChains") <<
         "Can't find rec hit collection with tag '" << rechits_[0] << "'" << std::endl;
      /* return; */
   }

   edm::Handle< edm::SortedCollection<HFRecHit> > hfhits;
   if (!event.getByLabel(rechits_[1], hfhits)) {
      edm::LogError("HcalCompareLegacyChains") <<
         "Can't find rec hit collection with tag '" << rechits_[1] << "'" << std::endl;
      /* return; */
   }

   setup.get<CaloGeometryRecord>().get(gen_geo_);

   auto isValid = [&](const auto& hit) {
      HcalDetId id(hit.id());
      auto s = status_->getValues(id);
      int level = comp_->getSeverityLevel(id, 0, s->getValue());
      return level <= max_severity_;
   };

   if (hits.isValid()) {
      for (auto& hit: *(hits.product())) {
         HcalDetId id(hit.id());
         if (not isValid(hit))
            continue;
         ev_rh_energy0_ += hit.eraw() / get_cosh(id);
         ev_rh_energy2_ += hit.energy() / get_cosh(id);
         ev_rh_energy3_ += hit.eaux() / get_cosh(id);

         auto tower_ids = tpd_geo.towerIds(id);
         for (auto& tower_id: tower_ids) {
            tower_id = HcalTrigTowerDetId(tower_id.ieta(), tower_id.iphi(), 1);
            rhits[tower_id].push_back(hit);
         }
      }
   }

   if (hfhits.isValid()) {
      for (auto& hit: *(hfhits.product())) {
         HcalDetId id(hit.id());
         if (not isValid(hit))
            continue;
         ev_rh_energy0_ += hit.energy() / get_cosh(id);
         ev_rh_energy2_ += hit.energy() / get_cosh(id);
         ev_rh_energy3_ += hit.energy() / get_cosh(id);

         auto tower_ids = tpd_geo.towerIds(id);
         for (auto& tower_id: tower_ids) {
            tower_id = HcalTrigTowerDetId(tower_id.ieta(), tower_id.iphi(), 1, tower_id.version());
            fhits[tower_id].push_back(hit);
         }
      }
   }

   // ESHandle<L1CaloHcalScale> hcal_scale;
   // setup.get<L1CaloHcalScaleRcd>().get(hcal_scale);
   // const L1CaloHcalScale* h = hcal_scale.product();

   ESHandle<CaloTPGTranscoder> decoder;
   setup.get<CaloTPGRecord>().get(decoder);

   for (const auto& digi: *digis) {
      HcalTrigTowerDetId id = digi.id();
      id = HcalTrigTowerDetId(id.ieta(), id.iphi(), 1, id.version());
      ev_tp_energy_ += decoder->hcaletValue(id, digi.t0());
      tpdigis[id].push_back(digi);

      tp_energy_ = decoder->hcaletValue(id, digi.t0());
      tp_ieta_ = id.ieta();
      tp_iphi_ = id.iphi();
      tp_soi_ = digi.SOI_compressedEt();

      tps_->Fill();
   }

   for (const auto& pair: tpdigis) {
      auto id = pair.first;

      auto new_id(id);
      if (swap_iphi_ and id.version() == 1 and id.ieta() > 28 and id.ieta() < 40) {
         if (id.iphi() % 4 == 1)
            new_id = HcalTrigTowerDetId(id.ieta(), (id.iphi() + 70) % 72, id.depth(), id.version());
         else
            new_id = HcalTrigTowerDetId(id.ieta(), (id.iphi() + 2) % 72 , id.depth(), id.version());
      }

      auto rh = rhits.find(new_id);
      auto fh = fhits.find(new_id);

      if (rh != rhits.end() and fh != fhits.end()) {
         assert(0);
      }

      mt_ieta_ = new_id.ieta();
      mt_iphi_ = new_id.iphi();
      mt_version_ = new_id.version();
      mt_tp_energy_ = 0;
      mt_tp_soi_ = 0;

      for (const auto& tp: pair.second) {
         mt_tp_energy_ += decoder->hcaletValue(new_id, tp.t0());
         mt_tp_soi_ = tp.SOI_compressedEt();
      }
      mt_rh_energy0_ = 0.;
      mt_rh_energy2_ = 0.;
      mt_rh_energy3_ = 0.;

      if (rh != rhits.end()) {
         for (const auto& hit: rh->second) {
            HcalDetId id(hit.id());
            auto tower_ids = tpd_geo.towerIds(id);
            auto count = std::count_if(std::begin(tower_ids), std::end(tower_ids),
                  [&](const auto& t) { return t.version() == new_id.version(); });
            mt_rh_energy0_ += hit.eraw() / get_cosh(id) / count;
            mt_rh_energy2_ += hit.energy() / get_cosh(id) / count;
            mt_rh_energy3_ += hit.eaux() / get_cosh(id) / count;
         }
         matches_->Fill();
         rhits.erase(rh);
      } else if (fh != fhits.end()) {
         for (const auto& hit: fh->second) {
            HcalDetId id(hit.id());
            auto tower_ids = tpd_geo.towerIds(id);
            auto count = std::count_if(std::begin(tower_ids), std::end(tower_ids),
                  [&](const auto& t) { return t.version() == new_id.version(); });
            mt_rh_energy0_ += hit.energy() / get_cosh(id) / count;
            mt_rh_energy2_ += hit.energy() / get_cosh(id) / count;
            mt_rh_energy3_ += hit.energy() / get_cosh(id) / count;
         }
         matches_->Fill();
         fhits.erase(fh);
      } else {
         ++ev_tp_unmatched_;
      }
   }

   for (const auto& pair: rhits) {
      ev_rh_unmatched_ += pair.second.size();
   }

   events_->Fill();
}

void
HcalCompareLegacyChains::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HcalCompareLegacyChains);
