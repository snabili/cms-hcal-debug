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
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <memory>
#include <math.h>
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

//PFJetCollection to get the  leading jet 
#include "DataFormats/JetReco/interface/PFJetCollection.h"


#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TTree.h"
#include <stdio.h>
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
      edm::EDGetTokenT<reco::PFJetCollection>      jetToken_;
      edm::EDGetTokenT<reco::VertexCollection>    pvToken_;

      edm::EDGetTokenT< reco::PFJetCollection > jetCollectionToken_;// to get the jet leading 

      bool swap_iphi_;

      TH2D *df_multiplicity_;
      TH2D *tp_multiplicity_;

      //====== 1D Histo Declaration =========
      TH1D *h_CHPVMPF_HS;
      TH1D *h_CHPVMPF_PU;
      TH1D *h_CHPVMPF;

      TH1D *h_PFRec_PU_D1;
      TH1D *h_PFRec_PU_D2;
      TH1D *h_PFRec_PU_D3;
      TH1D *h_PFRec_PU_D4;
      TH1D *h_PFRec_PU_D5;
      TH1D *h_PFRec_PU_D6;
      TH1D *h_PFRec_PU_D7;

      TH1D *h_PFRec_HS_D1;
      TH1D *h_PFRec_HS_D2;
      TH1D *h_PFRec_HS_D3;
      TH1D *h_PFRec_HS_D4;
      TH1D *h_PFRec_HS_D5;
      TH1D *h_PFRec_HS_D6;
      TH1D *h_PFRec_HS_D7;

      TH1D *h_RecMET;
      TH1D *h_L1SumEt;
      TH1D *h_L1MET;
      TH1D *h_PFMET;
      TH1D *h_PFRec_HS;
      TH1D *h_PFRec_PU;

      TH1D *h_pfjet;

      //====== 2D Histo Declaration =========
      /*TH2D *h_depth_RecHE;
      TH2D *h_depth_RecHE_HS;
      TH2D *h_depth_RecHE_PU;*/

      TTree *tps_;
      TTree *tpsplit_;
      TTree *events_;
      TTree *matches_;
      TTree *PFRecHit_PU;
      TTree *PFRecHit_HS;

      TTree *signalFile_;
      TTree *bckgndFile_;


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

      double PU_D1;
      double PU_D2;
      double PU_D3;
      double PU_D4;
      double PU_D5;
      double PU_D6;
      double PU_D7;

      double HS_D1;
      double HS_D2;
      double HS_D3;
      double HS_D4;
      double HS_D5;
      double HS_D6;
      double HS_D7;

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

   consumes<reco::PFJetCollection> (edm::InputTag("ak4PFJetsCHS"));

   sumToken_ = consumes<l1t::EtSumBxCollection>(config.getUntrackedParameter<edm::InputTag>("sumToken"));
   metToken_ = consumes<reco::PFMETCollection>(config.getUntrackedParameter<edm::InputTag>("metToken"));

   jetToken_ = consumes<reco::PFJetCollection>(config.getUntrackedParameter<edm::InputTag>("jetToken"));

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

   PFRecHit_PU = fs->make<TTree>("PFRecHit_PU","PFRecHit_PU depth");
   PFRecHit_PU->Branch("PFRecHit_PU_D1", &PU_D1);
   PFRecHit_PU->Branch("PFRecHit_PU_D1", &PU_D2);
   PFRecHit_PU->Branch("PFRecHit_PU_D1", &PU_D3);
   PFRecHit_PU->Branch("PFRecHit_PU_D1", &PU_D4);
   PFRecHit_PU->Branch("PFRecHit_PU_D1", &PU_D5);
   PFRecHit_PU->Branch("PFRecHit_PU_D1", &PU_D6);
   PFRecHit_PU->Branch("PFRecHit_PU_D1", &PU_D7);

   PFRecHit_HS = fs->make<TTree>("PFRecHit_HS","PFRecHit_HS depth");
   PFRecHit_HS->Branch("PFRecHit_HS_D1", &HS_D1);
   PFRecHit_HS->Branch("PFRecHit_HS_D1", &HS_D2);
   PFRecHit_HS->Branch("PFRecHit_HS_D1", &HS_D3);
   PFRecHit_HS->Branch("PFRecHit_HS_D1", &HS_D4);
   PFRecHit_HS->Branch("PFRecHit_HS_D1", &HS_D5);
   PFRecHit_HS->Branch("PFRecHit_HS_D1", &HS_D6);
   PFRecHit_HS->Branch("PFRecHit_HS_D1", &HS_D7);

   //====MET 1D Histograms=================================================================================
   // Charged Hadron Vertex 
   h_CHPVMPF_HS = fs->make<TH1D>("h_CHPVMPF_HS","Charged Hadrons Primary Vertex Z Position", 500,-15,15);
   h_CHPVMPF_PU = fs->make<TH1D>("h_CHPVMPF_PU" ,"Charged Hadrons PFCandidate Vertex Z Position", 500,-15,15);
   h_CHPVMPF = fs->make<TH1D>("h_CHPVMPF","Charged Hadrons Vertex Minus PFCandidate Z Position",500,-15,15);

   //====SumEt=================================================================================
   h_L1SumEt = fs->make<TH1D>("L1SumEt","L1SumEt",100,-10,2500);

   //====MET 2D Histograms=================================================================================
   /*h_PFRecSum = fs->make<TH1D>("h_PFRecSum","PFRecHit Sum in HE", 100,0,700);
   h_PFRecSum_HS = fs->make<TH1D>("h_PFRecSum_HS","PFRecHit Sum in HE for HS", 100,0,700);
   h_PFRecSum_PU = fs->make<TH1D>("h_PFRecSum_PU","PFRecHit Sum in HE for PU", 100,0,700);*/


   //====Depth Study=================================================================================
   h_PFRec_PU_D1 = fs->make<TH1D>("h_PFRec_PU_D1","",100,0,5);
   h_PFRec_PU_D2 = fs->make<TH1D>("h_PFRec_PU_D2","",100,0,5);
   h_PFRec_PU_D3 = fs->make<TH1D>("h_PFRec_PU_D3","",100,0,5);
   h_PFRec_PU_D4 = fs->make<TH1D>("h_PFRec_PU_D4","",100,0,5);
   h_PFRec_PU_D5 = fs->make<TH1D>("h_PFRec_PU_D5","",100,0,5);
   h_PFRec_PU_D6 = fs->make<TH1D>("h_PFRec_PU_D6","",100,0,5);
   h_PFRec_PU_D7 = fs->make<TH1D>("h_PFRec_PU_D7","",100,0,5);

   h_PFRec_HS_D1 = fs->make<TH1D>("h_PFRec_HS_D1","",100,0,5);
   h_PFRec_HS_D2 = fs->make<TH1D>("h_PFRec_HS_D2","",100,0,5);
   h_PFRec_HS_D3 = fs->make<TH1D>("h_PFRec_HS_D3","",100,0,5);
   h_PFRec_HS_D4 = fs->make<TH1D>("h_PFRec_HS_D4","",100,0,5);
   h_PFRec_HS_D5 = fs->make<TH1D>("h_PFRec_HS_D5","",100,0,5);
   h_PFRec_HS_D6 = fs->make<TH1D>("h_PFRec_HS_D6","",100,0,5);
   h_PFRec_HS_D7 = fs->make<TH1D>("h_PFRec_HS_D7","",100,0,5);



   h_RecMET = fs->make<TH1D>("h_RecMET","",500,0,300);
   h_L1MET = fs->make<TH1D>("h_L1MET","",500,0,300);
   h_PFMET = fs->make<TH1D>("h_PFMET","",500,0,300);
   h_PFRec_HS  = fs->make<TH1D>("h_PFRec_HS","",100,0.1,3);
   h_PFRec_PU  = fs->make<TH1D>("h_PFRec_PU","",100,0.1,3);
   h_pfjet  = fs->make<TH1D>("h_pfjet","",100,0,50);

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

   //============PFJet ===================================================================== 
   edm::Handle<reco::PFJetCollection> pfjetH;
   edm::Handle<reco::PFJetCollection> selectedJets_;
   event.getByToken(jetToken_, pfjetH);
   selectedJets_ = pfjetH;
   for ( reco::PFJetCollection::const_iterator jet = selectedJets_->begin(); jet != selectedJets_->end(); jet++ ) {
	double jeteta = jet->eta();
	if(1.305<=fabs(jeteta))h_pfjet->Fill(jet->pt());
   } 

   ievt++;
   //std::cout<<"Event number: "<<ievt<<" started"<<std::endl;

   //============Offline PFMET ===================================================================== 
   edm::Handle<reco::PFMETCollection> metLabel_;
   event.getByToken(metToken_, metLabel_);
   const reco::PFMETCollection *metCol = metLabel_.product();
   const reco::PFMET theMet = metCol->front();
   met_data->met     = theMet.et();
   met_data->metPhi  = theMet.phi();
   met_data->sumEt   = theMet.sumEt();
   h_PFMET->Fill(theMet.et());

   //============Online MET =====================================================================
   double L1MET = 0;
   //double L1SumEt = 0;
   edm::Handle<l1t::EtSumBxCollection> sums;
   event.getByToken(sumToken_, sums);
   if(sums.isValid()){
     for (int ibx = sums->getFirstBX(); ibx <= sums->getLastBX(); ++ibx) {
       for (l1t::EtSumBxCollection::const_iterator it=sums->begin(ibx); it!=sums->end(ibx); it++) {
         int type = static_cast<int>( it->getType() );
         if (type == 8 && ibx ==0){L1MET = it->et(); h_L1MET->Fill(L1MET);}
         //if (type == 2 && ibx ==0){double L1MET_HBE = it->et();}
         if (type == 0 && ibx ==0){h_L1SumEt->Fill(it->et());
	   //L1SumEt = it->et();
	 }
       }
     }
   }
   else{std::cout<<" invalid sums"<<std::endl;}

   //================ Primary Vertex Position =============================================
   edm::Handle<reco::VertexCollection> PVs;
   event.getByToken(pvToken_, PVs);
   math::XYZPoint leadPV(0,0,0);
   if (!PVs->empty() && !((*PVs)[0].isFake())) {//Offline Primary Vertices
     leadPV = math::XYZPoint((*PVs)[0].x(),(*PVs)[0].y(),(*PVs)[0].z());
   }

   double Eta = 0.0;
   //double Phi = 0.0;
   unsigned int ieta = 0;
   unsigned int iphi_he = 0;
   unsigned int depth = 0;
   unsigned int f = 0;
   unsigned int g = 0;
   unsigned int pfcandidate = 0;
   std::vector<unsigned int> Geometry;
   std::vector<double> PFRecHit;
   double PFRec[11][72][7] = {};
   double PFRec_Pt[11][72][7] = {};

   //================ PFRecHits from PFCandidates->PFBlocks->PFClusters->PFRecHitFractions =============================================
   //Process:	1-Loop over all Charged Hadrons PFCandidates. const_iterator ci later pointed to pfc
   //		2-Return elements in blocks by calling the function elementsInBlocks() from the vector ElementsInBlocks. Vector ElementsInBlocks is defined from the pair ElementInBlock. These all are defined in PFCandidate.h (https://github.com/cms-sw/cmssw/blob/CMSSW_10_1_X/DataFormats/ParticleFlowCandidate/interface/PFCandidate.h#L386-L387)
   //		3-Loop on the elements in block associated to the PFCandidate
   //		4-Loop on the first element of ElementInBlock pair: "blockRef" that contains TRACK, PS1, PS2, ECAL, HCAL, GSF, BREM, SC HO (https://github.com/cms-sw/cmssw/blob/CMSSW_10_1_X/DataFormats/ParticleFlowReco/interface/PFBlockElement.h#L33-L48). Type number represents element ID such as TRACK, PS1, PS2 ...
   //		5-Extract the PFCluster from blockelements components
   //		6-Loop over Hits within the PFCluster
   //================================================================================================================================================	


   edm::Handle<reco::PFCandidateCollection> pfCandidate;
   event.getByToken(inputTagPFCandidates_, pfCandidate);
   std::ofstream file;//to get the result of the big file and make ntuples
   file.open("file_4.txt");
   double posX,posY,posZ = 0;//Geometry Info.
   for( reco::PFCandidateCollection::const_iterator ci  = pfCandidate->begin(); ci!=pfCandidate->end(); ++ci)  { //PFCandidate Loop (1) [PFCandidates are the list of particle PF algorithm produces]
     const reco::PFCandidate& pfc = *ci;//(1)G
     if ( pfc.particleId() != 1 ) continue;// Charged Hadron Particles (1)
     //if(fabs(pfc.eta())<1.305 || 2.500<fabs(pfc.eta())) continue;
     PFRec[11][72][7] = {};
     math::XYZPoint vtx = pfc.vertex();//PFCandidate Vertices  
     double dz1 = vtx.z()-leadPV.z();//Cut on Hard Scatter and Pileup
     h_CHPVMPF->Fill(dz1);
     if(fabs(dz1)<0.1){h_CHPVMPF_HS->Fill(dz1);}
     else{h_CHPVMPF_PU->Fill(dz1);}
     pfcandidate++;//number of PFCandidates charged hadron particles
     const PFCandidate::ElementsInBlocks& theElements = pfc.elementsInBlocks();//returns elements in the block associated with the PFCandidate (2)
     for ( reco::PFCandidate::ElementsInBlocks::const_iterator Elements = theElements.begin(); Elements != theElements.end(); ++Elements ) {//Elements in block associated with PFCandidate (3)
	//std::cout<<"here the element in block loop starts"<<std::endl;
       const edm::OwnVector<reco::PFBlockElement>& PFBlockRef = Elements->first->elements();
       for ( edm::OwnVector<reco::PFBlockElement>::const_iterator blockComponent = PFBlockRef.begin(); blockComponent != PFBlockRef.end(); ++blockComponent ) {//Loop over block components such as Track, PS1, PS2 ... (4)
         if ( blockComponent->clusterRef().isNull()) continue;//If clusters refenced to pfBlockelement is not empty
         reco::PFClusterRef pfCluster = blockComponent->clusterRef();//get PFCluster from ECAL, HCAL, HO, HFEM, HFHAD (5)
	 if (pfCluster->layer() != PFLayer::HCAL_ENDCAP) continue;
         const std::vector<reco::PFRecHitFraction>& pfRecHitFractions = pfCluster->recHitFractions();
         for ( std::vector<reco::PFRecHitFraction>::const_iterator it = pfRecHitFractions.begin(); it != pfRecHitFractions.end(); ++it ) {//PFCluster Loop (6)
	   //file.open("file_4.txt");
	   bool repetetive = false;
           const reco::PFRecHitRef& pfRecHits = it->recHitRef();
           DetId idseed;//reference from:     https://github.com/cms-sw/cmssw/blob/master/RecoParticleFlow/PFClusterTools/src/PFPhotonClusters.cc#L27
           idseed = pfRecHits->detId();
           posX = pfRecHits->position().x();
           posY = pfRecHits->position().y();
           posZ = pfRecHits->position().z();
           math::XYZPoint pflowPos(posX,posY,posZ);
           Eta    = pflowPos.eta();
           //Phi    = pflowPos.phi();
           HcalDetId HEidSeed = HcalDetId(idseed.rawId());
           iphi_he = HEidSeed.iphi();
	   double PFRecHit_Et = pfRecHits->energy()/cosh(Eta);
	   double PFRecHit_Pt = pfc.pt()/cosh(Eta);
	   depth = pfRecHits->depth();// from depth =1 to depth =7
	     //ieta definition
	   if(1.305<=Eta && Eta<=2.500){//Restricting the area between 1.305<=|Eta|<=2.500 tracker region
	     if(2.322<Eta && Eta<=2.500){ieta = 1;}//ieta == 26
	     if(2.172<Eta && Eta<=2.322){ieta = 2;}//ieta == 25
	     if(2.043<Eta && Eta<=2.172){ieta = 3;}//ieta == 24
	     if(1.930<Eta && Eta<=2.043){ieta = 4;}//ieta == 23
	     if(1.830<Eta && Eta<=1.930){ieta = 5;}//ieta == 22
	     if(1.740<Eta && Eta<=1.830){ieta = 6;}//ieta == 21
	     if(1.653<Eta && Eta<=1.740){ieta = 7;}//ieta == 20
	     if(1.566<Eta && Eta<=1.653){ieta = 8;}//ieta == 19
	     if(1.479<Eta && Eta<=1.566){ieta = 9;}//ieta == 18
	     if(1.392<Eta && Eta<=1.479){ieta = 10;}//ieta == 17
	     if(1.305<=Eta && Eta<=1.392){ieta = 11;}//ieta == 16
	     Geometry.push_back(iphi_he);//(3*n)th element is iphi
             Geometry.push_back(ieta);//(3*n+1)th element is ieta
             Geometry.push_back(depth);//(3*n+2)th element is depth
             PFRecHit.push_back(pfRecHits->energy()/cosh(Eta));
	        //std::cout<<"HITS:           	"<<"    	ieta:   "<<ieta<<"   iphi:   "<<iphi_he<<"     	depth:  "<<depth<<std::endl;
             for(unsigned int i(0); i<(Geometry.size()/3-1); i++){
               if(iphi_he == Geometry[3*i] && ieta == Geometry[3*i+1] && depth == Geometry[3*i+2]){
	         repetetive = true;
	         f++;
	       }
	       else{repetetive = false;}
		//std::cout<<"GEOMETRY:		"<<"	ieta:	"<<Geometry[3*i+1]<<"	iphi:	"<<Geometry[3*i]<<"	depth:	"<<Geometry[3*i+2]<<std::endl;

	       if(repetetive == true) break;
             }
             if(repetetive == false){
               g++;
 	       PFRec[ieta-1][iphi_he-1][depth-1] += PFRecHit_Et;
	       PFRec_Pt[ieta-1][iphi_he-1][depth-1] += PFRecHit_Pt;
		//std::cout<<"PFRec:	*******		"<<PFRec[ieta-1][iphi_he-1][depth-1]<<std::endl;
	     }
	   }
         }//End loop over PFRechits or PFCluster Loop (6)
       }  //Loop over PFBlockElements
     }  //Loop over PFBlock
     for(unsigned int i=0; i<11; i++){
       for(unsigned int j=0; j<72; j++){
         for(unsigned int k=0; k<7; k++){
	   //std::cout<<3 <<"        "<<26-i<<"      "<<k+1<<"       "<<j+1<<"       "<<PFRec[i][j][k]<<"    "<<dz1<<std::endl;
	   if(fabs(dz1)<0.1 && PFRec[i][j][k] > 0.0){
	     file<<1 <<"	"<<26-i<<"	"<<k+1<<"	"<<j+1<<"	"<<PFRec[i][j][k]<<"	"<<dz1<<std::endl;
	     h_PFRec_HS->Fill(PFRec[i][j][k]);
	     h_PFRec_HS_D1->Fill(PFRec[i][j][0]);
             h_PFRec_HS_D2->Fill(PFRec[i][j][1]);
             h_PFRec_HS_D3->Fill(PFRec[i][j][2]);
             h_PFRec_HS_D4->Fill(PFRec[i][j][3]);
             h_PFRec_HS_D5->Fill(PFRec[i][j][4]);
             h_PFRec_HS_D6->Fill(PFRec[i][j][5]);
             h_PFRec_HS_D7->Fill(PFRec[i][j][6]);

	     HS_D1 = PFRec[i][j][0];
             HS_D2 = PFRec[i][j][1];
             HS_D3 = PFRec[i][j][2];
             HS_D4 = PFRec[i][j][3];
             HS_D5 = PFRec[i][j][4];
             HS_D6 = PFRec[i][j][5];
             HS_D7 = PFRec[i][j][6];
	     
	     std::cout<<1 <<"        "<<26-i<<"      "<<k+1<<"       "<<j+1<<"       "<<PFRec[i][j][k]<<"    "<<dz1<<"	"<<PFRec_Pt[i][j][k]<<std::endl;
	   }
	   if(fabs(dz1)>0.1 && PFRec[i][j][k] > 0.0){
	     file<<2 <<"	"<<26-i<<"	"<<k+1<<"	"<<j+1<<"	"<<PFRec[i][j][k]<<"	"<<dz1<<"	"<<PFRec_Pt[i][j][k]<<std::endl;
             h_PFRec_PU->Fill(PFRec[i][j][k]);
             h_PFRec_PU_D1->Fill(PFRec[i][j][0]);
             h_PFRec_PU_D2->Fill(PFRec[i][j][1]);
             h_PFRec_PU_D3->Fill(PFRec[i][j][2]);
             h_PFRec_PU_D4->Fill(PFRec[i][j][3]);
             h_PFRec_PU_D5->Fill(PFRec[i][j][4]);
             h_PFRec_PU_D6->Fill(PFRec[i][j][5]);
             h_PFRec_PU_D7->Fill(PFRec[i][j][6]);

             PU_D1 = PFRec[i][j][0];
             PU_D2 = PFRec[i][j][1];
             PU_D3 = PFRec[i][j][2];
             PU_D4 = PFRec[i][j][3];
             PU_D5 = PFRec[i][j][4];
             PU_D6 = PFRec[i][j][5];
             PU_D7 = PFRec[i][j][6];

             std::cout<<2 <<"        "<<26-i<<"      "<<k+1<<"       "<<j+1<<"       "<<PFRec[i][j][k]<<"    "<<dz1<<"       "<<PFRec_Pt[i][j][k]<<std::endl;
	   }
         }
       }
     }
   }  //Loop over PFCandidate
   //std::cout<<"pfcandidate:	"<<pfcandidate<<std::endl;
   
   file.close();
   //std::cout<<"true repetation is: "<<f<<"  false repetation is: "<<g<<std::endl;

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

   //==================================RecHit MET Computation ==================================
   std::ofstream myfile;//to get the result of the big file and make ntuples
   myfile.open("jetht3_3.txt");
   double EE_Ex, EE_Ey, EB_Ex, EB_Ey = 0;
   double HF_Ex, HF_Ey, HBE_Ex, HBE_Ey = 0 ;

   //============Computing HE Rechit Depth Segmentation Using Method Zero and MAHI ==================================

   Handle<HBHERecHitCollection> hRecHits;
   event.getByToken(hRhToken, hRecHits);

   //==================================Computing Rechit in HBHE ==================================

   double RHE_depth_[7];
   double RH_Dp_tot = 0;
   int nrhbe = 0;

   if (hits.isValid()) {
      for (auto& hit: *(hits.product())) {
         HcalDetId id(hit.id());
         if (not isValid(hit)) continue;
	 unsigned int depthHE = id.depth();
         RHE_depth_[depthHE] = hit.energy() / get_cosh(id);
         RH_Dp_tot += hit.energy() / get_cosh(id);

         //==================================  Rechit MET in HBHE  ==================================
         if(RHE_depth_[depthHE]>0){
           double phi_HBHE = id.iphi()*5*M_PI/180;
           HBE_Ex += hit.energy() / get_cosh(id) * cos(phi_HBHE) ;
           HBE_Ey += hit.energy() / get_cosh(id) * sin(phi_HBHE) ;
         }
         myfile<<1<<"  "<<HBE_Ex<<"   "<<HBE_Ey<<std::endl;

         ev_rh_energy0_ += hit.eraw() / get_cosh(id);
         ev_rh_energy2_ += hit.energy() / get_cosh(id);
         ev_rh_energy3_ += hit.eaux() / get_cosh(id);

         auto tower_ids = tpd_geo.towerIds(id);
         for (auto& tower_id: tower_ids) {
            tower_id = HcalTrigTowerDetId(tower_id.ieta(), tower_id.iphi(), 1);
            rhits[tower_id].push_back(hit);
         }
	 nrhbe++;
      }
   }

   //==================================  Computing Rechit in HF ==================================
    
   if (hfhits.isValid()) {
      for (auto& hit: *(hfhits.product())) {
         HcalDetId id(hit.id());
         if (not isValid(hit)) continue;
         ev_rh_energy0_ += hit.energy() / get_cosh(id);
         ev_rh_energy2_ += hit.energy() / get_cosh(id);
         ev_rh_energy3_ += hit.energy() / get_cosh(id);

         double RechitEt_M2_hf = hit.energy()/get_cosh(id);
         double phi_HF = id.iphi()*5*M_PI/180;
         HF_Ex += RechitEt_M2_hf * cos(phi_HF) ;
         HF_Ey += RechitEt_M2_hf * sin(phi_HF) ;
         auto tower_ids = tpd_geo.towerIds(id);
         for (auto& tower_id: tower_ids) {
            tower_id = HcalTrigTowerDetId(tower_id.ieta(), tower_id.iphi(), 1, tower_id.version());
            fhits[tower_id].push_back(hit);
         }
      }
   myfile<<2<<"  "<<HF_Ex<<"  "<<HF_Ey<<std::endl;
   }

     edm::Handle<EcalRecHitCollection> ecalhits;

     if (!event.getByToken(tok_EB_, ecalhits)) {
       edm::LogError("HcalCompareLegacyChains") <<
         "Can't find rec hit collection" << std::endl;
        /* return; */
     }
   //==================================  Computing Rechit in Ecal Barrel ==================================
   edm::Handle<EcalRecHitCollection> ecalEB;
   event.getByToken(tok_EB_,ecalEB);
   edm::ESHandle<CaloGeometry> EcalBarrelGeometry_;
   setup.get<CaloGeometryRecord>().get(EcalBarrelGeometry_);
   const CaloGeometry* geo_ecalb = EcalBarrelGeometry_.product();

   if (ecalhits.isValid()) {
     for (EcalRecHitCollection::const_iterator eItr=ecalEB->begin(); eItr!=ecalEB->end(); eItr++)  {
       const GlobalPoint& pos = geo_ecalb->getPosition(eItr->detid());
       double phihit = pos.phi();
       double etahit = pos.eta();
       double RechitEt_M2_eb = eItr->energy()/cosh(etahit);
       EB_Ex += RechitEt_M2_eb * cos(phihit) ;
       EB_Ey += RechitEt_M2_eb * sin(phihit) ;
     }
     myfile<<3<<"  "<<EB_Ex<<"  "<<EB_Ey<<std::endl;
   }

   //================================== Computing Rechit in Ecal Endcap ==================================
   edm::Handle<EcalRecHitCollection> ecalEE;
   event.getByToken(tok_EE_,ecalEE);
   edm::ESHandle<CaloGeometry> EcalEndcapGeometry_;
   setup.get<CaloGeometryRecord>().get(EcalEndcapGeometry_);
   const CaloGeometry* geo_ecale = EcalEndcapGeometry_.product();

   if (ecalhits.isValid()) {
         for (EcalRecHitCollection::const_iterator eItr=ecalEE->begin(); eItr!=ecalEE->end(); eItr++)  {
           const GlobalPoint& pos = geo_ecale->getPosition(eItr->detid());
           double phihit = pos.phi();
           double etahit = pos.eta();
           double RechitEt_M2_ee = eItr->energy()/cosh(etahit);
           EE_Ex += RechitEt_M2_ee * cos(phihit) ;
           EE_Ey += RechitEt_M2_ee * sin(phihit) ;
         }
   myfile<<4<<"  "<<EE_Ex<<"  "<<EE_Ey<<std::endl;
   }

   // ================================== MET Computation ==================================
   double Ex = HBE_Ex + HF_Ex + EB_Ex + EE_Ex;
   double Ey = HBE_Ey + HF_Ey + EB_Ey + EE_Ey;

   double RecMET = sqrt(Ex*Ex+Ey*Ey);

   h_RecMET->Fill(RecMET);
   myfile<<10<<"  "<<Ex<<"  "<<Ey<<"  "<<RecMET<<std::endl;
   myfile<<100<<"  "<<HBE_Ex<<"  "<<HF_Ex<<"  "<<EB_Ex<<"  "<<EE_Ex<<std::endl;
   myfile<<200<<"  "<<HBE_Ey<<"  "<<HF_Ey<<"  "<<EB_Ey<<"  "<<EE_Ey<<std::endl;
   //myfile.close();

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
   myfile.close();
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
