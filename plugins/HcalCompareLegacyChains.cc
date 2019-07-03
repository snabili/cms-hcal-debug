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
      TH1D *h_CHPVMPF_HS;
      TH1D *h_CHPVMPF_PU;
      TH1D *h_CHPVMPF;
      TH1D *h_PFRecHEdepth1;
      TH1D *h_PFRecHEdepth2;
      TH1D *h_PFRecHEdepth3;
      TH1D *h_PFRecHEdepth4;
      TH1D *h_PFRecHEdepth5;
      TH1D *h_PFRecHEdepth6;
      TH1D *h_PFRecHEdepth7;

      TH1D *h_PFRecHEdepth1_HS;
      TH1D *h_PFRecHEdepth2_HS;
      TH1D *h_PFRecHEdepth3_HS;
      TH1D *h_PFRecHEdepth4_HS;
      TH1D *h_PFRecHEdepth5_HS;
      TH1D *h_PFRecHEdepth6_HS;
      TH1D *h_PFRecHEdepth7_HS;

      TH1D *h_PFRecHEdepth1_PU;
      TH1D *h_PFRecHEdepth2_PU;
      TH1D *h_PFRecHEdepth3_PU;
      TH1D *h_PFRecHEdepth4_PU;
      TH1D *h_PFRecHEdepth5_PU;
      TH1D *h_PFRecHEdepth6_PU;
      TH1D *h_PFRecHEdepth7_PU;
      TH1D *h_PFRecSum;
      TH1D *h_PFRecSum_HS;
      TH1D *h_PFRecSum_PU;

      //====== 2D Histo Declaration =========
      TH2D *h_depth_RecHE;
      TH2D *h_depth_RecHE_HS;
      TH2D *h_depth_RecHE_PU;

      TTree *tps_;
      TTree *tpsplit_;
      TTree *events_;
      TTree *matches_;

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
   // Charged Hadron Vertex 
   h_CHPVMPF_HS = fs->make<TH1D>("h_CHPVMPF_HS","Charged Hadrons Primary Vertex Z Position", 500,-15,15);
   h_CHPVMPF_PU = fs->make<TH1D>("h_CHPVMPF_PU" ,"Charged Hadrons PFCandidate Vertex Z Position", 500,-15,15);
   h_CHPVMPF = fs->make<TH1D>("h_CHPVMPF","Charged Hadrons Vertex Minus PFCandidate Z Position",500,-15,15);

   //====MET 2D Histograms=================================================================================
   h_PFRecSum = fs->make<TH1D>("h_PFRecSum","PFRecHit Sum in HE", 100,0,700);
   h_PFRecSum_HS = fs->make<TH1D>("h_PFRecSum_HS","PFRecHit Sum in HE for HS", 100,0,700);
   h_PFRecSum_PU = fs->make<TH1D>("h_PFRecSum_PU","PFRecHit Sum in HE for PU", 100,0,700);


   //====Depth Study=================================================================================
   h_depth_RecHE = fs->make<TH2D>("h_depth_RecHE","RechitHE vs Depth",1000,0,8,1000,0,7.5);
   h_depth_RecHE_HS = fs->make<TH2D>("h_depth_RecHE_HS","RechitHE vs Depth for Dz<0.1",1000,0,8,1000,0,7.5);
   h_depth_RecHE_PU = fs->make<TH2D>("h_depth_RecHE_PU","RechitHE vs Depth for Dz>0.1",1000,0,8,1000,0,7.5);

   double arrayBins[1000] = {};
   for(int i(0);i<100;i++){arrayBins[i] = 0.02*i;}
   for(int i(100);i<600;i++){arrayBins[i] = 0.004*i+1.6;}
   for(int i(600);i<1000;i++){arrayBins[i] = 0.006*i+0.4;}


   h_PFRecHEdepth1 = fs->make<TH1D>("h_PFRecHEdepth1","PFRechit depth == 1",500,0,10);
   h_PFRecHEdepth2 = fs->make<TH1D>("h_PFRecHEdepth2","PFRechit depth == 2",500,0,10);
   h_PFRecHEdepth3 = fs->make<TH1D>("h_PFRecHEdepth3","PFRechit depth == 3",500,0,10);
   h_PFRecHEdepth4 = fs->make<TH1D>("h_PFRecHEdepth4","PFRechit depth == 4",500,0,10);
   h_PFRecHEdepth5 = fs->make<TH1D>("h_PFRecHEdepth5","PFRechit depth == 5",500,0,10);
   h_PFRecHEdepth6 = fs->make<TH1D>("h_PFRecHEdepth6","PFRechit depth == 6",500,0,10);
   h_PFRecHEdepth7 = fs->make<TH1D>("h_PFRecHEdepth7","PFRechit depth == 7",500,0,10);

   h_PFRecHEdepth1_HS = fs->make<TH1D>("h_PFRecHEdepth1_HS","PFRechit depth == 1 Hard Scatter",999,arrayBins);
   h_PFRecHEdepth2_HS = fs->make<TH1D>("h_PFRecHEdepth2_HS","PFRechit depth == 2 Hard Scatter",999,arrayBins);
   h_PFRecHEdepth3_HS = fs->make<TH1D>("h_PFRecHEdepth3_HS","PFRechit depth == 3 Hard Scatter",999,arrayBins);
   h_PFRecHEdepth4_HS = fs->make<TH1D>("h_PFRecHEdepth4_HS","PFRechit depth == 4 Hard Scatter",999,arrayBins);
   h_PFRecHEdepth5_HS = fs->make<TH1D>("h_PFRecHEdepth5_HS","PFRechit depth == 5 Hard Scatter",999,arrayBins);
   h_PFRecHEdepth6_HS = fs->make<TH1D>("h_PFRecHEdepth6_HS","PFRechit depth == 6 Hard Scatter",999,arrayBins);
   h_PFRecHEdepth7_HS = fs->make<TH1D>("h_PFRecHEdepth7_HS","PFRechit depth == 7 Hard Scatter",999,arrayBins);

   h_PFRecHEdepth1_PU = fs->make<TH1D>("h_PFRecHEdepth1_PU","PFRechit depth == 1 Pileup",999,arrayBins);
   h_PFRecHEdepth2_PU = fs->make<TH1D>("h_PFRecHEdepth2_PU","PFRechit depth == 2 Pileup",999,arrayBins);
   h_PFRecHEdepth3_PU = fs->make<TH1D>("h_PFRecHEdepth3_PU","PFRechit depth == 3 Pileup",999,arrayBins);
   h_PFRecHEdepth4_PU = fs->make<TH1D>("h_PFRecHEdepth4_PU","PFRechit depth == 4 Pileup",999,arrayBins);
   h_PFRecHEdepth5_PU = fs->make<TH1D>("h_PFRecHEdepth5_PU","PFRechit depth == 5 Pileup",999,arrayBins);
   h_PFRecHEdepth6_PU = fs->make<TH1D>("h_PFRecHEdepth6_PU","PFRechit depth == 6 Pileup",999,arrayBins);
   h_PFRecHEdepth7_PU = fs->make<TH1D>("h_PFRecHEdepth7_PU","PFRechit depth == 7 Pileup",999,arrayBins);

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
   if (!PVs->empty() && !((*PVs)[0].isFake())) {//Offline Primary Vertices
     leadPV = math::XYZPoint((*PVs)[0].x(),(*PVs)[0].y(),(*PVs)[0].z());
   }

   //================ PFRecHits from PFCandidates->PFBlocks->PFClusters->PFRecHitFractions =============================================
   //Process:	1-Loop over all Charged Hadrons PFCandidates. const_iterator ci later pointed to pfc
   //		2-Return elements in blocks by calling the function elementsInBlocks() from the vector ElementsInBlocks. Vector ElementsInBlocks is defined from the pair ElementInBlock. These all are defined in PFCandidate.h (https://github.com/cms-sw/cmssw/blob/CMSSW_10_1_X/DataFormats/ParticleFlowCandidate/interface/PFCandidate.h#L386-L387)
   //		3-Loop on the elements in block associated to the PFCandidate
   //		4-Loop on the first element of ElementInBlock pair: "blockRef" that contains TRACK, PS1, PS2, ECAL, HCAL, GSF, BREM, SC HO (https://github.com/cms-sw/cmssw/blob/CMSSW_10_1_X/DataFormats/ParticleFlowReco/interface/PFBlockElement.h#L33-L48). Type number represents element ID such as TRACK, PS1, PS2 ...
   //		5-Extract the PFCluster from blockelements components
   //		6-Loop over Hits within the PFCluster
   double PFSEt_T26 = 0;
   double PFSEt_T26_D1 = 0;
   double PFSEt_T26_D2 = 0;
   double PFSEt_T26_D3 = 0;
   double PFSEt_T26_D4 = 0;
   double PFSEt_T26_D5 = 0;
   double PFSEt_T26_D6 = 0;
   double PFSEt_T26_D7 = 0;
   unsigned int iphi = 0;
   edm::Handle<reco::PFCandidateCollection> pfCandidate;
   event.getByToken(inputTagPFCandidates_, pfCandidate);
   std::ofstream myfile;//to get the result of the big file and make ntuples
   myfile.open("tuple.txt");//write the results in a text file
   double posX,posY,posZ = 0;//Geometry Info.
   for( reco::PFCandidateCollection::const_iterator ci  = pfCandidate->begin(); ci!=pfCandidate->end(); ++ci)  { //PFCandidate Loop (1)
     const reco::PFCandidate& pfc = *ci;//(1)
     if ( pfc.particleId() != 1 ) continue;// Charged Hadron Particles (1)
     math::XYZPoint vtx = pfc.vertex();//PFCandidate Vertices  
     double dz = fabs(vtx.z()-leadPV.z());//Cut on Hard Scatter and Pileup
     h_CHPVMPF->Fill(dz);
     const PFCandidate::ElementsInBlocks& theElements = pfc.elementsInBlocks();//returns elements in the block associated with the PFCandidate (2)
     for ( reco::PFCandidate::ElementsInBlocks::const_iterator Elements = theElements.begin(); Elements != theElements.end(); ++Elements ) {//Elements in block associated with PFCandidate (3)
       //=== to cout the elements in the block:(https://github.com/cms-sw/cmssw/blob/CMSSW_10_1_X/RecoParticleFlow/Configuration/test/PFChargedHadronAnalyzer.cc#L188-L193) =========
       const PFCandidate::ElementsInBlocks& theElements = pfc.elementsInBlocks();
       const reco::PFBlockRef blockRef = theElements[0].first;// blockRef is the first object of the ElementInBlock pair (4)
       const edm::OwnVector<reco::PFBlockElement>& elements = blockRef->elements();
       for(unsigned iEle=0; iEle<elements.size(); iEle++) {//to give type of elements in blockRef: Track, PS1, PS2 ...
         std::cout<<"PFElementsInBlocks are: "<<elements[iEle].type()<<std::endl;
       }
       //================================================
       const edm::OwnVector<reco::PFBlockElement>& PFBlockRef = Elements->first->elements();
       for ( edm::OwnVector<reco::PFBlockElement>::const_iterator blockComponent = PFBlockRef.begin(); blockComponent != PFBlockRef.end(); ++blockComponent ) {//Loop over block components such as Track, PS1, PS2 ... (4)
         if ( blockComponent->clusterRef().isNonnull() ) {//If clusters refenced to pfBlockelement is valid
           reco::PFClusterRef pfCluster = blockComponent->clusterRef();//get PFCluster from ECAL, HCAL, HO, HFEM, HFHAD (5)
   	   if (pfCluster->layer() == PFLayer::HCAL_ENDCAP){//HCAL ENDCAP
             const std::vector<reco::PFRecHitFraction>& pfRecHitFractions = pfCluster->recHitFractions();
             for ( std::vector<reco::PFRecHitFraction>::const_iterator it = pfRecHitFractions.begin(); it != pfRecHitFractions.end(); ++it ) {//PFCluster Loop (6)
               const reco::PFRecHitRef& pfRecHits = it->recHitRef();
               posX = pfRecHits->position().x();
               posY = pfRecHits->position().y();
               posZ = pfRecHits->position().z();
               math::XYZPoint pflowPos(posX,posY,posZ);
               double Eta    = pflowPos.eta();
	       double Phi    = pflowPos.phi();
	       if(Phi>=0){iphi = Phi*180/5/M_PI;}
	       else{iphi = 72 + Phi*180/5/M_PI;}
	       std::cout<<"Phi is: "<<Phi<<" iphi is: "<<iphi<<std::endl;
	       //if(1.3<fabs(Eta) && fabs(Eta)<=2.5){//Tracker Area
	       if(2.322<=fabs(Eta) && fabs(Eta)<=2.500){//ieta == 26
		 PFSEt_T26 +=pfRecHits->energy()/cosh(Eta); 
	         unsigned int depth = pfRecHits->depth();// from depth =1 to depth =7
                 if(depth == 1){PFSEt_T26_D1 += pfRecHits->energy()/cosh(Eta);}
                 if(depth == 2){PFSEt_T26_D2 += pfRecHits->energy()/cosh(Eta);}
                 if(depth == 3){PFSEt_T26_D3 += pfRecHits->energy()/cosh(Eta);}
                 if(depth == 4){PFSEt_T26_D4 += pfRecHits->energy()/cosh(Eta);}
                 if(depth == 5){PFSEt_T26_D5 += pfRecHits->energy()/cosh(Eta);}
                 if(depth == 6){PFSEt_T26_D6 += pfRecHits->energy()/cosh(Eta);}
                 if(depth == 7){PFSEt_T26_D7 += pfRecHits->energy()/cosh(Eta);}

	       }//tower 26
             }//End loop over PFCluster
           }//Restrict HCAL ENDCAP
           if(fabs(dz)<0.1){myfile<<0<<" "<<PFSEt_T26<<" "<<PFSEt_T26_D1<<" "<<PFSEt_T26_D2<<" "<<PFSEt_T26_D3<<" "<<PFSEt_T26_D4<<" "<<PFSEt_T26_D5<<" "<<PFSEt_T26_D6<<" "<<PFSEt_T26_D7<<std::endl;}//HS
	   else{myfile<<1<<" "<<PFSEt_T26<<" "<<PFSEt_T26_D1<<" "<<PFSEt_T26_D2<<" "<<PFSEt_T26_D3<<" "<<PFSEt_T26_D4<<" "<<PFSEt_T26_D5<<" "<<PFSEt_T26_D6<<" "<<PFSEt_T26_D7<<std::endl;}//PU
         }//Pick PFClusters
       }  //Loop over PFBlockElements
     }  //Loop over PFBlock
   }  //Loop over PFCandidate
   myfile.close();
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
