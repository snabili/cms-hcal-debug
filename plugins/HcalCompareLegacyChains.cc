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

#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"

#include "CondFormats/DataRecord/interface/HcalChannelQualityRcd.h"
#include "CondFormats/DataRecord/interface/L1CaloGeometryRecord.h"
#include "CondFormats/HcalObjects/interface/HcalChannelQuality.h"
#include "CondFormats/L1TObjects/interface/L1CaloGeometry.h"

#include "CondFormats/L1TObjects/interface/L1RCTParameters.h"
#include "CondFormats/DataRecord/interface/L1RCTParametersRcd.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloHcalScaleRcd.h"

#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputer.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputerRcd.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TTree.h"
//
// class declaration
//

class HcalCompareLegacyChains : public edm::EDAnalyzer {
   public:
      explicit HcalCompareLegacyChains(const edm::ParameterSet&);
      ~HcalCompareLegacyChains();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);

      // ----------member data ---------------------------
      bool first_;

      std::vector<edm::InputTag> frames_;
      edm::InputTag digis_;
      edm::InputTag rechits_;

      TH2D *df_multiplicity_;
      TH2D *tp_multiplicity_;

      TTree *tps_;
      TTree *tpsplit_;
      TTree *events_;
      TTree *matches_;

      double tp_energy_;
      int tp_ieta_;
      int tp_iphi_;
      int tp_depth_max_;
      int tp_depth_start_;
      int tp_depth_end_;

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

      double ev_rh_energy_;
      double ev_tp_energy_;
      int ev_rh_unmatched_;
      int ev_tp_unmatched_;

      double mt_rh_energy_;
      double mt_tp_energy_;

      int mt_ieta_;
      int mt_iphi_;
};

HcalCompareLegacyChains::HcalCompareLegacyChains(const edm::ParameterSet& config) :
   edm::EDAnalyzer(),
   first_(true),
   frames_(config.getParameter<std::vector<edm::InputTag>>("DataFrames")),
   digis_(config.getParameter<edm::InputTag>("TriggerPrimitives")),
   rechits_(config.getParameter<edm::InputTag>("RecHits"))
{
   edm::Service<TFileService> fs;

   df_multiplicity_ = fs->make<TH2D>("df_multiplicity", "DataFrame multiplicity;ieta;iphi", 65, -32.5, 32.5, 72, 0.5, 72.5);
   tp_multiplicity_ = fs->make<TH2D>("tp_multiplicity", "TrigPrim multiplicity;ieta;iphi", 65, -32.5, 32.5, 72, 0.5, 72.5);

   tps_ = fs->make<TTree>("tps", "Trigger primitives");
   tps_->Branch("et", &tp_energy_);
   tps_->Branch("ieta", &tp_ieta_);
   tps_->Branch("iphi", &tp_iphi_);
   tps_->Branch("depth_max", &tp_depth_max_);
   tps_->Branch("depth_start", &tp_depth_start_);
   tps_->Branch("depth_end", &tp_depth_end_);

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
   events_->Branch("RH_energy", &ev_rh_energy_);
   events_->Branch("TP_energy", &ev_tp_energy_);
   events_->Branch("RH_unmatched", &ev_rh_unmatched_);
   events_->Branch("TP_unmatched", &ev_tp_unmatched_);

   matches_ = fs->make<TTree>("matches", "Matched RH and TP");
   matches_->Branch("RH_energy", &mt_rh_energy_);
   matches_->Branch("TP_energy", &mt_tp_energy_);
   matches_->Branch("ieta", &mt_ieta_);
   matches_->Branch("iphi", &mt_iphi_);
}

HcalCompareLegacyChains::~HcalCompareLegacyChains() {}

void
HcalCompareLegacyChains::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
   using namespace edm;

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

   ev_rh_energy_ = 0.;
   ev_rh_unmatched_ = 0.;
   ev_tp_energy_ = 0.;
   ev_tp_unmatched_ = 0.;

   std::map<HcalTrigTowerDetId, std::vector<HBHERecHit>> rhits;
   std::map<HcalTrigTowerDetId, std::vector<HcalTriggerPrimitiveDigi>> tpdigis;

   Handle<HcalTrigPrimDigiCollection> digis;
   if (!event.getByLabel(digis_, digis)) {
      LogError("HcalTrigPrimDigiCleaner") <<
         "Can't find hcal trigger primitive digi collection with tag '" <<
         digis_ << "'" << std::endl;
      return;
   }

   edm::Handle< edm::SortedCollection<HBHERecHit> > hits;
   if (!event.getByLabel(rechits_, hits)) {
      edm::LogError("HcalCompareLegacyChains") <<
         "Can't find rec hit collection with tag '" << rechits_ << "'" << std::endl;
      /* return; */
   }

   edm::ESHandle<CaloGeometry> gen_geo;
   setup.get<CaloGeometryRecord>().get(gen_geo);

   edm::ESHandle<HcalChannelQuality> h_status;
   // setup.get<HcalChannelQualityRcd>().get(h_status);
   // const HcalChannelQuality *status = h_status.product();

   // edm::ESHandle<HcalSeverityLevelComputer> h_comp;
   // setup.get<HcalSeverityLevelComputerRcd>().get(h_comp);
   // const HcalSeverityLevelComputer *comp = h_comp.product();

   if (hits.isValid()) {
      for (auto& hit: *(hits.product())) {
         HcalDetId id(hit.id());
         const auto *local_geo = gen_geo->getSubdetectorGeometry(id)->getGeometry(id);
         ev_rh_energy_ += hit.energy() / cosh(local_geo->getPosition().eta());

         auto tower_ids = tpd_geo.towerIds(id);
         for (auto& tower_id: tower_ids) {
            tower_id = HcalTrigTowerDetId(tower_id.ieta(), tower_id.iphi(), 1);
            rhits[tower_id].push_back(hit);
         }
      }
   }

   // ESHandle<L1CaloHcalScale> hcal_scale;
   // setup.get<L1CaloHcalScaleRcd>().get(hcal_scale);
   // const L1CaloHcalScale* h = hcal_scale.product();

   ESHandle<CaloTPGTranscoder> decoder;
   setup.get<CaloTPGRecord>().get(decoder);

   // edm::ESHandle<L1RCTParameters> rct;
   // setup.get<L1RCTParametersRcd>().get(rct);
   // const L1RCTParameters* r = rct.product();

   for (const auto& digi: *digis) {
      HcalTrigTowerDetId id = digi.id();
      ev_tp_energy_ += decoder->hcaletValue(id.ieta(), id.iphi(), digi.SOI_compressedEt());
      tpdigis[id].push_back(digi);

      tp_energy_ = decoder->hcaletValue(id.ieta(), id.iphi(), digi.SOI_compressedEt());
      tp_ieta_ = id.ieta();
      tp_iphi_ = id.iphi();

      tp_depth_start_ = -1;
      tp_depth_end_ = -1;
      tp_depth_max_ = -1;
      // int et_max = 0;
      // int et_sum = 0;
      // for (unsigned int i = 1; i < 6; ++i) {
      //    int depth = digi.SOI_depth_linear(i);
      //    if (depth > 0) {
      //       et_sum += depth;
      //       tp_depth_end_ = i;
      //       if (tp_depth_start_ < 0)
      //          tp_depth_start_ = i;
      //       if (depth > et_max) {
      //          tp_depth_max_ = i;
      //          et_max = depth;
      //       }
      //    }
      // }
      tps_->Fill();

      // if (et_sum > 0) {
      //    /* std::cout << "vvv" << std::endl; */
      //    for (unsigned int i = 1; i < 6; ++i) {
      //       int depth = digi.SOI_depth_linear(i);
      //       tpsplit_energy_ = tp_energy_ * float(depth) / et_sum;
      //       tpsplit_oot_ = tp_energy_ * float(digi.SOI_oot_linear(i)) / et_sum;
      //       /* std::cout << tpsplit_energy_ << std::endl; */
      //       tpsplit_ieta_ = tp_ieta_;
      //       tpsplit_iphi_ = tp_iphi_;
      //       tpsplit_depth_ = i;
      //       tpsplit_ettot_ = tp_energy_;
      //       tpsplit_rise_avg_ = digi.SOI_rising_avg(i);
      //       tpsplit_rise_rms_ = digi.SOI_rising_rms(i);
      //       tpsplit_fall_avg_ = digi.SOI_falling_avg(i);
      //       tpsplit_fall_rms_ = digi.SOI_falling_rms(i);
      //       tpsplit_->Fill();
      //    }
      // }
   }

   for (const auto& pair: tpdigis) {
      auto rh = rhits.find(pair.first);
      if (rh != rhits.end()) {
         mt_ieta_ = pair.first.ieta();
         mt_iphi_ = pair.first.iphi();
         mt_tp_energy_ = 0;
         for (const auto& tp: pair.second) {
            mt_tp_energy_ += decoder->hcaletValue(
                  pair.first.ieta(),
                  pair.first.iphi(),
                  tp.SOI_compressedEt());
         }
         mt_rh_energy_ = 0.;
         for (const auto& hit: rh->second) {
            HcalDetId id(hit.id());
            const auto *local_geo = gen_geo->getSubdetectorGeometry(id)->getGeometry(id);
            auto tower_ids = tpd_geo.towerIds(id);
            mt_rh_energy_ += hit.energy() / cosh(local_geo->getPosition().eta()) / tower_ids.size();
         }
         matches_->Fill();
         rhits.erase(rh);
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
