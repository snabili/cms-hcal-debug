// -*- C++ -*-
//
// Package:    HcalDebug
// Class:      AnalyzeTP
// 
/**\class AnalyzeTP AnalyzeTP.cc HcalDebug/CompareChans/src/AnalyzeTP.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  matthias wolf
//         Created:  Fri Nov 27 11:21:58 CET 2015
// $Id$
//
//


// system include files
#include <memory>
#include <unordered_map>

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
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TTree.h"
//
// class declaration
//

class AnalyzeTP : public edm::EDAnalyzer {
   public:
      explicit AnalyzeTP(const edm::ParameterSet&);
      ~AnalyzeTP();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);

      // ----------member data ---------------------------
      edm::InputTag digis_;
      bool detid_;
      double threshold_;

      int event_;

      TTree *match_;
      int m_ieta_;
      int m_iphi_;
      double old_et_;
      double new_et_;
      int new_count_;

      TTree *tps_;

      int tp_ieta_;
      int tp_iphi_;
      int tp_depth_;
      int tp_version_;
      int tp_soi_;
      double tp_et_;

      TTree *ev_;
      double ev_tp_v0_et_;
      double ev_tp_v1_et_;
};

AnalyzeTP::AnalyzeTP(const edm::ParameterSet& config) :
   edm::EDAnalyzer(),
   digis_(config.getParameter<edm::InputTag>("triggerPrimitives")),
   detid_(config.getUntrackedParameter<bool>("useDetIdForUncompression", true)),
   threshold_(config.getUntrackedParameter<double>("threshold", 0.))
{
   edm::Service<TFileService> fs;

   consumes<HcalTrigPrimDigiCollection>(digis_);

   tps_ = fs->make<TTree>("tps", "Trigger primitives");
   tps_->Branch("event", &event_);
   tps_->Branch("ieta", &tp_ieta_);
   tps_->Branch("iphi", &tp_iphi_);
   tps_->Branch("depth", &tp_depth_);
   tps_->Branch("version", &tp_version_);
   tps_->Branch("soi", &tp_soi_);
   tps_->Branch("et", &tp_et_);

   ev_ = fs->make<TTree>("evs", "Event quantities");
   ev_->Branch("event", &event_);
   ev_->Branch("tp_v0_et", &ev_tp_v0_et_);
   ev_->Branch("tp_v1_et", &ev_tp_v1_et_);

   match_ = fs->make<TTree>("ms", "TP matches");
   match_->Branch("event", &event_);
   match_->Branch("ieta", &m_ieta_);
   match_->Branch("iphi", &m_iphi_);
   match_->Branch("et1x1", &new_et_);
   match_->Branch("et2x3", &old_et_);
   match_->Branch("n1x1", &new_count_);
}

AnalyzeTP::~AnalyzeTP() {}

void
AnalyzeTP::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
   using namespace edm;

   event_ = event.id().event();

   Handle<HcalTrigPrimDigiCollection> digis;
   if (!event.getByLabel(digis_, digis)) {
      LogError("AnalyzeTP") <<
         "Can't find hcal trigger primitive digi collection with tag '" <<
         digis_ << "'" << std::endl;
      return;
   }

   ESHandle<CaloTPGTranscoder> decoder;
   setup.get<CaloTPGRecord>().get(decoder);

   std::unordered_map<int, std::unordered_map<int, double>> old_ets;
   std::unordered_map<int, std::unordered_map<int, double>> new_ets;
   std::unordered_map<int, std::unordered_map<int, int>> new_counts;

   ev_tp_v0_et_ = 0.;
   ev_tp_v1_et_ = 0.;

   ESHandle<HcalTrigTowerGeometry> tpd_geo;
   setup.get<CaloGeometryRecord>().get(tpd_geo);

   std::map<HcalTrigTowerDetId, HcalTriggerPrimitiveDigi> ttids;
   for (const auto& digi: *digis) {
      if (digi.id().version() == 1)
         ttids[digi.id()] = digi;
   }

   for (const auto& digi: *digis) {
      HcalTrigTowerDetId id = digi.id();

      tp_ieta_ = id.ieta();
      tp_iphi_ = id.iphi();
      tp_depth_ = id.depth();
      tp_version_ = id.version();
      tp_soi_ = digi.SOI_compressedEt();
      if (detid_)
         tp_et_ = decoder->hcaletValue(id, digi.t0());
      else
         tp_et_ = decoder->hcaletValue(tp_ieta_, tp_iphi_, tp_soi_);

      if (tp_et_ < threshold_)
         continue;

      tps_->Fill();

      if (tp_version_ == 0 and abs(tp_ieta_) >= 29) {
         ev_tp_v0_et_ += tp_et_;
      } else if (tp_version_ == 1) {
         ev_tp_v1_et_ += tp_et_;
      }

      if (abs(tp_ieta_) >= 29 and tp_version_ == 0) {
         std::set<HcalTrigTowerDetId> matches;
         for (const auto& detid: tpd_geo->detIds(id)) {
            for (const auto& ttid: tpd_geo->towerIds(detid)) {
               if (ttid.version() == 1)
                  matches.insert(ttid);
            }
         }

         m_ieta_ = tp_ieta_;
         m_iphi_ = tp_iphi_;
         new_et_ = 0;
         new_count_ = 0;
         old_et_ = tp_et_;
         for (const auto& m: matches) {
            new_et_ += decoder->hcaletValue(m, ttids[m].t0());
            ++new_count_;
         }
         match_->Fill();
      }
   }

   for (int i = -4; i <= 4; ++i) {
      if (i == 0)
         continue;
      for (int j = 0; j < 18; ++j) {
      }
   }

   ev_->Fill();
}

void
AnalyzeTP::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AnalyzeTP);
