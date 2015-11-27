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
   digis_(config.getParameter<edm::InputTag>("triggerPrimitives"))
{
   edm::Service<TFileService> fs;

   consumes<HcalTrigPrimDigiCollection>(digis_);

   tps_ = fs->make<TTree>("tps", "Trigger primitives");
   tps_->Branch("ieta", &tp_ieta_);
   tps_->Branch("iphi", &tp_iphi_);
   tps_->Branch("depth", &tp_depth_);
   tps_->Branch("version", &tp_version_);
   tps_->Branch("soi", &tp_soi_);
   tps_->Branch("et", &tp_et_);

   ev_ = fs->make<TTree>("evs", "Event quantities");
   ev_->Branch("tp_v0_et", &ev_tp_v0_et_);
   ev_->Branch("tp_v1_et", &ev_tp_v1_et_);
}

AnalyzeTP::~AnalyzeTP() {}

void
AnalyzeTP::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
   using namespace edm;

   Handle<HcalTrigPrimDigiCollection> digis;
   if (!event.getByLabel(digis_, digis)) {
      LogError("AnalyzeTP") <<
         "Can't find hcal trigger primitive digi collection with tag '" <<
         digis_ << "'" << std::endl;
      return;
   }

   ESHandle<CaloTPGTranscoder> decoder;
   setup.get<CaloTPGRecord>().get(decoder);

   ev_tp_v0_et_ = 0.;
   ev_tp_v1_et_ = 0.;

   for (const auto& digi: *digis) {
      HcalTrigTowerDetId id = digi.id();

      tp_ieta_ = id.ieta();
      tp_iphi_ = id.iphi();
      tp_depth_ = id.depth();
      tp_version_ = id.version();
      tp_soi_ = digi.SOI_compressedEt();
      tp_et_ = decoder->hcaletValue(id, digi.t0());

      if (tp_version_ == 0 and abs(tp_ieta_) >= 29) {
         ev_tp_v0_et_ += tp_et_;
      } else if (tp_version_ == 1) {
         ev_tp_v1_et_ += tp_et_;
      }

      tps_->Fill();
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
