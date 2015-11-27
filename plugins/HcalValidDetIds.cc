// -*- C++ -*-
//
// Package:    HcalDebug
// Class:      HcalValidDetIds
// 
/**\class HcalValidDetIds HcalValidDetIds.cc HcalDebug/CompareChans/src/HcalValidDetIds.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  matthias wolf
//         Created:  Fri Nov 27 09:45:23 CET 2015

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
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CondFormats/DataRecord/interface/HcalLutMetadataRcd.h"
#include "CondFormats/HcalObjects/interface/HcalLutMetadata.h"

#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

#include "Geometry/CaloTopology/interface/HcalTopology.h"

//
// class declaration
//

class HcalValidDetIds : public edm::EDAnalyzer {
   public:
      explicit HcalValidDetIds(const edm::ParameterSet&);
      ~HcalValidDetIds();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);

      // ----------member data ---------------------------
};

HcalValidDetIds::HcalValidDetIds(const edm::ParameterSet& config) : edm::EDAnalyzer()
{
}

HcalValidDetIds::~HcalValidDetIds() {}

void
HcalValidDetIds::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
   using namespace edm;

   edm::ESHandle<HcalLutMetadata> lutMetadata;
   setup.get<HcalLutMetadataRcd>().get(lutMetadata);

   const auto& topo = lutMetadata->topo();

   for (int version = 0; version <= 1; ++version) {
      for (int ieta = -41; ieta <= 41; ++ieta) {
         for (int iphi = 1; iphi <= 72; ++iphi) {
            for (int depth = 0; depth <= 3; ++depth) {
               HcalTrigTowerDetId id(ieta, iphi, depth, version);
               if (topo->validHT(id))
                  std::cout << "VALID: " << id << std::endl;
            }
         }
      }
   }
}

void
HcalValidDetIds::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HcalValidDetIds);
