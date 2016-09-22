#if 0
// -*- C++ -*-
//
// Package:    LinearizedTPProducer
// Class:      LinearizedTPProducer
// 
/**\class LinearizedTPProducer LinearizedTPProducer.cc Debug/LinearizedTPProducer/plugins/LinearizedTPProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Wolf
//         Created:  Thu, 12 Feb 2015 12:12:42 GMT
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "Debug/HcalDebug/interface/LinearizedTP.h"

//
// class declaration
//

class LinearizedTPProducer : public edm::EDProducer {
   public:
      explicit LinearizedTPProducer(const edm::ParameterSet&);
      ~LinearizedTPProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::InputTag digis_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
LinearizedTPProducer::LinearizedTPProducer(const edm::ParameterSet& config) :
   digis_(config.getParameter<edm::InputTag>("tps"))
{
   produces<std::vector<LinearizedTP>>();
}


LinearizedTPProducer::~LinearizedTPProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
LinearizedTPProducer::produce(edm::Event& event, const edm::EventSetup& setup)
{
   using namespace edm;

   std::auto_ptr<std::vector<LinearizedTP>> res(new std::vector<LinearizedTP>());

#if UP_THERE
   Handle<HcalUpgradeTrigPrimDigiCollection> digis;
   if (!event.getByLabel(digis_, digis)) {
      LogError("LinearizedTPProducer") <<
         "Can't find hcal trigger primitive digi collection with tag '" <<
         digis_ << "'" << std::endl;
      return;
   }

   ESHandle<CaloTPGTranscoder> decoder;
   setup.get<CaloTPGRecord>().get(decoder);
   decoder->setup(setup, CaloTPGTranscoder::HcalTPG);

   for (const auto& digi: *digis) {
      LinearizedTP lin(digi);
      auto soi = decoder->hcaletValue(
            lin.ieta,
            lin.iphi,
            digi.SOI_compressedEt());

      lin.soi_energy = soi;

      int total = 0;
      for (const auto& level: digi.getDepthData())
         total += level;
      for (const auto& level: digi.getDepthData())
         lin.summed_energies.push_back(total > 0 ? soi * level / float(total) : 0);

      lin.rising_times = digi.getRisingAvg();
      lin.falling_times = digi.getFallingAvg();

      res->push_back(lin);
   }
#endif

   event.put(res);
}

// ------------ method called once each job just before starting event loop  ------------
void 
LinearizedTPProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
LinearizedTPProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
LinearizedTPProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
LinearizedTPProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
LinearizedTPProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
LinearizedTPProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LinearizedTPProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LinearizedTPProducer);
#endif
