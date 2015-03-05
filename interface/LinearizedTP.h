#ifndef __Debug_LinearizedTP_h
#define __Debug_LinearizedTP_h

#include <vector>

#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "DataFormats/HcalDigi/interface/HcalUpgradeTriggerPrimitiveDigi.h"

/* #include "TObject.h" */

class LinearizedTP {
   public:
      LinearizedTP() {};
      LinearizedTP(const HcalUpgradeTriggerPrimitiveDigi& d);
      virtual ~LinearizedTP() {};

      int ieta;
      int iphi;

      double soi_energy;
      std::vector<double> summed_energies;

      std::vector<double> rising_times;
      std::vector<double> falling_times;

      /* ClassDef(LinearizedTP, 1); */
};

#endif
