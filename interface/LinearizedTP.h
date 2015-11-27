#define UP_THERE 0
#ifndef __Debug_LinearizedTP_h
#define __Debug_LinearizedTP_h

#include <vector>

#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#if UP_THERE
#include "DataFormats/HcalDigi/interface/HcalUpgradeTriggerPrimitiveDigi.h"
#endif

class LinearizedTP {
   public:
      LinearizedTP() {};
#if UP_THERE
      LinearizedTP(const HcalUpgradeTriggerPrimitiveDigi& d);
#endif
      virtual ~LinearizedTP() {};

      int ieta;
      int iphi;

      double soi_energy;
      std::vector<double> summed_energies;

      std::vector<double> rising_times;
      std::vector<double> falling_times;
};

#endif
