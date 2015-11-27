#include "Debug/HcalDebug/interface/LinearizedTP.h"

#if UP_THERE
LinearizedTP::LinearizedTP(const HcalUpgradeTriggerPrimitiveDigi& d) :
   ieta(d.id().ieta()),
   iphi(d.id().iphi())
{
}
#endif
