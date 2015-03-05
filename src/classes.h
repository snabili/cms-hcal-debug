#include "Debug/HcalCompareChains/interface/LinearizedTP.h"
#include "DataFormats/Common/interface/Wrapper.h"

namespace {
   struct dictionary {
      std::vector<LinearizedTP> vTP_;
      edm::Wrapper<std::vector<LinearizedTP> > anotherTP_;
   };
}
