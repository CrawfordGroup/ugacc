#ifndef PERTCC_H
#define PERTCC_H

#include <boost/shared_ptr.hpp>
#include "ccwavefunction.h"

namespace psi { namespace ugacc {

class PertCC {
  public:
    PertCC(boost::shared_ptr<CCWavefunction> cc, std::string op, double **pert); // do I need anything else?
    ~PertCC();

  protected:

};

}} // psi::ugacc

#endif // PERTCC_H
