#ifndef CCDENSITY_H
#define CCDENSITY_H

#include "ccrhwavefunction.h"
#include "cclhwavefunction.h"
#include <libmints/mints.h>
#include <boost/shared_ptr.hpp>

namespace psi { namespace ugacc {

class CCDensity {
public:
  CCDensity(boost::shared_ptr<CCRHWavefunction> CC, boost::shared_ptr<CCLHWavefunction> Lambda);
  virtual ~CCDensity();

protected:
  std::string wfn_;     // wfn type (CCSD, CCSD_T, etc.)
  bool ooc_;            // Use out-of-core algorithms?

  int no_;  // Number of active occupied MOs
  int nv_;  // Number of active virtual MOs

  // One-electron density components
  double **Doo_;
  double **Dvv_;
  double **Dov_;
  double **Dvo_;

  // Two-electron density components
  double ****Goooo_;
  double ****Gvvvv_;
  double ****Goovv_;
  double ****Govov_;
  double ****Gooov_;
  double ****Gvvvo_;

public:
  double onepdm();
  double twopdm();

}; // CCDensity

}} // psi::ugacc

#endif // CCDENSITY_H
