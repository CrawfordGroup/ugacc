#ifndef CCDENSITY_H
#define CCDENSITY_H

#include "ccwavefunction.h"
#include "cclambda.h"
#include <libmints/mints.h>
#include <boost/shared_ptr.hpp>

namespace psi { namespace ugacc {

class CCDensity {
public:
  CCDensity(boost::shared_ptr<CCWavefunction> CC, boost::shared_ptr<CCLambda> CCLambda);
  virtual ~CCDensity();

protected:
  int no_;
  int nv_;

  boost::shared_ptr<Hamiltonian> H_;
  boost::shared_ptr<HBAR> HBAR_;
  boost::shared_ptr<CCWavefunction> CC_;
  boost::shared_ptr<CCLambda> CCLambda_;

  double **t1_;
  double ****t2_;
  double ****tau_;
  double **l1_;
  double ****l2_;

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
