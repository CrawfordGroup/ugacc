#ifndef CCLHWAVEFUNCTION_H
#define CCLHWAVEFUNCTION_H

#include "hamiltonian.h"
#include <libmints/mints.h>
#include <boost/shared_ptr.hpp>

namespace psi { namespace ugacc {

class CCLHWavefunction {
public:
  CCLHWavefunction(boost::shared_ptr<CCRHWavefunction> CC, boost::shared_ptr<HBAR> HBAR)
  virtual ~CCLHWavefunction();

protected:
  int no_;  // Number of active occupied MOs
  int nv_;  // Number of active virtual MOs

  boost::shared_ptr<HBAR> HBAR_;

  // Energy denominators
  double **D1_;
  double ****D2_;

  // L-amplitude quantities
  double **l1_;      /* current l1 amplitudes */
  double **l1old_;   /* previous l1 amplitudes */
  double ****l2_;    /* current l2 amplitudes */
  double ****l2old_; /* previous l2 amplitudes */

  // Three-body intermediates
  double **Gvv_;
  double **Goo_;

  // In-core triples
  double ******l3_;

  // Extra inhomogeneous terms for Lambda from (T)-gradient
  double **s1_;
  double ****s2_;

public:
  void compute_lambda();

  void amp_save()
  double increment_amps()
  void init_lambda();
  void build_G();
  void build_l1();
  void build_l2();
  double pseudoenergy();

  friend class CCDensity;

}; // CCLHWavefunction

}} // psi::ugacc

#endif // CCLHWAVEFUNCTION_H
