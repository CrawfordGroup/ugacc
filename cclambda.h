#ifndef CCLAMBDA_H
#define CCLAMBDA_H

#include "hamiltonian.h"
#include "ccwfn.h"
#include "hbar.h"
#include "psi4/libmints/mintshelper.h"

using namespace std;

namespace psi { namespace ugacc {

class CCLambda {
public:
  CCLambda(shared_ptr<CCWfn>, shared_ptr<HBAR>);
  virtual ~CCLambda();

protected:
  int no_;  // Number of active occupied MOs
  int nv_;  // Number of active virtual MOs

  shared_ptr<Hamiltonian> H_;
  shared_ptr<CCWfn> CC_;
  shared_ptr<HBAR> HBAR_;

  // Energy denominators
  double **D1_;
  double ****D2_;

  // L-amplitude quantities
  double **l1_;      /* current l1 amplitudes */
  double **l1old_;   /* previous l1 amplitudes */
  double ****l2_;    /* current l2 amplitudes */
  double ****l2old_; /* previous l2 amplitudes */

  // DIIS-related vectors
  std::vector<double> l1diis_;
  std::vector<double> l2diis_;
  std::vector<double> l1err_;
  std::vector<double> l2err_;

  // Three-body intermediates
  double **Gvv_;
  double **Goo_;

  // In-core triples
  double ******l3_;

  // Extra inhomogeneous terms for Lambda from (T)-gradient
  double **s1_;
  double ****s2_;

  void amp_save();
  double increment_amps();
  void build_G();
  void build_l1();
  void build_l2();
  double pseudoenergy();
  void build_diis_error();
  void save_diis_vectors();

public:
  void compute_lambda();

  friend class CCDensity;
  friend class CCPert;

}; // CCLambda

}} // psi::ugacc

#endif // CCLAMBDA_H
