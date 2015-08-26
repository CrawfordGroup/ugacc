#ifndef CCPERT_H
#define CCPERT_H

#include <boost/shared_ptr.hpp>
#include <libmints/mints.h>

#include "hamiltonian.h"
#include "ccwavefunction.h"
#include "hbar.h"

namespace psi { namespace ugacc {

class CCPert {
public:
  CCPert(double **pert, double omega, boost::shared_ptr<CCWavefunction> CC, boost::shared_ptr<HBAR> HBAR);
  ~CCPert();
  void solve();

protected:
  int no_;
  int nv_;
  double **pert_;
  double omega_;

  boost::shared_ptr<CCWavefunction> CC_;
  boost::shared_ptr<Hamiltonian> H_;
  boost::shared_ptr<HBAR> HBAR_;

  // Energy denominators
  double **D1_;
  double ****D2_;

  // Similarity transformed perturbation operator components
  double **Xov_;
  double **Xoo_;
  double **Xvv_;
  double **Xvo_;
  double ****Xovoo_;
  double ****Xvvvo_;
  double ****Xvvoo_;

  // Perturbed wave function
  double **X1_;
  double ****X2_;
  double **X1old_;
  double ****X2old_;

  // DIIS-related vectors
  std::vector<double> X1diis_;
  std::vector<double> X2diis_;
  std::vector<double> X1err_;
  std::vector<double> X2err_;

  void xbar();
  void amp_save();
  double increment_amps();
  void build_X1();
  void build_X2();
  double pseudoresponse();
  void build_diis_error();
  void save_diis_vectors();
  void print_amps();

  friend class CCResp;
};

}} // psi::ugacc

#endif
