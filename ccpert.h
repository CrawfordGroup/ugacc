#ifndef CCPERT_H
#define CCPERT_H

#include <boost/shared_ptr.hpp>
#include <libmints/mints.h>

#include "hamiltonian.h"
#include "ccwfn.h"
#include "hbar.h"

namespace psi { namespace ugacc {

enum hand {left, right};

class CCPert {
public:
  CCPert(double **pert, double omega, boost::shared_ptr<CCWfn> CC, boost::shared_ptr<HBAR> HBAR);
  ~CCPert();
  void solve(enum hand);

protected:
  int no_;
  int nv_;
  double **pert_;
  double omega_;

  boost::shared_ptr<CCWfn> CC_;
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

  // Right-hand perturbed wave function
  double **X1_;
  double ****X2_;
  double **X1old_;
  double ****X2old_;

  // Left-hand perturbed wave function
  double **Y1_;
  double ****Y2_;
  double **Y1old_;
  double ****Y2old_;

  // DIIS-related vectors
  std::vector<double> X1diis_;
  std::vector<double> X2diis_;
  std::vector<double> X1err_;
  std::vector<double> X2err_;
  std::vector<double> Y1diis_;
  std::vector<double> Y2diis_;
  std::vector<double> Y1err_;
  std::vector<double> Y2err_;

  void pertbar();
  void amp_save(enum hand);
  double increment_amps(enum hand);
  void build_X1();
  void build_X2();
  void build_Y1();
  void build_Y2();
  double pseudoresponse(enum hand);
  void build_diis_error(enum hand);
  void save_diis_vectors(enum hand);
  void print_amps(enum hand);

  friend class CCResp;
};

}} // psi::ugacc

#endif
