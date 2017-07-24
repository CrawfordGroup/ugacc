#ifndef CCPERT_H
#define CCPERT_H

#include "psi4/libmints/mintshelper.h"

#include "hamiltonian.h"
#include "ccwfn.h"
#include "hbar.h"
#include "cclambda.h"

using namespace std;

namespace psi { namespace ugacc {

enum hand {left, right};

class CCPert {
public:
  CCPert(double **pert, double omega, shared_ptr<CCWfn> CC, shared_ptr<HBAR> HBAR, shared_ptr<CCLambda> CCLambda);
  ~CCPert();
  void solve(enum hand);

protected:
  int no_;
  int nv_;
  double **pert_;
  double omega_;

  shared_ptr<CCWfn> CC_;
  shared_ptr<Hamiltonian> H_;
  shared_ptr<HBAR> HBAR_;
  shared_ptr<CCLambda> CCLambda_;

  // Energy denominators
  double **D1_;
  double ****D2_;
  double **l1_;
  double ****l2_;

  // Similarity transformed perturbation operator components
  double **Aov_;
  double **Aoo_;
  double **Avv_;
  double **Avo_;
  double ****Aovoo_;
  double ****Avvvo_;
  double ****Avvoo_;

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

  // Three-body intermediates
  double **Gvv_;
  double **Goo_;

  void pertbar();
  void amp_save(enum hand);
  double increment_amps(enum hand);
  void build_G();
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
