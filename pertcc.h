#ifndef PERTCC_H
#define PERTCC_H

#include <boost/shared_ptr.hpp>
#include "ccwavefunction.h"

namespace psi { namespace ugacc {

class PertCC {

public:
  PertCC(boost::shared_ptr<CCWavefunction> cc, boost::shared_ptr<Perturbation> pert, double omega);
  ~PertCC();
  void solve();

protected:
  std::string wfn_;     // wfn type (CCSD, CCSD_T, etc.)
  double convergence_;  // conv. on RMS residual change between iterations
  int maxiter_;         // maximum number of iterations
  bool do_diis_;        // use DIIS algorithms?
  bool ooc_;            // Use out-of-core algorithms?

  boost::shared_ptr<Hamiltonian> H_;
  int nact_;
  int no_;              // Number of active occupied MOs
  int nv_;              // Number of active virtual MOs
  double **t1_;
  double **t2_;
  double **l1_;
  double ****l2_;
  double **Hoo_;
  double **Hvv_;
  double **Hov_;
  double ****Hoooo_;
  double ****Hvvvv_;
  double ****Hovov_;
  double ****Hovvo_;
  double ****Hvovv_;
  double ****Hooov_;
  double ****Hovoo_;
  double ****Hvvvo_;
};

}} // psi::ugacc

#endif // PERTCC_H
