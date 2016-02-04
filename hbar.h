#ifndef HBAR_H
#define HBAR_H

#include "hamiltonian.h"
#include "ccwfn.h"

namespace psi { namespace ugacc {

class HBAR {
public:
  HBAR(boost::shared_ptr<Hamiltonian> H, boost::shared_ptr<CCWfn> CC);
  ~HBAR();

protected:
  boost::shared_ptr<Hamiltonian> H_;
  boost::shared_ptr<CCWfn> CC_;
  int no_;
  int nv_;
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

  friend class CCWfn;
  friend class CCLambda;
  friend class CCPert;
};

}} // psi::ugacc

#endif // HBAR_H
