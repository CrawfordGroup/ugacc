#include "pertcc.h"

namespace psi { namespace ugacc {

PertCC::PertCC(boost::shared_ptr<CCWavefunction> cc, boost::shared_ptr<Perturbation> pert, double omega)
{
  // We're a friend of cc, so we have access to all its members, so we got that goin' for us...which is nice.
  wfn_ = cc->wfn_;
  maxiter_ = cc->maxiter_;
  convergence_ = cc->convergence_;
  do_diis_ = cc->do_diis_;
  ooc_ = cc->ooc_;
  H_ = cc->H_;
  no_ = cc->no_;
  nv_ = cc->nv_;
  t1_ = cc->t1_;
  t2_ = cc->t2_;
}

PertCC::~PertCC()
{

}

}} // psi::ugacc
