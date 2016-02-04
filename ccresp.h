#ifndef CCRESP_H
#define CCRESP_H

#include <boost/shared_ptr.hpp>
#include "ccpert.h"

namespace psi { namespace ugacc {

class CCResp {
public:
  CCResp(boost::shared_ptr<CCPert> Y, boost::shared_ptr<CCPert> X);
  CCResp(boost::shared_ptr<CCPert> X);
  ~CCResp();
  linresp();

protected:
  boost::shared_ptr<Hamiltonian> H_;
  boost::shared_ptr<CCWfn> CC_;
  boost::shared_ptr<HBAR> HBAR_;
  boost::shared_ptr<CCLambda> Lambda_;

  int no_;
  int nv_;
  boost::shared_ptr<CCPert> Y_;
  boost::shared_ptr<<CCPert> X_;
};

}} // psi::ugaccc

#endif
