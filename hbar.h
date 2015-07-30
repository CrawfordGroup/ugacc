#include "hamiltonian.h"
#include "ccrhwavefunction.h"
#include "cclhwavefunction.h"

namespace psi { namespace ugacc {

class HBAR {

public:
  HBAR(boost::shared_ptr<Hamiltonian> H, boost::shared_ptr<CCWavefunction> CC);
  ~HBAR();

protected:
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

  friend class CCRHWavefunction;
  friend class CCLHWavefunction;
};

}} // psi::ugacc
