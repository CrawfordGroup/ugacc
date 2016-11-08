#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "psi4/libmints/mintshelper.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.h"

using namespace std;

namespace psi { namespace ugacc {

class Hamiltonian {
public:
  Hamiltonian(shared_ptr<PSIO>, shared_ptr<Wavefunction>, std::vector<shared_ptr<MOSpace> >);
  virtual ~Hamiltonian();

protected:
  int nmo_;
  int nso_;
  int nact_;
  int nfzc_;
  int nfzv_;
  double efzc_;

  double **fock_;
  double ****ints_;
  double ****L_;

  friend class CCWfn;
  friend class HBAR;
  friend class CCLambda;
  friend class CCDensity;
  friend class CCPert;
}; // Hamiltonian

}} // psi::ugacc

#endif // HAMILTONIAN_H
