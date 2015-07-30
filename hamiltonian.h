#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <boost/shared_ptr.hpp>
#include <libmints/mints.h>
#include <libtrans/integraltransform.h>

// friends
class CCRHWavefunction;
class HBAR;

namespace psi { namespace ugacc {

class Hamiltonian {
public:
  Hamiltonian(boost::shared_ptr<PSIO>, boost::shared_ptr<Wavefunction>, std::vector<boost::shared_ptr<MOSpace> >);
  virtual ~Hamiltonian();

protected:
  int nmo_;
  int nso_;
  int nact_;
  int nfzc_;
  int nfzv_;

  double **fock_;
  double ****ints_;
  double ****L_;

  friend class CCRHWavefunction;
  friend class HBAR;

}; // Hamiltonian

}} // psi::ugacc

#endif // HAMILTONIAN_H
