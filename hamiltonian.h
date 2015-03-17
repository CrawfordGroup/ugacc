#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <boost/shared_ptr.hpp>
#include <libmints/mints.h>

namespace psi {

class Hamiltonian {
public:
  Hamiltonian(boost::shared_ptr<Wavefunction> reference);
  virtual ~Hamiltonian();
//  Hamiltonian(const boost::shared_ptr<Hamiltonian> &H);

  double ** fock_p() { return fock_; }
  double **** ints_p() { return ints_; }
  double **** L_p() { return L_; }

protected:
  int nmo_;
  int no_;
  int nact_;
  int nfzc_;
  int nfzv_;

  double **fock_;
  double ****ints_;
  double ****L_;

}; // Hamiltonian

} // psi

#endif // HAMILTONIAN_H
