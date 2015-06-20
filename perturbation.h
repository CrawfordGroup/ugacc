#ifndef PERTURBATION_H
#define PERTURBATION_H

#include <libmints/wavefunction.h>

// Allowed operators:
//
// Mu = electric dipole = -r
// P  = velocity-gauge electric dipole = -i Nabla
// P* = complex conjugate of P
// L  = magnetic dipole = -(1/2) r x P
// L* = complex conjugate of L
// Q  = traceless quadrupole
// 

namespace psi {

class Perturbation {
public:
  std::string operator_; // perturbation name
  Perturbation(std::string op, boost::shared_ptr<Wavefunction> ref);
  ~Perturbation();
  double **prop_p(int i) { return prop2_[i]; }
  double **prop_p(int i, int j) { return prop3_[i][j]; }
  void print(std::string out);
  void print(int i, std::string out);
  void print(int i, int j, std::string out);
  void print();
  void print(int i);
  void print(int i, int j);

protected:
  int nmo_;
  int nso_;
  int nact_;
  int nfzc_;
  int nfzv_;

  // Yes, a generalized structure would be nicer
  double ***prop2_; // dipolar quantities use this pointer
  double ****prop3_; // quadrupole quantities use this one

private:
  bool allowed(std::string op);
  bool onebody(std::string op);
  bool twobody(std::string op);
};

} // psi

#endif // PERTURBATION_H
