#ifndef CCWAVEFUNCTION_H
#define CCWAVEFUNCTION_H

#include "hamiltonian.h"
#include <libmints/mints.h>
#include <boost/shared_ptr.hpp>

namespace psi {

class CCWavefunction: public Wavefunction {
public:
  CCWavefunction(boost::shared_ptr<Wavefunction> reference,
                 boost::shared_ptr<Hamiltonian> H,
                 Options &options, boost::shared_ptr<PSIO> psio);
  virtual ~CCWavefunction();

protected:
  std::string wfn_;     // wfn type (CCSD, CCSD_T, etc.)
  double convergence_;  // conv. on RMS residual change between iterations
  int maxiter_;         // maximum number of iterations
  bool do_diis_;        // use DIIS algorithms?
  bool ooc_;            // Use out-of-core algorithms?

  int no_;  // Number of active occupied MOs
  int nv_;  // Number of active virtual MOs

  boost::shared_ptr<Hamiltonian> H_; // integrals and Fock matrix
  
  // Energy denominators
  double **D1_;
  double ****D2_;

  // Ground-state T amplitudes
  double **t1_;      // Current T1
  double **t1old_;   // Previous iteration T1
  double ****t2_;    // Current T2
  double ****t2old_; // Previous iteration T2

  // Effective doubles 
  double ****tau_;  // tau(ijab) = t2(ijab) + t1(ia) * t1(jb)
  double ****ttau_; // ttau(ijab) = t2(ijab) + (1/2) t1(ia) * t1(jb)
  
  // CCSD intermediates for amplitude equations (related to, but not the
  // same as corresponding HBAR quantities
  double **Fvv_;     
  double **Foo_;     
  double **Fov_;     
  double ****Woooo_; 
  double ****Wovvo_; 
  double ****Wovov_; 

public:
  int maxiter() { return maxiter_; }
  double convergence() { return convergence_; }
  std::string wfn() { return wfn_; } 
  bool do_diis() { return do_diis_; }
  bool ooc() { return ooc_; }

  double compute_energy();

  double energy();
  void build_tau();
  void amp_save();
  void build_F();
  void build_W();
  void build_t1();
  void build_t2();
  double t1norm();
  void diis(int iter);
  double increment_amps();

}; // CCWavefunction

} // psi

#endif // CCWAVEFUNCTION_H
