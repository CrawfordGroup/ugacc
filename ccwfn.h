#ifndef CCWFN_H
#define CCWFN_H

#include "hamiltonian.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/molecule.h"

using namespace std;

namespace psi { namespace ugacc {

class CCWfn: public Wavefunction {
public:
  CCWfn(shared_ptr<Wavefunction> reference, shared_ptr<Hamiltonian> H, Options &options);
  virtual ~CCWfn();

protected:
  std::string wfn_;     // wfn type (CCSD, CCSD_T, etc.)
  double convergence_;  // conv. on RMS residual change between iterations
  int maxiter_;         // maximum number of iterations
  bool do_diis_;        // use DIIS algorithms?
  bool ooc_;            // Use out-of-core algorithms?
  std::string dertype_; // Gradient level -- needed only for (T) gradients

  int no_;  // Number of active occupied MOs
  int nv_;  // Number of active virtual MOs

  shared_ptr<Hamiltonian> H_; // integrals and Fock matrix

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

  // Biorthogonal projection doubles
  double **t1s_;
  double ****t2s_;
  
  // CCSD intermediates for amplitude equations (related to, but not the
  // same as corresponding HBAR quantities)
  double **Fvv_;     
  double **Foo_;     
  double **Fov_;     
  double ****Woooo_; 
  double ****Wovvo_; 
  double ****Wovov_; 
  double ****Wovoo_; // Zmbij part of Wabef

  // Extra contributions for (T) gradients
  double ******t3_; // only for in-core code
  double ******l3_; // only for in-core code
  double **s1_;
  double ****s2_;
  double **Doo_;
  double **Dvv_;
  double **Dov_;
  double ****Gooov_;
  double ****Gvvvo_;
  double ****Goovv_;

  double energy();
  void build_tau();
  void amp_save();
  void build_F();
  void build_W();
  void build_t1();
  void build_t2();
  void build_diis_error(std::shared_ptr<Vector> t1err_, std::shared_ptr<Vector> t2_err_, std::shared_ptr<Vector> t1diis_, std::shared_ptr<Vector> t2diis_);
  void save_diis_vectors(std::shared_ptr<Vector> t1diis_, std::shared_ptr<Vector> t2diis_);
  double t1norm();
  double increment_amps();
  void build_tstar();
  void print_amps();

  double tcorr();
  double tcorr_ooc();
  double tcorr_ooc_TJL();
  void tgrad();
  void tgrad_ooc();

  void t3_ijk(double ***t3, int i, int j, int k, double ****t2, double **fock, double ****ints);
  void W3_ijk(double ***W3, int i, int j, int k, double ****t2, double ****ints);
  void t3_abc(double ***t3, int a, int b, int c, double ****t2, double **fock, double ****ints);
  void M3_ijk(double ***M3, int i, int j, int k, double ****t2, double **fock, double ****ints);
  void M3_abc(double ***M3, int a, int b, int c, double ****t2, double **fock, double ****ints);
  void N3_ijk(double ***N3, int i, int j, int k, double ****t2, double **t1, double **fock, double ****ints);
  void N3_abc(double ***N3, int a, int b, int c, double ****t2, double **t1, double **fock, double ****ints);

public:
  double compute_energy();

  friend class HBAR;
  friend class CCLambda;
  friend class CCDensity;
  friend class CCPert;
  friend class CCResp;
}; // CCWfn

}} // psi::ugacc

#endif // CCWFN_H
