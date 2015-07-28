#ifndef CCWAVEFUNCTION_H
#define CCWAVEFUNCTION_H

#include "hamiltonian.h"
#include <libmints/mints.h>
#include <boost/shared_ptr.hpp>

namespace psi { namespace ugacc {

struct onestack {
    double value;
    int i;
    int a;
};

struct twostack {
    double value;
    int i; int j;
    int a; int b;
};

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
  int dertype_;         // Gradient level

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

  // Biorthogonal projection doubles
  double **t1s_;
  double ****t2s_;
  
  // CCSD intermediates for amplitude equations (related to, but not the
  // same as corresponding HBAR quantities
  double **Fvv_;     
  double **Foo_;     
  double **Fov_;     
  double ****Woooo_; 
  double ****Wovvo_; 
  double ****Wovov_; 

  // L-amplitude quantities
  double **l1_;      /* current l1 amplitudes */
  double **l1old_;   /* previous l1 amplitudes */
  double ****l2_;    /* current l2 amplitudes */
  double ****l2old_; /* previous l2 amplitudes */

  // HBAR components
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

  // Three-body intermediates
  double **Gvv_;
  double **Goo_;

  // One-electron density components
  double **Doo_;
  double **Dvv_;
  double **Dov_;
  double **Dvo_;

  // Two-electron density components
  double ****Goooo_;
  double ****Gvvvv_;
  double ****Goovv_;
  double ****Govov_;
  double ****Gooov_;
  double ****Gvvvo_;

  // In-core triples
  double ******t3_;
  double ******l3_;

  // Extra inhomogeneous terms for Lambda from (T)-gradient
  double **s1_;
  double ****s2_;

public:
  int maxiter() { return maxiter_; }
  double convergence() { return convergence_; }
  std::string wfn() { return wfn_; } 
  bool do_diis() { return do_diis_; }
  bool ooc() { return ooc_; }
  int dertype() { return dertype_; }

  double compute_energy();
  void compute_lambda();

  double energy();
  void build_tau();
  void amp_save(std::string);
  void build_F();
  void build_W();
  void build_t1();
  void build_t2();
  double t1norm();
  void diis(int iter, std::string);
  double increment_amps(std::string);
  void build_tstar();

  void hbar();
  void init_lambda();
  void build_G();
  void build_l1();
  void build_l2();
  double pseudoenergy();

  void init_density();
  double onepdm();
  double twopdm();

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

  void amp_write(int, std::string);
  void onestack_insert(struct onestack *stack, double value, int i, int a, int level, int stacklen);
  void twostack_insert(struct twostack *stack, double value, int i, int j, int a, int b, int level, int stacklen);
  void amp_write_T1(double **T1, int no, int nv, int length, std::string label);
  void amp_write_T2(double ****T2, int no, int nv, int length, std::string label);

}; // CCWavefunction

}} // psi::ugacc

#endif // CCWAVEFUNCTION_H
