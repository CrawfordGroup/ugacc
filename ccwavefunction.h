#include <libmints/mints.h>

namespace psi {

class CCWavefunction: public Wavefunction {
public:
  CCWavefunction(boost::shared_ptr<Wavefunction> reference,
                 Options &options, boost::shared_ptr<PSIO> psio);
  virtual ~CCWavefunction();

protected:
  std::string wfn_;     // wfn type (CCSD, CCSD_T, etc.)
  double convergence_;  // conv. on RMS residual change between iterations
  unsigned int maxiter_;// maximum number of iterations
  bool do_diis_;        // use DIIS algorithms?
  bool ooc_;            // Use out-of-core algorithms?
  
  // Energy denominators
  double **D1_;
  double ****D2_;

  // Ground-state T amplitudes
  double **t1_;      // Current T1
  double **t1old_;   // Previous iteration T1
  double ****t2_;    // Current T2
  double ****t2old_; // Previous iteration T2

  // Amplitudes for biorthogonal (left-hand) projectors 
  double **t1s_;   // t1s(ia) = 2 t1(ia)
  double ****t2s_; // t2s(ijab) = 4 t2(ijab) - 2 t2(ijba)

  // Effective doubles 
  double ****tau_;  // tau(ijab) = t2(ijab) + t1(ia) * t1(jb)
  double ****ttau_; // ttau(ijab) = t2(ijab) + (1/2) t1(ia) * t1(jb)
  
  // Ground-state Lambda amplitudes
  double **l1_;      // Current L1
  double **l1old_;   // Previous iteration L1
  double ****l2_;    // Current L2
  double ****l2old_; // Previous iteration L2

  // CCSD intermediates for amplitude equations (related to, but not the
  // same as corresponding HBAR quantities
  double **Fae_;     
  double **Fmi_;     
  double **Fme_;     
  double ****Wmnij_; 
  double ****Wmbej_; 
  double ****Wmbje_; 

  // CCSD HBAR
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
  double **Gvv_;  // G(ae) = -t2(ijeb) * l2(ijab)
  double **Goo_;  // G(mi) = t2(mjab) * l2(ijab)

  // Inhomogeneous contributions to Lambda equations from (T) correction
  double **s1_;
  double ****s2_;

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

  // Triples (in-core algorithms only)
  double ******t3_;
  double ******l3_;

public:
  unsigned int maxiter() { return maxiter_; }
  double convergence() { return convergence_; }
  std::string wfn() { return wfn_; } 
  bool do_diis() { return do_diis_; }
  bool ooc() { return ooc_; }

  double compute_energy();

}; // CCWavefunction

} // psi
