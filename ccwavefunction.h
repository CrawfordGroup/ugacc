namespace psi {

class CCWavefunction: public Wavefunction {
public:
  CCWavefunction(boost::shared_ptr<Wavefunction> reference, Options &options);
  virtual ~CCWavefunction();

protected:
  // Energy denominators
  double **D1;
  double ****D2;

  // Ground-state T amplitudes
  double **t1;      // Current T1
  double **t1old;   // Previous iteration T1
  double ****t2;    // Current T2
  double ****t2old; // Previous iteration T2

  // Amplitudes for biorthogonal (left-hand) projectors 
  double **t1s;   // t1s(ia) = 2 t1(ia)
  double ****t2s; // t2s(ijab) = 4 t2(ijab) - 2 t2(ijba)

  // Effective doubles 
  double ****tau;  // tau(ijab) = t2(ijab) + t1(ia) * t1(jb)
  double ****ttau; // ttau(ijab) = t2(ijab) + (1/2) t1(ia) * t1(jb)
  
  // Ground-state Lambda amplitudes
  double **l1;      // Current L1
  double **l1old;   // Previous iteration L1
  double ****l2;    // Current L2
  double ****l2old; // Previous iteration L2

  // CCSD intermediates for amplitude equations (related to, but not the
  // same as corresponding HBAR quantities
  double **Fae;     
  double **Fmi;     
  double **Fme;     
  double ****Wmnij; 
  double ****Wmbej; 
  double ****Wmbje; 

  // CCSD HBAR
  double **Hoo;
  double **Hvv;
  double **Hov;
  double ****Hoooo;
  double ****Hvvvv;
  double ****Hovov;
  double ****Hovvo;
  double ****Hvovv;
  double ****Hooov;
  double ****Hovoo;
  double ****Hvvvo;

  // Three-body intermediates
  double **Gvv;  // G(ae) = -t2(ijeb) * l2(ijab)
  double **Goo;  // G(mi) = t2(mjab) * l2(ijab)

  // Inhomogeneous contributions to Lambda equations from (T) correction
  double **s1;
  double ****s2;

  // One-electron density components
  double **Doo;
  double **Dvv;
  double **Dov;
  double **Dvo;

  // Two-electron density components
  double ****Goooo;
  double ****Gvvvv;
  double ****Goovv;
  double ****Govov;
  double ****Gooov;
  double ****Gvvvo;

  // Triples (in-core algorithms only)
  double ******t3;
  double ******l3;

} // CCWavefunction

} // psi
