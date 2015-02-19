namespace psi {

class CCWavefunction: public Wavefunction {
public:
  CCWavefunction(boost::shared_ptr<Wavefunction> reference, Options &options);
  virtual ~CCWavefunction();

  

} // CCWavefunction

} // psi
