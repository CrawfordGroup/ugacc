class Perturbation {
public:
  std::string operator_; // perturbation name
  int nbody; // number of particles described by this operator
  Perturbation(std::string operator, boost::shared_ptr<Wavefunction> ref);

protected:
  int nmo_;
  int nso_;
  int nact_;
  int nfzc_;
  int nfzv_;
};
