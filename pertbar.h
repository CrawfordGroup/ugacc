#ifndef PERTBAR_H
#define PERTBAR_H

namespace psi {

class Pertbar {
public:
  Pertbar(boost::shared_ptr<Perturbation> pert, boost::shared_ptr<CCWavefunction> cc);
  virtual ~Pertbar();

protected:
  int no_;
  int nv_;
  int nact_;

  boost::shared_ptr<Perturbation> pert_;
  boost::shared_ptr<CCWavefunction> cc_;

  double ***oo_;
  double ***vv_;
  double ***ov_;
  double ***vo_;

  double *****vvvo_;
  double *****ooov_;
  double *****vvoo_;
};

} // psi

#endif // PERTBAR_H
