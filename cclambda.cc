#include "cclambda.h"

#include <boost/shared_ptr.hpp>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <cmath>
#include <libdiis/diismanager.h>

#include "array.h"

namespace psi { namespace ugacc {

CCLambda::CCLambda(boost::shared_ptr<CCWavefunction> CC, boost::shared_ptr<HBAR> HBAR)
{
  HBAR_ = HBAR;
  CC_ = CC;

  H_ = CC_->H_;

  D1_ = CC_->D1_;
  D2_ = CC_->D2_;

  no_ = CC_->no_;
  nv_ = CC_->nv_;

  int no = no_;
  int nv = nv_;

  l1_ = block_matrix(no,nv);
  l1old_ = block_matrix(no,nv);
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      l1_[i][a] = 2.0 * CC_->t1_[i][a];

  l2_ = init_4d_array(no,no,nv,nv);
  l2old_ = init_4d_array(no,no,nv,nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          l2_[i][j][a][b] = 2.0 * (2.0 * CC_->t2_[i][j][a][b] - CC_->t2_[i][j][b][a]);

  Gvv_ = block_matrix(nv, nv);
  Goo_ = block_matrix(no, no);

  // DIIS Vectors
  l1diis_.resize(no_*nv_);
  l2diis_.resize(no_*no_*nv_*nv_);
  l1err_.resize(no_*nv_);
  l2err_.resize(no_*no_*nv_*nv_);
}

CCLambda::~CCLambda()
{
  int no = no_;
  int nv = nv_;

  free_block(l1_);
  free_block(l1old_);
  free_4d_array(l2_, no, no, nv);
  free_4d_array(l2old_, no, no, nv);

  free_block(Gvv_);
  free_block(Goo_);
}

void CCLambda::amp_save()
{
  double ****l2tmp = l2_;
  l2_ = l2old_;
  l2old_ = l2tmp;

  double **l1tmp = l1_;
  l1_ = l1old_;
  l1old_ = l1tmp;
}

double CCLambda::increment_amps()
{
  int no = no_;
  int nv = nv_;

  double residual1 = 0.0;
  double residual2 = 0.0;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      residual1 += l1_[i][a] * l1_[i][a];     
      l1_[i][a] = l1old_[i][a] + l1_[i][a]/D1_[i][a];
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++) {
          residual2 += l2_[i][j][a][b] * l2_[i][j][a][b];
          l2_[i][j][a][b] = l2old_[i][j][a][b] + l2_[i][j][a][b]/D2_[i][j][a][b];
        }
    }

  return sqrt(residual1 + residual2);
}

void CCLambda::compute_lambda()
{
  outfile->Printf("\n\tThe Coupled-Cluster Lambda Iteration:\n");
  outfile->Printf(  "\t-------------------------------------\n");
  outfile->Printf(  "\t Iter   PseudoEnergy        RMS   \n");
  outfile->Printf(  "\t-------------------------------------\n");
  outfile->Printf(  "\t  %3d  %20.15f\n", 0, pseudoenergy());

  double rms = 0.0;
  boost::shared_ptr<DIISManager> diis(new DIISManager(8, "CCLambda DIIS",
    DIISManager::LargestError, DIISManager::InCore));
  diis->set_error_vector_size(2, DIISEntry::Pointer, no_*nv_, DIISEntry::Pointer, no_*no_*nv_*nv_);
  diis->set_vector_size(2, DIISEntry::Pointer, no_*nv_, DIISEntry::Pointer, no_*no_*nv_*nv_);
  for(int iter=1; iter <= CC_->maxiter_; iter++) {
    amp_save();    
    build_G();
    build_l1();
    build_l2();
    rms = increment_amps();
    if(rms < CC_->convergence_) break;
    if(CC_->do_diis_) {
      build_diis_error();
      diis->add_entry(4, l1err_.data(), l2err_.data(), l1diis_.data(), l2diis_.data());
      if(diis->subspace_size() > 2) diis->extrapolate(2, l1diis_.data(), l2diis_.data());
      save_diis_vectors();
    }
    outfile->Printf(  "\t  %3d  %20.15f  %5.3e\n",iter, pseudoenergy(), rms);
  }
  if(rms >= CC_->convergence_)
    throw PSIEXCEPTION("Computation has not converged.");
}

void CCLambda::build_G()
{
  int no = no_;
  int nv = nv_;
  double ****t2 = CC_->t2_;
  double ****l2 = l2old_;

  for(int m=0; m < no; m++)
    for(int i=0; i < no; i++) {
      double value = 0.0;
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          for(int j=0; j < no; j++)
            value += t2[m][j][a][b] * l2[i][j][a][b];
      Goo_[m][i] = value;
    }

  for(int a=0; a < nv; a++)
    for(int e=0; e < nv; e++) {
      double value = 0.0;
      for(int i=0; i < no; i++)
        for(int j=0; j < no; j++)
          for(int b=0; b < nv; b++)
            value -= t2[i][j][e][b] * l2[i][j][a][b];
      Gvv_[a][e] = value;
    }
}

void CCLambda::build_l1()
{
  int no = no_;
  int nv = nv_;
  double **l1 = l1old_;
  double ****l2 = l2old_;
  double **l1new = l1_;
  double **Gvv = Gvv_;
  double **Goo = Goo_;
  double **Hov = HBAR_->Hov_;
  double **Hvv = HBAR_->Hvv_;
  double **Hoo = HBAR_->Hoo_;
  double ****Hvvvo = HBAR_->Hvvvo_;
  double ****Hovoo = HBAR_->Hovoo_;
  double ****Hovvo = HBAR_->Hovvo_;
  double ****Hovov = HBAR_->Hovov_;
  double ****Hvovv = HBAR_->Hvovv_;
  double ****Hooov = HBAR_->Hooov_;

  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      double value = 2 * Hov[i][a];
      if(CC_->wfn_ == "CCSD_T" && CC_->dertype_ == "FIRST") value += CC_->s1_[i][a];

      for(int e=0; e < nv; e++)
        value += l1[i][e] * Hvv[e][a];

      for(int m=0; m < no; m++)
        value -= l1[m][a] * Hoo[i][m];

      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          value += l1[m][e] * (2*Hovvo[i][e][a][m] - Hovov[i][e][m][a]);

      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          for(int f=0; f < nv; f++)
            value += l2[i][m][e][f] * Hvvvo[e][f][a][m];

      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          for(int e=0; e < nv; e++)
            value -= l2[m][n][a][e] * Hovoo[i][e][m][n];

      for(int e=0; e < nv; e++)
        for(int f=0; f < nv; f++)
          value -= Gvv[e][f] * (2*Hvovv[e][i][f][a] - Hvovv[e][i][a][f]);

      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          value -= Goo[m][n] * (2*Hooov[m][i][n][a] - Hooov[i][m][n][a]);

      l1new[i][a] = value;
    }
}

void CCLambda::build_l2()
{
  int no = no_;
  int nv = nv_;
  double **l1 = l1old_;
  double ****l2 = l2old_;
  double ****l2new = l2_;
  double ****L = H_->L_;
  double **Gvv = Gvv_;
  double **Goo = Goo_;
  double **Hov = HBAR_->Hov_;
  double **Hvv = HBAR_->Hvv_;
  double **Hoo = HBAR_->Hoo_;

  double ****Hovvo = HBAR_->Hovvo_;
  double ****Hovov = HBAR_->Hovov_;
  double ****Hoooo = HBAR_->Hoooo_;
  double ****Hvvvv = HBAR_->Hvvvv_;
  double ****Hvovv = HBAR_->Hvovv_;
  double ****Hooov = HBAR_->Hooov_;

  double ****Z = init_4d_array(no, no, nv, nv);

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          double value = L[i][j][a+no][b+no];
          if(CC_->wfn_ == "CCSD_T" && CC_->dertype_ == "FIRST") value += 0.5 * CC_->s2_[i][j][a][b];

          value += 2.0*l1[i][a]*Hov[j][b] - l1[j][a]*Hov[i][b];

          for(int e=0; e < nv; e++)
            value += l2[i][j][e][b]*Hvv[e][a];

          for(int m=0; m < no; m++)
            value -= l2[m][j][a][b]*Hoo[i][m];

          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              value += 0.5 * l2[m][n][a][b] * Hoooo[i][j][m][n];

          for(int e=0; e < nv; e++)
            for(int f=0; f < nv; f++)
              value += 0.5 * l2[i][j][e][f] * Hvvvv[e][f][a][b];

          for(int e=0; e < nv; e++)
              value += l1[i][e]*(2*Hvovv[e][j][a][b] - Hvovv[e][j][b][a]);

          for(int m=0; m < no; m++)
              value -= l1[m][b]*(2*Hooov[j][i][m][a] - Hooov[i][j][m][a]);

          for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++) {
              value += (2*Hovvo[i][e][a][m] - Hovov[i][e][m][a])*l2[m][j][e][b];
              value -= Hovov[j][e][m][a]*l2[m][i][b][e];
              value -= Hovvo[j][e][a][m]*l2[m][i][e][b];
            }

          for(int e=0; e < nv; e++)
            value += Gvv[a][e]*L[i][j][e+no][b+no];
          for(int m=0; m < no; m++)
            value -= Goo[m][i]*L[m][j][a+no][b+no];

          Z[i][j][a][b] = value;
    }

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          l2new[i][j][a][b] = Z[i][j][a][b] + Z[j][i][b][a];

   free_4d_array(Z, no, no, nv);
}

/*
** pseudoenergy(): Evaluates an energy-like expression for the Lambda doubles amplitudes: 
**   E = <0|L2 H|0> = 1/2 <ab|ij> L(ij,ab)
** This expression is derived in the UGA formalism.
*/
double CCLambda::pseudoenergy()
{
  int no = no_;
  int nv = nv_;
  double ****ints = H_->ints_;
  double ****l2 = l2_;

  double energy=0.0;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          energy += 0.5*ints[i][j][a+no][b+no]*l2[i][j][a][b];

  return energy;
}

void CCLambda::build_diis_error()
{
  int no = no_;
  int nv = nv_;

  int t1len = 0;
  int t2len = 0;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      l1diis_[t1len] = l1_[i][a];
      l1err_[t1len++] = l1_[i][a] - l1old_[i][a];
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++) {
          l2diis_[t2len] = l2_[i][j][a][b];
          l2err_[t2len++] = l2_[i][j][a][b] - l2old_[i][j][a][b];
        }
    }
}

void CCLambda::save_diis_vectors()
{
  int no = no_;
  int nv = nv_;

  int t1len = 0;
  int t2len = 0;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      l1_[i][a] = l1diis_[t1len++];
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++) {
          l2_[i][j][a][b] = l2diis_[t2len++];
        }
    }
}

}} // psi::ugacc
