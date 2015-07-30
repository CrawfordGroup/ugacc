#include <boost/shared_ptr.hpp>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <cmath>
#include <libpsio/psio.h>

#include "hbar.h"
#include "ccrhwavefunction.h"
#include "cclhwavefunction.h"

#include "globals.h"

namespace psi { namespace ugacc {

CCLHWavefunction::CCLHWavefunction(boost::shared_ptr<CCRHWavefunction> CC, boost::shared_ptr<HBAR> HBAR)
{
  int no = no_;
  int nv = nv_;

  D1_ = CC->D1_;
  D2_ = CC->D2_;

  l1_ = block_matrix(no,nv);
  l1old_ = block_matrix(no,nv);
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      l1_[i][a] = 2.0 * CC->t1_[i][a];

  l2_ = init_4d_array(no,no,nv,nv);
  l2old_ = init_4d_array(no,no,nv,nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          l2_[i][j][a][b] = 2.0 * (2.0 * CC->t2_[i][j][a][b] - t2[i][j][b][a]);

  Gvv_ = block_matrix(nv, nv);
  Goo_ = block_matrix(no, no);
}

CCLHWavefunction::~CCLHWavefunction()
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

void CCLHWavefunction::amp_save()
{
  double ****t2tmp = l2_;
  l2_ = l2old_;
  l2old_ = t2tmp;

  double **t1tmp = l1_;
  l1_ = l1old_;
  l1old_ = t1tmp;
}

double CCLHWavefunction::increment_amps()
{
  int no = no_;
  int nv = nv_;
  double **D1 = D1_;
  double ****D2 = D2_;
  double **t1, **t1old, ****t2, ****t2old;

  t1 = l1_;
  t1old = l1old_;
  t2 = l2_;
  t2old = l2old_;

  double residual1 = 0.0;
  double residual2 = 0.0;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      residual1 += t1[i][a] * t1[i][a];     
      t1[i][a] = t1old[i][a] + t1[i][a]/D1[i][a];
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++) {
          residual2 += t2[i][j][a][b] * t2[i][j][a][b];
          t2[i][j][a][b] = t2old[i][j][a][b] + t2[i][j][a][b]/D2[i][j][a][b];
        }
    }

  return sqrt(residual1 + residual2);
}

void CCLHWavefunction::compute_lambda()
{
  outfile->Printf("\n\tThe Coupled-Cluster Lambda Iteration:\n");
  outfile->Printf(  "\t-------------------------------------\n");  outfile->Printf(  "\t Iter   Correlation Energy  RMS   \n");
  outfile->Printf(  "\t-------------------------------------\n");
  outfile->Printf(  "\t  %3d  %20.15f\n", 0, pseudoenergy());

  double rms = 0.0;
  for(int iter=1; iter <= maxiter(); iter++) {
    amp_save();    
    build_G();
    build_l1();
    build_l2();
    rms = increment_amps();
    if(rms < convergence()) break;
//    if(do_diis()) diis(iter, "L");
    outfile->Printf(  "\t  %3d  %20.15f  %5.3e\n",iter, pseudoenergy(), rms);
  }
  if(rms >= convergence())
    throw PSIEXCEPTION("Computation has not converged.");
}

void CCLHWavefunction::build_G()
{
  int no = no_;
  int nv = nv_;
  double ****t2 = t2_;
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

void CCLHWavefunction::build_l1()
{
  int no = no_;
  int nv = nv_;
  double **l1 = l1old_;
  double ****l2 = l2old_;
  double **l1new = l1_;
  double **Gvv = Gvv_;
  double **Goo = Goo_;
  double **Hov = Hov_;
  double **Hvv = Hvv_;
  double **Hoo = Hoo_;
  double ****Hvvvo = Hvvvo_;
  double ****Hovoo = Hovoo_;
  double ****Hovvo = Hovvo_;
  double ****Hovov = Hovov_;
  double ****Hvovv = Hvovv_;
  double ****Hooov = Hooov_;

  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      double value = 2 * Hov[i][a];

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

void CCLHWavefunction::build_l2()
{
  int no = no_;
  int nv = nv_;
  double **l1 = l1old_;
  double ****l2 = l2old_;
  double ****l2new = l2_;
  double ****L = H_->L_;
  double **Gvv = Gvv_;
  double **Goo = Goo_;
  double **Hov = Hov_;
  double **Hvv = Hvv_;
  double **Hoo = Hoo_;

  double ****Hovvo = Hovvo_;
  double ****Hovov = Hovov_;
  double ****Hoooo = Hoooo_;
  double ****Hvvvv = Hvvvv_;
  double ****Hvovv = Hvovv_;
  double ****Hooov = Hooov_;

  double ****Z = init_4d_array(no, no, nv, nv);

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          double value = L[i][j][a+no][b+no];

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
** pseudoenergy(): Evaluates an energy-like expression for the Lambda
doubles
** amplitudes: 
**   E = <0|L2 H|0> = 1/2 <ab|ij> L(ij,ab)
** This expression is derived in the UGA formalism.
*/
double CCLHWavefunction::pseudoenergy()
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

}} // psi::ugacc
