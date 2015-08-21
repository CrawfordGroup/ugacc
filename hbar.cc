#include "hbar.h"
#include "array.h"
#include <boost/shared_ptr.hpp>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>

namespace psi { namespace ugacc {

HBAR::HBAR(boost::shared_ptr<Hamiltonian> H, boost::shared_ptr<CCWavefunction> CC)
{
  CC_ = CC;
  H_ = H;

  no_ = CC_->no_;
  nv_ = CC_->nv_;
  double **fock = H_->fock_;
  double ****ints = H_->ints_;
  double ****L = H_->L_;
  double **t1 = CC_->t1_;
  double ****t2 = CC_->t2_;
  double ****tau = CC_->tau_;

  int no = no_;
  int nv = nv_;
  Hvv_ = block_matrix(nv, nv);
  Hoo_ = block_matrix(no, no);
  Hov_ = block_matrix(no, nv);
  Hoooo_ = init_4d_array(no, no, no, no);
  Hvvvv_ = init_4d_array(nv, nv, nv, nv);
  Hovov_ = init_4d_array(no, nv, no, nv);
  Hovvo_ = init_4d_array(no, nv, nv, no);
  Hvovv_ = init_4d_array(nv, no, nv, nv);
  Hooov_ = init_4d_array(no, no, no, nv);
  Hovoo_ = init_4d_array(no, nv, no, no);
  Hvvvo_ = init_4d_array(nv, nv, nv, no);

  for(int m=0; m < no; m++)
    for(int e=0; e < nv; e++) {
      double value = fock[m][e+no];
      for(int n=0; n < no; n++)
        for(int f=0; f < nv; f++)
          value += t1[n][f]*L[m][n][e+no][f+no];
      Hov_[m][e] = value;
    }

  for(int m=0; m < no; m++)
    for(int i=0; i < no; i++) {
      double value = fock[m][i];
      for(int e=0; e < nv; e++) {
        value += t1[i][e]*fock[m][e+no];
        for(int n=0; n < no; n++) {
          value += t1[n][e]*L[m][n][i][e+no];
          for(int f=0; f < nv; f++)
            value += tau[i][n][e][f]*L[m][n][e+no][f+no];
        }
      }
      Hoo_[m][i] = value;
    }

  for(int a=0; a < nv; a++)
    for(int e=0; e < nv; e++) {
      double value = fock[a+no][e+no];
      for(int m=0; m < no; m++) {
        value -= fock[m][e+no]*t1[m][a];
        for(int f=0; f < nv; f++) {
          value += t1[m][f]*L[a+no][m][e+no][f+no];
          for(int n=0; n < no; n++)
            value -= tau[m][n][f][a]*L[m][n][f+no][e+no];
        }
      }
      Hvv_[a][e] = value;
    }

  for(int m=0; m < no; m++)
    for(int n=0; n < no; n++)
      for(int i=0; i < no; i++)
        for(int j=0; j < no; j++) {
          double value = ints[m][n][i][j];
          for(int e=0; e < nv; e++) {
            value += 2.0*t1[j][e]*ints[m][n][i][e+no];
            for(int f=0; f < nv; f++)
              value += tau[i][j][e][f]*ints[m][n][e+no][f+no];
          }
          Hoooo_[m][n][i][j] = value;
        }

  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int e=0; e < nv; e++)
        for(int f=0; f < nv; f++) {
          double value = ints[a+no][b+no][e+no][f+no];
          for(int m=0; m < no; m++) {
            value -= 2.0*t1[m][b]*ints[a+no][m][e+no][f+no];
            for(int n=0; n < no; n++)
              value += tau[m][n][a][b]*ints[m][n][e+no][f+no];
          }
          Hvvvv_[a][b][e][f] = value;
        }

  for(int a=0; a < nv; a++)
    for(int m=0; m < no; m++)
      for(int e=0; e < nv; e++)
        for(int f=0; f < nv; f++) {
          double value = ints[a+no][m][e+no][f+no];
          for(int n=0; n < no; n++)
            value -= t1[n][a] * ints[n][m][e+no][f+no];
          Hvovv_[a][m][e][f] = value;
        }

  for(int m=0; m < no; m++)
    for(int n=0; n < no; n++)
      for(int i=0; i < no; i++)
        for(int e=0; e < nv; e++) {
          double value = ints[m][n][i][e+no];
          for(int f=0; f < nv; f++)
            value += t1[i][f] * ints[n][m][e+no][f+no];
          Hooov_[m][n][i][e] = value;
        }

  for(int m=0; m < no; m++)
    for(int b=0; b < nv; b++)
      for(int e=0; e < nv; e++)
        for(int j=0; j < no; j++) {
          double value = ints[m][b+no][e+no][j];
          for(int f=0; f < nv; f++)
            value += t1[j][f]*ints[m][b+no][e+no][f+no];
          for(int n=0; n < no; n++)
            value -= t1[n][b]*ints[m][n][e+no][j];
          for(int n=0; n < no; n++)
            for(int f=0; f < nv; f++)
              value -= ints[m][n][e+no][f+no]*tau[n][j][b][f];
          for(int n=0; n < no; n++)
            for(int f=0; f < nv; f++)
              value += t2[n][j][f][b]*L[m][n][e+no][f+no];
          Hovvo_[m][b][e][j] = value;
        }

  for(int m=0; m < no; m++)
    for(int b=0; b < nv; b++)
      for(int j=0; j < no; j++)
        for(int e=0; e < nv; e++) {
          double value = ints[m][b+no][j][e+no];
          for(int f=0; f < nv; f++)
            value += t1[j][f]*ints[b+no][m][e+no][f+no];
          for(int n=0; n < no; n++)
            value -= t1[n][b]*ints[m][n][j][e+no];
          for(int n=0; n < no; n++) {
            for(int f=0; f < nv; f++)
              value -= ints[n][m][e+no][f+no]*tau[j][n][f][b];
          }
          Hovov_[m][b][j][e] = value;
        }

  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int e=0; e < nv; e++)
        for(int i=0; i < no; i++) {
          double value = ints[a+no][b+no][e+no][i];
          for(int f=0; f < nv; f++)
            value += t1[i][f]*ints[a+no][b+no][e+no][f+no];
          for(int m=0; m < no; m++)
            value -= t1[m][b]*ints[a+no][m][e+no][i] +
                     t1[m][a]*ints[b+no][m][i][e+no];
          for(int m=0; m < no; m++)
            for(int f=0; f < nv; f++)
              value -= tau[i][m][f][a]*ints[m][b+no][e+no][f+no] +
                       tau[i][m][f][b]*ints[a+no][m][e+no][f+no];
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              value += tau[m][n][a][b]*ints[m][n][e+no][i];
          for(int m=0; m < no; m++)
            value -= fock[m][e+no]*t2[m][i][a][b];
          for(int m=0; m < no; m++)
            for(int f=0; f < nv; f++)
              value += t2[m][i][f][b]*L[a+no][m][e+no][f+no];
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              for(int f=0; f < nv; f++)
                value += ints[m][n][e+no][f+no]*(t1[i][f]*t2[m][n][a][b]+t1[m][a]*t2[i][n][f][b]+t1[n][b]*t2[m][i][a][f]);
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              for(int f=0; f < nv; f++)
                value -= (t1[m][f]*t2[n][i][a][b]+t1[n][a]*t2[m][i][f][b])*L[m][n][f+no][e+no];
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              for(int f=0; f < nv; f++)
                value += t1[i][f]*t1[m][a]*t1[n][b]*ints[m][n][e+no][f+no];

          Hvvvo_[a][b][e][i] = value;
        }

  for(int m=0; m < no; m++)
    for(int b=0; b < nv; b++)
      for(int i=0; i < no; i++)
        for(int j=0; j < no; j++) {
          double value = ints[m][b+no][i][j];
          for(int e=0; e < nv; e++)
            value += t1[j][e]*ints[m][b+no][i][e+no] +
                     t1[i][e]*ints[b+no][m][j][e+no];
          for(int n=0; n < no; n++)
            value -= t1[n][b]*ints[m][n][i][j];
          for(int n=0; n < no; n++)
            for(int e=0; e < nv; e++)
              value -= tau[i][n][e][b]*ints[n][m][j][e+no] +
                       tau[j][n][e][b]*ints[m][n][i][e+no];
          for(int e=0; e < nv; e++)
            for(int f=0; f < nv; f++)
              value += ints[m][b+no][e+no][f+no]*tau[i][j][e][f];
          for(int e=0; e < nv; e++)
            value += fock[m][e+no]*t2[i][j][e][b];
          for(int n=0; n < no; n++)
            for(int e=0; e < nv; e++)
              value += L[m][n][i][e+no]*t2[n][j][e][b];
          for(int e=0; e < nv; e++)
            for(int n=0; n < no; n++)
              for(int f=0; f < nv; f++)
                value -= ints[m][n][e+no][f+no]*(t1[j][f]*t2[i][n][e][b]+t1[i][e]*t2[j][n][f][b]+t1[n][b]*t2[i][j][e][f]);
          for(int e=0; e < nv; e++)
            for(int n=0; n < no; n++)
              for(int f=0; f < nv; f++)
                value += (t1[i][e]*t2[n][j][f][b]+t1[n][f]*t2[i][j][e][b])*L[m][n][e+no][f+no];
          for(int e=0; e < nv; e++)
            for(int n=0; n < no; n++)
              for(int f=0; f < nv; f++)
                value -= t1[j][f]*t1[i][e]*t1[n][b]*ints[m][n][e+no][f+no];

          Hovoo_[m][b][i][j] = value;
        }
}

HBAR::~HBAR()
{
  int no = no_;
  int nv = nv_;
  free_block(Hoo_);
  free_block(Hvv_);
  free_block(Hov_);
  free_4d_array(Hoooo_, no, no, no);
  free_4d_array(Hvvvv_, nv, nv, nv);
  free_4d_array(Hovov_, no, nv, no);
  free_4d_array(Hovvo_, no, nv, nv);
  free_4d_array(Hvovv_, nv, no, nv);
  free_4d_array(Hooov_, no, no, no);
  free_4d_array(Hovoo_, no, nv, no);
  free_4d_array(Hvvvo_, nv, nv, nv);
}

}} // psi::ugacc
