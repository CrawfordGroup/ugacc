#include <string>
#include <cstdio>
#include <psi4-dec.h>
#include <libciomr/libciomr.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void amp_write(int, double **, double ****, std::string);
void triples_gradient_viking(double ******);
void make_Z_amps(double **l1, double ****l2);


double triples(void)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  double **fock = moinfo.fock;
  double ****ints = moinfo.ints;
  double ****L = moinfo.L;
  double **t1 = moinfo.t1;
  double ****t2 = moinfo.t2;
  double ******t3;

  // Scandinavian expression for (T) correction
  t3 = init_6d_array(no, no, no, nv, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++) {

        for(int a=0; a < nv; a++)
          for(int b=0; b < nv; b++)
            for(int c=0; c < nv; c++) {
              double value = 0.0;
              for(int e=0; e < nv; e++) {
                value +=
                   + ints[i][e+no][a+no][b+no] * t2[k][j][c][e]
                   + ints[i][e+no][a+no][c+no] * t2[j][k][b][e]
                   + ints[k][e+no][c+no][a+no] * t2[j][i][b][e]
                   + ints[k][e+no][c+no][b+no] * t2[i][j][a][e]
                   + ints[j][e+no][b+no][c+no] * t2[i][k][a][e]
                   + ints[j][e+no][b+no][a+no] * t2[k][i][c][e];
              }
              for(int m=0; m < no; m++) {
                value -=
                   + ints[j][k][m][c+no] * t2[i][m][a][b]
                   + ints[k][j][m][b+no] * t2[i][m][a][c]
                   + ints[i][j][m][b+no] * t2[k][m][c][a]
                   + ints[j][i][m][a+no] * t2[k][m][c][b]
                   + ints[k][i][m][a+no] * t2[j][m][b][c]
                   + ints[i][k][m][c+no] * t2[j][m][b][a];
              }
              double denom = fock[i][i] + fock[j][j] + fock[k][k];
              denom -= fock[a+no][a+no] + fock[b+no][b+no] + fock[c+no][c+no];

              t3[i][j][k][a][b][c] = value/denom;
            }
      }
  moinfo.t3 = t3;

  double **X1 = block_matrix(no, nv);
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      for(int k=0; k < no; k++)
        for(int l=0; l < no; l++)
          for(int c=0; c < nv; c++)
            for(int d=0; d < nv; d++)
              X1[i][a] += (t3[i][k][l][a][c][d] - t3[l][k][i][a][c][d]) * L[k][l][c+no][d+no];

  double ****X2 = init_4d_array(no, no, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {

          for(int k=0; k < no; k++)
            for(int c=0; c < nv; c++)
              X2[i][j][a][b] += (t3[i][j][k][a][b][c] - t3[k][j][i][a][b][c]) * fock[k][c+no];

          for(int k=0; k < no; k++)
            for(int c=0; c < nv; c++)
              for(int l=0; l < no; l++) {
                X2[i][j][a][b] -= t3[i][k][l][a][b][c] * L[k][l][j][c+no];
                X2[i][j][a][b] += t3[l][k][i][a][b][c] * ints[k][l][j][c+no];
              }

          for(int k=0; k < no; k++)
            for(int c=0; c < nv; c++)
              for(int d=0; d < nv; d++) {
                X2[i][j][a][b] += t3[i][j][k][a][c][d] * L[b+no][k][c+no][d+no];
                X2[i][j][a][b] -= t3[k][j][i][a][c][d] * ints[b+no][k][c+no][d+no];
              }
        }
  

  double ET_UGA = 0.0;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      ET_UGA += t1[i][a] * X1[i][a];
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++)
          ET_UGA += (2.0*t2[i][j][a][b] - t2[i][j][b][a]) * X2[i][j][a][b];    
    }
  ET_UGA *= 2.0;

  if(params.dertype) triples_gradient_viking(t3);

  free_block(X1);
  free_4d_array(X2, no, no, nv);

  return ET_UGA;
}

void triples_gradient_viking(double ******t3)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  double **fock = moinfo.fock;
  double ****ints = moinfo.ints;
  double ****L = moinfo.L;
  double **t1 = moinfo.t1;
  double ****t2 = moinfo.t2;
  double **s1 = moinfo.s1;
  double ****s2 = moinfo.s2;

  // T3 --> Lambda1
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      for(int k=0; k < no; k++)
        for(int l=0; l < no; l++)
          for(int c=0; c < nv; c++)
            for(int d=0; d < nv; d++)
              s1[i][a] += 2.0 * (t3[i][k][l][a][c][d] - t3[k][i][l][a][c][d]) * L[k][l][c+no][d+no];

  // T3 --> Lambda2
  double ****Y2 = init_4d_array(no, no, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          double value = 0.0;         

          for(int k=0; k < no; k++)
            for(int c=0; c < nv; c++)
              value += (t3[i][j][k][a][b][c] - t3[k][j][i][a][b][c]) * fock[k][c+no];

          for(int k=0; k < no; k++)
            for(int c=0; c < nv; c++)
              for(int l=0; l < no; l++)
                value -= (2.0 * t3[j][k][l][b][a][c] - t3[j][l][k][b][a][c] - t3[l][k][j][b][a][c]) * ints[k][l][i][c+no];

          for(int k=0; k < no; k++)
            for(int c=0; c < nv; c++)
              for(int d=0; d < nv; d++)
                value += (2.0 * t3[j][i][k][b][c][d] - t3[j][k][i][b][c][d] - t3[k][i][j][b][c][d]) * ints[a+no][k][c+no][d+no];

          Y2[i][j][a][b] = value;
        }

  double ****X2 = init_4d_array(no, no, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
           X2[i][j][a][b] += Y2[i][j][a][b] + Y2[j][i][b][a];

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          s2[i][j][a][b] += 4.0 * X2[i][j][a][b] - 2.0 * X2[i][j][b][a];

  // Special left-projection T amplitudes
  double **t1s = block_matrix(no, nv);
  double ****t2s = init_4d_array(no, no, nv, nv);
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      t1s[i][a] = 2.0 * t1[i][a];

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          t2s[i][j][a][b] = 4.0 * t2[i][j][a][b] - 2.0 * t2[i][j][b][a];

  // Lambda3 amplitudes
  double ******l3 = init_6d_array(no, no, no, nv, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++) {

        for(int a=0; a < nv; a++)
          for(int b=0; b < nv; b++)
            for(int c=0; c < nv; c++) {
              double value = 0.0;
  
              value += fock[i][a+no] * t2s[j][k][b][c] - fock[j][a+no] * t2s[i][k][b][c];
              value += fock[i][a+no] * t2s[k][j][c][b] - fock[k][a+no] * t2s[i][j][c][b];
              value += fock[j][b+no] * t2s[i][k][a][c] - fock[i][b+no] * t2s[j][k][a][c];
              value += fock[j][b+no] * t2s[k][i][c][a] - fock[k][b+no] * t2s[j][i][c][a];
              value += fock[k][c+no] * t2s[i][j][a][b] - fock[i][c+no] * t2s[k][j][a][b];
              value += fock[k][c+no] * t2s[j][i][b][a] - fock[j][c+no] * t2s[k][i][b][a];

              value += L[i][j][a+no][b+no] * t1s[k][c] - L[i][j][a+no][c+no] * t1s[k][b];
              value += L[i][k][a+no][c+no] * t1s[j][b] - L[i][k][a+no][b+no] * t1s[j][c];
              value += L[j][i][b+no][a+no] * t1s[k][c] - L[j][i][b+no][c+no] * t1s[k][a];
              value += L[j][k][b+no][c+no] * t1s[i][a] - L[j][k][b+no][a+no] * t1s[i][c];
              value += L[k][i][c+no][a+no] * t1s[j][b] - L[k][i][c+no][b+no] * t1s[j][a];
              value += L[k][j][c+no][b+no] * t1s[i][a] - L[k][j][c+no][a+no] * t1s[i][b];

              for(int f=0; f < nv; f++) {
                value += L[f+no][j][a+no][b+no] * t2s[i][k][f][c];
                value += L[f+no][k][a+no][c+no] * t2s[i][j][f][b];
                value += L[f+no][i][b+no][a+no] * t2s[j][k][f][c];
                value += L[f+no][k][b+no][c+no] * t2s[j][i][f][a];
                value += L[f+no][i][c+no][a+no] * t2s[k][j][f][b];
                value += L[f+no][j][c+no][b+no] * t2s[k][i][f][a];

                value -= ints[k][f+no][a+no][b+no] * t2s[i][j][c][f];
                value -= ints[j][f+no][a+no][c+no] * t2s[i][k][b][f];
                value -= ints[k][f+no][b+no][a+no] * t2s[j][i][c][f];
                value -= ints[i][f+no][b+no][c+no] * t2s[j][k][a][f];
                value -= ints[j][f+no][c+no][a+no] * t2s[k][i][b][f];
                value -= ints[i][f+no][c+no][b+no] * t2s[k][j][a][f];
              }
              for(int n=0; n < no; n++) {
                value -= L[i][j][a+no][n] * t2s[n][k][b][c];
                value -= L[i][k][a+no][n] * t2s[n][j][c][b];
                value -= L[j][i][b+no][n] * t2s[n][k][a][c];
                value -= L[j][k][b+no][n] * t2s[n][i][c][a];
                value -= L[k][i][c+no][n] * t2s[n][j][a][b];
                value -= L[k][j][c+no][n] * t2s[n][i][b][a];

                value += ints[k][j][a+no][n] * t2s[n][i][b][c];
                value += ints[j][k][a+no][n] * t2s[n][i][c][b];
                value += ints[k][i][b+no][n] * t2s[n][j][a][c];
                value += ints[i][k][b+no][n] * t2s[n][j][c][a];
                value += ints[j][i][c+no][n] * t2s[n][k][a][b];
                value += ints[i][j][c+no][n] * t2s[n][k][b][a];
              }

              double denom = fock[i][i] + fock[j][j] + fock[k][k];
              denom -= fock[a+no][a+no] + fock[b+no][b+no] + fock[c+no][c+no];

              l3[i][j][k][a][b][c] = value/denom;
            }
      }
  moinfo.l3 = l3;

  // Lambda3 --> Lambda2
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          Y2[i][j][a][b] = 0.0;

          for(int k=0; k < no; k++)
            for(int c=0; c < nv; c++)
              for(int l=0; l < no; l++)
                Y2[i][j][a][b] -= l3[i][k][l][a][c][b] * ints[c+no][j][k][l];

          for(int k=0; k < no; k++)
            for(int c=0; c < nv; c++)
              for(int d=0; d < nv; d++)
                Y2[i][j][a][b] += l3[i][k][j][a][c][d] * ints[c+no][d+no][k][b+no];
        }

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          s2[i][j][a][b] += Y2[i][j][a][b] + Y2[j][i][b][a];

  free_block(t1s);
  free_4d_array(t2s, no, no, nv);
  free_4d_array(X2, no, no, nv);
  free_4d_array(Y2, no, no, nv);

  make_Z_amps(moinfo.s1, moinfo.s2);

  return;
}

}} // namespace psi::ugacc
