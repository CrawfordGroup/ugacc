#include <string>
#include <cstdio>
#include <psi4-dec.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void make_Z_amps(double **l1, double ****l2);

void tgrad(void)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  double **fock = moinfo.fock;
  double ****ints = moinfo.ints;
  double ****L = moinfo.L;
  double **t1 = moinfo.t1;
  double ****t2 = moinfo.t2;
  double ******t3 = moinfo.t3;
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

  // Lambda3 amplitudes
  double **t1s = moinfo.t1s;
  double ****t2s = moinfo.t2s;
  double ******l3_tmp = init_6d_array(no, no, no, nv, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++) {

        for(int a=0; a < nv; a++)
          for(int b=0; b < nv; b++)
            for(int c=0; c < nv; c++) {
              double value = 0.0;
  
              value += fock[i][a+no] * t2s[j][k][b][c] - fock[j][a+no] * t2s[i][k][b][c];
              value += L[i][j][a+no][b+no] * t1s[k][c] - L[i][j][a+no][c+no] * t1s[k][b];

              for(int f=0; f < nv; f++) {
                value += L[f+no][j][a+no][b+no] * t2s[i][k][f][c];
                value -= ints[k][f+no][a+no][b+no] * t2s[i][j][c][f];
              }
              for(int n=0; n < no; n++) {
                value -= L[i][j][a+no][n] * t2s[n][k][b][c];
                value += ints[k][j][a+no][n] * t2s[n][i][b][c];
              }

              double denom = fock[i][i] + fock[j][j] + fock[k][k];
              denom -= fock[a+no][a+no] + fock[b+no][b+no] + fock[c+no][c+no];

              l3_tmp[i][j][k][a][b][c] = value/denom;
            }
      }

  double ******l3 = init_6d_array(no, no, no, nv, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++)
        for(int a=0; a < nv; a++)
          for(int b=0; b < nv; b++)
            for(int c=0; c < nv; c++) 
              l3[i][j][k][a][b][c] = l3_tmp[i][j][k][a][b][c]
                                   + l3_tmp[i][k][j][a][c][b]
                                   + l3_tmp[j][i][k][b][a][c]
                                   + l3_tmp[j][k][i][b][c][a]
                                   + l3_tmp[k][i][j][c][a][b]
                                   + l3_tmp[k][j][i][c][b][a];

  moinfo.l3 = l3;
  free_6d_array(l3_tmp, no, no, no, nv, nv);

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

  free_4d_array(X2, no, no, nv);
  free_4d_array(Y2, no, no, nv);

  make_Z_amps(moinfo.s1, moinfo.s2);

  // (T) density contributions
  // OO
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++) {
      moinfo.Doo[i][j] = 0.0;
      for(int l=0; l < no; l++)
        for(int m=0; m < no; m++) 
          for(int d=0; d < nv; d++)
            for(int e=0; e < nv; e++)
              for(int f=0; f < nv; f++)
                moinfo.Doo[i][j] -= 0.5 * t3[i][l][m][d][e][f] * l3[j][l][m][d][e][f];
    }

  // VV
  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++) {
      moinfo.Dvv[a][b] = 0.0;
      for(int l=0; l < no; l++)
        for(int m=0; m < no; m++)
          for(int n=0; n < no; n++)
            for(int d=0; d < nv; d++)
              for(int e=0; e < nv; e++)
                moinfo.Dvv[a][b] += 0.5 * t3[l][m][n][b][d][e] * l3[l][m][n][a][d][e];
    }

  // OV
/*
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      moinfo.Dov[i][a] = 0.0;
      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          for(int e=0; e < nv; e++)
            for(int f=0; f < nv; f++)
              moinfo.Dov[i][a] += (t3[m][n][i][e][f][a] - t3[m][i][n][e][f][a]) * t2s[m][n][e][f];
    }
*/

  // OOOV
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++)
        for(int a=0; a < nv; a++) {
          moinfo.Gooov[i][j][k][a] = 0.0;
          for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++)
              for(int f=0; f < nv; f++)
                moinfo.Gooov[i][j][k][a] -= (4.0 * t2[k][m][e][f] - 2.0 * t2[k][m][f][e]) *
                        (2.0 * t3[j][i][m][a][e][f] - t3[i][j][m][a][e][f] - t3[m][i][j][a][e][f]);

          for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++)
              for(int f=0; f < nv; f++)
                moinfo.Gooov[i][j][k][a] -= t2[k][m][e][f] * l3[m][j][i][f][a][e];
        }

  // VVVO
  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int c=0; c < nv; c++)
        for(int i=0; i < no; i++) {
          moinfo.Gvvvo[a][b][c][i] = 0.0;
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              for(int e=0; e < nv; e++)
                moinfo.Gvvvo[a][b][c][i] += (4.0 * t2[m][n][e][c] - 2.0 * t2[m][n][c][e]) *
                      (2.0 * t3[n][i][m][a][b][e] - t3[n][i][m][b][a][e] - t3[n][i][m][a][e][b]);

          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              for(int e=0; e < nv; e++)
                moinfo.Gvvvo[a][b][c][i] += t2[m][n][e][c] * l3[n][i][m][a][b][e];
        }

  // OOVV
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          moinfo.Goovv[i][j][a][b] = 0.0;
          for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++)
              moinfo.Goovv[i][j][a][b] += 4.0 * t1[m][e] *
                  ( 2.0 * (t3[i][j][m][a][b][e] - t3[i][j][m][a][e][b]) - (t3[i][j][m][b][a][e] - t3[i][j][m][b][e][a]) );
        }

  return;
}

}} // namespace psi::ugacc
