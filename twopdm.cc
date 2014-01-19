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

void twopdm(void)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  double **t1 = moinfo.t1;
  double ****t2 = moinfo.t2;
  double ****tau = moinfo.tau;
  double **l1 = moinfo.l1;
  double ****l2 = moinfo.l2;
  double ****ints = moinfo.ints;

  double ****Goooo = init_4d_array(no, no, no, no);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++)
        for(int l=0; l < no; l++) {
          Goooo[i][j][k][l] = 0.0;
          for(int e=0; e < nv; e++)
            for(int f=0; f < nv; f++)
              Goooo[i][j][k][l] += tau[i][j][e][f] * l2[k][l][e][f];
        }
  moinfo.Goooo = Goooo;

  double Eoooo = 0.0;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++)
        for(int l=0; l < no; l++)
          Eoooo += 0.5 * ints[i][j][k][l] * Goooo[i][j][k][l];
  fprintf(outfile, "OOOO Energy = %20.14f\n", Eoooo);

  double ****Gvvvv = init_4d_array(nv, nv, nv, nv);
  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int c=0; c < nv; c++)
        for(int d=0; d < nv; d++) {
          Gvvvv[a][b][c][d] = 0.0;
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              Gvvvv[a][b][c][d] += tau[m][n][a][b] * l2[m][n][c][d];
        }
  moinfo.Gvvvv = Gvvvv;

  double Evvvv = 0.0;
  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int c=0; c < nv; c++)
        for(int d=0; d < nv; d++)
          Evvvv += 0.5 * ints[a+no][b+no][c+no][d+no] * Gvvvv[a][b][c][d];
  fprintf(outfile, "VVVV Energy = %20.14f\n", Evvvv);

  double ****Gooov = init_4d_array(no, no, no, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++)
        for(int a=0; a < nv; a++) {
          Gooov[i][j][k][a] = 0.0;
          for(int e=0; e < nv; e++)
            Gooov[i][j][k][a] -= l1[k][e] * (2.0 * tau[i][j][e][a] - tau[i][j][a][e]);
          for(int e=0; e < nv; e++)
            Gooov[i][j][k][a] -= t1[i][e] * l2[j][k][a][e];
          for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++)
              for(int f=0; f < nv; f++)
                Gooov[i][j][k][a] -= l2[k][m][e][f] *
                           (2.0 * (t1[j][a] * t2[i][m][e][f] + t1[i][e] * t2[j][m][a][f])
                                - (t1[i][a] * t2[j][m][e][f] + t1[j][e] * t2[i][m][a][f])
                                - (t1[m][a] * t2[i][j][e][f] + t1[i][e] * t2[m][j][a][f] + t1[j][f] * t2[i][m][e][a]));
          for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++)
              for(int f=0; f < nv; f++)
                Gooov[i][j][k][a] += l2[k][m][e][f] * t1[m][a] * t1[i][e] * t1[j][f];
        }
  moinfo.Gooov = Gooov;

  double Eooov = 0.0;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++)
        for(int a=0; a < nv; a++)
          Eooov += ints[i][j][k][a+no] * Gooov[i][j][k][a];
  fprintf(outfile, "OOOV Energy = %20.14f\n", Eooov);

  double ****Gvvvo = init_4d_array(nv, nv, nv, no);
  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int c=0; c < nv; c++)
        for(int i=0; i < no; i++) {
          Gvvvo[a][b][c][i] = 0.0;
          for(int m=0; m < no; m++)
            Gvvvo[a][b][c][i] += l1[m][c] * (2.0 * tau[m][i][a][b] - tau[i][m][a][b]);
          for(int m=0; m < no; m++)
            Gvvvo[a][b][c][i] += t1[m][a] * l2[i][m][b][c];
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              for(int e=0; e < nv; e++)
                Gvvvo[a][b][c][i] += l2[n][m][c][e] * 
                           (2.0 * (t1[n][a] * t2[i][m][b][e] + t1[i][b] * t2[n][m][a][e])
                                - (t1[n][b] * t2[i][m][a][e] + t1[i][a] * t2[n][m][b][e])
                                - (t1[m][b] * t2[n][i][a][e] + t1[n][a] * t2[m][i][b][e] + t1[i][e] * t2[n][m][a][b]));
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              for(int e=0; e < nv; e++)
                Gvvvo[a][b][c][i] -= l2[n][m][c][e] * t1[m][b] * t1[n][a] * t1[i][e];
        }
  moinfo.Gvvvo = Gvvvo;

  double Evvvo = 0.0;
  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int c=0; c < nv; c++)
        for(int i=0; i < no; i++)
          Evvvo += ints[a+no][b+no][c+no][i] * Gvvvo[a][b][c][i];
  fprintf(outfile, "VVVO Energy = %20.14f\n", Evvvo);

  double ****Govov = init_4d_array(no, nv, no, nv);
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++) {
          Govov[i][a][j][b] -= t1[i][a] * l1[j][b];
          for(int m=0; m < no; m++) 
            for(int e=0; e < nv; e++)
              Govov[i][a][j][b] -= (tau[m][i][b][e] * l2[j][m][e][a] + t2[i][m][b][e] * l2[m][j][e][a]);
        }
  moinfo.Govov = Govov;

  double Eovov = 0.0;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++)
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++)
          Eovov += ints[i][a+no][j][b+no] * Govov[i][a][j][b];
  fprintf(outfile, "OVOV Energy = %20.14f\n", Eovov);

  double ****Goovv = init_4d_array(no, no, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          Goovv[i][j][a][b] = 4.0 * t1[i][a] * l1[j][b];
          Goovv[i][j][a][b] += 4.0 * tau[i][j][a][b] - 2.0 * tau[i][j][b][a];
          Goovv[i][j][a][b] += l2[i][j][a][b];
          for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++)
              Goovv[i][j][a][b] += 2.0 * l1[m][e] * 
                         (2.0 * t1[i][a] * (2.0 * t2[j][m][b][e] - t2[j][m][e][b])
                              - t1[i][b] * (2.0 * t2[j][m][a][e] - t2[j][m][e][a])
                              - t1[i][e] * (2.0 * tau[j][m][b][a] - tau[j][m][a][b])
                              - t1[m][a] * (2.0 * t2[i][j][e][b] - t2[i][j][b][e]));
          for(int m=0; m < no; m++)
            for(int e=0; e < nv; e++)
              Goovv[i][j][a][b] += 4.0 * t2[i][m][b][e] * l2[m][j][e][a] - 2.0 * tau[m][j][b][e] * l2[i][m][a][e];
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              for(int e=0; e < nv; e++)
                for(int f=0; f < nv; f++)
                  Goovv[i][j][a][b] += 2.0 * l2[m][n][e][f] * 
                      (0.5 * (t2[m][n][a][b] * t2[i][j][e][f] + t2[m][i][a][e] * t2[n][j][b][f] + t2[n][j][a][e] * t2[i][m][f][b])
                     + 2.0 * (- tau[i][j][a][e] * t2[m][n][b][f] - tau[i][m][a][b] * t2[j][n][e][f] 
                              - tau[m][j][b][e] * t2[i][n][a][f] + t2[i][m][a][e] * t2[j][n][b][f])
                           - (- tau[i][j][b][e] * t2[m][n][a][f] - tau[i][m][b][a] * t2[j][n][e][f] 
                              - tau[m][j][a][e] * t2[i][n][b][f] + t2[i][m][b][e] * t2[j][n][a][f]));
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              for(int e=0; e < nv; e++)
                for(int f=0; f < nv; f++)
                  Goovv[i][j][a][b] += 2.0 * l2[m][n][e][f] *
                       (t1[m][a]*t1[n][b]*t2[i][j][e][f] + t1[m][a]*t1[i][e]*t2[n][j][b][f] + t1[n][a]*t1[j][e]*t2[i][m][f][b]);
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              for(int e=0; e < nv; e++)
                for(int f=0; f < nv; f++)
                  Goovv[i][j][a][b] += l2[m][n][e][f] * t1[m][a] * t1[n][b] * t1[i][e] * t1[j][f];
        }
  moinfo.Govov = Govov;

  double Eoovv = 0.0;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          Eoovv += 0.5 * ints[i][j][a+no][b+no] * Goovv[i][j][a][b];
  fprintf(outfile, "OOVV Energy = %20.14f\n", Eoovv);
  fprintf(outfile, "OVOV + OOVV = %20.14f\n", Eovov+Eoovv);

  return;
}

}} // namespace psi::ugacc
