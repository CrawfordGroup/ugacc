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

double onepdm(void)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  double **t1 = moinfo.t1;
  double ****t2 = moinfo.t2;
  double ******t3 = moinfo.t3;
  double ******l3 = moinfo.l3;
  double ****tau = moinfo.tau;
  double **l1 = moinfo.l1;
  double ****l2 = moinfo.l2;
  double **fock = moinfo.fock;

  double **Doo = moinfo.Doo;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++) {
      for(int e=0; e < nv; e++)
        Doo[i][j] -= t1[i][e] * l1[j][e];
      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          for(int f=0; f < nv; f++)
            Doo[i][j] -= t2[i][m][e][f] * l2[j][m][e][f];
    }

  double **Dvv = moinfo.Dvv;
  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++) {
      for(int m=0; m < no; m++)
        Dvv[a][b] += t1[m][b] * l1[m][a];
      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          for(int e=0; e < nv; e++)
            Dvv[a][b] += t2[m][n][b][e] * l2[m][n][a][e];
    }

  double **Dvo = moinfo.Dvo;
  for(int a=0; a < nv; a++)
    for(int i=0; i < no; i++) {
      Dvo[a][i] = l1[i][a];
    }

  double **Dov = moinfo.Dov;
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      Dov[i][a] += 2.0 * t1[i][a];
      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          Dov[i][a] += l1[m][e] * (2.0 * t2[i][m][a][e] - tau[m][i][a][e]);
      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          for(int e=0; e < nv; e++)
            for(int f=0; f < nv; f++)
              Dov[i][a] -= l2[m][n][e][f] * (t1[m][a] * t2[i][n][e][f] + t1[i][e] * t2[m][n][a][f]);
    }

  double energy = 0.0;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      energy += fock[i][j] * Doo[i][j];
  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
    energy += fock[a+no][b+no] * Dvv[a][b];

  return energy;
}

}} // namespace psi::ugacc
