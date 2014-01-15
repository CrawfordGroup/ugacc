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

void density(void)
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

  double **Doo = block_matrix(no, no);
  double **Dvv = block_matrix(nv, nv);
  double **Dov = block_matrix(no, nv);
  double **Dvo = block_matrix(nv, no);

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++) {
      Doo[i][j] = 0.0;
      for(int e=0; e < nv; e++)
        Doo[i][j] -= t1[i][e] * l1[j][e];
      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          for(int f=0; f < nv; f++)
            Doo[i][j] -= t2[i][m][e][f] * l2[j][m][e][f];
      if(params.wfn == "CCSD_T") {
        for(int l=0; l < no; l++)
          for(int m=0; m < no; m++)
            for(int d=0; d < nv; d++)
              for(int e=0; e < nv; e++)
                for(int f=0; f < nv; f++)
                  Doo[i][j] -= 0.5 * t3[i][l][m][d][e][f] * l3[j][l][m][d][e][f];
      }
      Doo[i][j] /= 2.0;
    }

  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++) {
      Dvv[a][b] = 0.0;
      for(int m=0; m < no; m++)
        Dvv[a][b] += t1[m][b] * l1[m][a];
      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          for(int e=0; e < nv; e++)
            Dvv[a][b] += t2[m][n][b][e] * l2[m][n][a][e];
      if(params.wfn == "CCSD_T") {
        for(int l=0; l < no; l++)
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              for(int d=0; d < nv; d++)
                for(int e=0; e < nv; e++)
                  Dvv[a][b] += 0.5 * t3[l][m][n][b][d][e] * l3[l][m][n][a][d][e];
      }
      Dvv[a][b] /= 2.0;
    }

  for(int a=0; a < nv; a++)
    for(int i=0; i < no; i++) {
      Dvo[a][i] = l1[i][a];
      Dvo[a][i] /= 2.0;
    }

  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      Dov[i][a] = 2.0 * t1[i][a];
      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          Dov[i][a] += l1[m][e] * (2.0 * t2[i][m][a][e] - tau[m][i][a][e]);
      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          for(int e=0; e < nv; e++)
            for(int f=0; f < nv; f++)
              Dov[i][a] -= l2[m][n][e][f] * (t1[m][a] * t2[i][n][e][f] + t1[i][e] * t2[m][n][a][f]);
/*
      if(params.wfn == "CCSD_T") {
        for(int m=0; m < no; m++)
          for(int n=0; n < no; n++)
            for(int e=0; e < nv; e++)
              for(int f=0; f < nv; f++)
                Dov[i][a] += (t3[m][n][i][e][f][a] - t3[m][i][n][e][f][a]) * (4.0 * t2[m][n][e][f] - 2.0 * t2[m][n][f][e]);
      }
*/
      Dov[i][a] /= 2.0;
    }

  moinfo.Doo = Doo;
  moinfo.Dvv = Dvv;
  moinfo.Dov = Dov;
  moinfo.Dvo = Dvo;

  fprintf(outfile, "\tDij Matrix:\n");
  mat_print(Doo, no, no, outfile);
  fprintf(outfile, "\tDab Matrix:\n");
  mat_print(Dvv, nv, nv, outfile);
  fprintf(outfile, "\tDai Matrix:\n");
  mat_print(Dvo, nv, no, outfile);
  fprintf(outfile, "\tDia Matrix:\n");
  mat_print(Dov, no, nv, outfile);

/*
  double energy = 0.0;
  for(int i=0; i < no; i++)
    energy += fock[i][i] * Doo[i][i];
  for(int a=0; a < nv; a++)
    energy += fock[a+no][a+no] * Dvv[a][a];
  fprintf(outfile, "\tOne-electron energy = %20.14f\n", energy);
*/

  return;
}

}} // namespace psi::ugacc
