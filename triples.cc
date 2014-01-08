#include <string>
#include <cstdio>
#include <psi4-dec.h>
#include <libciomr/libciomr.h>
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void amp_write(int, double **, double ****, std::string);

double triples(void)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  double **fock = moinfo.fock;
  double ****ints = moinfo.ints;
  double **t1 = moinfo.t1;
  double ****t2 = moinfo.t2;
  double **s1 = moinfo.s1;
  double ****s2 = moinfo.s2;

  double ***W = init_3d_array(nv, nv, nv);
  double ***V = init_3d_array(nv, nv, nv);
  double ***M = init_3d_array(nv, nv, nv);

  double ET = 0.0;
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

              W[a][b][c] = value;


              V[a][b][c] = ints[j][k][b+no][c+no] * t1[i][a]
                         + ints[i][k][a+no][c+no] * t1[j][b]
                         + ints[i][j][a+no][b+no] * t1[k][c];

            } /* abc loop */

        for(int a=0; a < nv; a++)
          for(int b=0; b < nv; b++)
            for(int c=0; c < nv; c++) {
              double denom = fock[i][i] + fock[j][j] + fock[k][k];
              denom -= fock[a+no][a+no] + fock[b+no][b+no] + fock[c+no][c+no];

              ET += (W[a][b][c] + V[a][b][c] - W[c][b][a] - V[c][b][a])*
                    (4.0*W[a][b][c] + W[b][c][a] + W[c][a][b])/(3.0*denom);
            } /* abc loop */

        for(int a=0; a < nv; a++)
          for(int b=0; b < nv; b++)
            for(int c=0; c < nv; c++) {

              double denom = fock[i][i] + fock[j][j] + fock[k][k];
              denom -= fock[a+no][a+no] + fock[b+no][b+no] + fock[c+no][c+no];

              s1[i][a] += ints[j][k][b+no][c+no] *
                 (4.0 * W[a][b][c] + W[b][c][a] + W[c][a][b]
                - 2.0 * W[c][b][a] - 2.0 * W[a][c][b]
                - 2.0 * W[b][a][c])/denom;

              M[a][b][c] = 6.0 * (
                  8.0 * W[a][b][c] + 2.0 * W[c][a][b] + 2.0 * W[b][c][a]
                - 4.0 * W[a][c][b] - 4.0 * W[b][a][c] - 4.0 * W[c][b][a]
                + 4.0 * V[a][b][c] + 1.0 * V[c][a][b] + 1.0 * V[b][c][a]
                - 2.0 * V[a][c][b] - 2.0 * V[b][a][c] - 2.0 * V[c][b][a]
                )/denom;

//              for(int d=0; d < nv; d++)
//                s2[k][j][c][d] += M[a][b][c] * ints[i][b+no][a+no][d+no];
//              for(int l=0; l < no; l++)
//                s2[k][l][c][b] -= M[a][b][c] * ints[i][j][a+no][l];

            } /* abc loop */
    
      } /* ijk loop */

  free_3d_array(W, nv, nv);
  free_3d_array(V, nv, nv);
  free_3d_array(M, nv, nv);

  amp_write(20, moinfo.s1, moinfo.s2, "S"); fprintf(outfile, "\n");

  // Also print non-UGA version of these amps for comparison to
  // spin-orbital code
  double **Z1 = block_matrix(no, nv);
  double ****Z2 = init_4d_array(no, no, nv, nv);
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      Z1[i][a] = 0.5 * moinfo.s1[i][a];
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++)
          Z2[i][j][a][b] = (1./3.)*moinfo.s2[i][j][a][b] + (1./6.)*moinfo.s2[i][j][b][a];
   }
  amp_write(20, Z1, Z2, "SZ"); fprintf(outfile, "\n");
  free_block(Z1);
  free_4d_array(Z2, no, no, nv);

  return ET;
}

}} // namespace psi::ugacc
