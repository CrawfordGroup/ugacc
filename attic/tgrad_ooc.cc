/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */
 
/*
  tgrad_ooc(): Computes triples contributes to the CCSD(T) gradient using a
  triples-driven algorithm.  Both T3 and L3 amplitudes are computed in VVV
  batches for a given combination of OOO indices.

  -TDC, 1/2014
*/

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

void M3_ijk(double ***, int, int, int, double ****, double **, double ****);
void N3_ijk(double ***, int, int, int, double ****, double **, double **, double ****);
void M3_abc(double ***, int, int, int, double ****, double **, double ****);
void N3_abc(double ***, int, int, int, double ****, double **, double **, double ****);

void tgrad_ooc(void)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  double **fock = moinfo.fock;
  double ****ints = moinfo.ints;
  double **t1 = moinfo.t1;
  double ****t2 = moinfo.t2;

  double **X1 = block_matrix(no, nv); // T3 --> L1
  double ****X2 = init_4d_array(no, no, nv, nv); // T3 & L3 --> L2

  double ***M3 = init_3d_array(nv, nv, nv);
  double ***N3 = init_3d_array(nv, nv, nv);
  double ***X3 = init_3d_array(nv, nv, nv);
  double ***Y3 = init_3d_array(nv, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++) {
        M3_ijk(M3, i, j, k, t2, fock, ints);
        N3_ijk(N3, i, j, k, t2, t1, fock, ints);

        for(int a=0; a < nv; a++) 
          for(int b=0; b < nv; b++)
            for(int c=0; c < nv; c++) {
              X3[a][b][c] = 8.0*M3[a][b][c]-4.0*M3[b][a][c]-4.0*M3[a][c][b]-4.0*M3[c][b][a]+2.0*M3[c][a][b]+2.0*M3[b][c][a];
              Y3[a][b][c] = 8.0*N3[a][b][c]-4.0*N3[b][a][c]-4.0*N3[a][c][b]-4.0*N3[c][b][a]+2.0*N3[c][a][b]+2.0*N3[b][c][a];
            }

        for(int a=0; a < nv; a++) 
          for(int b=0; b < nv; b++)
            for(int c=0; c < nv; c++) {

              moinfo.Dvv[a][a] += 0.5 * M3[a][b][c] * (X3[a][b][c] + Y3[a][b][c]);
              moinfo.s1[i][a] += (4.0*M3[a][b][c] - 2.0*M3[c][b][a] - 2.0*M3[a][c][b] + M3[b][c][a]) * ints[j][k][b+no][c+no];
              moinfo.Goovv[i][j][a][b] += 4.0 * t1[k][c] * (2.0*(M3[a][b][c] - M3[a][c][b]) - (M3[b][a][c] - M3[b][c][a]));

              for(int l=0; l < no; l++) {
                X2[i][l][a][b] -= (2.0 * X3[a][b][c] + Y3[a][b][b]) * ints[j][k][l][c+no];
                moinfo.Gooov[j][i][l][a] -= (2.0 * X3[a][b][c] + Y3[a][b][c]) * t2[l][k][b][c];
              }

              for(int d=0; d < nv; d++) {
                X2[i][j][a][d] += (2.0 * X3[a][b][c] + Y3[a][b][c]) * ints[d+no][k][b+no][c+no];
                moinfo.Gvvvo[a][b][d][j] += (2.0 * X3[a][b][c] + Y3[a][b][c]) * t2[k][i][c][d];
              }

            } // abc

      } // ijk
  free_3d_array(M3, nv, nv);
  free_3d_array(N3, nv, nv);
  free_3d_array(X3, nv, nv);
  free_3d_array(Y3, nv, nv);

  M3 = init_3d_array(no, no, no);
  N3 = init_3d_array(no, no, no);
  X3 = init_3d_array(no, no, no);
  Y3 = init_3d_array(no, no, no);
  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int c=0; c < nv; c++) {
        M3_abc(M3, a, b, c, t2, fock, ints);
        N3_abc(N3, a, b, c, t2, t1, fock, ints);

        for(int i=0; i < no; i++)
          for(int j=0; j < no; j++)
            for(int k=0; k < no; k++) {
              X3[i][j][k] = 8.0*M3[i][j][k]-4.0*M3[j][i][k]-4.0*M3[i][k][j]-4.0*M3[k][j][i]+2.0*M3[k][i][j]+2.0*M3[j][k][i];
              Y3[i][j][k] = 8.0*N3[i][j][k]-4.0*N3[j][i][k]-4.0*N3[i][k][j]-4.0*N3[k][j][i]+2.0*N3[k][i][j]+2.0*N3[j][k][i];
            }

        for(int i=0; i < no; i++)
          for(int j=0; j < no; j++)
            for(int k=0; k < no; k++)
                moinfo.Doo[i][i] -= 0.5 * M3[i][j][k] * (X3[i][j][k] + Y3[i][j][k]);
           
      } // abc
  free_3d_array(M3, no, no);
  free_3d_array(N3, no, no);
  free_3d_array(X3, no, no);
  free_3d_array(Y3, no, no);

  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) 
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++) {
          moinfo.s2[i][j][a][b] = X2[i][j][a][b] + X2[j][i][b][a];
        }

  free_block(X1);
  free_4d_array(X2, no, no, nv);

  return;
}

}} // namespace psi::ugacc
