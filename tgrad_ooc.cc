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

void t3_ijk(double ***, int, int, int, double ****, double **, double ****);
void t3_abc(double ***, int, int, int, double ****, double **, double ****);
void l3_ijk(double ***, int, int, int, double ****, double **, double **, double ****, double ****);
void l3_ijk_new(double ***, int, int, int, double ****, double **, double **, double ****);
void l3_abc(double ***, int, int, int, double ****, double **, double **, double ****, double ****);

void tgrad_ooc(void)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  double **fock = moinfo.fock;
  double ****ints = moinfo.ints;
  double ****L = moinfo.L;
  double **t1 = moinfo.t1;
  double ****t2 = moinfo.t2;
  double **t1s = moinfo.t1s;
  double ****t2s = moinfo.t2s;


  double **X1 = block_matrix(no, nv); // T3 --> L1
  double ****X2 = init_4d_array(no, no, nv, nv); // T3 --> L2
  double ****Z2 = init_4d_array(no, no, nv, nv); // L3 --> L2

  double ***t3 = init_3d_array(nv, nv, nv);
  double ***l3 = init_3d_array(nv, nv, nv);
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int k=0; k < no; k++) {
        t3_ijk(t3, i, j, k, t2, fock, ints);
//        l3_ijk(l3, i, j, k, t2s, t1s, fock, L, ints);
        l3_ijk_new(l3, i, j, k, t2, t1, fock, ints);

        for(int a=0; a < nv; a++) 
          for(int b=0; b < nv; b++)
            for(int c=0; c < nv; c++) {

//              X1[i][a] += (4.0*t3[a][b][c] - 2.0*t3[c][b][a] - 2.0*t3[a][c][b] + t3[b][c][a])*ints[j][k][b+no][c+no];
              moinfo.s1[i][a] += 2.0*(t3[a][b][c] - t3[c][b][a]) * L[j][k][b+no][c+no];
              X2[i][j][a][b] += (t3[a][b][c] - t3[c][b][a]) * fock[k][c+no];
              moinfo.Goovv[i][j][a][b] += 2.0 * t1s[k][c] * (2.0*(t3[a][b][c] - t3[a][c][b]) - (t3[b][a][c] - t3[b][c][a]));

              for(int l=0; l < no; l++) {
                X2[i][l][a][b] -= (2.0*t3[a][b][c] - t3[a][c][b] - t3[c][b][a]) * ints[j][k][l][c+no];
                Z2[i][l][a][b] -= l3[a][b][c] * ints[j][k][l][c+no];
                moinfo.Gooov[j][i][l][a] -= (2.0*t3[a][b][c] - t3[b][a][c] - t3[c][b][a]) * t2s[l][k][b][c] + l3[a][b][c] * t2[l][k][b][c];
              }

              for(int d=0; d < nv; d++) {
                X2[i][j][a][d] += (2.0*t3[a][b][c] - t3[a][c][b] - t3[c][b][a]) * ints[d+no][k][b+no][c+no];
                Z2[i][j][a][d] += l3[a][b][c] * ints[d+no][k][b+no][c+no];
                moinfo.Dvv[a][b] += 0.5 * t3[b][c][d] * l3[a][c][d];
                moinfo.Gvvvo[a][b][d][j] += (2.0*t3[a][b][c] - t3[b][a][c] - t3[a][c][b]) * t2s[k][i][c][d] + l3[a][b][c] * t2[k][i][c][d];
              }

            } // abc

      } // ijk
  free_3d_array(t3, nv, nv);
  free_3d_array(l3, nv, nv);

  t3 = init_3d_array(no, no, no);
  l3 = init_3d_array(no, no, no);
  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int c=0; c < nv; c++) {
        t3_abc(t3, a, b, c, t2, fock, ints);
        l3_abc(l3, a, b, c, t2s, t1s, fock, L, ints);

        for(int i=0; i < no; i++)
          for(int j=0; j < no; j++)
            for(int k=0; k < no; k++) {
              for(int l=0; l < no; l++)
                moinfo.Doo[i][j] -= 0.5 * t3[i][k][l] * l3[j][k][l];

            } // ijk
      } // abc
  free_3d_array(t3, no, no);
  free_3d_array(l3, no, no);

//  for(int i=0; i < no; i++)
//    for(int a=0; a < nv; a++) 
//      moinfo.s1[i][a] = X1[i][a];
  double ****Y2 = init_4d_array(no, no, nv, nv);
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) 
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++)
          Y2[i][j][a][b] = X2[i][j][a][b] + X2[j][i][b][a];
  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) 
      for(int j=0; j < no; j++)
        for(int b=0; b < nv; b++) {
          moinfo.s2[i][j][a][b] = 4.0 * Y2[i][j][a][b] - 2.0 * Y2[i][j][b][a];
          moinfo.s2[i][j][a][b] += Z2[i][j][a][b] + Z2[j][i][b][a];
        }
  free_4d_array(Y2, no, no, nv);

  free_block(X1);
  free_4d_array(X2, no, no, nv);
  free_4d_array(Z2, no, no, nv);

  return;
}

}} // namespace psi::ugacc
