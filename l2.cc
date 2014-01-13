#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void l2_build(void)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  double **l1 = moinfo.l1old;
  double ****l2 = moinfo.l2old;
  double ****l2new = moinfo.l2;
  double ****L = moinfo.L;
  double **Gvv = moinfo.Gvv;
  double **Goo = moinfo.Goo;
  double **Hov = moinfo.Hov;
  double **Hvv = moinfo.Hvv;
  double **Hoo = moinfo.Hoo;

  double ****Hovvo = moinfo.Hovvo;
  double ****Hovov = moinfo.Hovov;
  double ****Hoooo = moinfo.Hoooo;
  double ****Hvvvv = moinfo.Hvvvv;
  double ****Hvovv = moinfo.Hvovv;
  double ****Hooov = moinfo.Hooov;

  double ****Z = init_4d_array(no, no, nv, nv);

  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++) {
          double value = L[i][j][a+no][b+no];
          if(params.wfn == "CCSD_T") value += moinfo.s2[i][j][a][b];

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

}} // namespace psi::ugacc
