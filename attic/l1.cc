#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void l1_build(void)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  double **l1 = moinfo.l1old;
  double ****l2 = moinfo.l2old;
  double **l1new = moinfo.l1;
  double **Gvv = moinfo.Gvv;
  double **Goo = moinfo.Goo;
  double **Hov = moinfo.Hov;
  double **Hvv = moinfo.Hvv;
  double **Hoo = moinfo.Hoo;
  double ****Hvvvo = moinfo.Hvvvo;
  double ****Hovoo = moinfo.Hovoo;
  double ****Hovvo = moinfo.Hovvo;
  double ****Hovov = moinfo.Hovov;
  double ****Hvovv = moinfo.Hvovv;
  double ****Hooov = moinfo.Hooov;

  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      double value = 2 * Hov[i][a];
      if(params.wfn == "CCSD_T") value += moinfo.s1[i][a];

      for(int e=0; e < nv; e++)
        value += l1[i][e] * Hvv[e][a];

      for(int m=0; m < no; m++)
        value -= l1[m][a] * Hoo[i][m];

      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          value += l1[m][e] * (2*Hovvo[i][e][a][m] - Hovov[i][e][m][a]);

      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          for(int f=0; f < nv; f++)
            value += l2[i][m][e][f] * Hvvvo[e][f][a][m];

      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          for(int e=0; e < nv; e++)
            value -= l2[m][n][a][e] * Hovoo[i][e][m][n];

      for(int e=0; e < nv; e++)
        for(int f=0; f < nv; f++)
          value -= Gvv[e][f] * (2*Hvovv[e][i][f][a] - Hvovv[e][i][a][f]);

      for(int m=0; m < no; m++)
        for(int n=0; n < no; n++)
          value -= Goo[m][n] * (2*Hooov[m][i][n][a] - Hooov[i][m][n][a]);

      l1new[i][a] = value;
    }
}

}} // namespace psi::ugacc
