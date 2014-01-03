#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void t1_build(void)
{
  int no = moinfo.no;  
  int nv = moinfo.nv;
  double **t1new = moinfo.t1;
  double **t1 = moinfo.t1old;
  double **fock = moinfo.fock;
  double ****ints = moinfo.ints;
  double ****L = moinfo.L;
  double **Fae = moinfo.Fae;
  double **Fmi = moinfo.Fmi;
  double **Fme = moinfo.Fme;
  double ****t2 = moinfo.t2old;
  double value;

  for(int i=0; i < no; i++)
    for(int a=0; a < nv; a++) {
      value = fock[a+no][i];
      for(int e=0; e < nv; e++)
        value += t1[i][e]*Fae[a][e];
      for(int m=0; m < no; m++)
        value -= t1[m][a]*Fmi[m][i];
      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          value += (2*t2[i][m][a][e]-t2[i][m][e][a])*Fme[m][e];
      for(int n=0; n < no; n++)
        for(int f=0; f < nv; f++)
          value += t1[n][f]*L[n][a+no][f+no][i];
      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          for(int f=0; f < nv; f++)
            value += (2*t2[m][i][e][f]-t2[m][i][f][e])*ints[m][a+no][e+no][f+no];
      for(int m=0; m < no; m++)
        for(int e=0; e < nv; e++)
          for(int n=0; n < no; n++)
            value -= t2[m][n][a][e]*L[n][m][e+no][i];
      t1new[i][a] = value;
   }

  return;
}

}} // namespace psi::ugacc
