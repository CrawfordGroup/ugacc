#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void G_build(int iter)
{
  int no = moinfo.no;
  int nv = moinfo.nv;
  double ****t2 = moinfo.t2;
  double ****l2 = moinfo.l2old;

  if(iter == 1) {
    Gvv = moinfo.Gvv;
    Goo = moinfo.Goo;
  }

  for(int m=0; m < no; m++)
    for(int i=0; i < no; i++) {
      double value = 0.0;
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          for(int j=0; j < no; j++)
            value += t2[m][j][a][b] * l2[i][j][a][b];
      Goo[m][i] = value;
    }

  for(int a=0; a < nv; a++)
    for(int e=0; e < nv; e++) {
      double value = 0.0;
      for(int i=0; i < no; i++)
        for(int j=0; j < no; j++)
          for(int b=0; b < nv; b++)
            value -= t2[i][j][e][b] * l2[i][j][a][b];
      Gvv[a][e] = value;
    }

  return;
}

}} // namespace psi::ugacc
