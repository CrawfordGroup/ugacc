#include <stdio.h>
#include <strings.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

double pseudoenergy(void)
{
  int no = moinfo.no; 
  int nv = moinfo.nv;
  double ****L = moinfo.L;
  double ****l2 = moinfo.l2;
  double energy=0.0;

  double energy=0.0;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          energy += L[i][j][a+no][b+no]*l2[i][j][a][b];

  return energy;
}

}} // namespace psi::ugacc

