#include <stdio.h>
#include <strings.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

/*
** pseudoenergy(): Evaluates an energy-like expression for the Lambda doubles
** amplitudes: 
**   E = <0|L2 H|0> = 1/2 <ab|ij> L(ij,ab)
** This expression is derived in the UGA formalism.
*/

double pseudoenergy(void)
{
  int no = moinfo.no; 
  int nv = moinfo.nv;
  double ****ints = moinfo.ints;
  double ****l2 = moinfo.l2;

  double energy=0.0;
  for(int i=0; i < no; i++)
    for(int j=0; j < no; j++)
      for(int a=0; a < nv; a++)
        for(int b=0; b < nv; b++)
          energy += 0.5*ints[i][j][a+no][b+no]*l2[i][j][a][b];

  return energy;
}

}} // namespace psi::ugacc

