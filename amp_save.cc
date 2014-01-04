#include "MOInfo.h"
#define EXTERN
#include "globals.h"

/* 
** amp_save(): Save the current amplitudes for the next iteration.
*/
namespace psi { namespace ugacc {

void amp_save(double ****t1, double ****t1old, double ****t2, double ****t2old)
{
  double ****t2tmp = t2;
  t2 = t2old;
  t2old = t2tmp;

  double **t1tmp = t1;
  t1 = t1old;
  t1old = t1tmp;
}

}} // namespace psi::ugacc

