#include <libciomr/libciomr.h>
#include <cstring>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void init_density(void)
{
  int no = moinfo.no; 
  int nv = moinfo.nv;

  moinfo.Doo = block_matrix(no, no);
  moinfo.Dvv = block_matrix(nv, nv);
  moinfo.Dov = block_matrix(no, nv);
  moinfo.Dvo = block_matrix(nv, no);
  moinfo.Goooo = init_4d_array(no, no, no, no);
  moinfo.Gvvvv = init_4d_array(nv, nv, nv, nv);
  moinfo.Gooov = init_4d_array(no, no, no, nv);
  moinfo.Gvvvo = init_4d_array(nv, nv, nv, no);
  moinfo.Govov = init_4d_array(no, nv, no, nv);
  moinfo.Goovv = init_4d_array(no, no, nv, nv);

  return;
}

}} // namespace psi::ugacc
