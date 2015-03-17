#include <libciomr/libciomr.h>
#include <psi4-dec.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void F_build(void)
{
  int no = moinfo.no; 
  int nv = moinfo.nv;
  double **fock = moinfo.fock;
  double ****L = moinfo.L;
  double **t1 = moinfo.t1old;
  double ****ttau = moinfo.ttau;

  moinfo.Fae = block_matrix(nv,nv);
  for(int a=0; a < nv; a++)
    for(int e=0; e < nv; e++) {
      double value = fock[a+no][e+no]; 
      for(int m=0; m < no; m++) {
        value -= 0.5*fock[m][e+no]*t1[m][a];
        for(int f=0; f < nv; f++) {
          value += t1[m][f]*L[m][a+no][f+no][e+no];
          for(int n=0; n < no; n++) 
            value -= ttau[m][n][a][f]*L[m][n][e+no][f+no];
        }
      }
      moinfo.Fae[a][e] = value;
    }

  moinfo.Fmi = block_matrix(no,no);
  for(int m=0; m < no; m++)
    for(int i=0; i < no; i++) {
      double value = fock[m][i];
      for(int e=0; e < nv; e++) {
	value += 0.5*t1[i][e]*fock[m][e+no];
	for(int n=0; n < no; n++) {
	  value += t1[n][e]*L[m][n][i][e+no];
	  for(int f=0; f < nv; f++)
	    value += ttau[i][n][e][f]*L[m][n][e+no][f+no];
	}
      }
      moinfo.Fmi[m][i] = value;
    }

  moinfo.Fme = block_matrix(no,nv);
  for(int m=0; m < no; m++)
    for(int e=0; e < nv; e++) {
      double value = fock[m][e+no];
      for(int n=0; n < no; n++)
        for(int f=0; f < nv; f++)
          value += t1[n][f]*L[m][n][e+no][f+no];
      moinfo.Fme[m][e] = value;
    }
}

}} // namespace psi::ugacc
