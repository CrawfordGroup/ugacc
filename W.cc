#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void Wmnij_build(void)
{
  int no = moinfo.no; 
  int nv = moinfo.nv;
  double ****tau = moinfo.tau;
  double **t1 = moinfo.t1old;
  double ****ints = moinfo.ints;
  double ****Wmnij = moinfo.Wmnij;
  double value;

  for(int m=0; m < no; m++)
    for(int n=0; n < no; n++)
      for(int i=0; i < no; i++)
	for(int j=0; j < no; j++) {
	  value = ints[m][n][i][j];
	  for(int e=0; e < nv; e++) {
	    value += t1[j][e]*ints[m][n][i][e+no] +
                     t1[i][e]*ints[m][n][e+no][j];
	    for(int f=0; f < nv; f++)
	      value += tau[i][j][e][f]*ints[m][n][e+no][f+no];
	  }
	  Wmnij[m][n][i][j] = value;
	}
}

void Wmbje_build(void)
{
  int no = moinfo.no; 
  int nv = moinfo.nv;
  double ****t2 = moinfo.t2old;
  double **t1 = moinfo.t1old;
  double ****ints = moinfo.ints;
  double ****Wmbje = moinfo.Wmbje;
  double value;

  for(int m=0; m < no; m++)
    for(int b=0; b < nv; b++)
      for(int j=0; j < no; j++)
        for(int e=0; e < nv; e++) {
	  value = -ints[m][b+no][j][e+no];
	  for(int f=0; f < nv; f++)
	    value -= t1[j][f]*ints[m][b+no][f+no][e+no];
	  for(int n=0; n < no; n++) 
	    value += t1[n][b]*ints[m][n][j][e+no];
	  for(int n=0; n < no; n++) {
	    for(int f=0; f < nv; f++)
	      value += ints[m][n][f+no][e+no]*
		(0.5*t2[j][n][f][b] + t1[j][f]*t1[n][b]);
	  }
	  Wmbje[m][b][j][e] = value;
	}
}

void Wmbej_build(void)
{
  int no = moinfo.no; 
  int nv = moinfo.nv;
  double ****t2 = moinfo.t2old;
  double **t1 = moinfo.t1old;
  double ****ints = moinfo.ints;
  double ****L = moinfo.L;
  double ****Wmbej = moinfo.Wmbej;
  double value;

  for(int m=0; m < no; m++)
    for(int b=0; b < nv; b++)
      for(int e=0; e < nv; e++)
        for(int j=0; j < no; j++) {
          value = ints[m][b+no][e+no][j];
          for(int f=0; f < nv; f++)
            value += t1[j][f]*ints[m][b+no][e+no][f+no];
          for(int n=0; n < no; n++)
            value -= t1[n][b]*ints[m][n][e+no][j];
          for(int n=0; n < no; n++)
            for(int f=0; f < nv; f++)
              value -= ints[m][n][e+no][f+no]*
                (0.5*t2[j][n][f][b] + t1[j][f]*t1[n][b]);
          for(int n=0; n < no; n++) {
            for(int f=0; f < nv; f++)
              value += 0.5*L[m][n][e+no][f+no]*t2[n][j][f][b];
          }
          Wmbej[m][b][e][j] = value;
        }
}

}} // namespace psi::ugacc
