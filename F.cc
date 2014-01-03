#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void Fae_build(void)
{
  int no = moinfo.no; 
  int nv = moinfo.nv;
  double **fock = moinfo.fock;
  double ****L = moinfo.L;
  double **t1 = moinfo.t1old;
  double ****ttau = moinfo.ttau;
  double **Fae = moinfo.Fae;
  double value;

  for(int a=0; a < nv; a++)
    for(int e=0; e < nv; e++) {
      value = fock[a+no][e+no]; 
      for(int m=0; m < no; m++) {
        value -= 0.5*fock[m][e+no]*t1[m][a];
        for(int f=0; f < nv; f++) {
          value += t1[m][f]*L[m][a+no][f+no][e+no];
          for(int n=0; n < no; n++) 
            value -= ttau[m][n][a][f]*L[m][n][e+no][f+no];
        }
      }
      Fae[a][e] = value;
    }

/*
    fprintf(outfile, "\n\tFae Intermediate:\n");
    fprintf(outfile,   "\t-----------------\n");
    print_mat(Fae, nv, nv, outfile);
*/
}

void Fmi_build(void)
{
  int no = moinfo.no; 
  int nv = moinfo.nv;
  double **t1 = moinfo.t1old;
  double ****ttau = moinfo.ttau;
  double **fock = moinfo.fock;
  double ****L = moinfo.L;
  double **Fmi = moinfo.Fmi;
  double value;

  for(int m=0; m < no; m++)
    for(int i=0; i < no; i++) {
      value = fock[m][i];
      for(int e=0; e < nv; e++) {
	value += 0.5*t1[i][e]*fock[m][e+no];
	for(int n=0; n < no; n++) {
	  value += t1[n][e]*L[m][n][i][e+no];
	  for(int f=0; f < nv; f++)
	    value += ttau[i][n][e][f]*L[m][n][e+no][f+no];
	}
      }
      Fmi[m][i] = value;
    }
/*
    fprintf(outfile, "\n\tFmi Intermediate:\n");
    fprintf(outfile,   "\t-----------------\n");
    print_mat(Fmi, no, no, outfile);
*/
}

void Fme_build(void)
{
  int no = moinfo.no; 
  int nv = moinfo.nv;
  double **fock = moinfo.fock;
  double **t1 = moinfo.t1old;
  double ****L = moinfo.L;
  double **Fme = moinfo.Fme; 
  double value;

  for(int m=0; m < no; m++)
    for(int e=0; e < nv; e++) {
      value = fock[m][e+no];
      for(int n=0; n < no; n++)
        for(int f=0; f < nv; f++)
          value += t1[n][f]*L[m][n][e+no][f+no];
      Fme[m][e] = value;
    }
/*
  fprintf(outfile, "\n\tFme Intermediate:\n");
  fprintf(outfile,   "\t-----------------\n");
  print_mat(Fme, no, nv, outfile);
*/
}

}} // namespace psi::ugacc
