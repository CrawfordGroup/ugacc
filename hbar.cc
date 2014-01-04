#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void hbar(void)
{
  int no = moinfo.no;
  int nv = moinfo.nv;

  moinfo.Hvv = block_matrix(nv, nv);
  moinfo.Hoo = block_matrix(no, no);
  moinfo.Hov = block_matrix(no, nv);
  moinfo.Hoooo = init_4d_array(no, no, no, no);
  moinfo.Hvvvv = init_4d_array(nv, nv, nv, nv);
  moinfo.Hovov = init_4d_array(no, nv, no, nv);
  moinfo.Hovvo = init_4d_array(no, nv, nv, no);
  moinfo.Hvovv = init_4d_array(nv, no, nv, nv);
  moinfo.Hooov = init_4d_array(no, no, no, nv);
  moinfo.Hovoo = init_4d_array(no, nv, no, no);
  moinfo.Hvvvo = init_4d_array(nv, nv, nv, no);

  double **fock = moinfo.fock;
  double ****ints = moinfo.ints;
  double ****L = moinfo.L;
  double **t1 = moinfo.t1;
  double ****ttau = moinfo.ttau;
  double ****tau = moinfo.tau;
  double ****t2 = moinfo.t2;

  for(int a=0; a < nv; a++)
    for(int e=0; e < nv; e++) {
      double value = fock[a+no][e+no]; 
      for(int m=0; m < no; m++) {
        value -= 0.5*fock[m][e+no]*t1[m][a];
        value -= 0.5*Fme[m][e]*t1[m][a];
        for(int f=0; f < nv; f++) {
          value += t1[m][f]*L[m][a+no][f+no][e+no];
          for(int n=0; n < no; n++) 
            value -= ttau[m][n][a][f]*L[m][n][e+no][f+no];
        }
      }
      moinfo.Hvv[a][e] = value;
    }

  for(int m=0; m < no; m++)
    for(int i=0; i < no; i++) {
      double value = fock[m][i];
      for(int e=0; e < nv; e++) {
	value += 0.5*t1[i][e]*fock[m][e+no];
	value += 0.5*t1[i][e]*Fme[m][e];
	for(int n=0; n < no; n++) {
	  value += t1[n][e]*L[m][n][i][e+no];
	  for(int f=0; f < nv; f++)
	    value += ttau[i][n][e][f]*L[m][n][e+no][f+no];
	}
      }
      moinfo.Hoo[m][i] = value;
    }

  for(int m=0; m < no; m++)
    for(int e=0; e < nv; e++) {
      double value = fock[m][e+no];
      for(int n=0; n < no; n++)
        for(int f=0; f < nv; f++)
          value += t1[n][f]*L[m][n][e+no][f+no];
      moinfo.Hov[m][e] = value;
    }

  for(int m=0; m < no; m++)
    for(int n=0; n < no; n++)
      for(int i=0; i < no; i++)
	for(int j=0; j < no; j++) {
	  double value = ints[m][n][i][j];
	  for(int e=0; e < nv; e++) {
	    value += t1[j][e]*ints[m][n][i][e+no] +
                     t1[i][e]*ints[m][n][e+no][j];
	    for(int f=0; f < nv; f++)
	      value += tau[i][j][e][f]*ints[m][n][e+no][f+no];
	  }
	  moinfo.Hoooo[m][n][i][j] = value;
	}

  for(int m=0; m < no; m++)
    for(int b=0; b < nv; b++)
      for(int j=0; j < no; j++)
        for(int e=0; e < nv; e++) {
	  double value = -ints[m][b+no][j][e+no];
	  for(int f=0; f < nv; f++)
	    value -= t1[j][f]*ints[m][b+no][f+no][e+no];
	  for(int n=0; n < no; n++) 
	    value += t1[n][b]*ints[m][n][j][e+no];
	  for(int n=0; n < no; n++) {
	    for(int f=0; f < nv; f++)
	      value += ints[m][n][f+no][e+no]*
		(t2[j][n][f][b] + t1[j][f]*t1[n][b]);
	  }
	  moinfo.Hovov[m][b][j][e] = value;
	}

  for(int m=0; m < no; m++)
    for(int b=0; b < nv; b++)
      for(int e=0; e < nv; e++)
        for(int j=0; j < no; j++) {
          double value = ints[m][b+no][e+no][j];
          for(int f=0; f < nv; f++)
            value += t1[j][f]*ints[m][b+no][e+no][f+no];
          for(int n=0; n < no; n++)
            value -= t1[n][b]*ints[m][n][e+no][j];
          for(int n=0; n < no; n++)
            for(int f=0; f < nv; f++)
              value -= ints[m][n][e+no][f+no]*
                (t2[j][n][f][b] + t1[j][f]*t1[n][b]);
          for(int n=0; n < no; n++)
            for(int f=0; f < nv; f++)
              value += t2[n][j][f][b]*L[m][n][e+no][f+no];
          moinfo.Hovvo[m][b][e][j] = value;
        }

  for(int a=0; a < nv; a++)
    for(int m=0; m < no; m++)
      for(int e=0; e < nv; e++)
        for(int f=0; f < nv; f++) {
          double value = ints[a+no][m][e+no][f+no];
          for(int n=0; n < no; n++)
            value -= t1[n][a] * ints[n][m][e+no][f+no];
          moinfo.Hvovv[a][m][e][f] = value;
        }

  for(int m=0; m < no; m++)
    for(int n=0; n < no; n++)
      for(int i=0; i < no; i++)
        for(int e=0; e < nv; e++) {
          double value = ints[m][n][i][e+no];
          for(int f=0; f < nv; f++)
            value += t1[i][f] * ints[m][n][f+no][e+no];
          moinfo.Hooov[m][n][i][e] = value;
        }

  for(int m=0; m < no; m++)
    for(int b=0; b < nv; b++)
      for(int i=0; i < no; i++)
        for(int j=0; j < no; j++) {
          double value = ints[m][b+no][i][j];
          for(int e=0; e < nv; e++)
            value += Fme[m][e] * t2[i][j][e][b];
          for(int n=0; n < no; n++)
            value -= t1[n][b] * Wmnij[m][n][i][j];
          for(int e=0; e < nv; e++)
            for(int f=0; f < nv; f++)
              value += ints[m][b+no][e+no][f+no] * tau[i][j][e][f];
          for(int n=0; n < no; n++)
            for(int e=0; e < nv; e++)
              value += L[m][n][i][e+no] * t2[j][n][b][e] - ints[m][n][i][e+no] * t2[j][n][e][b] - ints[m][n][e+no][j] * t2[i][n][e][b];
          for(int e=0; e < nv; e++)
            value += t1[i][e] * ints[m][b+no][e+no][j]
                   + t1[j][e] * ints[m][b+no][i][e+no];
          for(int e=0; e < nv; e++)
            for(int n=0; n < no; n++)
              for(int f=0; f < nv; f++)
                value += t1[i][e] * ( t2[j][n][b][f] * L[m][n][e+no][f+no] - t2[j][n][f][b] * ints[m][n][e+no][f+no] ) - t1[j][e] * t2[n][i][b][f] * ints[m][n][f+no][e+no];

          moinfo.Hovoo[m][b][i][j] = value;
        }

  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int e=0; e < nv; e++)
        for(int i=0; i < no; i++) {
          /* I */
          double value = ints[a+no][b+no][e+no][i];
          /* II */
          for(int m=0; m < no; m++)
            value -= Fme[m][e] * t2[m][i][a][b];
          /* IIIa */
          for(int f=0; f < nv; f++)
            value += t1[i][f] * ints[a+no][b+no][e+no][f+no];
          /* IIIc + IIId + IV */
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              value += Wmnie[n][m][i][e] * tau[m][n][a][b];
          /* IIIb + V */
          for(int m=0; m < no; m++)
            for(int f=0; f < nv; f++)
              value += (2*t2[m][i][f][b] - t2[m][i][b][f] - t1[i][f] * t1[m][b]) * ints[m][a+no][f+no][e+no]  - t2[m][i][f][b] * ints[m][a+no][e+no][f+no] - tau[m][i][a][f] * ints[m][b+no][e+no][f+no];
          /* VI */
          for(int m=0; m < no; m++)
            value -= t1[m][a] * ints[m][b+no][e+no][i]
                   + t1[m][b] * ints[a+no][m][e+no][i];
          /* VII */
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              for(int f=0; f < nv; f++)
                value -= t1[m][a] * (t2[n][i][f][b] * L[m][n][e+no][f+no]
                                   - t2[n][i][b][f] * ints[m][n][e+no][f+no])
                      - t1[m][b] * t2[n][i][a][f] * ints[m][n][f+no][e+no];
          moinfo.Hvvvo[a][b][e][i] = value;
        }
}

}} // namespace psi::ugacc
