#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <stdlib.h>
#include "psi4-dec.h"
#include "libmints/mints.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

extern "C" int dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);

namespace psi { namespace ugacc {

/*!
 * Analyzes a set of singular values, printing out their sparsity pattern.
 * @param tensor_name - A title used for printing
 * @param S - The singular values
 * @param dim - The dimension of the singular values
 */
void analyze_svd(const char* tensor_name, double*S, int dim){
  int totcount = 0;
  int count2 = 0;
  int count3 = 0;
  int count4 = 0;
  int count5 = 0;
  fprintf(outfile, "\n%s singular values\n", tensor_name);
  for(int n = 0; n < dim; ++n){
    fprintf(outfile, " Singular values[%3d] = %16.10f\n", n, S[n]);
    if(fabs(S[n] > 1.0E-2)) count2++;
    if(fabs(S[n] > 1.0E-3)) count3++;
    if(fabs(S[n] > 1.0E-4)) count4++;
    if(fabs(S[n] > 1.0E-5)) count5++;
    totcount++;
  }
  fprintf(outfile, "\nSparsity report for %s\n", tensor_name);
  fprintf(outfile, "10^-2: %5d elements (%.2f\%)\n", count2, 100.0*(double)count2/(double)totcount);
  fprintf(outfile, "10^-3: %5d elements (%.2f\%)\n", count3, 100.0*(double)count3/(double)totcount);
  fprintf(outfile, "10^-4: %5d elements (%.2f\%)\n", count4, 100.0*(double)count4/(double)totcount);
  fprintf(outfile, "10^-5: %5d elements (%.2f\%)\n", count5, 100.0*(double)count5/(double)totcount);
  fprintf(outfile, "\n");
}


/*!
 * Transforms an integral tensor to the MO basis
 * @param Ints - On input, the AO basis integral.  On output, the MO basis integrals.
 */
void transform_integrals(double**Ints)
{
  /*
   * (pq|rs) = Ct Ct (mu nu | rho sigma) C C
   */
  boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
  double **pS = wfn->S()->pointer();
  double **pC = wfn->Ca()->pointer();
  int nmo = wfn->nmo();
  double**tmpPQRS = block_matrix(nmo*nmo, nmo*nmo);
  /*
   * 1st quarter transformation
   * (mu nu | rho s)  = (mu nu | rho sigma) C
   */
  for(int mu=0; mu < nmo; ++mu){
    for(int nu=0; nu < nmo; ++nu){
      for(int rho=0; rho < nmo; ++rho){
        for(int s=0; s < nmo; ++s){
          double val = 0.0;
          for(int sigma=0; sigma < nmo; ++sigma){
            val += Ints[mu*nmo + nu][rho*nmo + sigma] * pC[sigma][s];
          }
          tmpPQRS[mu*nmo+nu][rho*nmo+s] = val;
        }
      }
    }
  }
  /*
   * 2nd quarter transformation
   * (mu nu | r s)  = (mu nu | rho s) C
   */
  for(int mu=0; mu < nmo; ++mu){
    for(int nu=0; nu < nmo; ++nu){
      for(int r=0; r < nmo; ++r){
        for(int s=0; s < nmo; ++s){
          double val = 0.0;
          for(int rho=0; rho < nmo; ++rho){
            val += tmpPQRS[mu*nmo + nu][rho*nmo + s] * pC[rho][r];
          }
          Ints[mu*nmo + nu][r*nmo + s] = val;
        }
      }
    }
  }
  /*
   * 3rd quarter transformation
   * (mu q | r s)  = (mu nu | r s) C
   */
  for(int mu=0; mu < nmo; ++mu){
    for(int q=0; q < nmo; ++q){
      for(int r=0; r < nmo; ++r){
        for(int s=0; s < nmo; ++s){
          double val = 0.0;
          for(int nu=0; nu < nmo; ++nu){
            val += pC[nu][q] * Ints[mu*nmo + nu][r*nmo + s];
          }
          tmpPQRS[mu*nmo + q][r*nmo + s] = val;
        }
      }
    }
  }
  /*
   * 4th quarter transformation
   * (p q | r s)  = (mu q | r s) C
   */
  for(int p=0; p < nmo; ++p){
    for(int q=0; q < nmo; ++q){
      for(int r=0; r < nmo; ++r){
        for(int s=0; s < nmo; ++s){
          double val = 0.0;
          for(int mu=0; mu < nmo; ++mu){
            val += pC[mu][p] * tmpPQRS[mu*nmo + q][r*nmo + s];
          }
          Ints[p*nmo + q][r*nmo + s] = val;
        }
      }
    }
  }
  free_block(tmpPQRS);
}


/*!
 * Transforms an integral tensor back to the AO basis
 * @param Ints - On input, the MO basis integral.  On output, the AO basis integrals.
 */
void backtransform_integrals(double**Ints)
{
  /*
   * Ct S C = I
   *
   * => C^-1 = Ct S
   * => Ct^-1 = S C = (C^-1)t
   *
   * (pq|rs) = Ct Ct (mu nu | rho sigma) C C
   *
   * (mu nu | rho sigma) = Ct^-1 Ct^-1 (pq|rs) C^-1 C^-1
   */
  boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
  double **pS = wfn->S()->pointer();
  double **pC = wfn->Ca()->pointer();
  int nmo = wfn->nmo();
  double **CInv = block_matrix(nmo, nmo);
  for(int p=0; p < nmo; ++p){
    for(int q=0; q < nmo; ++q){
      double val = 0.0;
      for(int r=0; r < nmo; ++r){
        val  += pC[r][p] * pS[r][q];
      }
      CInv[p][q] = val;
    }
  }
  double**tmpPQRS = block_matrix(nmo*nmo, nmo*nmo);
  /*
   * 1st quarter transformation
   * (p q | r sigma) = (p q | r s) C^-1
   */
  for(int p=0; p < nmo; ++p){
    for(int q=0; q < nmo; ++q){
      for(int r=0; r < nmo; ++r){
        for(int sigma=0; sigma < nmo; ++sigma){
          double val = 0.0;
          for(int s=0; s < nmo; ++s){
            val += Ints[p*nmo + q][r*nmo + s] * CInv[s][sigma];
          }
          tmpPQRS[p*nmo+q][r*nmo+sigma] = val;
        }
      }
    }
  }
  /*
   * 2nd quarter transformation
   * (p q | rho sigma) = (p q | r sigma) C^-1
   */
  for(int p=0; p < nmo; ++p){
    for(int q=0; q < nmo; ++q){
      for(int rho=0; rho < nmo; ++rho){
        for(int sigma=0; sigma < nmo; ++sigma){
          double val = 0.0;
          for(int r=0; r < nmo; ++r){
            val += tmpPQRS[p*nmo + q][r*nmo + sigma] * CInv[r][rho];
          }
          Ints[p*nmo + q][rho*nmo + sigma] = val;
        }
      }
    }
  }
  /*
   * 3rd quarter transformation
   * (p nu | rho sigma) = Ct^-1 (p q | rho sigma)
   */
  for(int p=0; p < nmo; ++p){
    for(int nu=0; nu < nmo; ++nu){
      for(int rho=0; rho < nmo; ++rho){
        for(int sigma=0; sigma < nmo; ++sigma){
          double val = 0.0;
          for(int q=0; q < nmo; ++q){
            val += CInv[q][nu] * Ints[p*nmo + q][rho*nmo + sigma];
          }
          tmpPQRS[p*nmo + nu][rho*nmo + sigma] = val;
        }
      }
    }
  }
  /*
   * 4th quarter transformation
   * (mu nu | rho sigma) = Ct^-1 (p nu | rho sigma)
   */
  for(int mu=0; mu < nmo; ++mu){
    for(int nu=0; nu < nmo; ++nu){
      for(int rho=0; rho < nmo; ++rho){
        for(int sigma=0; sigma < nmo; ++sigma){
          double val = 0.0;
          for(int p=0; p < nmo; ++p){
            val += CInv[p][mu] * tmpPQRS[p*nmo + nu][rho*nmo + sigma];
          }
          Ints[mu*nmo + nu][rho*nmo + sigma] = val;
        }
      }
    }
  }
  free_block(tmpPQRS);
  free_block(CInv);
}




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
  double ****t2 = moinfo.t2;
  double ****tau = moinfo.tau;

  for(int m=0; m < no; m++)
    for(int e=0; e < nv; e++) {
      double value = fock[m][e+no];
      for(int n=0; n < no; n++)
        for(int f=0; f < nv; f++)
          value += t1[n][f]*L[m][n][e+no][f+no];
      moinfo.Hov[m][e] = value;
    }

  for(int m=0; m < no; m++)
    for(int i=0; i < no; i++) {
      double value = fock[m][i];
      for(int e=0; e < nv; e++) {
	value += t1[i][e]*fock[m][e+no];
	for(int n=0; n < no; n++) {
	  value += t1[n][e]*L[m][n][i][e+no];
	  for(int f=0; f < nv; f++)
	    value += tau[i][n][e][f]*L[m][n][e+no][f+no];
	}
      }
      moinfo.Hoo[m][i] = value;
    }

  for(int a=0; a < nv; a++)
    for(int e=0; e < nv; e++) {
      double value = fock[a+no][e+no]; 
      for(int m=0; m < no; m++) {
        value -= fock[m][e+no]*t1[m][a];
        for(int f=0; f < nv; f++) {
          value += t1[m][f]*L[a+no][m][e+no][f+no];
          for(int n=0; n < no; n++) 
            value -= tau[m][n][f][a]*L[m][n][f+no][e+no];
        }
      }
      moinfo.Hvv[a][e] = value;
    }

  for(int m=0; m < no; m++)
    for(int n=0; n < no; n++)
      for(int i=0; i < no; i++)
	for(int j=0; j < no; j++) {
	  double value = ints[m][n][i][j];
	  for(int e=0; e < nv; e++) {
	    value += 2.0*t1[j][e]*ints[m][n][i][e+no];
	    for(int f=0; f < nv; f++)
	      value += tau[i][j][e][f]*ints[m][n][e+no][f+no];
	  }
	  moinfo.Hoooo[m][n][i][j] = value;
	}

  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int e=0; e < nv; e++)
        for(int f=0; f < nv; f++) {
          double value = ints[a+no][b+no][e+no][f+no];
          for(int m=0; m < no; m++) {
            value -= 2.0*t1[m][b]*ints[a+no][m][e+no][f+no];
            for(int n=0; n < no; n++)
              value += tau[m][n][a][b]*ints[m][n][e+no][f+no];
          }
          moinfo.Hvvvv[a][b][e][f] = value;
        }
#define TEST_EVALS 0
#define SVD 1
#define PRINT 0

#if PRINT
    fprintf(outfile, "Wabcd\n");
    print_mat(Wvvvv, nv2, nv2, outfile);
#endif
#if TEST_EVALS
  double *work = new double[4*nv2];
  double *evalsR = new double[nv2];
  double *evalsI = new double[nv2];
  int info = C_DGEEV('n', 'n', nv*nv, Wvvvv[0], nv2, evalsR, evalsI, NULL, 1, NULL, 1, work, 4*nv2);
  fprintf(outfile, "Info: %d\n", info);
  for(int n = 0; n < nv2; ++n)
    fprintf(outfile, "Evals[%3d] = %16.10f + %16.10fi\n", n, evalsR[n], evalsI[n]);
#endif
#if SVD
  int nv2 = nv*nv;
  char jobu = 'A';
  char jobvt = 'A';
  int m = nv2;
  int n = nv2;
  int lda = nv2;
  double *S = new double[nv2];
  double **U = block_matrix(nv2, nv2);
  int ldu = nv2;
  double **Vt = block_matrix(nv2, nv2);
  int ldvt = nv2;
  int lwork = 8*nv2;
  double *work = new double[lwork];
  int info;
  double **Wvvvv = block_matrix(nv2, nv2);
  size_t totcount, count5, count4, count3, count2;
  // Chemists' notation W
  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int e=0; e < nv; e++)
        for(int f=0; f < nv; f++) {
          double value = ints[a+no][b+no][e+no][f+no];
          for(int m=0; m < no; m++) {
            value -= 2.0*t1[m][b]*ints[a+no][m][e+no][f+no];
            for(int n=0; n < no; n++)
              value += tau[m][n][a][b]*ints[m][n][e+no][f+no];
          }
          Wvvvv[a*nv+e][b*nv+f] = value;
        }
  dgesvd_(&jobu, &jobvt, &m, &n, Wvvvv[0], &lda, S, U[0], &ldu, Vt[0], &ldvt, work, &lwork, &info);
  analyze_svd("W (MO Basis, Mulliken notation)", S, nv2);
  //// Physicists' notation W
  //for(int a=0; a < nv; a++)
  //  for(int b=0; b < nv; b++)
  //    for(int e=0; e < nv; e++)
  //      for(int f=0; f < nv; f++) {
  //        double value = ints[a+no][b+no][e+no][f+no];
  //        for(int m=0; m < no; m++) {
  //          value -= 2.0*t1[m][b]*ints[a+no][m][e+no][f+no];
  //          for(int n=0; n < no; n++)
  //            value += tau[m][n][a][b]*ints[m][n][e+no][f+no];
  //        }
  //        Wvvvv[a*nv+b][e*nv+f] = value;
  //      }
  //dgesvd_(&jobu, &jobvt, &m, &n, Wvvvv[0], &lda, S, U[0], &ldu, Vt[0], &ldvt, work, &lwork, &info);
  //analyze_svd("W (MO Basis, Dirac notation)", S, nv2);
  // Chemists' notation H
  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int e=0; e < nv; e++)
        for(int f=0; f < nv; f++) {
          double value = ints[a+no][b+no][e+no][f+no];
          Wvvvv[a*nv+e][b*nv+f] = value;
        }
  dgesvd_(&jobu, &jobvt, &m, &n, Wvvvv[0], &lda, S, U[0], &ldu, Vt[0], &ldvt, work, &lwork, &info);
  analyze_svd("H (MO Basis, Mulliken notation)", S, nv2);
  // Physicists' notation H
  //for(int a=0; a < nv; a++)
  //  for(int b=0; b < nv; b++)
  //    for(int e=0; e < nv; e++)
  //      for(int f=0; f < nv; f++) {
  //        double value = ints[a+no][b+no][e+no][f+no];
  //        Wvvvv[a*nv+b][e*nv+f] = value;
  //      }
  //dgesvd_(&jobu, &jobvt, &m, &n, Wvvvv[0], &lda, S, U[0], &ldu, Vt[0], &ldvt, work, &lwork, &info);
  //analyze_svd("H (MO Basis, Dirac notation)", S, nv2);
  free_block(U);
  free_block(Vt);
  free_block(Wvvvv);
  delete [] S;
  delete [] work;

  /*
   * SVD of the AO basis quantities
   */
  int nmo = no+nv;
  int nmo2 = nmo*nmo;
  m = nmo2;
  n = nmo2;
  lda = nmo2;
  S = new double[nmo2];
  U = block_matrix(nmo2, nmo2);
  ldu = nmo2;
  Vt = block_matrix(nmo2, nmo2);
  ldvt = nmo2;
  lwork = 8*nmo2;
  work = new double[lwork];

  /*
   * Construct W in the MO basis.  For now, it's just Wabef
   */
  double **Ints = block_matrix(nmo2, nmo2);

  for(int a=0; a < nv; ++a)
    for(int b=0; b < nv; ++b)
      for(int e=0; e < nv; ++e)
        for(int f=0; f < nv; ++f)
          Ints[(a+no)*nmo + b+no][(e+no)*nmo + f+no] = moinfo.Hvvvv[a][e][b][f];
  backtransform_integrals(Ints);
#define TESTTRANS 0
#if TESTTRANS
  transform_integrals(Ints);
  for(int a=0; a < nv; ++a)
    for(int b=0; b < nv; ++b)
      for(int e=0; e < nv; ++e)
        for(int f=0; f < nv; ++f){
          double exact = moinfo.Hvvvv[a][e][b][f];
          double transformed = Ints[(a+no)*nmo + b+no][(e+no)*nmo + f+no];
          double diff = fabs(exact-transformed);
          fprintf(outfile, "Exact: %16.10f Trans: %16.10f Diff %16.10f\n", exact, transformed, diff);
        }
  exit(1);
#endif
  dgesvd_(&jobu, &jobvt, &m, &n, Ints[0], &lda, S, U[0], &ldu, Vt[0], &ldvt, work, &lwork, &info);
  analyze_svd("W (AO Basis, Mulliken notation)", S, nmo2);
  /*
   * Construct H in the MO basis.  For now, it's just Habef
   */
  for(int pq = 0; pq < nmo2; ++pq)
    for(int rs = 0; rs < nmo2; ++rs)
      Ints[pq][rs] = 0.0;

  for(int a=no; a < nv+no; ++a)
    for(int b=no; b < nv+no; ++b)
      for(int e=no; e < nv+no; ++e)
        for(int f=no; f < nv+no; ++f)
          Ints[a*nmo + b][e*nmo + f] = ints[a][e][b][f];
  backtransform_integrals(Ints);
  dgesvd_(&jobu, &jobvt, &m, &n, Ints[0], &lda, S, U[0], &ldu, Vt[0], &ldvt, work, &lwork, &info);
  analyze_svd("H (AO Basis, Mulliken notation)", S, nmo2);

#endif
  exit(1);
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
            value += t1[i][f] * ints[n][m][e+no][f+no];
          moinfo.Hooov[m][n][i][e] = value;
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
              value -= ints[m][n][e+no][f+no]*tau[n][j][b][f];
          for(int n=0; n < no; n++)
            for(int f=0; f < nv; f++)
              value += t2[n][j][f][b]*L[m][n][e+no][f+no];
          moinfo.Hovvo[m][b][e][j] = value;
        }

  for(int m=0; m < no; m++)
    for(int b=0; b < nv; b++)
      for(int j=0; j < no; j++)
        for(int e=0; e < nv; e++) {
	  double value = ints[m][b+no][j][e+no];
	  for(int f=0; f < nv; f++)
	    value += t1[j][f]*ints[b+no][m][e+no][f+no];
	  for(int n=0; n < no; n++) 
	    value -= t1[n][b]*ints[m][n][j][e+no];
	  for(int n=0; n < no; n++) {
	    for(int f=0; f < nv; f++)
	      value -= ints[n][m][e+no][f+no]*tau[j][n][f][b];
	  }
	  moinfo.Hovov[m][b][j][e] = value;
	}

  for(int a=0; a < nv; a++)
    for(int b=0; b < nv; b++)
      for(int e=0; e < nv; e++)
        for(int i=0; i < no; i++) {
          double value = ints[a+no][b+no][e+no][i];
          for(int f=0; f < nv; f++)
            value += t1[i][f]*ints[a+no][b+no][e+no][f+no];
          for(int m=0; m < no; m++)
            value -= t1[m][b]*ints[a+no][m][e+no][i] +
                     t1[m][a]*ints[b+no][m][i][e+no];
          for(int m=0; m < no; m++)
            for(int f=0; f < nv; f++)
              value -= tau[i][m][f][a]*ints[m][b+no][e+no][f+no] +
                       tau[i][m][f][b]*ints[a+no][m][e+no][f+no];
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              value += tau[m][n][a][b]*ints[m][n][e+no][i];
          for(int m=0; m < no; m++)
            value -= fock[m][e+no]*t2[m][i][a][b];
          for(int m=0; m < no; m++)
            for(int f=0; f < nv; f++)
              value += t2[m][i][f][b]*L[a+no][m][e+no][f+no];
          for(int m=0; m < no; m++)
            for(int n=0; n < no; n++)
              for(int f=0; f < nv; f++)
                value += ints[m][n][e+no][f+no]*(t1[i][f]*t2[m][n][a][b]+t1[m][a]*t2[i][n][f][b]+t1[n][b]*t2[m][i][a][f]);
          for(int m=0; m < no; m++)            
            for(int n=0; n < no; n++)
              for(int f=0; f < nv; f++)
                value -= (t1[m][f]*t2[n][i][a][b]+t1[n][a]*t2[m][i][f][b])*L[m][n][f+no][e+no];
          for(int m=0; m < no; m++)            
            for(int n=0; n < no; n++)
              for(int f=0; f < nv; f++)
                value += t1[i][f]*t1[m][a]*t1[n][b]*ints[m][n][e+no][f+no];

          moinfo.Hvvvo[a][b][e][i] = value;
        }

  for(int m=0; m < no; m++)
    for(int b=0; b < nv; b++)
      for(int i=0; i < no; i++)
        for(int j=0; j < no; j++) {
          double value = ints[m][b+no][i][j];
          for(int e=0; e < nv; e++)
            value += t1[j][e]*ints[m][b+no][i][e+no] +
                     t1[i][e]*ints[b+no][m][j][e+no];
          for(int n=0; n < no; n++)
            value -= t1[n][b]*ints[m][n][i][j];
          for(int n=0; n < no; n++)
            for(int e=0; e < nv; e++)
              value -= tau[i][n][e][b]*ints[n][m][j][e+no] +
                       tau[j][n][e][b]*ints[m][n][i][e+no];
          for(int e=0; e < nv; e++)
            for(int f=0; f < nv; f++)
              value += ints[m][b+no][e+no][f+no]*tau[i][j][e][f];
          for(int e=0; e < nv; e++)
            value += fock[m][e+no]*t2[i][j][e][b];
          for(int n=0; n < no; n++)
            for(int e=0; e < nv; e++)
              value += L[m][n][i][e+no]*t2[n][j][e][b];
          for(int e=0; e < nv; e++)
            for(int n=0; n < no; n++)
              for(int f=0; f < nv; f++)
                value -= ints[m][n][e+no][f+no]*(t1[j][f]*t2[i][n][e][b]+t1[i][e]*t2[j][n][f][b]+t1[n][b]*t2[i][j][e][f]);
          for(int e=0; e < nv; e++)
            for(int n=0; n < no; n++)
              for(int f=0; f < nv; f++)
                value += (t1[i][e]*t2[n][j][f][b]+t1[n][f]*t2[i][j][e][b])*L[m][n][e+no][f+no];
          for(int e=0; e < nv; e++)
            for(int n=0; n < no; n++)
              for(int f=0; f < nv; f++)
                value -= t1[j][f]*t1[i][e]*t1[n][b]*ints[m][n][e+no][f+no];

          moinfo.Hovoo[m][b][i][j] = value;
        }
}

}} // namespace psi::ugacc
