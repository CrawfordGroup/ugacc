#include <stdio.h>
#include <stdlib.h>

namespace psi {

double ***init_3d_array(int p, int q, int r)
{
  double ***A;
  int i,j,k;

  A = (double ***) malloc(p * sizeof(double **));
  for(i=0; i < p; i++) {
    A[i] = (double **) malloc(q * sizeof(double *));
    for(j=0; j < q; j++) {
      A[i][j] = (double *) malloc(r * sizeof(double));
      for(k=0; k < r; k++) {
	A[i][j][k] = 0.0;
      }
    }
  }

  return A;
}

void free_3d_array(double ***A, int p, int q)
{
  int i,j;

  for(i=0; i < p; i++)
    for(j=0; j < q; j++)
      free(A[i][j]);

  for(i=0; i < p; i++) free(A[i]);

  free(A);
}

double ****init_4d_array(int p, int q, int r, int s)
{
  double ****A;
  int i,j,k,l;

  A = (double ****) malloc(p * sizeof(double ***));
  for(i=0; i < p; i++) {
      A[i] = (double ***) malloc(q * sizeof(double **));
      for(j=0; j < q; j++) {
          A[i][j] = (double **) malloc(r * sizeof(double*));
          for(k=0; k < r; k++) {
              A[i][j][k] = (double *) malloc(s * sizeof(double));
              for(l=0; l < s; l++) {
                  A[i][j][k][l] = 0.0;
                }
            }
        }
    }

  return A;
}

void free_4d_array(double ****A, int p, int q, int r)
{
  int i,j,k;

  for(i=0; i < p; i++) {
      for(j=0; j < q; j++) {
          for(k=0; k < r; k++) {
              free(A[i][j][k]);
            }
        }
    }

  for(i=0; i < p; i++) {
      for(j=0; j < q; j++) {
          free(A[i][j]);
        }
    }

  for(i=0; i < p; i++) free(A[i]);
  free(A);

}

double ******init_6d_array(int p, int q, int r, int s, int t, int u)
{
  double ******A;
  int i,j,k,l,m,n;

  A = (double ******) malloc(p * sizeof(double *****));
  for(i=0; i < p; i++) {
    A[i] = (double *****) malloc(q * sizeof(double ****));
    for(j=0; j < q; j++) {
      A[i][j] = (double ****) malloc(r * sizeof(double ***));
      for(k=0; k < r; k++) {
	A[i][j][k] = (double ***) malloc(s * sizeof(double **));
	for(l=0; l < s; l++) {
	  A[i][j][k][l] = (double **) malloc(t * sizeof(double *));
	  for(m=0; m < t; m++) {
	    A[i][j][k][l][m] = (double *) malloc(u * sizeof(double));
	    for(n=0; n < u; n++) {
	      A[i][j][k][l][m][n] = 0.0;
	    }
	  }
	}
      }
    }
  }

  return A;
}

void free_6d_array(double ******A, int p, int q, int r, int s, int t)
{
  int i,j,k,l,m;

  for(i=0; i < p; i++) {
    for(j=0; j < q; j++) {
      for(k=0; k < r; k++) {
	for(l=0; l < s; l++) {
	  for(m=0; m < t; m++) {
	    free(A[i][j][k][l][m]);
	  }
	}
      }
    }
  }

  for(i=0; i < p; i++) {
    for(j=0; j < q; j++) {
      for(k=0; k < r; k++) {
	for(l=0; l < s; l++) {
	  free(A[i][j][k][l]);
	}
      }
    }
  }

  for(i=0; i < p; i++) {
    for(j=0; j < q; j++) {
      for(k=0; k < r; k++) {
	free(A[i][j][k]);
      }
    }
  }

  for(i=0; i < p; i++) {
    for(j=0; j < q; j++) {
      free(A[i][j]);
    }
  }

  for(i=0; i < p; i++) free(A[i]);
  free(A);
}

} // namespace psi

