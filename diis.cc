#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ugacc {

void diis(int error_file, int amp_file, int iter, double **t1, 
          double **t1old, double ****t2, double ****t2old)
{
  int nvector=8;  /* Number of error vectors to keep */
  int no,nv,word;
  int p,q,i,j,a,b;
  int diis_cycle;
  int vector_length;
  double **t1, **t1old;
  double ****t2, ****t2old;
  double *error;
  div_t fraction;
  double **B, *C, **vector;
  double product, determinant, maximum;
  psio_address start, end;

  no = moinfo.no; nv = moinfo.nv;

  /* Calculate the length of a single error vector */
  vector_length = no*nv + no*no*nv*nv;

  /* If we haven't already, open the vector files for reading/writing */
  if(iter == 1) { 
    psio_open(error_file, PSIO_OPEN_NEW);
    psio_open(amp_file, PSIO_OPEN_NEW);
  }

  /* Set the diis cycle value */
  fraction = div((iter-1),nvector);
  diis_cycle = fraction.rem;

  /* Build the current error vector and dump it to disk */
  error = init_array(vector_length);
  word=0;
  for(i=0; i < no; i++)
    for(a=0; a < nv; a++) {
      error[word++] = t1[i][a] - t1old[i][a];
    }

  for(i=0; i < no; i++)
    for(j=0; j < no; j++)
      for(a=0; a < nv; a++)
        for(b=0; b < nv; b++) {
          error[word++] = t2[i][j][a][b] - t2old[i][j][a][b];
        }

  start = psio_get_address(PSIO_ZERO, diis_cycle*vector_length*sizeof(double));
  psio_write(error_file, "DIIS Error Vectors", (char *) error, 
	     vector_length*sizeof(double), start, &end);

  /* Store the amplitudes, too */
  word=0;
  for(i=0; i < no; i++)
    for(a=0; a < nv; a++)
      error[word++] = t1[i][a];

  for(i=0; i < no; i++)
    for(j=0; j < no; j++)
      for(a=0; a < nv; a++)
        for(b=0; b < nv; b++)  {
          error[word++] = t2[i][j][a][b];
        }

  start = psio_get_address(PSIO_ZERO, diis_cycle*vector_length*sizeof(double));
  psio_write(amp_file, "DIIS Amplitude Vectors", (char *) error, 
	     vector_length*sizeof(double), start, &end);
  
  free(error);
    
  /* If we haven't run through enough iterations, set the correct dimensions
     for the extrapolation */
  if(!(iter >= (nvector))) {
    if(iter < 2) return; /* Leave if we can't extrapolate at all */
    nvector = iter;
  }

  /* Now grab the full set of error vectors from the file */
  vector = init_matrix(nvector, vector_length);
  for(p=0; p < nvector; p++) {
    start = psio_get_address(PSIO_ZERO, p*vector_length*sizeof(double));
    psio_read(error_file, "DIIS Error Vectors", (char *) vector[p], 
	      vector_length*sizeof(double), start, &end);
  }

  /* Build B matrix of error vector products */
  B = init_matrix(nvector+1,nvector+1);

  for(p=0; p < nvector; p++)
    for(q=0; q < nvector; q++) {
      dot_arr(vector[p], vector[q], vector_length, &product); 
      B[p][q] = product;
    }

  for(p=0; p < nvector; p++) {
    B[p][nvector] = -1;
    B[nvector][p] = -1;
  }

  B[nvector][nvector] = 0;

  /* Find the maximum value in B and scale all its elements */
  maximum = fabs(B[0][0]);
  for(p=0; p < nvector; p++)
    for(q=0; q < nvector; q++)
      if(fabs(B[p][q]) > maximum) maximum = fabs(B[p][q]);

  for(p=0; p < nvector; p++)
    for(q=0; q < nvector; q++)
      B[p][q] /= maximum; 

  /* Build the constant vector */
  C = init_array(nvector+1);
  C[nvector] = -1;

  /* Solve the linear equations */
  flin(B, C, nvector+1, 1, &determinant);

  /* Grab the old amplitude vectors */
  for(p=0; p < nvector; p++) {
    start = psio_get_address(PSIO_ZERO, p*vector_length*sizeof(double));
    psio_read(amp_file, "DIIS Amplitude Vectors", (char *) vector[p], 
	      vector_length*sizeof(double), start, &end);
  }
  
  /* Build the new amplitude vector from the old ones */
  word=0;
  for(i=0; i < no; i++)
    for(a=0; a < nv; a++) {
      t1[i][a] = 0.0;
      for(p=0; p < nvector; p++)
        t1[i][a] += C[p]*vector[p][word];
        word++;
    }

  for(i=0; i < no; i++)
    for(j=0; j < no; j++)
      for(a=0; a < nv; a++)
        for(b=0; b < nv; b++) {
          t2[i][j][a][b] = 0.0;
          for(p=0; p < nvector; p++)
            t2[i][j][a][b] += C[p]*vector[p][word];
            word++;
        }

  /* Release memory and return */
  free_matrix(vector, nvector);
  free_matrix(B, nvector+1);
  free(C);

  return;
}
}} // namespace devel::rhfccenergy

