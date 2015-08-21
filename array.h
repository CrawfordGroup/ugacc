#ifndef ARRAY_H
#define ARRAY_H

namespace psi {

double ***init_3d_array(int, int, int);
void free_3d_array(double ***, int, int);

double ****init_4d_array(int, int, int, int);
void free_4d_array(double ****, int, int, int);

double ******init_6d_array(int, int, int, int, int, int);
void free_6d_array(double ******, int, int, int, int, int);

} // namespace psi

#endif // ARRAY_H
