
#include <R.h>
#include <math.h>

double *create_real_vector(int size);

void destroy_real_vector(double *vector, int size);

double **create_real_matrix(int rows, int cols);

double **create_and_copy_real_matrix(int rows, int cols, double *data);

void real_matrix_copy(double **mat, int d1, int d2, double *data);

void destroy_real_matrix(double **matrix, int rows, int cols);

double **create_real_array3d(int d1, int d2, int d3);

void real_array3d_copy(double ***array, int d1, int d2, int d3, double *data);

void destroy_real_array3d(double ***array, int d1, int d2, int d3);

void chol(double **L, double **M, int n); 

void chol_up(double **L_new, double **L_old, double *x, int n, double a, double b, double *work);

void fwdinv(double **L_inv, double **L, int n);

void matvec(double *y, double **A, double *x, int n);
