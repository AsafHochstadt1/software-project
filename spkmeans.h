#include <math.h>
#ifndef SPKMEANS_H
#define SPKMEANS_H
void error_message();
double calc_distance(double *v1 , double *v2, int dim);
double** build_adjacency_matrix(double **vectors, int n, int dim);
double** build_dd_matrix(int n, double **am);
void build_lm_matrix(int n, double **am, double **ddm);
int* find_maximal_entry(double **lm, int n);
double** build_p_matrix(double **lm,int n, int *pos);
void mult_ptap(double **am, double **pm,int n, int *pos);
double** matrix_mult(double **p1, double **p2, int n);
double of_f(double **m, int n);
double* jacobi_algo(double **lm, int n);
double **copy_matrix(double **sm, int n);
void free_2d_double_array(double **arr, int len);
#endif