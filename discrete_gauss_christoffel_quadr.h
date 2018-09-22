#ifndef DISCRETE_GAUSS_CHRISTOFFEL_QUADR_H
#define DISCRETE_GAUSS_CHRISTOFFEL_QUADR_H

#ifdef __cplusplus
extern "C" {            /* Assume C declarations for C++ */
#endif /* __cplusplus */



/* perform discrete gauss christoffel quadrature with final diagonalization. The output is a downsampled function at L points whose nodes are stored in nodes and the corresponding weights i weights

 *  N - length of arrays w, x
 *  w - array with function values
 *  x - array with the arguments at which function is evaluated
 *  L - number of downsampled points
 *
 *  nodes   - the downsampled nodes at which the function is evaluated, dimension L
 *  weights - the downsampled weights corresponding to the nodes, dimension L
 *
 */
void discrete_gauss_christoffel_quadr(const size_t N, const double *w, const double *x, const double total_weight, const size_t L, double * nodes, double *weights);



/* perform discrete gauss christoffel quadrature without final diagonalization. The output is the tridiagonal Jacobi matrix whose main diagonal elements are stored in alpha_vec_out and its secondary diagonal elements are stored in beta_vec_out

 *  N - length of arrays w, x
 *  w - array with function values
 *  x - array with the arguments at which function is evaluated
 *  L - number of downsampled points
 *
 *  alpha_vec_out   - the main diagonal elements of the output matrix, dimension L
 *  beta_vec_out    - the secondary diagonal elements of the output matrix, dimension L-1
 *
 */
void discrete_gauss_christoffel_quadr_pre_diag(const size_t N, const double *w, const double *x, const size_t L, double *alpha_vec_out, double *beta_vec_out);



/* Arnoldi algorithm to generate effective opeators and an orthogonal basis in the Krylov subspace {state, H*state, ..., H^(dim-1)*state} whose first basis vector is the normalized input state
 this algorithm is for diagonal, real-valued operators H

 *  N - length of diag(H)
 *  H - operator
 *  state0 - input/initial state
 *  dim - the dimension of the Krylov subspace
 *
 *  alpha_vec_out - main diagonal elements of the Krylov Hamiltonian
 *  beta_vec_out  - secondary diagonal elements of the Krylov Hamiltonian 

*/
void arnoldi_ddiag(const size_t N, const double *H, const double *state0, const size_t dim, double *alpha_vec_out, double *beta_vec_out, double *basis_subspace_out);



/* reorthogonalization of a vector w.r.t to a given set of already mutually orthonormal vectors
algorithm is taken from 'Matrix Algorithms - Volume 1' by G W Stewart and can be found in Chapter 4, algorithm 1.13
nput parameters:
This routine is for real-valued data in double precision
 
 * N_row    - number of rows of Q
 * N_col    - number of cols of Q
 * K_col    - number of accessed cols of Q
 * Q        - matrix storing a set of mutually orthonormal vectors in its columns
 * x        - vector which is to be made orthonormal to the columns of Q
 *
 * q_out    - vector which is orthogonal to Q up to machine precision
 * r_out    - vector containing overlap of x_orth with all column vectors of Q from all iteration steps
 * rho_out  - norm of x_orth after final iteration step 
 
*/
void gsreorthog_d(const size_t N_row, const size_t N_col, const size_t K_col, const double *Q, const double *x, double *q_out, double *r_out, double *rho_out);



/* function declarations called from BLAS */
extern double   dnrm2_(const int *N, const double *X, const int *incX);
extern double   dasum_(const int *N, const double *X, const int *incX);
extern void     dcopy_(const int *N, const double *X, const int *incX, double *Y, const int *incY);
extern void     daxpy_(const int *N, const double *alpha, const double *X, const int *incX, double *Y, const int *incY);
extern void     dscal_(const int *N, const double *alpha, double *X, const int *incX);

extern void     dgemv_(const char *TRANS, const int *M, const int *N,
                 const double *alpha, const double *A, const int *lda,
                 const double *X, const int *incX, const double *beta,
                 double *Y, const int *incY);
extern void     dgbmv_(const char *TRANS, const int *M, const int *N,
                 const int *KL, const int *KU, const double *alpha,
                 const double *A, const int *lda, const double *X,
                 const int *incX, const double *beta, double *Y, const int *incY);

/* function declarations called from LAPACK */
extern void dpteqr_(const char* compz, const int* n, double* d, double* e, double* z, const int* ldz, double* work, int *info );

#ifdef __cplusplus
}
#endif

#endif /* DISCRETE_GAUSS_CHRISTOFFEL_QUADR_H */
