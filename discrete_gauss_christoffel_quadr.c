#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "discrete_gauss_christoffel_quadr.h"



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
void discrete_gauss_christoffel_quadr(const size_t N, const double *w, const double *x, const double total_weight, const size_t L, double * nodes, double *weights){

    const int i_one = 1, i_zero = 0;
    const double d_zero = 0., d_one = 1.;
    const char opt = 'N';

    // scale function and arguments
    double *w_scal = (double *) malloc(N*sizeof(double));
    double *x_scal = (double *) malloc(N*sizeof(double));

    double min_x; double max_x;

    for(int k=0;k<N;k++){
        w_scal[k] = sqrt(w[k]);
        if(k==0){
            min_x = x[k];
            max_x = x[k];
        }
        else{
            if(x[k] < min_x)
                min_x = x[k];
            if(x[k] > max_x)
                max_x = x[k];
        }
    }

    double alpha = 0.9/(max_x - min_x);
    double beta = 0.1 - min_x*alpha;

    for(int k = 0;k < N; k++)
        x_scal[k] = alpha*x[k]+beta;

    double *beta_vec = (double *) malloc((L-1)*sizeof(double));

    // perform discretization and obtain tridiagonal Jacobi matrix
    discrete_gauss_christoffel_quadr_pre_diag(N, w_scal, x_scal, L, nodes, beta_vec);

    const char compz = 'I';
    double *ev = (double *) malloc(L*L*sizeof(double));
    double *work = (double *) malloc(4*L*sizeof(double));
    int info;

    // diagonalize tridiagonal Jacobi matrix
    dpteqr_(&compz, (const int *) &L, nodes, beta_vec, ev, (const int *) &L, work, &info);

    // determine weights from the eigenvectors stored in ev + 
    // scale and shift nodes back to the original support
    for(int k = 0;k < L; k++){
        nodes[k] = (nodes[k]-beta)/alpha;
        weights[k] = ev[k*L]*ev[k*L]*total_weight;
    }

    free(w_scal); free(x_scal);
    free(beta_vec);
    free(ev); free(work);
}



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
void discrete_gauss_christoffel_quadr_pre_diag(const size_t N, const double *w, const double *x, const size_t L, double *alpha_vec_out, double *beta_vec_out){

    double *basis_subspace_out = (double *) malloc(N*L*sizeof(double ));

    arnoldi_ddiag(N, x, w, L, alpha_vec_out, beta_vec_out, basis_subspace_out);

    free(basis_subspace_out);

}



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
void arnoldi_ddiag(const size_t N, const double *H, const double *state0, const size_t dim, double *alpha_vec_out, double *beta_vec_out, double *basis_subspace_out){

    const int i_one = 1, i_zero = 0;
    const double d_one = 1., d_zero = 0.; 
    const char opt = 'N';

    double norm_state = dnrm2_((const int *) &N, state0, &i_one);
    double *state = (double *) malloc(N*sizeof(double));
    double *v = (double *) malloc(N*sizeof(double));
    double *h = (double *) malloc(N*sizeof(double));
    memset(v, 0, N*sizeof(double));
    memset(h, 0, N*sizeof(double));

    dcopy_((const int *) &N, state0, &i_one, state, &i_one);

    const double a1 = 1./norm_state;
    dscal_((const int *) &N, &a1 , state, &i_one);

    memset(basis_subspace_out, 0, N*dim*sizeof(double));

    memset(alpha_vec_out, 0, dim*sizeof(double));
    memset(beta_vec_out, 0, (dim-1)*sizeof(double));

    dcopy_((const int *) &N, state, &i_one, basis_subspace_out, &i_one);

    double beta = 0;

    for(int j=1; j<=dim; j++){

        dgbmv_(&opt, (const int *) &N, (const int *) &N , &i_zero, &i_zero, &d_one, H, &i_one, (basis_subspace_out + (j -1)*N), &i_one, &d_zero, v, &i_one);

        gsreorthog_d(N, dim, j, basis_subspace_out, v, v, h, &beta);
        
        alpha_vec_out[j-1] = h[j-1];
        if(j>1){
            beta_vec_out[j-2] += h[j-2];
            beta_vec_out[j-2] /= 2.;
        }

        if(j<dim){
            beta_vec_out[j-1] = beta;

            dcopy_((const int *) &N, v, &i_one, (basis_subspace_out + j*N), &i_one);
        }
    }
    free(state); free(v); free(h);
}



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
void gsreorthog_d(const size_t N_row, const size_t N_col, const size_t K_col, const double *Q, const double *x, double *q_out, double *r_out, double *rho_out){

    const int i_one = 1, i_zero = 0;
    const double d_one = 1., d_zero = 0., d_negone = -1.;
    const char optT = 'T', optN = 'N';

    size_t k, ndx;
    double tmp1, tmp2;

    // alpha is a critical value for the ratio of norms of the reorthogonalized vector x_orth to the input vector x
    // if ||x_orth||/||x|| > alpha than good, otherwise stay in the iteration loop
    // alpha = 0.5 is suggested by G W Stewart in his book

    double alpha = 0.5;
    double sigma = dnrm2_((const int *) &N_row, x, &i_one);
    double nu = sigma;
    
    double *x_orth = (double *) malloc(N_row*sizeof(double));
    dcopy_((const int *) &N_row, x, &i_one, x_orth, &i_one);

    double *s = (double *) malloc(K_col*sizeof(double));
    memset(s, 0, K_col*sizeof(double));

    dgemv_(&optT, (const int *) &N_row, (const int *) &K_col, &d_one, Q, (const int *) &N_row, x_orth, &i_one, &d_zero, s, &i_one);

    memset(r_out, 0, N_col*sizeof(double));
    daxpy_((const int *) &K_col, &d_one, s, &i_one, r_out, &i_one);

    dgemv_(&optN, (const int *) &N_row, (const int *) &K_col, &d_negone, Q, (const int *) &N_row, s, &i_one, &d_one, x_orth, &i_one);

    *rho_out = dnrm2_((const int *) &N_row, x_orth, &i_one);

    while(*rho_out/sigma < alpha){

        if(*rho_out > 0.1*nu*DBL_EPSILON){

            sigma = *rho_out;

        }
        else{
            
            sigma = 0.1*sigma*DBL_EPSILON;
            nu = sigma;

            ndx = 0;
            tmp1 = 0.; tmp2 = 0.;
            for(k=0; k< N_row; k++){

                tmp2 = dasum_((const int *) &K_col, (Q+k), (const int *) &N_row);
                if(k==0){
                    tmp1 = tmp2;
                }
                else{
                    if(tmp2 <= tmp1){
                        tmp1 = tmp2;
                        ndx = k;
                    }
                }
            }

            x_orth[ndx] = sigma;

        }

        dgemv_(&optT, (const int *) &N_row, (const int *) &K_col, &d_one, Q, (const int *) &N_row, x_orth, &i_one, &d_zero, s, &i_one);

        daxpy_((const int *) &K_col, &d_one, s, &i_one, r_out, &i_one);

        dgemv_(&optN, (const int *) &N_row, (const int *) &K_col, &d_negone, Q, (const int *) &N_row, s, &i_one, &d_one, x_orth, &i_one);

        *rho_out = dnrm2_((const int *) &N_row, x_orth, &i_one);
    }

    dcopy_((const int *) &N_row, x_orth, &i_one, q_out, &i_one);

    const double a1 = 1./(*rho_out);
    dscal_((const int *) &N_row, &a1, q_out, &i_one);

    free(x_orth);
    free(s);
}
