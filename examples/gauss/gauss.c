#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "discrete_gauss_christoffel_quadr.h"


int main(){

    int NN = 1000;
    double sigma2 = 2;
    double x0 = 0;
    double x[NN+1], w[NN+1];
    double psum = 0.;
    for(int k=0; k<=NN; k++){
        x[k] = 5.*((double)k/NN - 0.5);
        w[k] = exp(-(x[k]-x0)*(x[k]-x0)/2./sigma2)/sqrt(2*M_PI*sigma2);
    }

    FILE *F;



    F = fopen("gauss_sample_1000.txt","w+");
    for(int k=0; k<=NN; k++)
        fprintf(F, "%f\t%f\n", x[k], w[k]);
    fclose(F);



    int dim = 5;
    double nodes[dim], weights[dim];
    discrete_gauss_christoffel_quadr(NN+1, w, x, 1., dim, nodes, weights);

    F = fopen("gauss_downsample_N5.txt","w+");
    for(int k=0; k<dim; k++)
        fprintf(F, "%f\t%f\n", nodes[k], weights[k]);
    fclose(F);



    dim = 30;
    double nodes2[dim], weights2[dim];
    discrete_gauss_christoffel_quadr(NN+1, w, x, 1., dim, nodes2, weights2);

    F = fopen("gauss_downsample_N100.txt","w+");
    for(int k=0; k<dim; k++){
        fprintf(F, "%f\t%f\n", nodes2[k], weights2[k]);
        psum += weights2[k];
    }
    printf("\n check total weight of 30 node sampling: %f\n", psum);
    fclose(F);



    return 0;
}
