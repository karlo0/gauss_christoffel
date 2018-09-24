#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gauss_christoffel.h"


int main(){

    int NN = 2000;
    double gamma = 1;
    double x1 = -1.5;
    double x2 = 1.5;
    double x[NN+1], w[NN+1];
    double psum = 0.;
    for(int k=0; k<=NN; k++){
        x[k] = 10.*((double)k/NN - 0.5);
        w[k] = gamma/((x[k]-x1)*(x[k]-x1) + gamma*gamma)/M_PI;
        w[k] += gamma/((x[k]-x2)*(x[k]-x2) + gamma*gamma)/M_PI;
    }

    FILE *F;



    F = fopen("double_cauchy_sample_1000.txt","w+");
    for(int k=0; k<=NN; k++)
        fprintf(F, "%f\t%f\n", x[k], w[k]);
    fclose(F);



    int dim = 5;
    double nodes[dim], weights[dim];
    gauss_christoffel(NN+1, w, x, 2., dim, nodes, weights);

    F = fopen("double_cauchy_downsample_N5.txt","w+");
    for(int k=0; k<dim; k++)
        fprintf(F, "%f\t%f\n", nodes[k], weights[k]);
    fclose(F);



    dim = 30;
    double nodes2[dim], weights2[dim];
    gauss_christoffel(NN+1, w, x, 2., dim, nodes2, weights2);

    F = fopen("double_cauchy_downsample_N100.txt","w+");
    for(int k=0; k<dim; k++){
        fprintf(F, "%f\t%f\n", nodes2[k], weights2[k]);
        psum += weights2[k];
    }
    printf("\n check total weight of 30 node sampling: %f\n", psum);
    fclose(F);



    return 0;
}
