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
        x[k] = 2*M_PI*(double)k/NN;
        w[k] = sin(x[k])*sin(x[k]);
    }

    FILE *F;



    F = fopen("sin_squared_sample_1000.txt","w+");
    for(int k=0; k<=NN; k++)
        fprintf(F, "%f\t%f\n", x[k], w[k]);
    fclose(F);



    int dim = 5;
    double nodes[dim], weights[dim];
    gauss_christoffel(NN+1, w, x, M_PI, dim, nodes, weights);

    F = fopen("sin_squared_downsample_N5.txt","w+");
    for(int k=0; k<dim; k++)
        fprintf(F, "%f\t%f\n", nodes[k], weights[k]);
    fclose(F);



    dim = 30;
    double nodes2[dim], weights2[dim];
    gauss_christoffel(NN+1, w, x, M_PI, dim, nodes2, weights2);

    F = fopen("sin_squared_downsample_N100.txt","w+");
    for(int k=0; k<dim; k++){
        fprintf(F, "%f\t%f\n", nodes2[k], weights2[k]);
        psum += weights2[k];
    }
    printf("\n check total weight of 30 node sampling: %f\n", psum);
    fclose(F);



    return 0;
}
