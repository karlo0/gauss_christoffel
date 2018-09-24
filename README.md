# Gauß-Christoffel quadrature

The function 'discrete\_gauss\_christoffel\_quadr(N, w, x, total\_weight, L, nodes, weights)' performs a Gauß-Christoffel quadrature rule
to calculate for an arbitrary weight function w(x) [w(x) >= 0 for all x] a number L of nodes (nodes`[0]`, ..., nodes`[L-1]`) and 
weights (weights`[0]`, ..., weights`[L-1]`). The calculated weights conserve the first N moments of the given weight function w(x).
For more informations about Gauß-Christoffel quadrature rule see (https://dlmf.nist.gov/3.5#v).
<br />
<br />
The Gauß-Christoffel quadrature rule can be used to: <br />
(1) to obtain a downsampling of any positive function that needs to conserve the first
N moments of the function <br />
( zeroth moment == total weight: integral w(x), <br />
  first moment = expectation value: integral x\*w(x), <br />
  ..., <br />
  Nth moment: integral x^N\*w(x) <br />
)
<br />
<br />
(2) to calculate integrals of the form: integral w(x)\*f(x),  with w(x) the weight function and f(x) an arbitray function.
The nodes x\_i and weights w\_i from the Gauß-Christoffel quadrature are used to calculate sum\_{i=1}^{N} w\_i f(x\_i) that approximates
the integral.

## Usage
The following example performs a Gauß-Christoffel quadrature for a Gaussian with variance 2 centered around x = 0.
An extension of this example can be found under /examples/gauss/gauss.c
<br />
<br />
```
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gauss_christoffel.h"


int main(){

    int NN = 1000;                                                          // length of the weight function w
    double sigma2 = 2;                                                      // the variance of the gaussian
    double x0 = 0;                                                          // the mean value of the gaussian
    double x[NN+1], w[NN+1];                                                // allocate arrays for the arguments and the weight function
    for(int k=0; k<=NN; k++){                                               
        x[k] = 5.*((double)k/NN - 0.5);                                     // assign argument values
        w[k] = exp(-(x[k]-x0)*(x[k]-x0)/2./sigma2)/sqrt(2*M_PI*sigma2);     // assign function values for a gaussian
    }

    FILE *F;

    F = fopen("gauss_sample_1000.txt","w+");                                // output the gaussian in the file "gauss_sample_1000.txt"
    for(int k=0; k<=NN; k++)                                                
        fprintf(F, "%f\t%f\n", x[k], w[k]);                                 // write in the first column the arguments and in the second column the values of the weight function
    fclose(F);


    // perform Gauß-Christoffel quadrature with 5 sampling points
    int dim = 5;
    double nodes[dim], weights[dim];

    // call the the Gauß-Christoffel quadrature rule function
    gauss_christoffel(NN+1, w, x, 1., dim, nodes, weights);

    F = fopen("gauss_downsample_N5.txt","w+");                              // output the gaussian in the file "gauss_downsample_N5.txt"
    for(int k=0; k<dim; k++)
        fprintf(F, "%f\t%f\n", nodes[k], weights[k]);                       // write in the first column the nodes and in the second column the weights of the downsampled weight function
    fclose(F);

    return 0;
}
```
<br />
<br />
Compiling requires linking against a blas and lapack
