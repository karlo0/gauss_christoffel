# Gauß-Christoffel quadrature

The function `gauss_christoffel(N, w, x, total_weight, L, nodes, weights)` performs a Gauß-Christoffel quadrature rule
to calculate for an arbitrary weight function *w\(x\)* \(*w(x) >= 0, &forall; x*\) *L* nodes
<span style="font-size:1.2em;">\(</span>
  `nodes[0]` \(*=x<sub>1</sub>*\), ..., `nodes[L-1]` \(*=x<sub>L</sub>*\)
  <span style="font-size:1.2em;">\)</span>
  and weights
  <span style="font-size:1.2em;">\(</span>
`weights[0]` \(*=w<sub>1</sub>*\), ..., `weights[L-1]` \(*=w<sub>L</sub>*\)
<span style="font-size:1.2em;">\)</span>.
The calculated weights conserve the first N moments of the given weight function *w\(x\)*.
The correct scaling of all the moments is ensured by giving to the function as an input parameter the `total_weight` \(*=&#x222b; w\(x\) dx*\). For more informations about the Gauß-Christoffel quadrature rule, see <https://dlmf.nist.gov/3.5#v>.
<br />

The Gauß-Christoffel quadrature rule can be used to: <br />

(1) obtain a downsampling of any positive function that needs to conserve the first
L moments of the function:<br />
<span style="color:white;opacity: 0.0;">0.</span>0th moment: *&Sigma;<sub>i=1,...,N</sub> w<sub>i</sub> = &#x222b; w\(x\) dx*,<br />
<span style="color:white;opacity: 0.0;">1.</span>1st moment: *&Sigma;<sub>i=1,...,N</sub> w<sub>i</sub> x<sub>i</sub> = &#x222b; w\(x\)x dx*,<br />
<span style="color:white;opacity: 0.0;">2.</span>&#8942;<br />
<span style="color:white;opacity: 0.0;">3.</span>Lth moment: *&Sigma;<sub>i=1,...,N</sub> w<sub>i</sub> x<sub>i</sub><sup>L</sup> = &#x222b; w\(x\)x<sup>L</sup> dx*,<br />
(2) calculate integrals of the form: *&#x222b; w\(x\) f\(x\) dx*,  with *w\(x\)* the weight function and *f\(x\)* an arbitray function.
The nodes *x<sub>i</sub>* and weights *w<sub>i</sub>* obtained from the Gauß-Christoffel quadrature are then used to calculate *&Sigma;<sub>i=1,...,N</sub> w<sub>i</sub> f(x<sub>i</sub>)* that is an numerical approximation to the integral.

---
## Usage
The following example performs a Gauß-Christoffel quadrature for a Gaussian with variance 2 centered around x = 0.
An extended version of this example is [gauss.c](/examples/gauss/gauss.c)<br />
Compiling requires linking against a blas and lapack.
<br />
```C
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
