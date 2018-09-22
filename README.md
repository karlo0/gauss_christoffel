# gauss_christoffel

The function 'discrete_gauss_christoffel_quadr(N, w, x, total_weight, L, nodes, weights)' performs a Gauß-Christoffel quadrature rule
to calculate for an arbitrary weight function w(x) [w(x) >= 0 for all x] a number L of nodes x_i
and weights w_i. The calculated weights conserve the total weight of the weight function [sum_{i=1}^{L} wi == integral w(x)].

The Gauß-Christoffel quadrature rule can be used:
(1) to obtain a downsampling of any positive function that needs to conserve the first
N moments of the function 
( zeroth moment == total weight: integral w(x), 
  first moment = expectation value: integral x*w(x),
  ...,
  Nth moment: integral x^N*w(x)
)

(2) to calculate integrals of the form: integral w(x)*f(x),  with w(x) the weight function and f(x) an arbitray function.
The nodes x_i and weights w_i from the Gauß-Christoffel quadrature are used to calculate sum_{i=1}^{N} w_i f(x_i) that approximates
the integral.

Usage:
