# Approximation Kernels

This code produces the coefficients of kernels which may be used to approximate n-variate functions over the hypercube. 
The software is based on this paper (insert arxiv link).

We also include the coefficients of kernels obtained by running the program for different values of $n$ and $r$. The folder SigmaKernels contains coefficients of kernels with minimal $\sigma_r^2$ the folder SigmaBarKernels contains coefficients of kernels with minimal $\bar \sigma_r^2$. The order of the kernels is lexicographic (see https://en.wikipedia.org/wiki/Lexicographic_order).

To run the programm set $n$ to be the number of variables, $r$ to the degree, bar = true if you want to minimize $\bar \sigma_r^2$ and bar = false if you want to minimize $\sigma_r^2$. Set silent = true if you want to suppress the solver output. 

Plan: Also add a PDF explaining in details the theoretical details of the symmetry adapted basis construction.
