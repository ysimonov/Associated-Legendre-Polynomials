The program generates

1) Normalized Associated Legendre Polynomials of size [L+1,L+1]

2) Derivatives of Normalized Associated Legendre Polynomials of size [L+1,L+1]

Normalization Constant: https://latex.codecogs.com/png.latex?%5CLARGE%20%5Chat%7BP%7D_%7BN%7D%5E%7B%7CM%7C%7D%28%5Ccos%5Ctheta%29%3D%5Csqrt%7B%5Cfrac%7B2N&plus;1%7D%7B2%7D%20%5Cfrac%7B%28N-%7CM%7C%29%21%7D%7B%28N&plus;%7CM%7C%29%21%7D%7D%20P_%7BN%7D%5E%7B%7CM%7C%7D%28%5Ccos%5Ctheta%29

Phase Factor (-1)^M is omitted.

References:

Limpanuparb, T. and Milthorpe, J., 2014. Associated Legendre Polynomials and Spherical Harmonics Computation for Chemistry Applications. arXiv preprint arXiv:1410.1748.

Li, P. and Jiang, L.J., 2012. The far field transformation for the antenna modeling based on spherical electric field measurements. Progress In Electromagnetics Research, 123, pp.243-261.

Wysin, G.M., 2011. Associated legendre functions and dipole transition matrix elements.

Associated_Legendre_Decimal.py is the extended precision version of Associated_Legendre.py. While it is more accurate, in can be significantly slower compared to floating point precision arithmetic:

Timing Comparison (in seconds): 

L = 30 : 0.2833747863769531 | 1.2349169254302979

L = 70 : 1.4793822765350342 | 5.806833505630493

L = 100 : 2.87738299369812 | 11.526304721832275

L = 200 : 12.345782041549683 | 46.637011766433716 

L = 300 : 25.131556034088135 | 103.22485566139221

L = 1000 : 302.46400713920593 | ????

Test machine: Intel(R) Core(TM) i7-2620M CPU @ 2.70GHz 4096 KB 4 CPUs 16Gb RAM


