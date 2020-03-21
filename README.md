# Associated-Legendre-Polynomials
Derivatives of Normalized Associated Legendre Polynomials

So, I've spend a HUGE amount of time trying to calculate derivatives of normalised associated Legendre polynomials, 
and finally did it! This was required for calculation of far-field electric fields from FEKO simulation by expanding them
in terms of normalized vector spherical harmonics, commonly known as multipole expansion. 

For those people out there who are struggling to find a working version of normalized associated Legendre polynomials and their
derivatives, can grab this code and use it freely.

This program allows extended precision for calculation of very high order normalized Legendre functions. 
This was acomplished with module decimal, extending beyond standard python precision. 

The associated Legendre polynomials programmed here are currently valid only for m>=0 and theta = [0,pi/2].
The phase factor (-1)^m was omitted.

The script is quite raw, but fully operational. 

#For more details regarding computation of fully normalized Associated Legendre polynomials and 
#Spherical harmonics, refer to the article below:
#Limpanuparb, T. and Milthorpe, J., 2014. 
#Associated Legendre Polynomials and Spherical Harmonics Computation for Chemistry Applications. 
#arXiv preprint arXiv:1410.1748.

I also used the recurrence relation after computing Pnm(cos(theta))/sin(theta) to evaluate dP/d(theta) (Just ask Wolframalpha
for derivative of LegendreP[n,m,cos(theta)] w.r.t. theta recurrence relation).

I used this reference for computation of P(cos(theta))/sin(theta) term:
#For more information regarding P/sin term, read
#Li, P. and Jiang, L.J., 2012. 
#The far field transformation for the antenna modeling based on spherical electric field measurements. 
#Progress In Electromagnetics Research, 123, pp.243-261.

Furthermore, I found out that the floor(n/2) has to be replaced with ceil(n/2) to give correct results! 
The limit of Legendre/sin at theta->0+ has to be negative of what is described in the above paper.

For more information, contact me on simonov.yevgeniy@gmail.com
