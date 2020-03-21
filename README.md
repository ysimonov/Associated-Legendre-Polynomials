# Associated-Legendre-Polynomials
Derivatives of Normalized Associated Legendre Polynomials

So, I've spend a HUGE amount of time trying to calculate derivatives of normalised associated Legendre polynomials, 
and finally did it! This was required for calculation of far-field electric fields from FEKO simulation by expanding them
in terms of normalized vector spherical harmonics, commonly known as multipole expansion. 

For those people out there who are struggling to find a working version of normalized associated Legendre polynomials and their
derivatives, can grab this code and use it freely.

This program allows extended precision for calculation of very high order normalized Legendre functions. 
This was acomplished with module decimal, extending beyond standard python precision. 

The associated Legendre polynomials programmed here are currently valid only for m>=0.
The phase factor (-1)^m was omitted.

The script is quite raw, but fully operational. 
