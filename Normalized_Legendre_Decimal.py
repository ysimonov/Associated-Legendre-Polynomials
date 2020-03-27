import numpy as np
import math
import matplotlib.pyplot as plt 
import decimal
from decimal import Decimal
from decimal import getcontext
import time

decimal.getcontext().prec = 50

decimal_type = np.dtype(Decimal)

def deg2rad(angle_deg):
    """
    multiplicative conversion factor from radians to degrees
    angle_deg - angle in degrees
    """
    return angle_deg * pi() / Decimal(180)

def pi():
    """Compute Pi to the current precision.

    >>> print pi()
    3.141592653589793238462643383

    """
    getcontext().prec += 2  # extra digits for intermediate steps
    three = Decimal(3)      # substitute "three=3.0" for regular floats
    lasts, t, s, n, na, d, da = 0, three, 3, 1, 0, 0, 24
    while s != lasts:
        lasts = s
        n, na = n+na, na+8
        d, da = d+da, da+32
        t = (t * n) / d
        s += t
    getcontext().prec -= 2
    return +s               # unary plus applies the new precision

def cos(x):
    """Return the cosine of x as measured in radians.

    >>> print cos(Decimal('0.5'))
    0.8775825618903727161162815826
    >>> print cos(0.5)
    0.87758256189
    >>> print cos(0.5+0j)
    (0.87758256189+0j)

    """
    getcontext().prec += 2
    i, lasts, s, fact, num, sign = 0, 0, 1, 1, 1, 1
    while s != lasts:
        lasts = s
        i += 2
        fact *= i * (i-1)
        num *= x * x
        sign *= -1
        s += num / fact * sign
    getcontext().prec -= 2
    return +s

def sin(x):
    """Return the sine of x as measured in radians.

    >>> print sin(Decimal('0.5'))
    0.4794255386042030002732879352
    >>> print sin(0.5)
    0.479425538604
    >>> print sin(0.5+0j)
    (0.479425538604+0j)

    """
    getcontext().prec += 2
    i, lasts, s, fact, num, sign = 1, 0, x, 1, x, 1
    while s != lasts:
        lasts = s
        i += 2
        fact *= i * (i-1)
        num *= x * x
        sign *= -1
        s += num / fact * sign
    getcontext().prec -= 2
    return +s

def computeP(L,x,y,A,B):
    """
    computes Normalized Associated Legendre Polynomials
    A, B - arrays used in recurrence
    x - cos(theta)
    y - sin(theta)
    L - degree of P

    Returns P[0:L+1,0:L+1] for all l = 0..L, |m| = 0..L
    """
   P = np.zeros((L+1,L+1),dtype=decimal_type)
   one = Decimal(1)
   two = Decimal(2)
   three = Decimal(3)
   temp = one / two.sqrt() #1/sqrt(2)
   xd = Decimal(x)
   yd = Decimal(y)
   P[0,0] = temp
   if(L>0):
      SQRT3 = three.sqrt() #sqrt(3)
      SQRT3DIV2 = -SQRT3 * temp #-sqrt(3/2)
      P[1,0] = -xd*SQRT3DIV2
      temp = SQRT3DIV2*Decimal(y)*temp
      P[1,1]=temp
      for l in range(2,L+1):
         ld = Decimal(l)
         for m in range(0,l-1):
            P[l,m] = A[l,m]*(xd * P[l-1,m] + B[l,m] * P[l-2,m]) #(14)
         P[l,l-1]=(two*(ld-one)+three).sqrt()*xd*temp #(13)
         temp=-(one+one/(two*ld)).sqrt()*yd*temp
         P[l,l]=temp
   return P

def generateP(L,theta_min,theta_max,theta_inc):
    """
    generates Normalized Associated Legendre Polynomials and
    their derivateves for all n,m = 0 .. L
    
    Input: L - integer, order of LegendreP
           theta_min - minimum value of theta, in degrees
           theta_max - maximum value of theta, in degrees
           theta_inc - incremental step, in degrees

    Output: LegendreP = sqrt( ( 2 * l + 1 ) / 2 * ( l - |m| )! / ( l + |m| )! ) * Legendre[l,|m|,cos(theta)]
            LegendreD = d/dtheta { LegendreP[l,|m|,cos(theta)] }

    There is no phase factor (-1)^M in this implementation

    Legendre Polynomials = 0 for all |m| > l

    Note, LegendreP / sin(theta) -> Infinity for m=0 at theta = 0 and theta = pi
    To prevent NaNs, I assumed these terms to be 0 instead of Infinity

    """

    xmax = int((theta_max + theta_inc - theta_min) / theta_inc)
    theta = np.zeros((xmax),dtype=decimal_type)

    dtheta = Decimal(theta_inc) * deg2rad(1)
    theta[0] = Decimal(theta_min) * deg2rad(1)

    eps = 1e-14

    for i in range(1,xmax):
        theta[i] = theta[i-1] + dtheta

    theta[abs(theta)<eps] = Decimal(0)

    one = Decimal(1)
    two = Decimal(2)
    three = Decimal(3)
    four = Decimal(4)

    #construct arrays that will be used in Legendre recurrence
    A = np.zeros((L+2,L+2),dtype=decimal_type)
    B = np.zeros((L+2,L+2),dtype=decimal_type)
    for l in range(2,L+2):
        ld = Decimal(l)
        ls = ld * ld
        lm1s = (ld - one) * (ld - one)
        ls1 = four * ls - one
        lm1s1 = four * lm1s - one
        for m in range(0,l-1):
            md = Decimal(m)
            ms = md * md
            A[l,m] = Decimal(ls1/(ls-ms)).sqrt()
            B[l,m] = -Decimal((lm1s-ms)/lm1s1).sqrt()

    #construct cos(theta) and sin(theta)
    x = np.zeros((xmax),dtype=decimal_type)
    y = np.zeros((xmax),dtype=decimal_type)
    for i in range(0,xmax):
        x[i] = cos(theta[i])
        y[i] = sin(theta[i])
        # remove zeros to prevent underflow
        if(abs(x[i])<eps): x[i] = Decimal(0)
        if(abs(y[i])<eps): y[i] = Decimal(0)

    #compute Normalized Associated Legendre Polynomials
    LegendreP = np.zeros((L+2,L+2,xmax),dtype=decimal_type)
    for i in range(0,xmax):
        LegendreP[:,:,i] = computeP(L+1,x[i],y[i],A,B)
    LegendreP[abs(LegendreP) < 1e-15] = Decimal(0)

    #find asymptotes for computation of Legendre / sin term
    mask1 = (abs(theta - Decimal(0)) < eps)
    mask2 = (abs(theta - pi()) < eps)
    mask = np.ma.mask_or(mask1,mask2)

    asymptote = np.where(mask)[0]
    regular = np.where(~mask)[0]

    #construct Legendre / sin
    LegendreS = np.zeros((L+2,L+2,xmax),dtype=decimal_type)
    LegendreS[:,:,regular] = LegendreP[:,:,regular] / y[regular]

    #consider special case when m = 1
    m = 1
    p5 = Decimal(0.5)
    for l in range(1,L+2):
        ld = Decimal(l)
        lp1 = ld + one
        LegendreS[l,m,asymptote[0]] = -p5 * Decimal(ld*(two*ld+one)*lp1*p5).sqrt()
        LegendreS[l,m,asymptote[1]] = -one ** lp1 * LegendreS[l,m,asymptote[0]]

    #compute Normalized Associated Legendre derivatives using recurrence relation
    LegendreD = np.zeros((L+1,L+1,xmax),dtype=decimal_type)
    for l in range(0,L+1):
        ld = Decimal(l)
        lp1 = ld + one
        lp3 = two * ld + three
        for m in range(0,l):
            md = Decimal(m)
            lmp1 = ld - md + one
            LegendreD[l,m,:] = -lp1 * x[:] * LegendreS[l,m,:] + lmp1 * LegendreS[l+1,m,:] * \
                                Decimal((two*ld+one)*(ld+md+one)/lp3/lmp1).sqrt()

    return np.float64(LegendreP[0:L+1,0:L+1]), np.float64(LegendreD)

if __name__ == "__main__":

    L = 100 #maximum order of associated Legendre polynomials

    theta_decimal_max = Decimal(180)

    #grid parameters in degrees
    theta_min = 0 #minimum value, in degrees
    theta_max = 180 #maximum value, in degrees
    theta_inc = 0.5 #increment, in degrees

    theta = np.arange(theta_min,theta_max+theta_inc,theta_inc,dtype=np.float64)

    start_time = time.time()

    #call functions to generate Legendre Polynomials in Decimal Precision and return float64 arrays
    LegendreP, LegendreD = generateP(L,theta_min,theta_max,theta_inc)

    stop_time = time.time()
    print('Time elapsed = ', stop_time-start_time)

    #plot results
    m = L-5
    for l in range(m,L+1):
        plt.plot(theta,LegendreD[l,m,:],linewidth=0.5)

    plt.title('Derivatives of Normalized Associated Polynomials up to l=' + str(L) + ' for m='+str(m))
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'$\partial\{\hat{P}_{l}^{m}(\cos\theta)\}/\partial\theta$')
    plt.xlim(theta_min,theta_max)
    plt.grid()
    plt.savefig('Lmax='+str(L)+'_M='+str(m)+'.png',dpi=800)

