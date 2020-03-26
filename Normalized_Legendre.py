import numpy as np
import math
import matplotlib.pyplot as plt 

def computeP(L,x,y,A,B):
   P = np.zeros((L+1,L+1),dtype=np.float64)
   temp = 1.0/math.sqrt(2.0)
   P[0,0] = temp
   if(L>0):
      SQRT3 = math.sqrt(3.0)
      P[1,0] = x*SQRT3*temp
      SQRT3DIV2 = -SQRT3*temp
      temp = SQRT3DIV2*y*temp
      P[1,1]=temp
      for l in range(2,L+1):
         for m in range(0,l-1):
            P[l,m]=A[l,m]*(x*P[l-1,m] + B[l,m]*P[l-2,m]) #(14)
         P[l,l-1]=x*math.sqrt(2.0*(l-1.0)+3.0)*temp #(13)
         temp=-math.sqrt(1.0+0.5/l)*y*temp
         P[l,l]=temp
   return P

def generateP(L,theta):
    """
    generates Normalized Associated Legendre Polynomials and
    their derivateves for all n,m = 0 .. L
    
    Input: L - integer, order of LegendreP
           theta - angular grid in radians

    Output: LegendreP = sqrt( ( 2 * l + 1 ) / 2 * ( l - |m| )! / ( l + |m| )! ) * Legendre[l,|m|,cos(theta)]
            LegendreD = d/dtheta { LegendreP[l,|m|,cos(theta)] }

    There is no phase factor (-1)^M in this implementation

    Legendre Polynomials = 0 for all |m| > l

    Note, LegendreP / sin(theta) -> Infinity for m=0 at theta = 0 and theta = pi
    To prevent NaNs, I assumed these terms to be 0 instead of Infinity

    """
    A = np.zeros((L+2,L+2),dtype=np.float64)
    B = np.zeros((L+2,L+2),dtype=np.float64)
    for l in range(2,L+2):
        ls = l**2
        lm1s = (l-1)**2
        for m in range(0,l-1):
            ms = m**2
            A[l,m] = math.sqrt((4*ls-1)/(ls-ms))
            B[l,m] = -math.sqrt((lm1s-ms)/(4*lm1s-1))

    xmax = np.size(theta)
    x = np.cos(theta)
    y = np.sin(theta)

    LegendreP = np.zeros((L+2,L+2,xmax),dtype=np.float64)
    for i in range(0,xmax):
        LegendreP[:,:,i] = computeP(L+1,x[i],y[i],A,B)

    LegendreS = np.zeros((L+2,L+2,xmax),dtype=np.float64)

    #find asymptotes
    eps = 1e-13
    mask1 = (abs(theta-0.0)<eps)
    mask2 = (abs(theta-math.pi)<eps)
    mask = np.ma.mask_or(mask1,mask2)
    asymptote = np.where(mask)[0]
    regular = np.where(~mask)[0]

    LegendreS[:,:,regular] = LegendreP[:,:,regular] / y[regular]

    m = 1
    for l in range(1,L+2):
        LegendreS[l,m,asymptote[0]] = -0.5 * math.sqrt(l*(2.0*l+1.0)*(l+1.0)/2.0)
        LegendreS[l,m,asymptote[1]] = (-1) ** (l+1) * LegendreS[l,m,asymptote[0]]

    LegendreD = np.zeros((L+1,L+1,xmax),dtype=np.float64)
    for l in range(0,L+1):
        for m in range(0,l):
            LegendreD[l,m,:] = -(l+1) * x[:] * LegendreS[l,m,:] + \
                                    (l-m+1) * LegendreS[l+1,m,:] * \
                                    math.sqrt((2*l+1)*(l+m+1)/(2*l+3)/(l-m+1))

    return LegendreP[0:L+1,0:L+1], LegendreD

if __name__ == "__main__":

    L = 100 #maximum order of associated Legendre polynomials

    dtheta = 0.5 #angular spacing
    theta = np.arange(0,180+dtheta,dtheta,dtype=np.float64) * math.pi / 180
    theta_min = theta[0]
    theta_max = theta[-1]

    LegendreP, LegendreD = generateP(L,theta)

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

