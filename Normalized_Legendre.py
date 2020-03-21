import numpy as np
import math
import matplotlib.pyplot as plt 
from scipy.special import poch,factorial
from decimal import *
import decimal
getcontext().prec = 100 #this allows computation of L/sin up to 161-th order of n

def pochs(z,m):
   prod = Decimal(z) 
   if(m>0):
      j = Decimal(z)
      while(j < z + m - 1):
         prod *= Decimal(j + 1)
         j += Decimal(1)
   return prod
   
def fact(k):
   if(k<0):
      return None
   elif(k==1 or k==0):
      return Decimal(1)
   else:
      prod = Decimal(1)
      j = Decimal(1)
      while(j <= k):
         prod = prod * Decimal(j)
         j += Decimal(1)
      return prod
      
def PT(l,m):
   return int((m)+((l)*((l)+1))/2)

def computeP(L,idx,A,B,P,x):
   sintheta = math.sqrt(1.-x*x)
   temp = math.sqrt(0.5)
   P[idx[0,0]] = 0.7071067811865475244
   if(L>0):
      SQRT3 = 1.7320508075688772935
      P[idx[1,0]] = x*SQRT3*temp
      SQRT3DIV2 = -1.2247448713915890491
      temp = SQRT3DIV2*sintheta*temp
      P[idx[1,1]]=temp
      for l in range(2,L+1):
         for m in range(0,l-1):
            P[idx[l,m]]=A[idx[l,m]]*(x*P[idx[l-1,m]] + \
                       B[idx[l,m]]*P[idx[l-2,m]])
         P[idx[l,l-1]]=x*math.sqrt(2*(l-1)+3)*temp
         temp=-math.sqrt(1.0+0.5/l)*sintheta*temp
         P[idx[l,l]]=temp
         
   return P

LL = 200 #Order and degree of legendre polynomials
size = PT(LL+1,LL+1)
idx = np.zeros((LL+1,LL+1),np.int_)
for l in range(0,LL+1):
   for m in range(0,LL+1):
      idx[l,m] = PT(l,m) #indexes of Legendre polynomials
      
A = np.zeros((size),dtype=np.float_)
B = np.zeros((size),dtype=np.float_)
P = np.zeros((size),dtype=np.float_)

for l in range(2,LL+1):
   ls = l*l
   lm1s = (l-1)*(l-1)
   for m in range(0,l-1):
      ms = m*m
      A[idx[l,m]] = math.sqrt((4*ls-1.)/(ls-ms))
      B[idx[l,m]] = -math.sqrt((lm1s-ms)/(4*lm1s-1.))

degs_per_pixel = 0.5
deg2rad = math.pi / 180
NPTS = int(90 / degs_per_pixel) + 1
theta = np.arange(0,NPTS)*degs_per_pixel*deg2rad 

x = np.cos(theta)
y = np.sin(theta)

Leg = np.zeros((size,NPTS),np.float_) # L,M,Theta

for i in range(0,len(x)):
   Leg[:,i] = computeP(LL,idx,A,B,P,x[i])

#compute P/sin functions
Leg_sin = np.zeros((size,NPTS),np.float_)
Leg_sin[:,1:] = Leg[:,1:] / y[1:]

#consider the case when theta=0:
m=1
for n in range(1,LL+1):
   summ = Decimal(0)
   for k in range(0,math.ceil(n/2)):
      summ += Decimal(-1) ** Decimal(k) * ( pochs(n-k+1,n-k) / fact(k) / fact(n-2*k-1) ) * Decimal(1)**Decimal(n-2*k-1)
              
   Leg_sin[idx[n,m],0] = -float(summ / Decimal(2) ** Decimal(n) * (Decimal(2*n+1)/Decimal(2)* \
                          fact(n-m)/fact(n+m)) ** Decimal(0.5))
                          
   print(n,1,Leg_sin[idx[n,1],0])

#Evaluate derivatives of normalized Associated Legendre polynomials using recurrence relation 
   Leg_deriv = np.zeros((size,NPTS),np.float_)
   for n in range(0,LL):
      for m in range(0,LL+1):
         Leg_deriv[idx[n,m],:] = -(n+1) * x[:] * Leg_sin[idx[n,m],:] + (n-m+1) * Leg_sin[idx[n+1,m],:] 

fig = plt.figure()
ax = plt.gca()
ax.set_yscale('symlog')

for n in range(195,200):
   for m in range(195,200):
      ax.scatter(theta[:],Leg[idx[n,m],:], label = 'n = ' + str(n) + ',' + 'm = ' + str(m))
plt.xlabel('theta')
plt.ylabel('L(cos(theta))')
plt.title('Normalised Associated Legendre Polynomialsfor n=1 to 200 and m=1')
plt.grid()
plt.legend()
plt.show()

for n in range(195,200):
   for m in range(195,200):
      ax.scatter(theta[:],Leg_sin[idx[n,m],:], label = 'n = ' + str(n) + ',' + 'm = ' + str(m))
plt.xlabel('theta')
plt.ylabel('L(cos(theta))/sin(theta)')
plt.title('Normalised Associated Legendre Polynomials over sin(theta) for n=1 to 200 and m=1')
plt.grid()
plt.legend()
plt.show()

for n in range(1,200):
   for m in range(195,200):
      ax.scatter(theta[:],Leg_deriv[idx[n,m],:], label = 'n = ' + str(n) + ',' + 'm = ' + str(m))
plt.xlabel('theta')
plt.ylabel('dL(cos(theta))/d(theta)')
plt.title('Derivatives of Normalised Associated Legendre Polynomials for n=1 to 200 and m=1')
plt.grid()
plt.legend()
plt.show()
