import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import scipy
from numpy import random

def vibration(w,t,m,m1f,m1r,J,ktf,ktr,k1,k2,c1,c2,a,b,u):
	z = w[0]
	zp = w[1]
	z1 = w[2]
	z1p = w[3]
	z2 = w[4]
	z2p = w[5]
	theta = w[6]
	thetap = w[7]
	qr = np.sin(t)
	qf = np.sin(t-(a+b)/u)
	#mu = 10
	#sigma = 5
	#qf = np.sin(t)
	#qr = qf
	print(t)
	z1pp = (ktf*qf-c1*(z1p+a*thetap-zp)-(k1+ktf)*z1-k1*(a*theta-z))/m1f
	z2pp = (ktr*qr-c2*(z2p-b*thetap-zp)-(k2+ktr)*z2+k2*(b*theta+z))/m1r
	zpp = (c1*z1p+c2*z2p-(c1+c2)*zp-(b*c2-a*c1)*thetap+k1*z1+k2*z2-(k1+k2)*z-(b*k2-a*k1)*theta)/m
	thetapp = (-a*c1*z1p+b*c2*z2p-(b*c2-a*c1)*zp-(a*a*c1+b*b*c2)*thetap-a*k1*z1+b*k2*z2-(b*k2-a*k1)*z-(a*a*k1+b*b*k2)*theta)/J
	#return [z,zp,zpp,z1,z1p,z1pp,z2,z2p,z2pp,theta,thetap,thetapp]
	return [zp,zpp,z1p,z1pp,z2p,z2pp,thetap,thetapp]
	#return [z,zp,z1,z1p,z2,z2p,theta,thetap]
	
def integer_z(x):
	global p
	global t
	#print((p[:x],t[:x]))
	return 	scipy.integrate.trapz(p[:x],t[:x])


t = np.arange(0,1000,1)
track1 = odeint(vibration,[0,0,0,0,0,0,0,0],t,args=(690,40,45,1220,200.1,200,20,24,1.5,1.51,1.3,1.52,20.1))
p = track1[:,0]
interger_z_np = np.frompyfunc(integer_z,1,1)
x = np.arange(0,1000,1)
z = interger_z_np(x)
#plt.plot(t,track1[:,0],t,track1[:,1])#,t,track1[:,0],t,track1[:,1])
plt.plot(t,z,t,track1[:,0],t,track1[:,1])



plt.show()
