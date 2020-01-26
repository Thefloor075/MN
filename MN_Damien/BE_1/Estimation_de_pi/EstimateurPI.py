from pylab import figure


import random
import math
import matplotlib.pyplot as plt



pi = math.pi #pi
s = 0

A = 10
D = 10 
B = 10
P = 5

Th = []
L = []
Err = []
piT = []
Approx_Pi = []

for i in range(D,1000000):
	x = random.random() 
	y = random.random()
	if (x*x + y*y <= 1):
		s+=1
	approx_pi = 4*s/i
	piT.append(pi)
	Approx_Pi.append(approx_pi)
	L.append(i) 
	Err.append(abs(pi - approx_pi))
	Th.append(i**(-0.5))

biais = sum(Approx_Pi)/len(Approx_Pi)
print(biais)
std = sum(list((element - biais)**2 for element in Approx_Pi))/len(Approx_Pi) 
print(std)
	



