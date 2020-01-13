import math
from random import random
import matplotlib.pyplot as plt
import numpy as np

def MonteCarlo():
	return 0



def main():
	print("Hello world")



def trace(N):
	L = [random() for i in range(N)]
	Lp = L[1:len(L)]
	L = L[0:len(L)-1]
	for i in range(len(L)-1):
		plt.scatter(Lp[i],L[i],s=1,c=np.random.rand(3,))
	plt.show()



if __name__ == "__main__":
	C = 10000
	main()
	trace(C)

