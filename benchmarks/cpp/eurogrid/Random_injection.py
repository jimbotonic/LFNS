import numpy as np
import random

def main():

	num_nodes=6021
	P=[random.uniform(-0.1, 0.1) for i in range(0,num_nodes)]
	P=np.array(P)
	P=P-P.sum()/len(P)
	print "somma=",P.sum()
	np.savetxt("P.dat",P, fmt=("%.10f"))
	
if __name__ == '__main__':
    main()
