import os
import sys
from os import listdir
from os.path import isfile, join
import shutil
import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def plot_arrows(ax,theta):
	length=0.6
	Nx=49
	for i in range(0,len(theta)):
		xc,yc=i%Nx,-(i-i%Nx)/Nx #arrow center
		xs,ys=xc-0.5*length*np.cos(theta[i]),yc-0.5*length*np.sin(theta[i])
		ax.arrow(xs, ys, length*np.cos(theta[i]), length*np.sin(theta[i]), head_width=0.3, head_length=0.3, fc='k', ec='k')

if __name__ == '__main__':
	mypath = sys.argv[1] 
	fn = [f for f in listdir(mypath) if isfile(join(mypath, f))]
	for i in range(0,len(fn)):
		theta=np.loadtxt(mypath + "/" + fn[i])
		theta+=np.pi/4
		### PLOTTING PART 
		fig, ax = plt.subplots()
		fig.subplots_adjust(right=0.95, left=0.05,top=0.95,bottom=0.05)
		ax.set_aspect('equal')
		ax.set_xlim(0.0,49.0)
		ax.set_ylim(0.0,-49.0)
		plot_arrows(ax,theta)
		positions=np.array([n for n in range(0,50)])
		ax.xaxis.set_ticks(positions)
		ax.yaxis.set_ticks(-positions)
		fig.set_size_inches(18.5, 18.5, forward=True)
		ax.set_title(fn[i])
		#plt.savefig(fn[i]+".pdf",format='pdf',dvi=700)
		plt.savefig(fn[i]+".jpg",format='jpg',dvi=700)
		ax.grid()
	


