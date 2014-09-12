

def naca4_sym(c,t,n):
	"""
	NACA4_SYMMETRIC evaluates y(x) for a NACA symmetric 4-digit airfoil.
	C the chord length.
	T the maximum relative thickness.
	N the number of sample points.
	
	python wrapper for the code by  John Burkardt
	
	by Sukhbinder Singh
	
	"""
	import numpy as np
	import nacalib

	x = np.linspace(0.0,c,n)
	y=np.zeros(n)
	y=nacalib.naca4_symmetric(t,c,x,n)
	xy=np.zeros((2,2*n))
	xy[0,0:n]=x
	xy[0,n:]=x[n::-1]
	xy[1,0:n]=-y
	xy[1,n:]=y[n::-1]
	return xy

def naca4_cam(m,p,t,c,n):
	"""
	Generates (xu,yu), (xl,yl) for a NACA cambered 4-digit airfoil.
	
	M, the maximum camber.
	P, the location of maximum camber.
	T, the maximum relative thickness.
	C, the chord length.
	N, the number of sample points.
	
	python wrapper for the code by  John Burkardt
	
	by Sukhbinder Singh
	
	"""
	import numpy as np
	import nacalib
	
	xc=np.linspace(0.0,c,n)
	xu,yu,xl,yl = nacalib.naca4_cambered ( m, p, t, c,xc,n)
	xy=np.zeros((2,2*n))
	xy[0,0:n]=xl
	xy[1,0:n]=yl
	xy[0,n:]=xu[n::-1]
	xy[1,n:]=yu[n::-1]
	return xy
