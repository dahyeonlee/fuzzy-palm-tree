import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.special import yn, hankel1
from mpl_toolkits.mplot3d import Axes3D

def makesquare(L=5., N = 50):
	"""Makes a square of side L, 2N points per side. Total 8N points"""
	s = L/(2.*N)	# spacing between points
	x = []
	y = []
	for i in range(N+1):
		x.append(L/2.)
		y.append(s*i)
	for i in range(1,2*N+1):
		x.append(L/2. - s*i)
		y.append(L/2.)
	for i in range(1,2*N+1):
		x.append(-L/2.)
		y.append(L/2.-s*i)
	for i in range(1,2*N+1):
		x.append(-L/2.+s*i)
		y.append(-L/2.)
	for i in range(1,N+1):
		x.append(L/2.)
		y.append(-L/2.+s*i)
	x = np.array(x)
	y = np.array(y)
	return x, y

def cos_ij(xi,eta,j):
	if x[j]==xi and y[j]==eta:
		return 0.
	else:
		nx = (y[j+1]-y[j])/s[j]	# normal vector
		ny = (x[j]-x[j+1])/s[j]
		dx = x[j]-xi
		dy = y[j]-eta
		r = np.sqrt(dx**2 + dy**2)
		return (nx*dx + ny*dy)/r

def compute_c(k,xi,eta,j):
	"""Computes matrix elements of C used in Method III"""
	if x[j]==xi and y[j]==eta:
		c = 1. + 0.5j*k*s[j]
	else:
		c = 0.5j*k*s[j]*cos_ij(xi,eta,j)*hankel1(1,k*np.sqrt((x[j]-xi)**2+(y[j]-eta)**2))
	return c

def construct_cij(k,dim):
	"""Constructs Cij of dimension dim"""
	C = [[0. for col in range(dim)] for row in range(dim)]
	for i in range(dim):
		for j in range(dim):
			C[i][j] = compute_c(k,x[i],y[i],j)
	return np.array(C)

def find_min(a):
	"""Returns a filter for minima. Can be improved by fancy slicing"""
	indices=[]
	for i in range(len(a)):
		if a[i]<a[i-1] and a[i]<a[i+1]:
			indices.append(i)
	return indices

def null(A, eps=0.01):
    u, s, vh = sp.linalg.svd(A)
    null_mask = (s <= eps)
    null_space = sp.compress(null_mask, vh, axis=0)
    return sp.transpose(null_space[0])

def find_psi(k,xi,eta,dim):
	sum = 0
	for j in range(dim):
		sum -= p[j]*compute_c(k,xi,eta,j)
	return sum

vfind_psi = np.vectorize(find_psi)

# Start with a square
N = 40
x, y = makesquare(1.,N)
dx = x[1:]-x[:-1]
dy = y[1:]-y[:-1]
s = np.sqrt(dx**2+dy**2)

# This code shows the shape of the well
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.scatter(x,y)
# ax.set_aspect('equal')
# plt.show()

# This code plots abs(det(C)) vs k
# k = np.linspace(8.9643,8.9652,10)
# C = np.array([construct_cij(w,8*N) for w in k])
# C = np.log10(abs(np.linalg.det(C)))
# # print k[find_min(C)]
# plt.scatter(k,C)
# plt.show()

# k starting from 4.5019 (ground state) for N=20
#                 4.47247 for N=40
#				  8.9648

k = 4.47247
Cij = construct_cij(k,8*N)
print np.linalg.det(Cij)
p = null(Cij)
print Cij.dot(p)

u = np.linspace(-0.49,0.49,30)
v = np.linspace(-0.40,0.49,30)
uu, vv = np.meshgrid(u,v)
psi = vfind_psi(k,uu,vv,8*N)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(uu,vv,psi)
plt.show()
