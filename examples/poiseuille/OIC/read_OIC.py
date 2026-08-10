import sys, os
import matplotlib.pyplot as plt
import scipy.io
import pymech as pm

data = scipy.io.loadmat('OIC_mat.mat')

y = data['ycoord']

u = data['uinvar']
v = data['vinvar']
w = data['winvar']

fig, ax = plt.subplots(1,3)

ax[0].plot(u.real, y, label='Re(u)')
ax[0].plot(u.imag, y, label='Im(u)')
ax[1].plot(v.real, y, label='Re(v)')
ax[1].plot(v.imag, y, label='Im(v)')
ax[2].plot(w.real, y, label='Re(w)')
ax[2].plot(w.imag, y, label='Im(w)')

plt.show()
