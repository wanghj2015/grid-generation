import numpy as np
import matplotlib.pyplot as plt

from scipy.io import FortranFile

import sys
np.set_printoptions(threshold=sys.maxsize)



# ragged array:
# https://tonysyu.github.io/ragged-arrays.html#.XwTreJNKjlw

def stack_ragged(array_list, axis=0):
    lengths = [np.shape(a)[axis] for a in array_list]
    idx = np.cumsum(lengths[:-1])
    stacked = np.concatenate(array_list, axis=axis)
    return stacked, idx

def save_stacked_array(fname, array_list, axis=0):
    stacked, idx = stack_ragged(array_list, axis=axis)
    np.savez(fname, stacked_array=stacked, stacked_index=idx)

def load_stacked_arrays(fname, axis=0):
    npzfile = np.load(fname)
    idx = npzfile['stacked_index']
    stacked = npzfile['stacked_array']
    return np.split(stacked, idx, axis=axis)


#nlp = 45; nmp = 1; nion = 1

#print ('nlp, nmp, nion =', nlp, nmp, nion)


# read in grid

f = FortranFile('../run/igrf/general/fort.20', 'r')

f = FortranFile('../run/igrf/general/fort.20', 'r')

[nlp,nmp,npts] = f.read_ints(np.int32)

print (nlp,nmp,npts)

# skip some
f.read_ints(np.int32)

galt = []; glat = []; glon = []
dist = []; bmag = []
for m in range(0, nmp):
    for l in range(0, nlp):
        x = f.read_reals(dtype='float64')
        galt.append(x)
        x = f.read_reals(dtype='float64')
        glat.append(x)
        x = f.read_reals(dtype='float64')
        glon.append(x)
        x = f.read_reals(dtype='float64')
        bmag.append(x)
        x = f.read_reals(dtype='float64')
        dist.append(x)


f.close()


stacked, idx = stack_ragged(galt)
r = np.array(stacked) * 1.0e-3

stacked, idx = stack_ragged(glat)
theta = np.array(stacked)

print (r.shape)


x = (r-0)*np.sin(theta)
y = (r-0)*np.cos(theta)

#print (x)
#print (y)


fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_aspect('equal')


#plt.subplot(2, 1, 1)

#plt.plot(x, y, 'o-')

plt.scatter(x, y, s=0.2)

#plt.title('A tale of 2 subplots')
#plt.ylabel('Damped oscillation')


#thetaxx = np.linspace(0, np.pi, 100)
thetaxx = np.linspace(np.pi/4, np.pi*3/4, 100)


x1 = 6371.2*np.sin(thetaxx)
y1 = 6371.2*np.cos(thetaxx)

#plt.plot(x1, y1)

x2 = (6371.2+90)*np.sin(thetaxx)
y2 = (6371.2+90)*np.cos(thetaxx)

#plt.plot(x2, y2)


#plt.xlim(-2500, 500)
#plt.ylim(6000+00, 6371+1000)


#plt.subplot(2, 1, 2)
#plt.plot(x2, y2, '.-')
#plt.xlabel('time (s)')
#plt.ylabel('Undamped')

#plt.xlim(-3000, 0000)
#plt.ylim(6000+00, 6371+600)

#plt.xlim(-6000, 4000)
#plt.ylim(4000+00, 6371+6000)

#plt.xlim(-20000, 20000)
#plt.ylim(0, 40000)

#plt.xlim(-50000, 40000)
#plt.ylim(-4000+00, 6371+37000)

#plt.xlim(-500000, 400000)
#plt.ylim(-4000+00, 6371+370000)

#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))

#ax.grid()


fig.savefig("grid.png")


plt.show()



