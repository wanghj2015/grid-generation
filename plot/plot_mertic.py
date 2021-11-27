import numpy as np
import matplotlib.pyplot as plt

from scipy.io import FortranFile

import sys
np.set_printoptions(threshold=sys.maxsize)


# read in 'ell'

f = FortranFile('../run/dipole/fort.51', 'r')

ell = f.read_reals(dtype=np.float64)

f.close()

#print (ell)


# read in 'dipole' (dp) metric

f = FortranFile('../run/dipole/fort.20', 'r')

[nlp,nmp,npts] = f.read_ints(np.int32)

print (nlp,nmp,npts)

# skip some
f.read_ints(np.int32)

# read in again
r = []; theta = []; b = []; s = []; s1 = []
for m in range(0, nmp):
    for l in range(0, nlp):
        x = f.read_reals(dtype=np.float64)
        r.append(x)
        x = f.read_reals(dtype=np.float64)
        theta.append(x)
        x = f.read_reals(dtype=np.float64)
        b.append(x)
        x = f.read_reals(dtype=np.float64)
        s.append(x)
        s1.append(x[len(x)-1])

f.close()

#print (r[0])

#print (type(r[0]))

#print (s)

r_dp = r[0]
theta_dp = theta[0]
b_dp = b[0]
s_dp = s[0]
ell_dp = np.array(s1)


f = FortranFile('../run/dipole/fort.21', 'r')

h1 = []; h2 = []; h3 = []
for m in range(0, nmp):
    for l in range(0, nlp):
        x = f.read_reals(dtype='float64')
        h1.append(x)
        x = f.read_reals(dtype='float64')
        h2.append(x)
        x = f.read_reals(dtype='float64')
        h3.append(x)

f.close()


h1_dp = h1[0]
h2_dp = h2[0]
h3_dp = h3[0]



# read in 'general coordinates' (gc) metric

f = FortranFile('../run/igrf/dipole/fort.20', 'r')

[nlp,nmp,npts] = f.read_ints(np.int32)

print (nlp,nmp,npts)

# skip some
f.read_ints(np.int32)

# read in again
r = []; theta = []; b = []; s = []; s1 = []
for m in range(0, nmp):
    for l in range(0, nlp):
        x = f.read_reals(dtype='float64')
        r.append(x)
        x = f.read_reals(dtype='float64')
        theta.append(x)
        x = f.read_reals(dtype='float64')
        # skip
        x = f.read_reals(dtype='float64')
        b.append(x)
        x = f.read_reals(dtype='float64')
        s.append(x)
        s1.append(x[len(x)-1])

f.close()

#print (r[0])

#print (type(r[0]))

#print (s)

r_gc = r[0]
theta_gc = theta[0]
b_gc = b[0]
s_gc = s[0]
ell_gc = np.array(s1)


f = FortranFile('../run/igrf/dipole/fort.21', 'r')

h1 = []; h2 = []; h3 = []
for m in range(0, nmp):
    for l in range(0, nlp):
        x = f.read_reals(dtype='float64')
        h1.append(x)
        x = f.read_reals(dtype='float64')
        h2.append(x)
        x = f.read_reals(dtype='float64')
        h3.append(x)

f.close()


h1_gc = h1[0]
h2_gc = h2[0]
h3_gc = h3[0]



# process data and plot


#fig, ax = plt.subplots()

fig = plt.figure()
gs = fig.add_gridspec(6, hspace=0.5)
#ax = gs.subplots(sharex=True, sharey=True)
ax = gs.subplots(sharex=True, sharey=False)

#fig.suptitle('Relative errors')


r_earth = 6371.2

s_dp = s_dp / 1000


print (npts//2)

s_dp = s_dp - s_dp[npts//2]


gm = 29376.299999999999
gm = 29469.014999999999
gm = 29619.400000000001
gm = 29619.4

h1_gc = h1_gc * (gm * r_earth)

#ax.plot(s_dp, h1_gc / h1_dp)
#ax.plot(s_dp, h1_gc - h1_dp)
#ax[3].plot(s_dp, (h1_gc - h1_dp) / h1_dp)

yy = (h1_gc - h1_dp) / h1_dp
print (' h1 = ', '{:.4e}'.format(np.min(yy)), '{:.4e}'.format(np.max(yy)))

ax[3].plot(s_dp, yy)


h2_gc = h2_gc / r_earth

#ax.plot(h2_gc / h2_dp)

#ax[4].plot(s_dp, (h2_gc - h2_dp) / h2_dp * 1.0e-3)
ax[4].plot(s_dp, (h2_gc - h2_dp) / h2_dp)

yy = (h2_gc - h2_dp) / h2_dp
print (' h2 = ', '{:.4e}'.format(np.min(yy)), '{:.4e}'.format(np.max(yy)))

#ax.plot(h3_gc / h3_dp)
#ax.plot(s_dp, h3_gc - h3_dp)
ax[5].plot(s_dp, (h3_gc - h3_dp) / h3_dp)

yy = (h3_gc - h3_dp) / h3_dp
print (' h3 = ', '{:.4e}'.format(np.min(yy)), '{:.4e}'.format(np.max(yy)))

#ax.plot(s_dp, r_gc / r_dp)
#ax.plot(s_dp, r_gc - r_dp)
ax[0].plot(s_dp, (r_gc - r_dp) / r_dp)

yy = (r_gc - r_dp) / r_dp
print (' r = ', '{:.4e}'.format(np.min(yy)), '{:.4e}'.format(np.max(yy)))
#print (' r = ', np.min((r_gc - r_dp) / r_dp), np.max((r_gc - r_dp) / r_dp))

theta_dp = np.pi/2 - theta_dp
theta_gc = np.pi/2 - theta_gc

#print ((theta_gc - theta_dp)*180/np.pi)

#ax.plot(s_dp, theta_gc / theta_dp)
#ax.plot(s_dp, theta_gc - theta_dp)
#ax.plot(s_dp, (theta_gc - theta_dp) / theta_dp)

def safe_div(x,y):
    if y==0: return 0
    return x/y

#npts = len(theta_dp)
y = np.zeros(npts)
for i in range(0, npts):
    y[i] = safe_div(theta_gc[i] - theta_dp[i], theta_dp[i])

ax[1].plot(s_dp, y)

yy = y
print (' theta = ', '{:.4e}'.format(np.min(yy)), '{:.4e}'.format(np.max(yy)))
#print (' theta = ', np.min(y), np.max(y))

#ax.plot(s_dp, b_gc / b_dp)
ax[2].plot(s_dp, (b_gc - b_dp) / b_dp)

yy = (b_gc - b_dp) / b_dp
print (' B = ', '{:.4e}'.format(np.min(yy)), '{:.4e}'.format(np.max(yy)))
#print (' B = ', np.min((b_gc - b_dp) / b_dp), np.max((b_gc - b_dp) / b_dp))

#ax.plot(ell_dp / ell)
#ax.plot((ell_dp - ell) / ell)

#ax.plot(ell_gc / ell_dp)
#ax.plot((ell_gc - ell_dp) / ell_dp)
#ax.plot((ell_gc - ell) / ell)


print ('ell_gc = ', '{:.4f}'.format(ell_gc[0]/1000))
print ('ell_dc = ', '{:.4f}'.format(ell_dp[0]/1000))
print ('ell    = ', '{:.4f}'.format(ell   [0]/1000))

print ('r_gc/dc  = ', '{:.4e}'.format((ell_gc[0] - ell_dp[0]) / ell_dp[0]))
print ('r_gc/ell = ', '{:.4e}'.format((ell_gc[0] - ell[0]) / ell[0]))
print ('r_dc/ell = ', '{:.4e}'.format((ell_dp[0] - ell[0]) / ell[0]))


#ax.set(xlabel='arc distance (km)', ylabel='relative error',
#       title='Relative errors')

# add grid
for i in range(0, 6):
    ax[i].grid()


#ax[0].text(0.97, 0.36, 'a', transform=ax[0].transAxes, fontsize=10, va='top')
#ax[1].text(0.97, 0.36, 'b', transform=ax[1].transAxes, fontsize=10, va='top')
#ax[2].text(0.97, 0.36, 'c', transform=ax[2].transAxes, fontsize=10, va='top')
#ax[3].text(0.97, 0.36, 'd', transform=ax[3].transAxes, fontsize=10, va='top')
#ax[4].text(0.97, 0.36, 'e', transform=ax[4].transAxes, fontsize=10, va='top')
#ax[5].text(0.97, 0.36, 'f', transform=ax[5].transAxes, fontsize=10, va='top')


#plt.rcParams["text.usetex"] =True


ax[0].text(0.965, 0.36, '$r$', transform=ax[0].transAxes, fontsize=10, va='top')
#ax[1].text(0.97, 0.36, '$\varphi$', transform=ax[1].transAxes, fontsize=10, va='top')
ax[1].text(0.965, 0.36, '$φ$', transform=ax[1].transAxes, fontsize=10, va='top')
ax[2].text(0.965, 0.36, '$B$', transform=ax[2].transAxes, fontsize=10, va='top')
ax[3].text(0.965, 0.36, '$h_1$', transform=ax[3].transAxes, fontsize=10, va='top')
ax[4].text(0.965, 0.36, '$h_2$', transform=ax[4].transAxes, fontsize=10, va='top')
ax[5].text(0.965, 0.36, '$h_3$', transform=ax[5].transAxes, fontsize=10, va='top')

#fig.text(0.5, 0.12, "Arc length along a dipole field line [10^3$ km]", fontsize=14, ha='center', va='center')


#ax[5].set(xlabel='Arc length along a dipole field line [$10^3$ km]')
ax[5].set(xlabel='Arc length along a dipole field line [km]')



ax[0].ticklabel_format(style='sci', axis='y', scilimits=(0,0))#, useMathText=True)
ax[1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))#, useMathText=True)
ax[2].ticklabel_format(style='sci', axis='y', scilimits=(0,0))#, useMathText=True)
ax[3].ticklabel_format(style='sci', axis='y', scilimits=(0,0))#, useMathText=True)
ax[4].ticklabel_format(style='sci', axis='y', scilimits=(0,0))#, useMathText=True)
ax[5].ticklabel_format(style='sci', axis='y', scilimits=(0,0))#, useMathText=True)


#ax[0].yaxis.offsetText.set_fontsize(9)
#ax[1].yaxis.offsetText.set_fontsize(9)
#ax[2].yaxis.offsetText.set_fontsize(9)
#ax[3].yaxis.offsetText.set_fontsize(9)
#ax[4].yaxis.offsetText.set_fontsize(9)
#ax[5].yaxis.offsetText.set_fontsize(9)



fig.savefig("relative_errors.png")
plt.show()


