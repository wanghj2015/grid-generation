import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri


from mpl_toolkits import axes_grid1
from matplotlib import ticker
from matplotlib import colors

from matplotlib.ticker import ScalarFormatter

import matplotlib.ticker

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


def add_colorbar(im, aspect=15, pad_fraction=0.5, **kwargs):
    """Add a vertical color bar to an image plot."""
    divider = axes_grid1.make_axes_locatable(im.axes)
    width = axes_grid1.axes_size.AxesY(im.axes, aspect=1./aspect)
    pad = axes_grid1.axes_size.Fraction(pad_fraction, width)
    current_ax = plt.gca()
    cax = divider.append_axes("right", size=width, pad=pad)
    plt.sca(current_ax)
    cbar = im.axes.figure.colorbar(im, cax=cax, **kwargs)
    cbar.ax.yaxis.set_major_locator(ticker.AutoLocator())
    cbar.ax.yaxis.set_minor_locator(ticker.AutoLocator())
    cbar.ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True, useOffset=True))
    cbar.ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    cbar.ax.ticklabel_format(style='sci', scilimits=(0, 0))
    return cbar


# read in 'general coordinates' (gc) basis vectors

f = FortranFile('../run/igrf/general/fort.20', 'r')

[nlp,nmp,npts] = f.read_ints(np.int32)

print (nlp,nmp,npts)


# skip: could get klm, lpts here
#f.read_ints(np.int32)

klm  = np.zeros((nlp,nmp), dtype=np.int32)
lpts = np.zeros((npts,nmp), dtype=np.int32)
x = f.read_ints(np.int32)

klm [:,0] = x[0:nlp]
lpts[:,0] = x[nlp:nlp+npts]

#print (' klm  = ', klm)
#print (' lpts = ', lpts)


# read in again
r = []; theta = []; b = []; s = []
for m in range(0, nmp):
    for l in range(0, nlp):
        x = f.read_reals(dtype='float64')
        r.append(x)
        x = f.read_reals(dtype='float64')
        theta.append(x)
        x = f.read_reals(dtype='float64')
        # skip 'phi'
        x = f.read_reals(dtype='float64')
        b.append(x)
        x = f.read_reals(dtype='float64')
        s.append(x)

f.close()

#print (s)

stacked, idx = stack_ragged(r)
r = np.array(stacked) * 1.0e-3

stacked, idx = stack_ragged(theta)
theta = np.array(stacked)

#print (r.shape)

#x = (r-0)*np.sin(theta)
#y = (r-0)*np.cos(theta)


f = FortranFile('../run/igrf/general/fort.25', 'r')

ang12 = []; ang13 = []; ang23  = []; ratio = []
for m in range(0, nmp):
    for l in range(0, nlp):
        x = f.read_reals(dtype='float64')
        ang12.append(x)
        x = f.read_reals(dtype='float64')
        ang13.append(x)
        x = f.read_reals(dtype='float64')
        ang23.append(x)
        x = f.read_reals(dtype='float64')
        ratio.append(x)

f.close()

stacked, idx = stack_ragged(ang12)
ang12 = np.array(stacked)

stacked, idx = stack_ragged(ang13)
ang13 = np.array(stacked)

stacked, idx = stack_ragged(ang23)
ang23 = np.array(stacked)

stacked, idx = stack_ragged(ratio)
ratio = np.array(stacked)


#idx = np.array(idx)

#print (' idx = ', np.shape(idx), idx)

# plot basis angles

fig, ax = plt.subplots()

ax.set_aspect('equal')

r_earth = 6371.2

x = (r-0)*np.sin(theta)
y = (r-0)*np.cos(theta)


triang = tri.Triangulation(x, y)


# https://stackoverflow.com/questions/42426095/matplotlib-contour-contourf-of-concave-non-gridded-data

def apply_mask(triang, alpha=0.4):
    # Mask triangles with sidelength bigger some alpha
    triangles = triang.triangles
    # Mask off unwanted triangles.
    xtri = x[triangles] - np.roll(x[triangles], 1, axis=1)
    ytri = y[triangles] - np.roll(y[triangles], 1, axis=1)
    maxi = np.max(np.sqrt(xtri**2 + ytri**2), axis=1)
    # apply masking
    triang.set_mask(maxi > alpha)

apply_mask(triang, alpha=1.e3)


#plt.rcParams.update({'font.size': 5})
#plt.tick_params(labelsize=5)
#plt.tight_layout()

#plt.rc('xtick',labelsize=5)
#plt.rc('ytick',labelsize=5)


print (' ratio = ', min(ratio), max(ratio))


minval = min(ratio)
minval = 0.94
maxval = max(ratio)
maxval = 1.08
step   = (maxval - minval) / 10
step   = 0.02

levels = np.arange(minval,maxval,step)


im = ax.tricontourf(triang, ratio, levels=levels)#, cmap="jet")#, extend="both")

#cb = fig.colorbar(im, ax=ax, orientation='vertical')
#cb.formatter.set_powerlimits((0, 0))
#cb.update_ticks()


cbar = add_colorbar(im)


# plot grid points too
#plt.scatter(x, y, s=0.3)


# plot apex points too?


# x & y labeling:
#x: distance from earth's center
#y: distance from south to north, centered at the earth's center


ax.xaxis.offsetText.set_visible(False)
ax.yaxis.offsetText.set_visible(False)

plt.xlabel("Distance from earth's center [$10^3$ km]")
plt.ylabel("Distance from earth's center [$10^3$ km]")


plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useMathText=True)


class OOMFormatter(matplotlib.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_order_of_magnitude(self):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin=None, vmax=None):
        self.format = self.fformat
        if self._useMathText:
            self.format = r'$\mathdefault{%s}$' % self.format

ax.xaxis.set_major_formatter(OOMFormatter(3, "%1.0f"))


thetaxx = np.linspace(0, np.pi, 100)
thetaxx = np.linspace(np.pi/4, np.pi*3/4, 100)

x1 = 6371.2*np.sin(thetaxx)
y1 = 6371.2*np.cos(thetaxx)

#plt.plot(x1, y1)

x2 = (6371.2+90)*np.sin(thetaxx)
y2 = (6371.2+90)*np.cos(thetaxx)

#plt.plot(x2, y2)




#for l in range(0, nlp-2, 45):
#for l in range(16, nlp-2, 6):
for l in [12, 22, 28, 32, 36, 40]:
    rx = r[idx[l]:idx[l+1]]
    thetax = theta[idx[l]:idx[l+1]]
    xx = rx * np.sin(thetax)
    yy = rx * np.cos(thetax)
    plt.plot(xx, yy)

#print (r[0:idx[0]])


fig.savefig("euler_ratio.png", dpi=150, bbox_inches='tight')


plt.show()


