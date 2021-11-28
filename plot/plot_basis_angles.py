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
    cbar.ax.yaxis.set_offset_position('left')
    return cbar



# read in 'general coordinates' (gc) basis vectors

f = FortranFile('../run/igrf/general/fort.20', 'r')

[nlp,nmp,npts] = f.read_ints(np.int32)

print (nlp,nmp,npts)

# skip: could get klm, lpts here
f.read_ints(np.int32)

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


# plot basis angles

fig, ax = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(10, 5))

#ax.set_aspect('equal')

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


ang12 = ang12 - 90
ang13 = ang13 - 90
ang23 = ang23 - 90


print (' beta12 = ', '{:.4e}'.format(min(ang12)), '{:.4e}'.format(max(ang12)))
print (' beta13 = ', '{:.4e}'.format(min(ang13)), '{:.4e}'.format(max(ang13)))
print (' beta23 = ', '{:.4f}'.format(min(ang23)), '{:.4f}'.format(max(ang23)))


minval = min(ang12)
minval = -4.5e-2
minval = -4e-2
maxval = max(ang12)
maxval = 3.5e-2
maxval = 3e-2
step   = (maxval - minval) / 10
step   = 1.0e-2
step   = 1e-2
levels = np.arange(minval,maxval,step)

im = ax[0].tricontourf(triang, ang12, levels=levels)#, cmap="jet")#, extend="both")

#cb = plt.colorbar(im, ax=ax[0], fraction=0.05, pad=0.05, orientation='vertical')
#cb.formatter.set_powerlimits((0, 0))
#cb.update_ticks()

cb = add_colorbar(im)

ax[0].set_aspect('equal')
ax[0].ticklabel_format(style='sci', axis='both', scilimits=(0,0), useMathText=True)

ax[0].yaxis.set_offset_position('left')

#plt.ylabel("distance from earth's center")

#ax[0].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

#plt.scatter(x, y, s=0.2)

#plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))

#ax[0].set_yticklabels([])


minval = min(ang13)
minval = -1.0e-3
maxval = max(ang13)
maxval = 2.5e-3
step   = (maxval - minval) / 10
step   = 0.5e-3
levels = np.arange(minval,maxval,step)

im = ax[1].tricontourf(triang, ang13, levels=levels)#, cmap="jet")#, extend="both")

#cb = plt.colorbar(im, ax=ax[1], fraction=0.05, pad=0.05, orientation='vertical')
#cb.formatter.set_powerlimits((0, 0))
#cb.update_ticks()

cb = add_colorbar(im)

ax[1].set_aspect('equal')
#ax[1].set_yticklabels([])
ax[1].ticklabel_format(style='sci', axis='x', scilimits=(0,0), useMathText=True)

#ax[1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))


# Extend ScalarFormatter
#class MyScalarFormatter(ScalarFormatter):
#    # Override '_set_format' with your own
#    def _set_format(self):
#        #self.format = '%.2f'  # Show 2 decimals
#        self.format = '%.1f'  # Show 2 decimals

#custom_formatter = MyScalarFormatter(useMathText=True)
#ax[1].xaxis.set_major_formatter(custom_formatter)


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

ax[1].xaxis.set_major_formatter(OOMFormatter(3, "%1.0f"))




minval = min(ang23)
minval = -6 
maxval = max(ang23)
maxval = 14
step   = (maxval - minval) / 10
step   = 2
levels = np.arange(minval,maxval,step)

im = ax[2].tricontourf(triang, ang23, levels=levels)#, cmap="jet")#, extend="both")

cb = plt.colorbar(im, ax=ax[2], fraction=0.05, pad=0.05, orientation='vertical')

#cb.formatter.set_powerlimits((0, 0))
#cb.update_ticks()

ax[2].set_aspect('equal')
#ax[2].set_yticklabels([])


# plot grid points too
#plt.scatter(x, y, s=0.2)



ax[2].ticklabel_format(style='sci', axis='x', scilimits=(0,0), useMathText=True)



thetaxx = np.linspace(0, np.pi, 100)
thetaxx = np.linspace(np.pi/4, np.pi*3/4, 100)

x1 = 6371.2*np.sin(thetaxx)
y1 = 6371.2*np.cos(thetaxx)

#plt.plot(x1, y1)

x2 = (6371.2+90)*np.sin(thetaxx)
y2 = (6371.2+90)*np.cos(thetaxx)

#plt.plot(x2, y2)


plt.tight_layout()



#fig.text(0.5, 0.04, "distance from earth's center [km]", ha='center')
#fig.text(0.04, 0.5, "distance from earth's center [km]", va='center', rotation='vertical')

#plt.xlabel("distance from earth's center [km]")
#plt.ylabel("distance from earth's center [km]")

# Set common labels

ax[0].xaxis.offsetText.set_visible(False)
ax[1].xaxis.offsetText.set_visible(False)
ax[2].xaxis.offsetText.set_visible(False)

ax[0].yaxis.offsetText.set_visible(False)
ax[1].yaxis.offsetText.set_visible(False)
ax[2].yaxis.offsetText.set_visible(False)

#plt.xlabel("Distance from earth's center [$10^4$ km]")
#plt.ylabel("Distance from earth's center [$10^3$ km]")

fig.text(0.5, 0.12, "Distance from earth's center [$10^3$ km]", fontsize=14, ha='center', va='center')
fig.text(0.01, 0.5, "Distance from earth's center [$10^3$ km]", fontsize=12, ha='center', va='center', rotation='vertical')


#ax[0].set_title("$\mathbf{e}^1 x \mathbf{e}^2$")
#ax[1].set_title("$\mathbf{e}^1 x \mathbf{e}^3$")
#ax[2].set_title("$\mathbf{e}^2 x \mathbf{e}^3$")


#ax[0].legend(loc="upper right")
#ax[1].legend(loc="upper right")
#ax[2].legend(loc="upper right")


ax[0].text(0.90, 0.96, 'a', transform=ax[0].transAxes, fontsize=12, va='top')
ax[1].text(0.90, 0.96, 'b', transform=ax[1].transAxes, fontsize=12, va='top')
ax[2].text(0.90, 0.96, 'c', transform=ax[2].transAxes, fontsize=12, va='top')


# increase font size
#plt.rcParams.update({'font.size': 14})


#fig.savefig("basis_angles.png", dpi=150, bbox_inches='tight')
fig.savefig("basis_angles.png", dpi=fig.dpi, bbox_inches='tight')


plt.show()


