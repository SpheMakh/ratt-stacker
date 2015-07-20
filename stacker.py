#!/usr/bin/env python

## Image plane continuum stacking script
## Sphesihle Makhathini <sphemakh@gmail.com>

import pyfits
from astLib.astWCS import WCS
import math
import numpy
import Tigger
from Tigger.Tools import gaussfitter2
import sys
import pylab
from scipy.ndimage import filters


def subregion(data, centre, radius):
    
    lx, ly = data.shape
    xstart = centre[1] - radius
    xend = centre[1] + radius
    ystart = centre[0] - radius
    yend = centre[0] + radius

    imslice = [slice(xstart if xstart>=0 else 0, xend if xend<= lx else lx),
               slice(ystart if ystart>=0 else 0, yend if yend<= ly else ly)]
    
    return data[imslice]
    
    

def stackem(imagename, positions, width, pixels=False, savefig=None, showfig=False):
    """ Perfom continuun stacking """
    
    with pyfits.open(imagename) as hdu:
        hdr = hdu[0].header
        _data = hdu[0].data

    ndim = hdr["NAXIS"]
    # I only want the image data in the first plane
    # May deal with multiple channels/polarizations later
    imslice = [slice(None)]*ndim
    imslice[:-2] = [0]*(ndim-2)

    data = _data[imslice]

    if not pixels:
        wcs = WCS(hdr, mode="pyfits")
        positions = [wcs.wcs2pix(*pos) for pos in positions]

    if width%2==1:
        width += 1
    stacked = numpy.zeros([width,width])

    for pos in positions:
        stacked += subregion(data, pos, width/2)

    # Get average stacked signal
    #stacked /= len(positions)
    # Fit gaussian to stacked emission

    # First create a threshhold mask
    std = stacked.std()
    smooth = filters.gaussian_filter(stacked,[2,2])
    mask = smooth>2*std
    #p = gaussfitter2.gaussfit(stacked)
    #peak = p[1]
    #cell = numpy.deg2rad(abs(hdr["CDELT1"]))
    #tot = peak*(2*math.pi)*abs(p[4]*p[5])

    if savefig or showfig:
        pylab.imshow(stacked)
        pylab.colorbar()
        #text = "Peak=%.4g  : Tot=%s"%(peak, tot)
        #pylab.title(text)
        x,y = [numpy.linspace(0,width,width)]*2
        xx,yy = numpy.meshgrid(x,y)
        #pylab.contour( gaussfitter2.twodgaussian(p,0,1,1)(xx,yy))

        if savefig:
            pylab.savefig(savefig)
        if showfig:
            pylab.show()
        else:
            pylab.clf()

    return stacked.max() 


if __name__ == "__main__":
    imagename = sys.argv[1]
    lsmname = sys.argv[2]
    
    import Tigger
    model = Tigger.load(lsmname)
    positions = [ map(numpy.rad2deg, [src.pos.ra, src.pos.dec]) for src in model.sources]

    stackem(imagename, positions, 50, showfig=True)
