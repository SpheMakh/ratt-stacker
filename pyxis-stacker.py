## Continuum stacking pipeline

import Pyxis
from Pyxis.ModSupport import *
import ms
import mqt
import lsm
import stacker

import numpy
import math
import Tigger
import glob
import scipy.optimize as optimize
import pylab


def stackem(start=1, stop=100):
    
    v.MS = MS or II("${LSM:BASE}.MS")
    
    DO_STEP = lambda step: step>=float(start) and step<=float(stop)

    if not exists(MS) or MS_REDO:
        ms.create_empty_ms(tel=OBSERVATORY, pos=ANTENNAS, direction=DIRECTION, synthesis=SYNTHESIS, dtime=DTIME, 
                           freq0=FREQ0, dfreq=DFREQ)
        reload(ms)
        start = 1

    if DO_STEP(1):
        simsky()
        im.make_image()
     
    if DO_STEP(2):
        model = Tigger.load(LSM)
        radec = []
        flux = 0
        for src in model.sources:
            ra = src.pos.ra
            dec = src.pos.dec
            radec.append( map(numpy.rad2deg, [ra,dec]))
            flux += src.flux.I

        flux /=len(model.sources)
        
        sflux = stacker.stackem(im.DIRTY_IMAGE, radec, width=WIDTH)
        _add("%.4g %.4g"%(flux, sflux), STACKFILE, delimeter="\n")


def stackem_all(indir, lsmext="lsm.html", prefix="qaz-ratt-stacker", **kw):

    v.LSM_List = glob.glob("%s/*%s"%(indir, lsmext))
   
    #x.sh("rm -f $STACKFILE") 
    #pper("LSM", lambda: stackem(**kw))
    data = numpy.genfromtxt(STACKFILE)
    xx,yy = data[:,0],data[:,1]

    line = lambda x, m, c: m*x +c
    m = (yy[1]-yy[0])/(xx[1]-xx[0])
    parms, cov = optimize.curve_fit(line, xx, yy, [m,0])
    err = numpy.diagonal(cov)**2
    m,c = parms

    pylab.scatter(xx, yy)
    pylab.plot(xx, line(xx,m,c), "r-")

    text = "m=%.3g $\pm$  %.3g \n c=%.3g $\pm$ %.3g"%(m, err[0], c, err[1])
    pylab.text(xx.min(), yy.max(), text,size=15)
    pylab.savefig(II("${OUTDIR>/}${prefix}.png"))
            
    
def get_random_positions(imagename, npos, border=0.1):

    with pyfits.open(imagename) as hdu:
        hdr = hdu[0].header
    
    nx,ny = hdr["NAXIS1"], hdr["NAXIS2"]
    ra0,dec0 = hdr["CRVAL1"], hdr["CRVAL2"]
    cell = abs(hdr["CDELT1"])
    
    dra = nx*cell*(1-border)
    ddec = ny*cell*(1-border)
    ra = numpy.random.random(npos)*dra -(ra0 - dra/2) 
    dec = numpy.random.random(npos)*ddec -(dec0 - ddec/2) 

    return map(list, [ra,dec])

def simsky(msname="$MS", lsmname="$LSM", column="$COLUMN",
           tdlconf="$TDLCONF", tdlsec="$SIMSEC", addnoise=True,
           noise=0, sefd=0, recenter=True, options={} ,args=[],**kw):
    """ 
    Simulates visibilities into a MS.
    msname : MS name
    lsmname : LSM name
    column : Column to simulate visibilities into
    tdlconf : Meqtrees TDL configuration profiles file (required to run MeqTrees pipeliner) 
    tdlsec : Section to execute in tdlconf
    noise : Visibility noise to add to simulation.
    args, kw : extra arguments to pass to the MeqTrees pipeliner
    """
    msname,lsmname,column,tdlsec,tdlconf = interpolate_locals('msname lsmname column tdlsec tdlconf')

    # recenter LSM if required    
    if recenter:
        x.sh('tigger-convert --recenter=$DIRECTION $LSM -f')

    args = ["${ms.MS_TDL} ${lsm.LSM_TDL}"] + list(args)

    options['ms_sel.output_column'] = column

    if addnoise:
        sefd = sefd or SEFD
        options['noise_stddev'] = noise or compute_vis_noise(sefd)
    options.update(kw) # extra keyword args get preference
    mqt.run(TURBO_SIM,job='_tdl_job_1_simulate_MS', config=tdlconf,section=tdlsec,options=options,args=args)


def compute_vis_noise (sefd):
    """Computes nominal per-visibility noise"""
    #sefd = sefd or SEFD
    tab = ms.ms()
    spwtab = ms.ms(subtable="SPECTRAL_WINDOW")

    freq0 = spwtab.getcol("CHAN_FREQ")[ms.SPWID,0]
    wavelength = 300e+6/freq0
    bw = spwtab.getcol("CHAN_WIDTH")[ms.SPWID,0]
    dt = tab.getcol("EXPOSURE",0,1)[0]
    dtf = (tab.getcol("TIME",tab.nrows()-1,1)-tab.getcol("TIME",0,1))[0]

    # close tables properly, else the calls below will hang waiting for a lock...
    tab.close()
    spwtab.close()

    info(">>> $MS freq %.2f MHz (lambda=%.2fm), bandwidth %.2g kHz, %.2fs integrations, %.2fh synthesis"%(freq0*1e-6,wavelength,bw*1e-3,dt,dtf/3600))
    noise = sefd/math.sqrt(2*bw*dt)
    info(">>> SEFD of %.2f Jy gives per-visibility noise of %.2f mJy"%(sefd,noise*1000))

    return noise


def _add(addfile,filename, delimeter=","):
    """ Keeps track of Files when running multiple threads.
        The files names are stored into a file which can be specified 
        via filename. The files can then be rertieved as a python list using the function get_list().
    """

    filename = interpolate_locals('filename')

    
    if exists(filename):
        try :
            ms_std = open(filename,"r")
            if delimeter == "\n":
                lines = ms_std.readlines()
            else:
                lines = ms_std.readline().split(delimeter)
        except IOError:
            lines = []
        ms_std.close()
    else:
        lines = []
    
    for line in lines:
        if line:
            pass
        else:
            lines.remove(line)

    file_std = open(filename,"w")
    lines.append(addfile)
    file_std.write(delimeter.join(lines))
    file_std.close()
    info("New list : $lines")


def get_list(filename,delimeter=","):
    """ See help for _add"""

    filename = interpolate_locals('filename')
    if not exists(filename):
        return False

  
    with open(filename) as ms_std:
        ms_std = open(filename)
        if delimeter == "\n":
            oulist = ms_std.redalines()
        else:
            outlist = ms_std.readline().split(delimeter)

    info('Found %d files.'%(len(outlist)))

    return outlist

