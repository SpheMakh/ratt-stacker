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
        
        sflux = stacker.stackem(im.DIRTY_IMAGE, radec, width=WIDTH, showfig=True)
        _add("%.4g %.4g"%(flux, sflux), STACKFILE, delimeter="\n")


def stackem_all(indir, lsmext="lsm.html"):

    v.LSM_List = glob.glob("%s/*%s"%(indir, lsmext))
    
    pper("LSM", stackem)
    data = numpy.genfromtxt(get_list(STACKEM))
    print data
    
            
    
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
    mqt.run(TURBO_SIM,job='_tdl_job_1_simulate_MS',config=tdlconf,section=tdlsec,options=options,args=args)


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

  addfile,filename = interpolate_locals('addfile filename')

  try :
    ms_std = open(filename,"r")
    line = ms_std.readline().split("\n")[0]
    ms_std.close()
  except IOError: line = ''

  file_std = open(filename,"w")
  line += "%s%s"%(delimeter, addfile) if len(line.split()) > 0 else addfile
  file_std.write(line)
  file_std.close()
  info("New list : $line")


def get_list(filename):
    """ See help for _add"""

    filename = interpolate_locals('filename')
    if not exists(filename):
        return False

    with open(filename) as ms_std:
      ms_std = open(filename)
      mslist = ms_std.readline().split('\n')[0].split(',')
    info('Found %d files.'%(len(mslist)))

    return mslist

