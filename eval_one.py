#!/usr/bin/env python

'''
Synopsis:

Completely process a single proposal producing all of the outputs
needed to determine whether visits in the proposal should be bad
actors

Command line usage (if any):

    usage: eval_one.py foo.apt

Description:

    These routines parse an apt file to determine what exposures
    are contained in the file, and then retrieve data from
    IPAC that is intended to helpful in estimating the amount
    of persistence that the proposed observations will likely
    induce.

    If the apt file is not in the working directory, the rouitnes
    attempt to retrieve the apt file from the public repositiory
    of HST apt files.

    The routine creates an apt file for viewing the results, e.g.
    foo.html  (Subdiary files will be found in a directory created
    by the program, e.g. foo_dir.

Primary routines:

    doit

Notes:
    Apt files can be complicated and their detailed structure
    is not documented so it is possible that the routine will
    fail. 

History:

150805 ksl Coding begun

'''
import argparse
import os
import sys
import subprocess

from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table
from collections import OrderedDict
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import pylab

import read_apt
import actor
import persist_2mass

##### Relevant level for saturation
SATLEVEL = 75000

def read_table(filename='foo.txt'):
    '''
    Read a file using astropy.io.ascii and
    return this

    Notes:

    History:


    '''
    try:
        data=ascii.read(filename)
    except IOError:
        print('Error: file %s does not appear to exist' % filename)
        return

    print('Here are the column names:')

    print(data.colnames)

    return data

def image_plot(ra=66.76957,dec=26.10453,xfilter='F110W',exp=100.,rad=3,outroot='test_image', label=None):
    '''
    Create 2Mass images of the source region and analyze them for persistence.

    Notes

    Insofar as possilbe I will try to use Gabe's routines

    History

    150806  ksl Started work using Gabes persist_plot routine as a functional template
    '''

    print('Starting')

    x=1

    xxfilter=persist_2mass.get_2mass_filter(xfilter)

    print('Two mass filter',xxfilter)

    fits_file=outroot+'.fits'

    if (not os.path.exists(fits_file)) :
        print('Fetch 2MASS... %s' % fits_file)
        persist_2mass.montage_cutout(ra=ra,dec=dec,output=fits_file,survey='2MASS', band=xxfilter, verbose=False, clean=False)
    else:
        print('2MASS file %s already exits' % fits_file)


    im=fits.open(fits_file)

    im[0].data -= np.percentile(im[0].data, 5)
    plate_scale = np.abs(im[0].header['CDELT1']*3600)

    #### Scaled to WFC3/IR VEGAmag ZP
    scaled = persist_2mass.scale_image(fits_file, xfilter)
    print(xfilter)
    print("Scale for calc saturated pixels: ", scaled)

    time_to_saturate = SATLEVEL/scaled

    # Next step replaces negative values
    time_to_saturate[time_to_saturate < 0] = 1e8

    # Calcualate the total number of saturated pixels
    total_fluence = exp*scaled

    sat_levels = [1,5,10]
    NSATPIX = OrderedDict()
    for level in sat_levels:
        saturated = total_fluence > SATLEVEL*level
        NSATPIX[level] = saturated.sum() * plate_scale**2 / 0.128254**2

    print('NSATPIX',NSATPIX)

    #### Make figure
    fig=pylab.figure(2,(12,5.5))
    fig.clf()

    #### Image, saturation diagnostic, and color bar
    gs = gridspec.GridSpec(1, 3, width_ratios=[1,1,0.05])

    #### image itself
    ax = fig.add_subplot(gs[0])

    max = np.percentile(im[0].data[np.isfinite(im[0].data)], 95)
    im[0].data[~np.isfinite(im[0].data)] = 0
    ii=ax.imshow(persist_2mass.clipLog(im[0].data/max, scale=[-0.1, 10]), vmin=0, vmax=1, interpolation='Nearest', cmap='gray', origin='lower')
    pylab.colorbar(ii,fraction=0.046,pad=0.04)
    title=r'2MASS Image $-\ %s$' %(xxfilter.upper())

    ### Add 1' circle
    platescale = np.abs(im[0].header['CDELT1'])*3600
    R = int(np.round(1 * 60/platescale))
    sh = im[0].data.shape
    x0, y0 = sh[1]/2., sh[0]/2.
    xarr = np.arange(-R, R+1)
    yarr = np.sqrt(R**2-xarr**2)
    ax.plot(x0+xarr, y0+yarr, color='red', linewidth=2, alpha=0.5)
    ax.plot(x0+xarr, y0-yarr, color='red', linewidth=2, alpha=0.5)
    ax.set_xlim(0,sh[1]); ax.set_ylim(0,sh[1])

    ax.set_title(title)
    ax.set_xticklabels([]); ax.set_yticklabels([])
    
    #### Saturated imagge
    ax = fig.add_subplot(gs[1])

    # This is a depature from Gabe's routine
    saturated=total_fluence/SATLEVEL
    
    ### Discrete colormap
    ### from http://stackoverflow.com/questions/14777066/matplotlib-discrete-colorbar
    cmap = plt.cm.gist_earth_r
    cmaplist = [cmap(i) for i in range(cmap.N)]

    # create the new map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

    # define the bins and normalize
    bounds = np.linspace(-0.25,5.25,12)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    ii = ax.imshow(saturated, vmin=0, vmax=5, interpolation='Nearest', cmap='gist_earth_r', origin='lower', norm=norm)
    ax.plot(x0+xarr, y0+yarr, color='red', linewidth=2, alpha=0.5)
    ax.plot(x0+xarr, y0-yarr, color='red', linewidth=2, alpha=0.5)
    ax.set_xlim(0,sh[1]); ax.set_ylim(0,sh[1])

    pylab.colorbar(ii,fraction=0.046,pad=0.04, ticks=bounds+0.25, boundaries=bounds)
    title='Predicted Saturation'
    ax.set_title(title)

    if label is not None:
        ax.text(0.05, 0.95, label, ha='left', va='top', transform=ax.transAxes)

    ax.set_xticklabels([]); ax.set_yticklabels([])
    fig.tight_layout(pad=0.2)
    fig.savefig(outroot+'_image.png')

    return NSATPIX


def eval_one_main(filename,outroot):
    '''
    This is the main routine for carrying out
    a complete generation of products for
    estimating perisistence in the images

    Outroot is used a in creating the file names

    History

    150804  ksl Added to the routine
    150805  ksl Moved a version of the routine from actor to this
            routine in an attempt to add the IR images from
            Gabes routines
    150901  ksl Combined the old routine do_proposal into doit
            to make it easier to control where files are written
    150902  ksl Updated names of working files so that when identical
            exposures occur we do not retrieve data from IPAC
    151112  ksl	Updated so that the number of saturated stars and pixels
    		were logged in a way that could be compared with the badactor
		routine which is part of the persistence software

    '''

    if outroot=='':
        root=filename[0:filename.rindex('.')]
        html_file=root+'.html'
    else:
        root=outroot
        html_file=root+'.html'


    # Create an output directory for storing some fraction of the
    # output files

    directory=root+'_dir'

    if os.path.isdir(directory)==False:
        subprocess.call('mkdir %s' % directory, shell=True)
        print('Creating directory %s ' % directory)
    else:
        print('Directory %s already exists' % directory)


    table='%s/%s.txt' % (directory,root)

    # First see if the apt file is in the local directory
    if os.path.isfile(filename)==True:
        print('Using apt file %s in local directory' % filename)
    else:
        print('Fetching apt file %s' % filename)
        x=read_apt.fetch_apt(filename)
        if x==True:
            print('Fetch complete')
        else:
            print('Could not fetch apt file so quitting')
            return

    x=read_apt.read_apt_main(filename,table)

    if x==True:
        print('Analysis of apt file is complete')
    else:
        print('Analysis of apt file failed')
        return

    # Next line reads back the table created by read_apt
    xtable=read_table(table)
    root=table[0:table.rindex('.')]

    # At this point, root should be the dir and the name of the apt file
    logfile=root+'.log'
    g=open(logfile,'w')
    g.write('Processing %s\n' % table)
    g.write('Table %s\n' % table)

    # This is a list in which to store numbers of saturated pixels and stars
    summary=[]

    for one in xtable:
        print(one)
        saturated=[0,0,0,0,0,0]
        string='Header Visit %s Exp %03d %s' % (one['Visit'],one['ExpNo'],one['Target'])
        g.write('%s\n' % string)

        print('Test Config', one['Config'].count('IR'))
        if one['Config'].count('IR')==0:
            g.write('This is a UVIS exposure, so persistence is not an issue.\n')
        elif one['RA']!=-99.:
            name='%s_%s_%010.6f_%010.6f_%s_%04.0f' %(root,one['Target'],one['RA'],one['Dec'],one['Filter'],one['Exptime'])

            num_1x,num_5x,num_10x=actor.actor_main(ra=one['RA'], dec=one['Dec'], xfilter=one['Filter'], exp=one['Exptime'], rad=3, outroot=name)
            saturated[3]=num_1x
            saturated[4]=num_5x
            saturated[5]=num_10x

            exposure_log = 'This is an exposure at RA Dec = %f %f with the %s filter lasting %.2f seconds.' % (one['RA'],one['Dec'],one['Filter'],one['Exptime'])

            if (one['ScanRate'] is not None) & (one['ScanRate'] != 'None'):
                print('SCAN!!!!',  float(one['ScanRate']), one['Exptime'])
                scan_time = read_apt.get_read_time(one['SAMP-SEQ'], one['NSAMP'], APER=one['Aper'], scan_rate=None)
                exposure_log += '  <b> Note: It is a <span style="color:red">spatial scan</span> with trails %.1f pixels long (%.1f sec x %s"/sec).</b>' %(float(one['ScanRate'])*scan_time/0.128, scan_time, one['ScanRate'])

            print('XXX')
            exit

            g.write(exposure_log+'\n')
            string='Image %s.png' % name
            g.write('%s\n' % string)


            x=image_plot(ra=one['RA'],dec=one['Dec'],xfilter=one['Filter'],exp=one['Exptime'],rad=3,outroot=name, label='%s / %s / %.0f s' %(one['Target'], one['Filter'], one['Exptime']))
            string='Image %s_image.png' % name
            g.write('%s\n' % string)

            string='Saturated pixels in image: 1x  %d  5x %d  10x %d' % (x[1],x[5],x[10])
            saturated[0]=x[1]
            saturated[1]=x[5]
            saturated[2]=x[10]
            g.write('%s\n' % string)
            string='Number saturated stars: 1x  %d  5x %d 10x %d' % (num_1x,num_5x,num_10x)
            g.write('%s\n' % string)

        else:
            g.write('This exposure did not involve a target with a fixed position observed with WFC3/IR\n')
    summary.append(saturated)
    g.close()

    # Now add the associated with the nunber of saturated pixels to the table
    summary=np.array(summary)
    summary=np.transpose(summary)

    xtable['image_1x']=summary[0]
    xtable['image_5x']=summary[1]
    xtable['image_10x']=summary[2]
    xtable['star_1x']=summary[3]
    xtable['star_5x']=summary[4]
    xtable['star_10x']=summary[5]

    xtable['image_1x'].format='5.1f'
    xtable['image_5x'].format='5.1f'
    xtable['image_10x'].format='5.1f'
    xtable['star_1x'].format='5.1f'
    xtable['star_5x'].format='5.1f'
    xtable['star_10x'].format='5.1f'

    xtable.write(table,format='ascii.fixed_width_two_line',overwrite=True)

    print('The name of the output html file is ',html_file)
    actor.make_html(logfile,outfile=html_file)

    return

# -----------------------------------------------------------------------------
# For command line execution
# -----------------------------------------------------------------------------

def parse_args():
    """Parses command line arguments.
    Returns
    -------
        args : argparse.Namespace object
            An argparse object containing all of the added arguments.
    """

    # Create help string
    filename_help = 'Apt file name. The default is 15101.apt.'
    outroot_help = 'The outroot name. The default is '' to set it to the apt root.'

    # Add time arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('-filename --f',
        dest='filename',
        action='store',
        type=str,
        required=False,
        help=filename_help,
        default="15101.apt")
    parser.add_argument('-outroot --out',
        dest='outroot',
        action='store',
        type=str,
        required=False,
        help=outroot_help, 
        default='')

    # Parse args
    args = parser.parse_args()

    return args
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":

    args = parse_args()

    eval_one_main(args.filename, args.outroot)
