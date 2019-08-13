#!/usr/bin/env python 

'''
Make diagnostic persistence plots based on 2MASS cutouts for WFC3/IR exposures
in an APT file.
'''

import os

import astropy.table
import astropy.io.fits as pyfits
from collections import OrderedDict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import montage_wrapper as montage

import parse_apt

SATLEVEL = 75000
def make_plots(apt_file='gbb_snap.apt', program=None):
    '''    
    Make plots for all IR exposures in an APT file
    '''
    pid = parse_apt.parse_apt_main(apt_file=apt_file, program=program)

    t = astropy.table.Table()
    report = t.read('%s.IR.report' %(pid), format='ascii.commented_header')

    idx = np.arange(len(report))
    for i in idx:
        persist_plot(report[i])

def persist_plot(table_line=None, force_2mass=False, *args, **kwargs):
    '''
    Make diagnostic plot for a given exposure in an APT file.

    `table_line` is a line of the the table produced by the "parse_apt" script

    If `force_2mass` is set, force generate the 2MASS thumb.  Otherwise won't be updated
    if the FITS filename is unchanged (PID_Visit_Target_filter.fits).

    Also accepts the following keywords, with the indicated default values:

    PID=14095; Visit='78'; Exp=1; Target='NGC1068'; RA='02:42:40.7700'; Dec='-00:00:47.80'
    POSX=0.0; POSY=0.0; Config='WFC3/IR'; Aper='IR-UVIS'; Filter='F110W';
    SAMPSEQ='STEP50'; NSAMP=8
    '''
    default = {'PID': 14095, 'Visit': '78', 'Exp': 1, 'Target': 'NGC1068', 'RA': '02:42:40.7700', 'Dec': '-00:00:47.80', 'POSX': 0.0, 'POSY': 0.0, 'Config': 'WFC3/IR', 'Aper': 'IR-UVIS', 'Filter': 'F110W', 'SAMPSEQ': 'STEP50', 'NSAMP': 8}

    if table_line is None:
        report_line = default
    else:
        report_line = OrderedDict()
        for col in table_line.colnames:
            report_line[col] = table_line[col]

    for key in list(kwargs.keys()):
        report_line[key] = kwargs[key]

    twomass_filter = get_2mass_filter(report_line['Filter'])

    fits_file = '%d_%s_%s_%s.fits' %(report_line['PID'], str(report_line['Visit']), report_line['Target'], twomass_filter)

    if (not os.path.exists(fits_file)) | (force_2mass):
        print('Fetch 2MASS...')
        montage_cutout(ra=report_line['RA'], dec=report_line['Dec'], output=fits_file, survey='2MASS', band=twomass_filter, verbose=False, clean=True)

    #### 2MASS image
    im = pyfits.open(fits_file)
    im[0].data -= np.percentile(im[0].data, 5)
    plate_scale = np.abs(im[0].header['CDELT1']*3600)

    #### Scaled to WFC3/IR VEGAmag ZP
    scaled = scale_image(fits_file, report_line['Filter'])

    time_to_saturate = SATLEVEL/scaled
    time_to_saturate[time_to_saturate < 0] = 1e8

    samples = get_read_times(SAMPSEQ=report_line['SAMPSEQ'], NSAMP=report_line['NSAMP'])
    EXPTIME = samples[-1]

    #NSAMP = report_line['NSAMP']
    sat_read = time_to_saturate*0+10000
    for i in range(report_line['NSAMP'])[::-1]:
        sat = time_to_saturate < samples[i]
        sat_read[sat] = i

    ### Calculate number of saturated pixels
    total_fluence = EXPTIME*scaled

    sat_levels = [1,5,10]
    NSATPIX = OrderedDict()
    for level in sat_levels:
        saturated = total_fluence > SATLEVEL*level
        NSATPIX[level] = saturated.sum() * plate_scale**2 / 0.128254**2

    #### Make figure
    xs = 6
    fig = plt.figure(figsize=[xs,xs/2.*1.05])

    #### Image, saturation diagnostic, and color bar
    gs = gridspec.GridSpec(1, 3, width_ratios=[1,1,0.05])

    #### image itself
    ax = fig.add_subplot(gs[0])

    max = np.percentile(im[0].data[np.isfinite(im[0].data)], 95)
    im[0].data[~np.isfinite(im[0].data)] = 0
    ax.imshow(clipLog(im[0].data/max, scale=[-0.1, 10]), vmin=0, vmax=1, interpolation='Nearest', cmap='gray', origin='lower')
    ax.set_title('2MASS, %s' %(twomass_filter.upper()))
    ax.text(0.5, 0.98, report_line['Target'], fontsize=10, ha='center', va='top', backgroundcolor='white', transform=ax.transAxes)

    #### Read saturation
    ax = fig.add_subplot(gs[1])
    ii = ax.imshow(sat_read, vmin=0, vmax=report_line['NSAMP'], interpolation='Nearest', cmap='spectral', origin='lower')

    title = '%d visit: %s exp: %02d' %(report_line['PID'], str(report_line['Visit']), report_line['Exp'])
    ax.set_title(title)
    ax.text(0.5, 0.95, '%s  %s /%2d (%d s)' %(report_line['Filter'], report_line['SAMPSEQ'], report_line['NSAMP'], EXPTIME), ha='center', va='top', fontsize=10, transform=ax.transAxes)

    di = 0.05
    for ilev, level in enumerate(sat_levels):
        ax.text(0.95, 0.1+ilev*di, '%3dx: %6d (%5.1f' %(level, NSATPIX[level], NSATPIX[level]/1014**2*100) + '%)', ha='right', va='bottom', transform=ax.transAxes, fontsize=8)

    for ax in fig.axes:
        ax.set_xticklabels([]); ax.set_yticklabels([])

    #### Colorbar
    cax = fig.add_subplot(gs[2])

    cb = fig.colorbar(ii, cax=cax)
    cb.set_ticks(np.arange(report_line['NSAMP'], dtype=int)); cb.set_ticklabels(np.arange(report_line['NSAMP'], dtype=int)+1)
    cb.set_label(report_line['SAMPSEQ'])
    cax.tick_params(axis='both', which='major', labelsize=7)

    fig.tight_layout(pad=0.2)

    fig_file = '%d_%s_%02d_%s_%s.png' %(report_line['PID'], str(report_line['Visit']), report_line['Exp'], report_line['Filter'], report_line['Target'])

    fig.savefig(fig_file)
    plt.close()

    fp = open(fig_file.replace('.png', '.dat'), 'w')
    sat_str = ' '.join(['%d' %(val) for val in list(NSATPIX.values())])
    fp.write('%s %s %s %02d %s %s %d %s\n' %(fig_file, report_line['PID'], str(report_line['Visit']), report_line['Exp'], report_line['Filter'], report_line['SAMPSEQ'], report_line['NSAMP'], sat_str))
    fp.close()

    print(fits_file, fig_file)

def compare_flt_2mass(flt_file='../RAW/ib6wr9b3q_flt.fits'):
    '''
    Compare FLT to scaled 2MASS prediction
    '''
    import scipy.ndimage as nd

    im = pyfits.open(flt_file)
    h = im[0].header
    h1 = im[1].header

    twomass_filter = get_2mass_filter(h['FILTER'])

    print('Fetch 2MASS...')
    montage_cutout(ra=h1['CRVAL1'], dec=h1['CRVAL2'], output='tmp_2mass.fits', survey='2MASS', band=twomass_filter, verbose=False, clean=False, width=2.2, rotation=270-h['PA_V3']+44.6)
    scaled = scale_image('tmp_2mass.fits', h['FILTER'],twomass_filter)
    xs = 6
    fig = plt.figure(figsize=[xs,xs/2.*1.13])

    #### 2MASS image
    ax = fig.add_subplot(121)

    max = np.percentile(scaled[np.isfinite(scaled)], 85)
    scaled[~np.isfinite(scaled)] = 0

    ax.imshow(clipLog(scaled/max, scale=[-0.1, 10]), vmin=0, vmax=1, interpolation='Nearest', cmap='gray', origin='lower')
    ax.set_title('2MASS, %s' %(twomass_filter.upper()))
    
    #### FLT image
    ax = fig.add_subplot(122)
    flt_data = im['SCI'].data*1.
    mask = im['DQ'].data == 0
    flt_data -= np.percentile(flt_data[mask], 5)

    sm = nd.gaussian_filter(flt_data[::2,::2], 0.5)
    ax.imshow(clipLog(sm/max, scale=[-0.1, 10]), vmin=0, vmax=1, interpolation='Nearest', cmap='gray', origin='lower')
    ax.set_title(os.path.basename(flt_file))

    for ax in fig.axes:
        ax.set_xticklabels([]); ax.set_yticklabels([])

    fig.tight_layout(pad=0.1)

    out_file = os.path.basename(flt_file).replace('.fits', '_%s.png' %(twomass_filter))
    print(out_file)
    fig.savefig(out_file)
    plt.close()

def get_2mass_filter(filter='F110W'):
    '''
    Choose 2MASS J/H filters for comparision to WFC3
    '''
    j_filters = ['F105W', 'F110W', 'F125W', 'F098M', 'F127M', 'F139M', 'F126N', 'F128N', 'F130N', 'F132N']
    if filter in j_filters:
        twomass_filter = 'j'
    else:
        twomass_filter = 'h'

    return twomass_filter

def get_read_times(SAMPSEQ='STEP50', NSAMP=8):

    samples = {'RAPID': 2.932+np.arange(15)*2.932,
               'STEP25': np.array([2.93, 5.89, 8.80, 11.73, 24.23, 49.23, 74.23, 99.23, 124.23, 149.23, 174.23, 199.23, 224.23, 249.23, 274.23]),
               'STEP50': np.array([2.93, 5.87, 8.80, 11.73, 24.23, 49.23, 99.23, 149.23, 199.23, 249.23, 299.23, 349.23, 399.23, 449.23, 499.23]),
               'STEP100': np.array([2.93, 5.8, 8.8, 11.7, 24.23, 49, 99, 199, 299, 399, 499, 599, 699, 799, 899]),
               'STEP200': np.array([2.93, 5.8, 8.8, 11.7, 24.23, 49, 99, 199, 399, 599, 799, 999, 1199, 1399, 1599]),
               'STEP400': np.array([2.93, 5.8, 8.8, 11.7, 24.23, 49, 99, 199, 399, 799, 1199, 1599, 1999, 2399, 2799]),
               'SPARS5': 2.93+np.arange(15)*5,
               'SPARS10': 2.93+np.arange(15)*10,
               'SPARS25': 2.93+np.arange(15)*25,
               'SPARS50': 2.93+np.arange(15)*50,
               'SPARS100': 2.93+np.arange(15)*100,
               'SPARS200': 2.93+np.arange(15)*200}

    return samples[SAMPSEQ][:NSAMP]

def scale_image(img='twomass-j.fits', filter='F110W',twomass_filter='j'):
    '''
    Scale a 2MASS cutout to a WFC3 filter
    '''
    ### ABMAG
    #WFC3_ZPS = {'F110W':26.822}

    ### Vegamag, 2MASS MAG ZP in Vegamag
    #WFC3_ZPS = {'F110W':26.0628, 'F128N':21.9355, 'F105W':25.6236, 'F127M':23.6799, 'F139M':23.400}
    ### VEGAmag zeropoints, 2MASS MAGZP given in VEGAmag
    WFC3_ZPS = {'F105W': 25.6236,
                'F110W': 26.0628,
                'F125W': 25.329,
                'F140W': 25.376,
                'F160W': 24.695,
                'F098M': 25.106,
                'F127M': 23.678,
                'F139M': 23.401,
                'F153M': 23.210,
                'F126N': 21.940,
                'F128N': 21.936,
                'F130N': 22.014,
                'F132N': 21.945,
                'F164N': 21.534,
                'F167N': 21.595,
                'G102' : 19.6,
                'G141' : 19.4}

    im = pyfits.open(img)
    print(im[0].header)

    #### Scale factor includes the ZPs and the square of the plate scales (WFC3/IR hard-coded)
    try:
        ZP = im[0].header['MAGZP'] #+ 0.885
    except:
        print('Error: There is no MAGZP value in the header. Using hard-coded value.')
        h = im[0].header
        c = h.cards
        for i in range(len(c)):
            if 'MAGZP' in ' '.join(np.cast[str](c[i])):
                break
        if twomass_filter == 'j':
            ZP = 20.905 #average of the jband read_2 from 2MASS
        elif twomass_filter == 'h':
            ZP = 20.3243041848514 

    plate_scale = np.abs(im[0].header['CDELT1']*3600)
    try:
        scl = 10**(-0.4*(ZP-WFC3_ZPS[filter]))/(plate_scale/0.128254)**2
    except KeyError:
        print('Error: scale_image: No filter %s' % filter)

    #### With simple background subtraction
    im[0].data = (im[0].data - np.percentile(im[0].data, 5))*scl

    return im[0].data

def clipLog(im, lexp=1000, cmap=[-1.4914, 0.6273], scale=[-0.1,10]):
    '''
    DS9-like log scaling
    '''
    contrast, bias = cmap
    clip = (np.clip(im, scale[0], scale[1])-scale[0])/(scale[1]-scale[0])
    clip_log = np.clip((np.log10(lexp*clip+1)/np.log10(lexp)-bias)*contrast+0.5, 0, 1)

    return clip_log

def montage_cutout(ra=172.1305000, dec=+58.5615833, target=None, pix_size=None, width=3, survey='2MASS', band='j', output='montage.fits', verbose=True, clean=True, rotation=None):
    '''
    Use Montage to automatically generate cutouts. https://github.com/Caltech-IPAC/Montage

    Propagate MAGZP keyword from the raw image header.
    '''
    if verbose:
        print('Make header')

    if target is None:
        if isinstance(ra, str):
            target = '%s %s' %(ra, dec)
        else:
            target = '%.6f %.6f' %(ra, dec)
    
    montage.mHdr(target, width/60., 'montage.hdr', system=None, equinox=None, height=None, pix_size=pix_size, rotation=rotation)
    #montage.mHdr(target, width/60., 'montage.hdr', system=None, equinox=None, height=None, pix_size=pix_size, rotation=rotation, twomass_band=band)

    if os.path.exists('montage-tmp'):
        os.system('rm -rf montage-tmp')

    if verbose:
        print('Generate %s/%s cutout:\n %s, width=%f arcmin' %(survey, band, target, width))

    # status = montage.mExec(survey, band, raw_dir=None, n_tile_x=None, n_tile_y=None, level_only=False, keep=True, corners=False, output_image=output, debug_level=1, region_header='montage.hdr', header=None, workspace_dir='montage-tmp')
    status = montage.mExec(survey, band, raw_dir=None, n_tile_x=None, n_tile_y=None, level_only=False, keep=True, remove=False, output_image=output, debug_level=1, region_header='montage.hdr', header=None, workspace_dir='montage-tmp')

    ### Get ZP
    file_ix = -1
    scale_ix = -3
    if survey == 'SDSS':
        file_ix = -3
        scale_ix = -2

    line = open('montage-tmp/remote.tbl').readlines()[3].split()
    print(line)
    scale = float(line[scale_ix])

    raw_img = pyfits.open('montage-tmp/raw/%s' %(line[file_ix]))

    im = pyfits.open(output, mode='update')

    if survey == '2MASS':
        if band == 'j':
            MAGZP = 20.905 #average of the jband read_2 from 2MASS
        elif band == 'h':
            MAGZP = 20.3243041848514 #average of the 2 H band images from this region https://irsa.ipac.caltech.edu/workspace/TMP_M7TTAa_20769/IM/inv_qxqRD8/found.html 
        #ZP = raw_img[0].header['MAGZP']-2.5*np.log10(scale)
        ZP = MAGZP-2.5*np.log10(scale)

    if survey == 'SDSS':
        if 'FLUX20' not in list(raw_img[0].header.keys()):
            ZP = 27
        else:
            ZP = 20+2.5*np.log10(raw_img[0].header['FLUX20'])

    #im[0].header['MAGZP'] = ZP
    im.flush()


