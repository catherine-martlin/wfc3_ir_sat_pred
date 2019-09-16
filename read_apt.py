#!/usr/bin/env python 

'''
Read an APT file and extract a summary of information from
it in order to enabble an estimat of the likely persistence.

The results are written out as a astropy.ascii table


Command line usage:

    usage: read_apt.py filename

Description:  

    The routine reads the APT file and abstracts certain
    information from it.  The information is written
    to an ascii table where the default name is
    'root'_sum.txt where root is the root name of the 
    apt file.

Primary routines:

    doit

Notes:
                                       
History:

150723 ksl Coding begun

'''

import sys

from astropy.io import ascii
from astropy.table import Table,join
import numpy as np
import urllib.request, urllib.error, urllib.parse
import xml.etree.ElementTree as ET

def read_apt_main(filename='test.apt',outfile=''):
    '''
    Read a standard apt file and create a summary of certain
    elements of the file.  Produce an output file that contains
    one line per exposure.

    Notes:
    The routine reads the APT file using element tree.  In ET
    one can burrow down into an XML file as if it were
    a directory structure. Wildcards are allowed.  There are
    two main ways to burrow down.  (This is useful, because
    differnt HST proposals have different depths, that 
    may or may or may not be of intereste, e.g the example
    proposal I looked at had an ExpposureGroup, becasue it
    contained a Non-Interuptable Sequence, which may or may 
    not always occur.)
    
    If your are lookiong for a single "container" use "find".  
    If you want to find all containers of a given type use 
    findall instead. This produces something like a list of 
    visits, and you can then read the information in it.  
    Once you are at the right level, and are looking for a 
    value or a keyword, then you use "get".

    The root level for an HST proposal is HSTProposl, so one
    does not use this in burrowing down.

    History:

    150723  ksl Coded
    150805  ksl Added check on NSAMP to see whether it is defined
            If not, which occurs for UVIS observations, set the
            exposure time to -99
    150901  ksl Modified so the outputfile can be specified
    150901  ksl Add a comment field, and put information about the 
            type of vist there.

    '''
    # Determine a default file name
    if outfile=='':
        try:
            outfile=filename[0:filename.rindex('.')]+'_sum.txt'
        except ValueError:
            outfile=filename+'_sum.txt'

    # Read the XML file into an Element tree
    try:
        tree=ET.parse(filename)
    except IOError:
        print('Error: %s does not appear to exist' % filename)
        return  False

    apt=tree.getroot()    # The entire apt file

    print('Proposal:',apt.get('Phase2ID'))

    # Get some basic information about the proposal

    pi=apt.find('ProposalInformation/PrincipalInvestigator')
    print('PI:',pi.get('LastName'))

    # Get information about the targets
    # First create lists for the targets to store the infomration we 
    # want
    xname=[]
    xpos=[]

    targets=apt.findall('Targets/FixedTarget')

    for one_target in targets:
        position=one_target.find('EquatorialPosition')

        xname.append(one_target.get('Name'))
        
        ### Older APT had Coords in different format (e.g. 11712)
        pos = position.get('Value')
        if pos is None:
            RA = position.find('RA')
            DEC = position.find('DEC')
            if RA is None:
                xpos.append(pos)
            else:
                rr = RA.attrib
                dd = DEC.attrib
                pos = '%s %s %s %s %s %s' %(rr['Hrs'], rr['Mins'], rr['Secs'], dd['Degrees'], dd['Arcmin'], dd['Arcsec'])
                
        xpos.append(pos)

    # Create a table of the target information
    target_table=Table([xname,xpos],names=['Target','Pos'])
    print('This is the Target Table')
    print(target_table)
    print('Now look at any patterns')

    pattern_type=[]
    pattern_number=[]

    patterns=apt.findall('Patterns/Pattern')
    for pattern in patterns:
        ptype=pattern.find('PatternType')
        ptype=ptype.text
        number=pattern.find('Number')
        number=number.text
        print('Test', ptype,number)
        pattern_type.append(ptype)
        pattern_number.append(number)

    print('Now process the visits')

    # Now get information about the visits, First create lists for the infomation we want to keep
    xvisit=[]
    xtarg=[]
    xsamp_seq=[]
    xnsamp=[]
    xiter=[]
    xno=[]
    xfilt=[]
    xaper=[]
    xconfig=[]

    xscan = []
    
    vcomment=[]
    ecomment=[]

    visits=apt.findall('Visits/Visit')
    for one in visits:
        vcomment.append('Undithered')
    mosaics=apt.findall('Visits/HstQuadPattern')
    for one in mosaics:
        vcomment.append('Mosaic Visit')

    visits=visits+mosaics

    i=0
    for visit in visits:
        exposures =visit.findall('.//Exposure') # This locates all of the exposures
        print('The number of expsosures is %d\n' %  len(exposures))
        # Now get infomation about the indidual exposures
        # These may be buried in other containers, Right now we will ignore these
        # And just collect the exposure information
        this_visit=visit.get('Number')
        print('There are %d exposures (ignoring repeats) in visit %s' % (len(exposures),this_visit))

        for exposure in exposures:
            ipars=exposure.find('InstrumentParameters')
            xno.append(exposure.get('Number'))
            xvisit.append(visit.get('Number'))
            xtarg.append(exposure.get('TargetName'))
            xfilt.append(exposure.get('SpElement'))
            xaper.append(exposure.get('Aperture'))
            xconfig.append(exposure.get('Config'))
            xsamp_seq.append(ipars.get('SAMP-SEQ'))
            xnsamp.append(ipars.get('NSAMP'))
            xiter.append(exposure.get('NumberOfIterations'))
            ecomment.append(vcomment[i])  # Add the visit commetn to this
            
            ### Spatial scans
            scan = exposure.find('SpatialScan')
            if scan is None:
                xscan.append(None)
            else:
                xscan.append(float(scan.get('Rate')))
                
        # Now locate the exposures which have patterns and associate this with an exposure number
        egroups=visit.findall('.//Exposure/..[@Pattern]')
        for group in egroups:
            this_pattern=group.get('Pattern')
            # Now find the associated exposures
            xexposures=group.findall('.//Exposure')
            for one in xexposures:
                this_exposure=one.get('Number')
                print('whoopee',this_exposure,this_pattern)
                # Now we need to locate this visit and exposure number in our list
                jj=0
                while jj<len(ecomment):
                    if this_visit==xvisit[jj] and this_exposure==xno[jj]:
                        kk=0
                        coment='Dither'
                        while kk<len(pattern_number):
                            if pattern_number[kk]==this_pattern:
                                comment=pattern_type[kk]
                                break
                            kk=kk+1
                        if ecomment[jj]=='Undithered':
                            ecomment[jj]=comment
                        else:
                            ecomment.append('; '+comment)
                        print('Ta Da')
                        break
                    jj=jj+1
        i=i+1

    exposure_table=Table([xvisit,xno,xtarg,xconfig,xaper,xfilt,xsamp_seq,xnsamp,xiter, ecomment, xscan],names=['Visit','ExpNo','Target','Config','Aper','Filter','SAMP-SEQ','NSAMP','Repeats','Comment', 'ScanRate'])

    print('This is the Exposure Table')
    print(exposure_table)

    # There is apparently some kind of bug that sometimes occurs in np sorting so I need to add a row number or the table before joinint them
    n=np.arange(1,len(exposure_table)+1)
    exposure_table['line']=n

    # Check that both tables have rows, and return an error if this is not the case
    if len(exposure_table)==0 or len(target_table)==0:
        print('There are %d and %d rows in the exposere and target tables' % (len(exposure_table),len(target_table)))
        print('Error: read_apt: There should be rows in both exposure and target tables')
        return False

    # Now join the two tables
    test=join(exposure_table,target_table,join_type='left',keys='Target')

    # Now get this back into visit and exposure order
    q=np.argsort(test['line'])
    test=test[q]

    # Now get rid of the extra column
    test.remove_column('line')
    print('This is the combined table')

    # Now convert the positions to degees
    ra=[]
    dec=[]

    for one in test:
        xpos=one['Pos']
        xra,xdec=pos2deg(xpos)
        ra.append(xra)
        dec.append(xdec)


    test['RA']=ra
    test['Dec']=dec

    test['RA'].format='12.7f'
    test['Dec'].format='12.7f'

    # Now add the exposure times
    time=[]
    for one in test:
        print(one['SAMP-SEQ'],one['NSAMP'],type(one['NSAMP']), one['Aper'])
        xtime=get_read_time(one['SAMP-SEQ'], one['NSAMP'], APER=one['Aper'], scan_rate=one['ScanRate'])
        # 150809: ksl - deleted this checks because they were not working,
        # however a check may not be needed if we are sure to include ony
        # IR observations.  If they are need in future, it might be better
        # to check based on the value of NSAMP rather than the type
        # if type(one['NSAMP']) is str or type(one['NSAMP']) is np.str:
        #   xtime=get_read_time(one['SAMP-SEQ'],one['NSAMP'])
        # else:
        #   xtime=-99.
        print(xtime)
        time.append(xtime)
    
    test['Exptime']=time

    # Now reorder the rows of the table:
    xtest=test['Visit','ExpNo','Target','RA','Dec','Config','Aper','Filter','SAMP-SEQ','NSAMP','Exptime','Repeats','Comment', 'ScanRate']
    xtest.write(outfile,format='ascii.fixed_width_two_line',overwrite=True)
    print(xtest)

    return True

def pos2deg(pos='13 26 46.2800 -47 28 44.60'):
    '''
    Convert a position to degrees

    History

    150805  ksl Instead of eval, used float
            to convert strings, because
            a number like 08 is not converted
            by eval.
    '''
    try:
        z=pos.split()
    except AttributeError:
        return -99,-99

    # Convert strings to numbers
    zz=[]
    for one in z:
        # zz.append(eval(one))
        zz.append(float(one))
    z=zz

    ra=z[0]+z[1]/60.+z[2]/3600.
    ra=15.*ra

    dec=z[3]
    q=z[4]/60.+z[5]/3600.
    if dec<0:
        dec=dec-q
    else:
        dec=dec+q

    return ra,dec
    

def get_read_time(SAMPSEQ='STEP50', NSAMP=8, APER='IR', scan_rate=None):
    '''
    Get the exposure time from the SAMSEQ and NSAMP

    History

    150730  ksl Adapted from Gabe's routine get_times
    150804  ksl Debugged. The problem was that that the 
            NSAMP starts with 1, whereas the dictionarys
            are 0 based.
    150805  ksl Checked that NSAMP is actually defined and return
            -99 if is not
    '''

    try:
        NSAMP=int(NSAMP)-1
    except TypeError:
        return -99.

    
    samples = {'RAPID': np.arange(15)*2.932+2.932,
            'STEP25': np.array([2.9, 5.89, 8.80, 11.73, 24.23, 49.23, 74.23, 99.23, 124.23, 149.23, 174.23, 199.23, 224.23, 249.23, 274.23]),
               'STEP50': np.array([2.9, 5.87, 8.80, 11.73, 24.23, 49.23, 99.23, 149.23, 199.23, 249.23, 299.23, 349.23, 399.23, 449.23, 499.23]),
               'STEP100': np.array([2.9, 5.8, 8.8, 11.7, 24.23, 49, 99, 199, 299, 399, 499, 599, 699, 799, 899]),
               'STEP200': np.array([2.9, 5.8, 8.8, 11.7, 24.23, 49, 99, 199, 399, 599, 799, 999, 1199, 1399, 1599]),
               'STEP400': np.array([2.9, 5.8, 8.8, 11.7, 24.23, 49, 99, 199, 399, 799, 1199, 1599, 1999, 2399, 2799]),
                'SPARS5': 2.93+np.arange(15)*5,
               'SPARS10': 2.93+np.arange(15)*10,
               'SPARS25': 2.93+np.arange(15)*25,
               'SPARS50': 2.93+np.arange(15)*50,
               'SPARS100': 2.93+np.arange(15)*100,
               'SPARS200': 2.93+np.arange(15)*200}
    
    if '512' in APER:
        samples = {'RAPID': 0.853*(np.arange(15)+1),
                   'STEP25': np.array([0.853, 1.71, 2.56, 3.412, 13.8, 37.8, 59.7, 82.6, 105.5, 128.4, 151.4, 174.3, 197.2, 220.1, 243.0]),
                    'SPARS5': 0.853+np.arange(15)*2.92,
                   'SPARS25': 0.853+np.arange(15)*22.92}        
    
    if '256' in APER:
        samples = {'RAPID': 0.278*(np.arange(15)+1),
                    'SPARS5': 0.28+np.arange(15)*2.35,
                   'SPARS10': 0.28+np.arange(15)*7.34,
                   'SPARS25': 0.28+np.arange(15)*22.346}        

    if '128' in APER:
        samples = {'RAPID': 0.113*(np.arange(15)+1),
                   'SPARS10': 0.113+np.arange(15)*7.18}
    
    if '64' in APER:
        samples = {'RAPID': 0.061*(np.arange(15)+1)}
            
    try:
        xtime=samples[SAMPSEQ][NSAMP]
    except KeyError:
        print('Error: unknown SAMPSEQ %s with %d reads' % (SAMPSEQ,NSAMP))
        return -99.
    
    if scan_rate is not None:
        print('Spatial scan! %.1f' %(scan_rate))
        xtime = 2*0.13/float(scan_rate)
        
    return xtime


def fetch_apt(file='14183.apt'):
    '''
    Retrieve an apt file from the from the public location

    Notes;  

    Here is a typical link

    http://www.stsci.edu/hst/phase2-public/14258.apt

    History:

    150805  ksl Coded

    '''

    filename='http://www.stsci.edu/hst/phase2-public/'+file
    try:
        response = urllib.request.urlopen(filename)
        html = response.read()
    except urllib.error.HTTPError:
        print('Error: Could not retrieve %s' % file)
        return False
    # print html

    x=open(file,'wb')
    x.write(html)
    x.close()
    return True


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        # doit(int(sys.argv[1]))
        read_apt_main(sys.argv[1])
    else:
        print('usage: read_apt.py filename')
