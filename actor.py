#!/usr/bin/env python 

'''
Retrieve 2 mass objects near a given RA and DEC and estimate
the number of objects that are likely to cause persistence 
problems.


Command line usage (if any):

    usage: actor.py rootname  ra dec filter exposure 

    rootname - rootname of the ouput files
    ra - right ascension of the center of the field in degrees
    dec - dec of the center of the field in degress
    filter - filter name, e.g F11O@
    exposure - exposure time

Description:  
    The routine retrieveds 2Mass data near the center of the field
    and uses the J-band magnitude to estimate the persistence

    The saturation level is estimated from a lookup table which 
    contains the filter name and a rate.

    The retrieved 2Mass data is in the file rootname.2Mass.txt and is
    in the so-called IPAC table format, which astropy can read

    Information about the targets within the FOV of the IR channel
    are contained in the file rootname.per.txt. This is a more readable
    astropy table format

    To see visually how likely a filed is to cause persistnece a 
    plot file, rootname.png, is created.  The size of the points indicate
    how much saturation a particular object will cause.

Notes:
    At present this routine is just a prototype, and one of the items
    that really needs to be checked is whether the numbers from the 
    ETC are really correct, or want.  We also only make use of the j-band
    2 Mass data and in principle one might be able to do beter
                                       
History:
150716 ksl Coding begun

'''
import math
import os
import subprocess
import sys

from astropy.io import ascii
from astropy.table import Table, Column
import numpy as np
import pylab

import read_apt
import persist_2mass
import xhtml as html
import ql_html


RADIAN=57.29578
INFINITY=1.e32
SATLEVEL=70000

def radec2deg(ra='05:13:06.2',dec='-10:13:14.2'):
    ''' 
    Convert an ra/dec string to degrees.  The string can already
    be in degrees in which case all that happens is a conversion to
    a float

    If what is transferred is a float, the routine assumes it has been
    given ra and dec in degrees and just returns ra,dec
    
    '''

    # print 'Before',ra,dec
    try:
        r=ra.split(':')
        d=dec.split(':')
    except AttributeError:
        return ra,dec

    # print 'After',ra,dec
    rr=float(r[0])
    if len(r)>1:
        rr=rr+float(r[1])/60.
    if len(r)>2:
        rr=rr+float(r[2])/3600.
    if len(r)>1:
        rr=15.*rr  # Since we assume ra was in hms
    
    dd=float(d[0])
    x=0
    if len(d)>1:
        x=x+float(d[1])/60.
    if len(d)>2:
        x=x+float(d[2])/3600.
    
    if dd > 0:
        dd=dd+x
    else:
        dd=dd-x
    
    return rr,dd  


def radec2hms(ra='225.2',dec='-17.35'):
    '''
    Convert an ra dec in degress to hms.  The input formats may be
    strings or floats.

    '''

    try: 
        ra=eval(ra)
    except TypeError:
        ra=float(ra)

    try:
        dec=eval(dec)
    except TypeError:
        dec=float(dec)

    xra=ra/15.
    ra=int(xra)
    xra=xra-ra
    xra=xra*60
    min=int(xra)
    xra=xra-min
    xra=xra*60.
    ra_string='%02d:%02d:%06.3f' % (ra,min,xra)

    xdec=math.fabs(dec)
    deg=int(xdec)
    xdec=xdec-deg
    xdec=xdec*60
    min=int(xdec)
    xdec=xdec-min
    xdec=xdec*60.
    if dec<0:
        deg=deg*-1
    dec_string='%3d:%02d:%05.2f' % (deg,min,xdec)
    return ra_string,dec_string

def distance(r1,d1,r2,d2):
    '''
    distance(r1,d1,r2,d2)
    Return the angular offset between two ra,dec positions
    All variables are expected to be in degrees.
    Output is in degrees

    Note - This routine could easily be made more general
    This version supports np arrays
    '''
#   print 'distance',r1,d1,r2,d2
    r1=r1/RADIAN
    d1=d1/RADIAN
    r2=r2/RADIAN
    d2=d2/RADIAN

    xlambda=np.sin(d1)*np.sin(d2)+np.cos(d1)*np.cos(d2)*np.cos(r1-r2)
    xlambda=np.arccos(xlambda)
    xlambda=np.array(xlambda*RADIAN)
    return xlambda


def read_file(filename,char=''):
    '''
    Read a file and split it into words, eliminating comments
    
    char is an optional parameter used as the delimiter for
    splitting lines into words.  Otherwise white space is
    assumed.

    History:
    
    110729  ksl Added optional delimiters
    141209  ksl Reinstalled in my standard startup
            script so there was flexibility to
            read any ascii file
    '''

    try:
        f=open(filename,'r')
        xlines=f.readlines()
        f.close()
    except IOError :
        print("The file %s does not exist" % filename)
        return []   
    
    lines=[]
    
    i=0
    while i<len(xlines):
        z=xlines[i].strip()
        if char=='':
            z=z.split()
        else:
            z=z.split(char)
        if len(z)>0:
            if z[0][0]!='#':
                lines=lines+[z]
        i=i+1
    return lines


def read_table(filename='foo.txt'):
    '''
    Read a file using astropy.io.ascii and return it 
    '''
    try:
        data=ascii.read(filename)
    except IOError:
        print('Error: file %s does not appear to exist' % filename)
        return
    print('Here are the column names:')
    print(data.colnames)

    return data

def gen_circle(rad=1.):
    '''
    Generate a set of data points that describe a circle
    '''
    theta=np.linspace(0,360.1000)
    theta=theta/RADIAN
    x=rad*np.cos(theta)
    y=rad*np.sin(theta)
    return x,y

def get_counts(mag=13,exp=500,xfilt='F105W'):
    '''
    Get the number of counts in the central pixel for a source
    with a magnitude of mag in an exposure time of exp through
    xfilt

    Notes:

        Rates are calculated from the ETC and correspond to the
        counts in 1 sec from the peak pixel for a 15th mag
        source. See for example

        WFC3IR.im.727729 

        The second column in xdata is calculated for Vega mags
        in J, the second is for K mags

    (Why k mag here ???)
    
    History

    150808  ksl Updated so that rates were taken from the H
            band catalog if the WFC3 filter was close to
            that range.  
    150909  ksl  Added gratings, but it is not clear that these 
            are accurate

    '''
    xdata=[
            ['F098M',1961,1357],
            ['F105W',3781,1395],
            ['F110W',7134,2634],
            ['F125W',4663,1722],
            ['F140W',6856,2431],
            ['F126N', 256,78.6],
            ['F127M',1114,411],
            ['F128N', 242,89],
            ['F130N', 250,92],
            ['F132N', 245,90.5],
            ['F139M',1069,395],
            ['F153M',1198,442],
            ['F160W',4696,1734],
            ['F164N', 307,114],
            ['F167N', 327,121],
            ['G102' ,34, 5.5],
        ['G141', 86 , 13.3]
        ]
    
    two_mass_filter=persist_2mass.get_2mass_filter(xfilt)

    ifilter=1
    if two_mass_filter=='h':
        ifilter=2

    rate=10000.
    ok=False
    for one in xdata:
        if one[0]==xfilt:
            rate=one[ifilter]
            ok=True
            break
    if ok==False:
        print('Error: filter %s not found' % xfilt)
    
    x=exp*rate*10**(0.4*(15-mag))

    return x


def actor_main(ra=66.76957,dec=26.10453,xfilter='F110W',exp=100.,rad=3,outroot='out'):
    '''
    Retrieve the 2mass catalog for this position in the sky, and estimate
    the number of stars which are saturated

    ra,dec - of center of the field in degrees 
    xfilter - filter used
    exp - exposure time in seconds
    rad - radius in arcmin of the region from which stars are retrieved

    Notes:
    Normally, when called by other routines, outroot will be the name of the 
    working directory and the root name of the apt file

    History:
    15088   ksl Commented better. Still need to select which band is used
    150902  ksl Modified so that files are not regenerated
    151105  ksl Modified so returns the number of stars with 1x, 5x, and 10x
            saturationw wigthin 91 arscec of center position

    '''
    radius=rad/60.

    two_mass_table=outroot+'.2Mass.txt'

    if os.path.isfile(two_mass_table)==False:
        command_string='''curl -o %s "https://irsa.ipac.caltech.edu:443/TAP/sync?QUERY=SELECT+*+FROM+fp_psc+WHERE+CONTAINS(POINT('J2000',ra,dec),CIRCLE('J2000',%f,%f,%f))=1&format=IPAC_TABLE"'''  % (two_mass_table,ra,dec,radius)
        proc=subprocess.Popen(command_string,shell=True,stdout=subprocess.PIPE)
        x=proc.communicate()[0]
    else:
        print('2Mass objects have already been retrieved to file %s' % two_mass_table)


    data=ascii.read(two_mass_table,format='ipac')

    print('Retrieved %d 2Mass objects' % len(data))

    z=60.*distance(ra,dec,data['ra'],data['dec'])

    data['distance']=z

    # Determine which two mass magnitude to use and then estimate the rates
    twomass_filter=persist_2mass.get_2mass_filter(xfilter)
    q=data['j_m']
    if twomass_filter=='h':
        q=data['h_m']
    print('Using ',twomass_filter)

    if exp<0:
        print('Error: Exposure time %f is less than 0, which is impossible' % exp)

    counts=get_counts(q,exp,xfilter)

    data['counts']=counts
    data['saturation']=counts/SATLEVEL

    xdata=data['distance','j_m','h_m','counts','saturation']
    xdata=xdata.group_by('distance')

    # We want to create an output table that only has some rows of the input table
    # Select the rows that we want
    xx=xdata['distance','j_m','h_m','counts','saturation']

    i=0
    nrows=0
    rows_to_keep=[]
    while i<len(xdata):
        one=xdata[i]
        if one['saturation']>0.5 and one['distance']<(91.78/60):
            rows_to_keep.append(i)
        i=i+1
    xxx=xx[rows_to_keep]

    xxx['distance'].format='6.2f'
    xxx['j_m'].format='6.2f'
    xxx['h_m'].format='6.2f'
    xxx['counts'].format='8.2e'
    xxx['saturation'].format='6.1f'

    # print 'Here are the results'
    if len(xxx)>0:
        xxx.write(outroot+'.stars.txt',format='ascii.fixed_width_two_line',overwrite=True)
    else:
        print('This field has no bright stars of concern')

    # Now make a plot
    data['x']=60*(data['ra']-ra)*math.cos(dec/RADIAN)
    data['y']=60*(data['dec']-dec)

    pylab.figure(1,(14,6))
    pylab.clf()
    pylab.subplot(121)
    pylab.scatter(data['x'],data['y'],s=0.1*data['saturation'],alpha=0.3)

    x,y=gen_circle(64.6/60.)
    pylab.plot(x,y,'r-')
    x,y=gen_circle(91.78/60.)
    pylab.plot(x,y,'g-')


    pylab.axis((rad,-rad,-rad,rad))
    pylab.xlabel('Delta RA (arcmin)')
    pylab.ylabel('Delta Dec (arcmin)')
    pylab.draw()

    pylab.subplot(122)

    dist=91.78/60

    x,y=return_one_histogram(xdata['distance'],xdata['saturation'],1)
    pylab.plot(x,y,ls='steps-post',label='1x')
    i=0
    while i<len(x) and x[i]<dist:
        i=i+1
    num_1x=i

    x,y=return_one_histogram(xdata['distance'],xdata['saturation'],5)
    pylab.plot(x,y,ls='steps-post',label='5x')
    i=0
    while i<len(x) and x[i]<dist:
        i=i+1
    num_5x=i

    x,y=return_one_histogram(xdata['distance'],xdata['saturation'],10)
    pylab.plot(x,y,ls='steps-post',label='10x')
    i=0
    while i<len(x) and x[i]<dist:
        i=i+1
    num_10x=i

    pylab.legend(loc='best')
    pylab.xlabel('Distance (arcmin)')
    pylab.ylabel('Number of stars < Distance')
    a=pylab.axis()
    pylab.axis((0,2,0,a[3]))

    # Draw some vertical lines on the plot
    x1=64.6/60
    x2=91.78/60
    pylab.plot([x1,x1],[0,10e5],'r-')
    pylab.plot([x2,x2],[0,10e5],'g-')

    pylab.draw()
    pylab.savefig(outroot+'.png')

    return  num_1x,num_5x,num_10x

def return_one_histogram(xdis,xsat,saturation=5):
    '''
    From a col of a table containing distances from
    the central position and a column of a table
    containing saturation levels return two arrays
    the distance and number of objects within that 
    distance
    '''
    xdis=np.array(xdis)
    xsat=np.array(xsat)

    z=np.select([xsat>saturation],[True],default=False)

    elements=[]
    i=0
    while i<len(z):
        if z[i]==True:
            elements.append(i)
        i=i+1

    xx=xdis[elements]
    nobj=np.linspace(1,len(xx),len(xx))
    return xx,nobj

def select_rows(x,condition='x>1'):
    '''
    Return the indices that satisfy a condition
    '''
    x=np.array(x)
    z=np.select([eval(condition)],[True],default=False)

    elements=[]
    i=0
    while i<len(z):
        if z[i]==True:
            elements.append(i)
        i=i+1
    return elements

def do_proposal(table='test_sum.txt'):
    '''
    Process a proposal from the summary file generated by read_apt

    History

    150804  ksl Added to the routine
    '''
    x=read_table(table)

    root=table[0:table.rindex('.')]

    logfile=root+'.log'

    g=open(logfile,'w')
    g.write('Processing %s\n' % table)

    for one in x:
        print(one)
        string='Header Visit%s ExpNo%03d %s' % (one['Visit'],one['ExpNo'],one['Target'])
        g.write('%s\n' % string)

        print('Test Config', one['Config'].count('IR'))
        if one['Config'].count('IR')==0:
            g.write('This is a UVIS exposure, so persistence is not an issue.\n')
        elif one['RA']!=-99.:
            name='Visit%s_ExpNo%03d_%s' % (one['Visit'],one['ExpNo'],one['Target'])
            doit(ra=one['RA'],dec=one['Dec'],xfilter=one['Filter'],exp=one['Exptime'],rad=3,outroot=name)
            g.write('This is an exposure at %f %f with the %s filter of %.2f seconds.\n' % (one['RA'],one['Dec'],one['Filter'],one['Exptime']))
            string='Image %s.png' % name
            g.write('%s\n' % string)
        else:
            g.write('This exposure did not involve a target with a fixed position observed with WFC3/IR\n')
    g.close()

    make_html(logfile)
    make_ql_html(logfile)

    return

def make_html(logfile='test_sum.log',outfile=''):
    '''
    Create an html file using information in the log file

    History

    150901 ksl  fixed so that the htmlfile could be specified
    151112 ksl  Modified so that number of saturating stars could be printed out
            and so that only some of the columns in the table created by actor.actor_main
        were included in the html table that is part of the output.
    '''
    xsize=800
    ysize=400

    try:
        root=logfile[0:logfile.rindex('_sum')]
    except ValueError:
        root=logfile[0:logfile.rindex('.')]

    hstring=html.begin('Summary for %s' % root)

    x=open(logfile,'r')
    x=x.readlines()
    expcount = 0
    for line in x:
        xline=line.strip()
        words=xline.split()
        if len(words)>0:
            if words[0]=='Header':
                expcount += 1
                xline="<a name=\"EXP%04d\">%s</a>" %(expcount, xline.replace('Header','')) 
                hstring=hstring+html.h2(xline)
            elif words[0]=='Table':
                xxx=read_table(words[1])
                # Drop some columns, and reorder sligthly
                xxx=xxx['Visit','ExpNo','Target','RA','Dec','Config','Aper','Filter','SAMP-SEQ','NSAMP','Exptime','Repeats','ScanRate','Comment']
                ### make links from exposure number
                links = []

                for i in range(len(xxx)):
                    link = '<a href=#EXP%04d>%d</a>' %(i+1, xxx['ExpNo'][i])
                    links.append(link)
                
                xxx.add_column(Column(data=links, name='Exp'), 1)
                xxx.remove_column('ExpNo')
                
                # Note: ksl.  I had troble getting the column names to form a 
                # sincle record.  This is what worked in the end
                my_table=[xxx.colnames[:]]
                                
                for one in xxx:
                    my_table.append(one)

                hstring=hstring+html.table(my_table)
            elif words[0]=='Image':
                ximage=words[1]
                hstring=hstring+html.image(ximage,'Bright stars image',xsize,ysize)
            else: 
                if xline.startswith('Saturated'):
                    print(xline)
                    data = xline.split()[-5::2]
                    colors = ['green','blue','red']
                    for i in range(3):
                        if int(data[i]) == 0:
                            data[i] = "<span style=\"color:grey\">%s</span>" %(data[i])
                        else:
                            data[i] = "<span style=\"color:%s; font-weight:bold;\">%s</span>" %(colors[i] ,data[i])
                            
                    hstring += html.paragraph('Saturated pixels in image:') + html.table([xline.split()[-6::2], data], width="200px")
                elif xline.startswith('Number saturated'):
                    data = xline.split()[-5::2]
                    colors = ['green','blue','red']
                    for i in range(3):
                        if int(data[i]) == 0:
                            data[i] = "<span style=\"color:grey\">%s</span>" %(data[i])
                        else:
                            data[i] = "<span style=\"color:%s; font-weight:bold;\">%s</span>" %(colors[i] ,data[i])
                            
                    hstring += html.paragraph('Saturated stars in field:') + html.table([xline.split()[-6::2], data], width="200px")
                else:
                    hstring=hstring+html.paragraph(xline)            
    hstring=hstring+html.end()

    # Create the html file
    if outfile=='':
        html_file=root+'.html'
    elif outfile.count('html')==0:
        html_file=outfile+'.html'
    else:
        html_file=outfile
        
    print('test ',html_file)
    g=open(html_file,'w')
    g.write(hstring)
    g.close()

def make_ql_html(logfile='test_sum.log',outfile=''):
    '''
    Create an html file using information in the log file that's readable by QuickLook website code
    History
    130717 cam  edited make_html to make this to run on logfile but make different img 
    150901 ksl  fixed so that the htmlfile could be specified
    '''


    xsize=800
    ysize=400

    try:
        root=logfile[0:logfile.rindex('_sum')]
    except ValueError:
        root=logfile[0:logfile.rindex('.')]

    hstring=ql_html.begin('Summary for %s' % root)

    x=open(logfile,'r')
    x=x.readlines()
    expcount = 0
    for line in x:
        xline=line.strip()
        words=xline.split()
        if len(words)>0:
            if words[0]=='Header':
                expcount += 1
                xline="<a name=\"EXP%04d\">%s</a>" %(expcount, xline.replace('Header','')) 
                hstring=hstring+ql_html.h2(xline)
            elif words[0]=='Table':
                xxx=read_table(words[1])
                     
                ### make links from exposure number
                #print xxx, len(xxx)
                links = []
                for i in range(len(xxx)):
                    link = '<a href=#EXP%04d>%d</a>' %(i+1, xxx['ExpNo'][i])
                    links.append(link)
                
                xxx.add_column(Column(data=links, name='Exp'), 1)
                xxx.remove_column('ExpNo')
                
                # Note: ksl.  I had troble getting the column names to form a 
                # sincle record.  This is what worked in the end
                my_table=[xxx.colnames[:]]
                                
                for one in xxx:
                    my_table.append(one)

                hstring=hstring+ql_html.table(my_table)
            elif words[0]=='Image':
                ximage=words[1]
                hstring=hstring+ql_html.image(ximage,'Bright stars image',xsize,ysize)
            else: 
                if xline.startswith('Saturated'):
                    data = xline.split()[-5::2]
                    colors = ['green','blue','red']
                    for i in range(3):
                        if int(data[i]) == 0:
                            data[i] = "<span style=\"color:grey\">%s</span>" %(data[i])
                        else:
                            data[i] = "<span style=\"color:%s; font-weight:bold;\">%s</span>" %(colors[i] ,data[i])
                            
                    hstring += ql_html.paragraph('Saturated pixels in image:') + ql_html.table([xline.split()[-6::2], data], width="200px")
                else:
                    hstring=hstring+ql_html.paragraph(xline)
                #print xline
                
    hstring=hstring+ql_html.end()


    # Create the html file

    if outfile=='':
        ql_html_file=root+'_ql.html'
    elif outfile.count('_ql.html')==0:
        ql_html_file=outfile+'_ql.html'
    else:
        ql_html_file=outfile
        
    line = 'test ' + ql_html_file
    g=open(ql_html_file,'w')
    g.write(hstring)
    g.close()

# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        print('OK: lets try it')
        rad=3
        outroot=sys.argv[1]
        ra=eval(sys.argv[2])
        dec=eval(sys.argv[3])
        xfilter=sys.argv[4]
        exp=eval(sys.argv[5])
        actor_main(ra,dec,xfilter,exp,rad,outroot)
    else:
        print('usage: actor.py outroot ra(deg) dec(deg)  filter  expossre(s)')
        print(__doc__)
