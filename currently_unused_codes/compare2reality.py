#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

This routine compares results of the persistence prediction tool
to actual observations (as calculated with bad_actor.py)


Command line usage (if any):

    usage: compare2reality.py 

    When run from the command line the program assumes standard
    names for inputs and outputs


Description:  

Primary routines:

    doit

Notes:
                                       
History:

151110 ksl Coding begun

'''

import sys
from astropy.io import ascii
from astropy.table import Table,join
import numpy
import pylab


def read_file(filename,char=''):
    '''
    Read a file and split it into words, eliminating comments
    
    char is an optional parameter used as the delimiter for
    splitting lines into words.  Otherwise white space is
    assumed.

    History:
    
    110729    ksl    Added optional delimiters
    141209    ksl    Reinstalled in my standard startup
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
    Read a file using astropy.io.ascii and 
    return this 

    Description:

    Notes:

    History:


    '''
    try:
        data=ascii.read(filename,format='fixed_width_two_line',fill_values=['--',-99])
    except IOError:
        print('Error: file %s does not appear to exist' % filename)
        return

    print('Here are the column names:')
    
    print(data.colnames)

    return data

def plot_results(summary_file='compare.txt'):
    '''
    Compare the results of the prediction program to acual data
    '''

    data=read_table(summary_file)

    rootname=summary_file[0:summary_file.rindex('.')]

    keep=[]

    print('test',data.mask)
    
    i=0
    while i<len(data):
        if data['RA'][i]>0 and data['image_1x'].mask[i]==False and data['x_1'].mask[i]==False:
            keep.append(i)
        i=i+1
    
    xdata=data[keep]

    print(xdata)

    pylab.figure(1,(6,6))
    pylab.clf()
    pylab.plot(xdata['x_1'],100*xdata['image_1x']/1e6,'o',label='1x')
    pylab.plot(xdata['x_5'],100*xdata['image_5x']/1e6,'o',label='5x')
    pylab.plot(xdata['x_10'],100*xdata['image_10x']/1e6,'o',label='10x')
    pylab.xlabel('Saturated Pixels Observed (%)')
    pylab.ylabel('Saturated Pixels Predicted(%)')
    pylab.axis((0,100,0,100))
    pylab.legend(loc='best')
    pylab.draw()
    pylab.savefig(rootname+'_pix.png')



    pylab.figure(2,(6,6))
    pylab.clf()
    pylab.plot(xdata['x_1'],xdata['star_1x'],'o',label='1x')
    pylab.plot(xdata['x_5'],xdata['star_5x'],'o',label='5x')
    pylab.plot(xdata['x_10'],xdata['star_10x'],'o',label='10x')
    pylab.xlabel('Saturaded Pixels Observed (%)')
    pylab.ylabel('Saturated Stars Predicted')
    pylab.legend(loc='best')
    pylab.draw()
    pylab.savefig(rootname+'_star.png')




def doit(prediction_summary='eval_summary.txt',observation_summary='bad_actor_table.txt',outputfile='compare.txt'):
    '''
    Compre predicted levels of saturation  to actual levels of saturation, where
        prediction_summary is file produced by eval_summary.phy and
        observation_file  is a file produced by bad_actor.py
    
    The program produces a table with the join between the two tables and some plots

    Description:

        The program simply associates lines in the two tables based on ProgID and LineID.  

    Notes:
        There may be more than one line associated with each prediction 

    History:

    151112    Modified to accept astropy table inputs


    '''
    
    


    prediction=read_table(prediction_summary)
    reality=read_table(observation_summary)


    # There are some formatting changes required in the 
    # columns that are dealt with below, because of the
    # ways astropy interprets various fields

    # x=numpy.array(reality['ProgID'],dtype='str')
    # reality.remove_column('ProgID')
    # reality['ProgID']=x

    x=numpy.array(reality['LineID'],dtype='str')
    i=0
    while i<len(x):
        x[i]='%6s'% x[i]
        x[i]=x[i].replace(' ','0')
        i=i+1
    reality.remove_column('LineID')
    reality['LineID']=x


    line_id=[]
    for one in prediction:
        xline_id='%2s.%3s' % (one['Visit'],one['ExpNo'])
        xline_id=xline_id.replace(' ','0')
        line_id.append(xline_id)
    
    prediction['LineID']=line_id

    

    x=join(prediction,reality,keys=['ProgID','LineID'])
    # x=join(prediction,reality,keys=['ProgID'])

    print(x)




    


    # This format is the easy to read back automatically
    ascii.write(x,outputfile,format='fixed_width_two_line',overwrite=True)

    plot_results(outputfile)

    return


    





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    doit()
