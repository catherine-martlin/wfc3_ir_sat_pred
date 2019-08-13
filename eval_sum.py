#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

A routine which simply concatenates all of the ascii tables produced by 
various runs of eval_one so that one can identify those programs likely 
to be bad actors. 

The program also creates data  to allow a comparison with actual
bad actors in executed programs


Command line usage (if any):

    usage: eval_sum.py [filename]

Description:  

    The program simply uses glob to find all the tables from individual 
    runs of eval_one, and then uses vstack to put everything into one file

    If a filename is not provided then the resutls will be in eval_summary.txt


Primary routines:

    doit

Notes:

    The routine needs to be run in the direcorty where eval_one has
    been run previously.

    The program currently used to make the comparison between observations
    and reality is called compare2reality
                                       

History:

151110 ksl  Coding begun
170809 ksl  Adapted to run in Python3

'''

import sys
from astropy.io import ascii
from astropy.table import Table,vstack
import numpy
import glob
import astropy


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
        data=ascii.read(filename)
    except IOError:
        print('Error: file %s does not appear to exist' % filename)
        return False


    return data


def doit(outputfile='eval_summary.txt'):
    '''
    This routine simply looks for all of the directories with a valid
    summary file created by eval_one and concatenates them together
    into a single tabble, adding the program ID

    Description:

    Notes:

    History:

    151110    ksl    Coded


    '''

    directories=glob.glob('*_dir')

    xmaster=False

    for one_dir in directories:
        # Get the root program ID from the directory name
        print(one_dir)
        prog_id=one_dir.replace('_dir','')
        xtable=read_table('%s/%s.txt' % (one_dir,prog_id))

        # eliminate UVIS
        i=0
        keep=[]
        while i <len(xtable):
            if xtable['NSAMP'][i]!='None':
                keep.append(i)
            i=i+1
        xtable=xtable[keep]


        xtable['ProgID']=prog_id
        print('test',xtable.colnames)
        xtable.remove_column('NSAMP')

        if xtable!=False:

            # Put ProgID at the beginning

            z=['ProgID']
            q=xtable.colnames
            for one in q:
                if one!='ProgID':
                    z.append(one)
            xtable=xtable[z]


            if xmaster==False:
                master=xtable.copy()
                xmaster=True
            else:
                try:
                    master=vstack([master,xtable])
                    print('merged',one_dir)
                except astropy.table.np_utils.TableMergeError:
                    print('skip',one_dir)
    
    print(master)
    master.write(outputfile,format='ascii.fixed_width_two_line',overwrite=True)
    



    return


    





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        doit(sys.argv[1])
    else:
        print('No output file name supplied, using default : eval_summary.txt')
        doit()
