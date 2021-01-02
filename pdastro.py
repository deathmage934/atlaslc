#!/usr/bin/env python
'''
wrapper around pandas with convenience functions to ease handling of tables
A. Rest
'''
import sys,os,re,types,copy
import numpy as np
from astropy.time import Time
import astropy.io.fits as fits
import astropy
#import scipy
import pandas as pd
from astropy.nddata import bitmask

from scipy.interpolate import interp1d

def makepath(path,raiseError=True):
    if path == '':
        return(0)
    if not os.path.isdir(path):
        os.makedirs(path)
        if not os.path.isdir(path):
            if raiseError:
                raise RuntimeError('ERROR: Cannot create directory %s' % path)
            else:
                return(1)
    return(0)

def makepath4file(filename,raiseError=True):
    path = os.path.dirname(filename)
    if not os.path.isdir(path):
        return(makepath(path,raiseError=raiseError))
    else:
        return(0)
    
#https://numpy.org/doc/stable/reference/routines.set.html
def AorB(A,B):
    if len(A) == 0:
        return(B)
    if len(B) == 0:
        return(A)
    return(np.union1d(A,B))

def AandB(A,B,assume_unique=False):
    return(np.intersect1d(A,B,assume_unique=assume_unique))

def AnotB(A,B):
    return(np.setdiff1d(A,B))

def not_AandB(A,B):
    return(np.setxor1d(A,B))

class pdastroclass:
    def __init__(self,**kwargs):
        self.t = pd.DataFrame(**kwargs)
   
        self.verbose = 0
        
        # if self.auto_convert_dtypes==True, then before .to_string() is called in self.write, self.t.convert_dtypes() is run 
        # if self.default_dtypeMapping != None, then the dtype mapping is applied to the table .to_string() is called in self.write
        # This makes sure that formatters don't throw errors if the type of a column got changed to float or object during
        # one of the table operations
        self.auto_convert_dtypes = True
        self.default_dtypeMapping = None
        # example:
        # self.default_dtypeMapping = {'counter':np.int64}
        
        # default_formatters are passed to to_string() in self.write
        self.default_formatters = None
        # example:
        # self.default_formatters = {'MJD':'{:.6f}'.format,'counter':'{:05d}'.format}
       
        # dictionary for the splines. arguments are the y columns of the spline
        self.spline={}


    def load_spacesep(self,filename,test4commentedheader=True,namesMapping=None,roundingMapping=None,
                      na_values=['None','-','--'],**kwargs):
        
        kwargs['delim_whitespace']=True

        #also test for commented header to make it compatible to old format.
        self.load(filename,na_values=na_values,test4commentedheader=test4commentedheader,
                  namesMapping=namesMapping,roundingMapping=roundingMapping,**kwargs)

        return(0)

    def load(self,filename,raiseError=True,test4commentedheader=False,namesMapping=None,roundingMapping=None,**kwargs):
        #self.t = ascii.read(filename,format='commented_header',delimiter='\s',fill_values=[('-',0),('--',0)])
        try:
            self.t = pd.read_table(filename,**kwargs)
        except Exception as e:
            print('ERROR: could not read %s!' % filename)
            if raiseError:
                raise RuntimeError(str(e))
            return(1)
        
        if test4commentedheader:
            # This is to make it compatible to my old-style commented header files!
            # commented header: make sure it doesn't count the '#' as a column!
            if self.t.columns[0]=='#':
                renamemapping = {}
                for i in range(len(self.t.columns)-1):
                    renamemapping[self.t.columns[i]]=self.t.columns[i+1]
                renamemapping[self.t.columns[-1]]='__delme'
                self.t = self.t.rename(columns=renamemapping)
                self.t = self.t.drop(columns=['__delme'])

            # make sure the # is not kept in column name!
            if self.t.columns[0][0]=='#':
                self.t = self.t.rename(columns={self.t.columns[0]:self.t.columns[0][1:]})
            
        self.formattable(namesMapping=namesMapping,roundingMapping=roundingMapping)
        
        
        return(0)

    def write(self,filename=None,indices=None,columns=None,formatters=None,raiseError=True,overwrite=True,verbose=False, 
              index=False, makepathFlag=True,**kwargs):

        # make sure indices are converted into a valid list
        indices=self.getindices(indices)
        
        # make sure columns are converted into a valid list
        columns=self.getcolnames(columns)

        # make the path to the file if necessary
        if not (filename is None):
            if makepathFlag:
                if self.verbose: print('Clobbering %s' % filename)
                if makepath4file(filename,raiseError=False):
                    errorstring='ERROR: could not make directory for %s' % filename
                    if raiseError:
                        raise RuntimeError(errorstring)
                    #print(errorstring)
                    return(1)

            # if overwrite, then remove the old file first
            if os.path.isfile(filename):
                if overwrite:
                    os.remove(filename)
                    if os.path.isfile(filename):
                        errorstring='ERROR: could not clobber %s' % filename
                        if raiseError:
                            raise RuntimeError(errorstring)
                        #print(errorstring)
                        return(2)
                else:
                    print('Warning: file exists, not deleting it, skipping! if you want to overwrite, use overwrite option!')
                    return(0)
        
        # Fix the dtypes if wanted
        if self.auto_convert_dtypes:
            self.t=self.t.convert_dtypes()
        if not(self.default_dtypeMapping is None):
            self.formattable(dtypeMapping=self.default_dtypeMapping)            
        
        # if both formatters and self.defaultformatters are None, then no formatting. formatters supersedes self.defaultformatters
        if formatters is None:
            formatters = self.default_formatters
        
        if indices is None:
            # no indices are passed
            if self.verbose and not(filename is None): print('Saving %d rows into %s' % (len(self.t),filename))
            if len(self.t)==0:
                # just save the header
                if filename is None:
                    print(' '.join(columns)+'\n')
                else:
                    open(filename,'w').writelines(' '.join(columns)+'\n')
            else:
                if filename is None:
                    pass
                    print(self.t.to_string(index=index, columns=columns, formatters=formatters, **kwargs))
                else:
                    self.t.to_string(filename, index=index, columns=columns, formatters=formatters, **kwargs)
        else:
            if self.verbose and not(filename is None): print('Saving %d rows into %s' % (len(indices),filename))
            if len(indices)==0:
                # just save the header
                if filename is None:
                    print(' '.join(columns)+'\n')
                else:
                    if columns is None:
                        columns = []
                    open(filename,'w').writelines(' '.join(columns)+'\n')
            else:
                if filename is None:
                    print(self.t.loc[indices].to_string(index=index, columns=columns, formatters=formatters, **kwargs))
                else:
                    self.t.loc[indices].to_string(filename, index=index, columns=columns, formatters=formatters, **kwargs)

        if not (filename is None):
            # some extra error checking...
            if not os.path.isfile(filename):
                errorstring='ERROR: could not save %s' % filename
                if raiseError:
                    raise RuntimeError(errorstring)
                #print(errorstring)
                return(3)
                
        return(0)
        
    def formattable(self,namesMapping=None,roundingMapping=None,dtypeMapping=None):
             
        if not(namesMapping is None):
            self.t = self.t.rename(columns=namesMapping)
            
        if not(roundingMapping is None):
            self.t = self.t.round(roundingMapping)

        if not(dtypeMapping is None):
            for col in dtypeMapping:
                self.t[col] = self.t[col].astype(dtypeMapping[col])
            
        return(0)
    
    def getindices(self,indices=None):
        """make indices conform (input can be None,([list],), int, str, or list). The output is a list """
        
        #If indices is tuple, return the first entry which is the array list of indicies
        if indices is None: 
            return(self.t.index.values)
        
        # If indices=([indiceslist],), then it needs to be changed to indices=indiceslist
        if isinstance(indices,tuple):
            if len(indices)==0:
                #return an empty list instead of empty tuple
                return([])
            # teh first entry is a list, then these are the relevant indices!
            if isinstance(indices[0],list) or isinstance(indices[0],np.ndarray):
                return(indices[0])
            else:
                return(list(indices))
        
        # if the passed value is an integer or str, make it a list!
        if isinstance(indices,int) or isinstance(indices,str) or isinstance(indices,float):
            return([indices])           
        
        indices=np.array(indices)
        
        #if not (isinstance(indices,list) or isinstance(indices,np.ndarray)):
        #    raise RuntimeError("Can't convert this to an indices list!",type(indices),indices)
            
        return(indices)
    
    def getcolnamesDELME(self,colnames=None):
        """Return a list of all colnames of colnames=None. If colnames=string, return a list"""
        if (colnames is None) or colnames.lower()=='all':
            colnames = self.t.columns
        else:
            if isinstance(colnames,str):
                colnames=[colnames]
        return(colnames)
            

    def getcolnames(self,colnames=None):
        """Return a list of all colnames of colnames=None. If colnames=string, return a list"""
        if (colnames is None):
             colnames = self.t.columns[:] 
        elif isinstance(colnames,str):
            if colnames.lower()=='all':
                colnames = self.t.columns[:]
            else:
                colnames=[colnames]
        return(colnames)
            

    def ix_remove_null(self,colnames=None,indices=None):
        # get the indices based on input.
        indices=self.getindices(indices)
        
        # get the column names over which to iterate
        colnames=self.getcolnames(colnames)
        
        for colname in colnames:
            #print('XXX',indices)
            (notnull,) = np.where(pd.notnull(self.t.loc[indices,colname]))
            indices = indices[notnull]
            #print('YYY',notnull)
        return(indices)

    def ix_equal(self,colnames,val,indices=None):
        # get the indices based on input.
        indices=self.getindices(indices)
        
        # get the column names over which to iterate
        colnames=self.getcolnames(colnames)
        for colname in colnames:
            (keep,) = np.where(self.t.loc[indices,colname].eq(val))
            indices = indices[keep]
            
        return(indices)
        
    def ix_inrange(self,colnames=None,lowlim=None,uplim=None,indices=None,
                   exclude_lowlim=False,exclude_uplim=False):

        # get the indices based on input.
        indices=self.getindices(indices)
        
        # get the column names over which to iterate
        colnames=self.getcolnames(colnames)
        #print(colnames)
        for colname in colnames:
            if not(lowlim is None):
                if exclude_lowlim:
                    (keep,) = np.where(self.t.loc[indices,colname].gt(lowlim))
                else:
                    (keep,) = np.where(self.t.loc[indices,colname].ge(lowlim))
                indices = indices[keep]
                #print('lowlim cut:',keep)

            if not(uplim is None):
                if exclude_uplim:
                    (keep,) = np.where(self.t.loc[indices,colname].lt(uplim))
                else:
                    (keep,) = np.where(self.t.loc[indices,colname].le(uplim))
                indices = indices[keep]
                #print('uplim cut:',keep)
        return(indices)
    
    def ix_outrange(self,colnames=None,lowlim=None,uplim=None,indices=None,
                    exclude_lowlim=False,exclude_uplim=False):

        # get the indices based on input.
        indices=self.getindices(indices)        
        
        # get the column names over which to iterate
        colnames=self.getcolnames(colnames)
        
        #print('BBB',indices)
        for colname in colnames:
            if not(lowlim is None):
                if exclude_lowlim:
                    (keeplow,) = np.where(self.t.loc[indices,colname].lt(lowlim))
                else:
                    (keeplow,) = np.where(self.t.loc[indices,colname].le(lowlim))
                #print('lowlim cut:',keeplow)
            else:
                keeplow=[]

            if not(uplim is None):
                if exclude_uplim:
                    (keepup,) = np.where(self.t.loc[indices,colname].gt(uplim))
                else:
                    (keepup,) = np.where(self.t.loc[indices,colname].ge(uplim))
                #print('uplim cut:',keepup)
            else:
                keepup=[]
                
            indices = indices[AorB(keeplow,keepup)]
            
        return(indices)
    
    def ix_unmasked(self,maskcol,maskval=None,indices=None):

        # get the indices based on input.
        indices=self.getindices(indices)  
        
        if maskval is None:
            (keep,) = np.where(self.t.loc[indices,maskcol].eq(0))
        else:
            (keep,) = np.where(bitmask.bitfield_to_boolean_mask(self.t.loc[indices,maskcol].astype('int'),ignore_flags=~maskval,good_mask_value=True))
        indices = indices[keep]
        return(indices)           
    
    def ix_masked(self,maskcol,maskval=None,indices=None):

        # get the indices based on input.
        indices=self.getindices(indices)  
        
        if maskval is None:
            (keep,) = np.where(self.t.loc[indices,maskcol].ne(0))
        else:
            (keep,) = np.where(bitmask.bitfield_to_boolean_mask(self.t.loc[indices,maskcol].astype('int'),ignore_flags=~maskval))
        indices = indices[keep]
        return(indices)    
    
    def ix_sort_by_cols(self,cols,indices=None):

        # get the indices based on input.
        indices=self.getindices(indices)  
        
        # get the column names (makes sure that it is a list)
        cols=self.getcolnames(cols)

        ix_sorted = self.t.loc[indices].sort_values(cols).index.values

        return(ix_sorted)

    def newrow(self,dicti=None):
        self.t = self.t.append(dicti,ignore_index=True)
        index = self.t.index.values[-1]
        return(index)
        
    def fitsheader2table(self,fitsfilecolname,indices=None,requiredfitskeys=None,optionalfitskey=None,raiseError=True,skipcolname=None,headercol=None):

        indices = self.getindices(indices)        

        # initialize columns if necessary
        if requiredfitskeys!=None:
            for fitskey in requiredfitskeys:
                if not (fitskey in self.t.columns):
                    self.t[fitskey]=None
        if optionalfitskey!=None:
            for fitskey in optionalfitskey:
                if not (fitskey in self.t.columns):
                    self.t[fitskey]=None

        if headercol!=None and (not (headercol in self.t.columns)):
            self.t[headercol]=None

        # loop through the images
        for index in indices:
            header = fits.getheader(self.t[fitsfilecolname][index])
            if headercol!=None:
                self.t[headercol]=header
                
            if requiredfitskeys!=None:
                for fitskey in requiredfitskeys:
                    if fitskey in header:
                        self.t[fitskey][index]=header[fitskey]
                    else:
                        if raiseError:
                            raise RuntimeError("fits key %s does not exist in file %s" % (fitskey,self.t[fitsfilecolname][index]))
                        else:
                            self.t[fitskey][index]=None
                            if skipcolname!=None:
                                 self.t[skipcolname][index]=True
                                 
            if optionalfitskey!=None:
                for fitskey in optionalfitskey:
                    if fitskey in header:
                        self.t[fitskey][index]=header[fitskey]
                    else:
                        self.t[fitskey][index]=None

    def dateobs2mjd(self,dateobscol,mjdcol,timeobscol=None):
        if not (mjdcol in self.t.columns):
            self.t[mjdcol]=None

        if timeobscol!=None:
            dateobslist = list(self.t[dateobscol]+'T'+self.t[timeobscol])
        else:
            dateobslist = list(self.t[dateobscol])

        dateobjects = Time(dateobslist,  format='isot', scale='utc')
        mjds = dateobjects.mjd

        self.t[mjdcol]=mjds
        
    def initspline(self,xcol,ycol,indices=None,
                   kind='cubic',bounds_error=False,fill_value='extrapolate', 
                   **kwargs):
        if not (xcol in self.t.columns):
            raise RuntimeError("spline: x column %s does not exist in table!" % xcol)
        if not (ycol in self.t.columns):
            raise RuntimeError("spline: y column %s does not exist in table!" % ycol)

        # make sure there are no nan values
        indices = self.ix_remove_null(colnames=[xcol,ycol])        

        # initialize the spline and save it in self.spline with the key ycol
        self.spline[ycol]= interp1d(self.t.loc[indices,xcol],self.t.loc[indices,ycol],
                                    kind=kind,bounds_error=bounds_error,fill_value=fill_value,**kwargs)
        
    def getspline(self,xval,ycol):
        if not(ycol in self.spline):
            raise RuntimeError('Spline for column %s is not defined!' % ycol)
        return(self.spline[ycol](xval))
    
class pdastrostatsclass(pdastroclass):
    def __init__(self,**kwargs):
        pdastroclass.__init__(self,**kwargs)
        self.reset()
        self.set_statstring_format()
        self.c4_smalln = [0.0, 0.0, 0.7978845608028654, 0.8862269254527579, 0.9213177319235613, 0.9399856029866251, 0.9515328619481445]

    def c4(self,n):
        #http://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation
        if n<=6:
            return(self.c4_smalln[n])
        else:
            return(1.0 - 1.0/(4.0*n) - 7.0/(32.0*n*n) - 19.0/(128.0*n*n*n))
        
    def reset(self):
        self.statparams = {}
        for k  in ['mean','mean_err','stdev','stdev_err','X2norm','ix_good','ix_clip']:
            self.statparams[k]=None
        for k  in ['Ngood','Nclip','Nchanged','Nmask','Nnan','converged','i']:
            self.statparams[k]=0
        self.calc_stdev_X2_flag = True

    def set_statstring_format(self,format_floats='{:f}',format_ints='{:d}',format_none='{}', format_X2norm='{:.2f}'):
        self.str_format1 = "mean:%s(%s) stdev:%s(%s) X2norm:%s Nchanged:%s Ngood:%s Nclip:%s" % (format_floats,format_floats,format_floats,format_floats,format_X2norm,format_ints,format_ints,format_ints)
        self.str_format_none = "mean:%s(%s) stdev:%s(%s) X2norm:%s Nchanged:%s Ngood:%s Nclip:%s" % (format_none,format_none,format_none,format_none,format_none,format_none,format_none,format_none)
        #self.str_format2 = "mean:%s(%s) stdev:%s X2norm:%s Nchanged:%s Ngood:%s Nclip:%s" % (format_floats,format_floats,format_floats,format_floats,format_ints,format_ints,format_ints)

    def statstring(self):
        if self.statparams['mean'] is None or self.statparams['stdev'] is None:
            formatstring = "WARNING! i:{:02d} "+ self.str_format_none
        else:            
            formatstring = "i:{:02d} "+ self.str_format1

        s = formatstring.format(self.statparams['i'],self.statparams['mean'],self.statparams['mean_err'],
                                self.statparams['stdev'],self.statparams['stdev_err'],self.statparams['X2norm'],
                                self.statparams['Nchanged'],self.statparams['Ngood'],self.statparams['Nclip'])
        return(s)
        
        

    def calcaverage_errorcut(self,datacol, noisecol, indices=None, 
                             mean=None,Nsigma=None,medianflag=False,
                             return_ix=False, verbose=0):

        # get the indices based on input.
        indices=self.getindices(indices)
            
        # If N-sigma cut and second iteration (i.e. we have a stdev from the first iteration), skip bad measurements.
        if not(Nsigma is None) and not(mean is None):
            ix_good_bkp = copy.deepcopy(self.statparams['ix_good'])  
            (keep,) = np.where(np.absolute(self.t.loc[indices,datacol]-mean)<=Nsigma*self.t.loc[indices,noisecol])
            ix_good  = indices[keep]
        else:
            ix_good_bkp = None
            ix_good  = indices
 
        
        if verbose>3 and not(ix_good_bkp is None) and not(ix_good is None):
            print('{} good data after sigma clipping, {} clipped'.format(len(ix_good),len(indices)-len(ix_good)))
            self.write(indices=ix_good)
            #self.write(columns=['objID','filter','psfFlux','psfFluxErr','psfMag','psfMagErr','psfMagErr_tot'],indices=ix_good)
 
        Ngood = len(ix_good)      
        if Ngood>1:
            if medianflag:
                mean = self.t.loc[ix_good,datacol].median()
                if verbose>1: print('median: {:f}'.format(mean))
                stdev =  np.sqrt(1.0/(Ngood-1.0)*np.sum(np.square(self.t.loc[ix_good,datacol] - mean)))/self.c4(Ngood)
                mean_err = stdev/np.sqrt(Ngood-1)
                
            else:
                c1 = np.sum(1.0*self.t.loc[ix_good,datacol]/np.square(self.t.loc[ix_good,noisecol]))
                c2 = np.sum(1.0/np.square(self.t.loc[ix_good,noisecol]))
                mean = c1/c2
                mean_err = np.sqrt(1.0/c2)
                stdev = self.t.loc[ix_good,datacol].std()
                
            stdev_err = 1.0*stdev/np.sqrt(2.0*Ngood)
            X2norm = 1.0/(Ngood-1.0)*np.sum(np.square((self.t.loc[ix_good,datacol] - mean)/self.t.loc[ix_good,noisecol]))                
                
        else:
            if Ngood==1:
                mean = self.t.loc[ix_good[0],datacol]*1.0
                mean_err = self.t.loc[ix_good[0],noisecol]*1.0
            else:
                mean = None
                mean_err = None
                
            X2norm   = None
            stdev     = None
            stdev_err = None
            
        self.statparams['ix_good']=ix_good
        self.statparams['Ngood']=Ngood
        self.statparams['ix_clip']=AnotB(indices,ix_good)
        self.statparams['Nclip']=len(indices) - Ngood
        if not(ix_good_bkp is None):
            self.statparams['Nchanged'] = len(not_AandB(ix_good_bkp,ix_good))            
        else:
            self.statparams['Nchanged'] = 0
        
        self.statparams['mean']      = mean    
        self.statparams['stdev']     = stdev    
        self.statparams['mean_err']  = mean_err
        self.statparams['stdev_err'] = stdev_err
        self.statparams['X2norm']    = X2norm

        if Ngood<1:
            return(1)
        return(0)

    def calcaverage_sigmacut(self,datacol, noisecol=None, indices=None, 
                             mean=None,stdev=None,Nsigma=None,
                             percentile_cut=None,percentile_Nmin=3,
                             medianflag=False,
                             return_ix=False, verbose=0,rescol='__tmp_residuals'):

        # get the indices based on input.
        indices=self.getindices(indices)
        if len(indices)==0:
            print('Warning!! no data passed for sigma cut!')
            self.reset()
            return(2)
        
        ix_good_bkp = None
        if (percentile_cut is None) or (len(indices)<=percentile_Nmin):
            # If N-sigma cut and second iteration (i.e. we have a stdev from the first iteration), skip bad measurements.
            if not(Nsigma is None) and not(stdev is None) and not(mean is None):
                ix_good_bkp = copy.deepcopy(self.statparams['ix_good'])    
                (keep,) = np.where(np.absolute(self.t.loc[indices,datacol]-mean)<=Nsigma*stdev)
                ix_good  = indices[keep]
                if verbose>3:
                    print('good data after sigma clipping:')
                    self.write(indices=ix_good)
                    #self.write(columns=['objID','filter','psfFlux','psfFluxErr','psfMag','psfMagErr'],indices=ix_good)
            else:
                if verbose>3:
                    print('No sigma clipping yet...')
                ix_good  = indices        
        else:
            # percentile clipping!!!
            if mean is None:
                if medianflag:
                    median = self.t.loc[indices,datacol].median()
                    if verbose>1: print('median: {:f}'.format(median))
                    mean = median
                else:
                    mean = self.t.loc[indices,datacol].mean()
                    if verbose>1: print('mean: {:f}'.format(mean))
            
            self.t[rescol]=np.absolute(self.t.loc[indices,datacol]-mean)
            max_residual = np.percentile(np.absolute(self.t.loc[indices,datacol]-mean),percentile_cut)
            if verbose: 
                print('%f percentile cut: max residual for cut: %f' % (percentile_cut,max_residual))
            ix_good  = self.ix_inrange(rescol,None,max_residual,exclude_uplim=True)
            

            #print('all')
            #self.write(columns=['objID','psfMag','psfMagErr',rescol],indices=indices)
            #print('good')
            #self.write(columns=['objID','psfMag','psfMagErr',rescol],indices=ix_good)

            # make sure the percentile clipping does not cut too much away for small numbers: always keep the best percentile_Nmin
            if len(ix_good)<percentile_Nmin:
                print('Warning: %d<%d made it through the percentile cut, too few, taking the %d with the least residuals' % (len(ix_good),percentile_Nmin,percentile_Nmin))
                residuals = np.sort(self.t.loc[indices,rescol])
                max_residual = residuals[percentile_Nmin-1]
                ix_good  = self.ix_inrange(rescol,None,max_residual,exclude_uplim=False)
                if len(ix_good)<percentile_Nmin:
                    raise RuntimeError('%d<%d in percentile cut!' % (len(ix_good),percentile_Nmin))

            if verbose>3:
                print('good data after percentile clipping:')
                self.write(indices=ix_good)
                #self.write(columns=['objID','filter','psfFlux','psfFluxErr','psfMag','psfMagErr',rescol],indices=ix_good)


            #sys.exit(0)
            
            #residuals = np.sort(np.absolute(self.t.loc[indices,datacol]-mean))
            #print(residuals)
            #max_residual = np.percentile(np.absolute(self.t.loc[indices,datacol]-mean),percentile_cut)

            #sys.exit(0)
            
 
        Ngood = len(ix_good)      
        if Ngood>1:
            if medianflag:
                median = self.t.loc[ix_good,datacol].median()
                #mean = scipy.median(self.t.loc[ix_good,datacol])
                if verbose>1: print('median: {:f}'.format(median))
                stdev =  np.sqrt(1.0/(Ngood-1.0)*np.sum(np.square(self.t.loc[ix_good,datacol] - median)))/self.c4(Ngood)
                mean = median
            else:
                mean = self.t.loc[ix_good,datacol].mean()
                if verbose>1: print('mean: {:f}'.format(mean))
                stdev = self.t.loc[ix_good,datacol].std()
                
            mean_err = stdev/np.sqrt(Ngood-1)
            stdev_err = 1.0*stdev/np.sqrt(2.0*Ngood)
            if noisecol is None:
                X2norm = 1.0/(Ngood-1.0)*np.sum(np.square((self.t.loc[ix_good,datacol] - mean)/stdev))                
            else:
                X2norm = 1.0/(Ngood-1.0)*np.sum(np.square((self.t.loc[ix_good,datacol] - mean)/self.t.loc[ix_good,noisecol]))     
        else:
            if Ngood==1:
                mean = self.t.loc[ix_good[0],datacol]*1.0
                mean_err = self.t.loc[ix_good[0],noisecol]*1.0
            else:
                mean = None
                mean_err = None
                
            X2norm   = None
            stdev     = None
            stdev_err = None
           
        self.statparams['ix_good']=ix_good
        self.statparams['Ngood']=Ngood
        self.statparams['ix_clip']=AnotB(indices,ix_good)
        self.statparams['Nclip']=len(indices) - Ngood
        if not(ix_good_bkp is None):
            self.statparams['Nchanged'] = len(not_AandB(ix_good_bkp,ix_good))            
        else:
            self.statparams['Nchanged'] = 0
        
        self.statparams['mean']      = mean    
        self.statparams['stdev']     = stdev    
        self.statparams['mean_err']  = mean_err
        self.statparams['stdev_err'] = stdev_err
        self.statparams['X2norm']    = X2norm

        if Ngood<1:
            return(1)
        return(0)

    def calcaverage_sigmacutloop(self,datacol, indices=None, noisecol=None, maskcol=None, maskval=None, 
                                 removeNaNs = True,
                                 Nsigma=3.0,Nitmax=10,verbose=0,
                                 percentile_cut_firstiteration=None,
                                 median_firstiteration=True):
        """
        mask must have same dimensions than data. If mask[x]=True, then data[x] is not used.
        noise must have same dimensions than data. If noise != None, then the error weighted mean is calculated.
        if saveused, then self.use contains array of datapoints used, and self.clipped the array of datapoints clipped
        median_firstiteration: in the first iteration, use the median instead the mean. This is more robust if there is a population of bad measurements
        """

        # get the indices based on input.
        indices=self.getindices(indices)
        
        # remove null values if wanted
        if removeNaNs:
            print('REMOVING NANS!!!')
            colnames = [datacol]
            if not(noisecol is None): colnames.append(noisecol)
            if not(maskcol is None): colnames.append(maskcol)
            Ntot = len(indices)
            indices = self.ix_remove_null(colnames,indices=indices)
            self.statparams['Nnan']= Ntot-len(indices)
            if verbose>1: print('Keeping {:d} out of {:d}, skippin {:d} because of null values in columns {:s}'.format(len(indices),Ntot,Ntot-len(indices),",".join(colnames)))
        else:
            self.statparams['Nnan']= 0
            

        #self.write(columns=['objID','filter','psfFlux','psfFluxErr','psfMag','psfMagErr','mask'],indices=indices)
        # exclude data if wanted
        if maskcol!=None:
            Ntot = len(indices)
            print('maskval:',maskval)
            indices = self.ix_unmasked(maskcol,maskval=maskval,indices=indices)
            self.statparams['Nmask']= Ntot-len(indices)
            if verbose>1: print('Keeping {:d} out of {:d}, skippin {:d} because of masking in column {} (maskval={})'.format(len(indices),Ntot,Ntot-len(indices),maskcol,maskval))
        else:
            self.statparams['Nmask']= 0

        self.reset()
        while ((self.statparams['i']<Nitmax) or (Nitmax==0)) and (not self.statparams['converged']):
            # median only in first iteration and if wanted
            medianflag = median_firstiteration and (self.statparams['i']==0) and (Nsigma!=None)
            # percentile_cut only in first interation and if wanted
            percentile_cut = None
            if (self.statparams['i']==0):
                percentile_cut = percentile_cut_firstiteration

            if noisecol is None:
                errorflag = self.calcaverage_sigmacut(datacol, indices=indices, 
                                                      mean = self.statparams['mean'], stdev = self.statparams['stdev'],
                                                      Nsigma=Nsigma, 
                                                      medianflag=medianflag, percentile_cut=percentile_cut,
                                                      verbose=verbose)
            else:
                errorflag = self.calcaverage_errorcut(datacol, noisecol, indices=indices, 
                                                      mean = self.statparams['mean'],
                                                      Nsigma=Nsigma, medianflag=medianflag, verbose=verbose)
            #if self.statparams['i']==0:
            #    self.statparams['stdev']=0.05

            if verbose>2:
                print(self.statstring())
                
                
            # Not converged???
            if errorflag or self.statparams['stdev']==None or self.statparams['stdev']==0.0 or self.statparams['mean']==None:
                self.statparams['converged']=False
                break
            # Only do a sigma cut if wanted
            if Nsigma == None or Nsigma == 0.0:
                self.statparams['converged']=True
                break
            # No changes anymore? If yes converged!!!
            if (self.statparams['i']>0) and (self.statparams['Nchanged']==0):
                self.statparams['converged']=True
                break
            self.statparams['i']+=1
        
        if not(self.statparams['converged']):
            if self.verbose>1:
                print('WARNING! no convergence!')

        return(not self.statparams['converged'])
    
    def statresults2table(self,desttable,destindex=None,colmapping={},prefix='',suffix='',skipcolumns=[]):
        
        skipcolumns = self.getcolnames(skipcolumns)
        
        cols=[]
        vals=[]
        resultdict={}
        for k  in ['mean','mean_err','stdev','stdev_err','X2norm','Ngood','Nclip','Nmask','Nnan','converged','i']:
            if k in skipcolumns:
                continue
            
            # set outcol basename to
            outcol=k
            if k in colmapping:
                outcol=colmapping[k]
                
            outcol='{}{}{}'.format(prefix,outcol,suffix)
            if destindex is None:
                resultdict[outcol]:self.statparams[k]
            else:
                cols.append(outcol)
                vals.append(self.statparams[k])
        if destindex is None:
            outindex = self.newrow(resultdict)
        else:
            outindex = destindex
            print('Cols:',cols)
            print('Vals:',vals)
            desttable.loc[destindex,cols]=vals
            