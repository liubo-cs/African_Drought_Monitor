import scipy.stats as ss
cimport numpy as np
import time
from cython.parallel import prange

#Script to reproduce percentileofscore
def percentileofscore(scores_all,scores_array):
 cdef int nbelow
 nrows = scores_all.shape[0]
 ncols = scores_all.shape[1]
 pct = []
 #for i in prange(nrows,nogil=True):#range(nrows):
 for i in range(nrows):
  #find number of scores below this score
  nbelow = 0
  for j in prange(ncols,nogil=True):
   if scores_all[i,j] <= scores_array[i]:
    nbelow = nbelow + 1
  #nbelow = scores_all[i,scores_all[i,:] <= scores_array[i]].size
  #find percentile
  pct.append(100*float(nbelow)/float(ncols))
 #nall = scores.size
 ##find number of scores below this score
 #nbelow = scores[scores <= score].size
 ##find percentile
 #pct = 100*float(nbelow)/float(nall)
 
 return pct

#Script to reproduce scoreofpercentile
def scoreofpercentile():
 return
