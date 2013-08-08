import scipy.stats as ss
import numpy as np
import time
import library_f90
import library
import random

ncols = 100
nrows = 10000
scores_all = np.random.rand(nrows,ncols)
score_array = np.random.rand(nrows)
#percentile of score
a = []
b = []
c = []
tic = time.clock()
for i in range(nrows):
 #provide a uniform random number
 a.append(library_f90.percentileofscore(scores_all[i,:],score_array[i],random.random()))
toc = time.clock()
time1 = toc - tic
#tic = time.clock()
#for i in range(nrows):
# b.append(library.percentileofscore(scores_all[i,:],score_array[i]))
#toc = time.clock()
time2 = toc - tic
tic = time.clock()
for i in range(nrows):
 c.append(ss.percentileofscore(scores_all[i,:],score_array[i],kind='mean'))
toc = time.clock()
time3 = toc - tic
print "new fortran: %f" % time1
print np.mean(np.array(a) - np.array(c))
#print "new python: %f" % time2
#print np.mean(np.array(b) - np.array(c))
print "old: %f" % time3
'''
#Score of percentile
a = []
b = []
c = []
pct_array = 100*np.random.rand(nrows)
tic = time.clock()
for i in range(nrows):
 a.append(library_f90.scoreatpercentile(scores_all[i,:],pct_array[i]))
toc = time.clock()
time1 = toc - tic
tic = time.clock()
for i in range(nrows):
 b.append(library.scoreatpercentile(scores_all[i,:],pct_array[i]))
toc = time.clock()
time2 = toc - tic
tic = time.clock()
for i in range(nrows):
 c.append(ss.scoreatpercentile(scores_all[i,:],pct_array[i]))
toc = time.clock()
time3 = toc - tic
print "new fortran: %f" % time1
print np.std((np.array(a) - np.array(c))**2)**0.5
print "new python: %f" % time2
print np.std((np.array(b) - np.array(c))**2)**0.5
print "old: %f" % time3
'''
