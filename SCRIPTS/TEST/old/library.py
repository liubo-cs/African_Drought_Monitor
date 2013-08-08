import numpy as np

#Script to reproduce percentileofscore
def percentileofscore(scores,score):
 pct = []
 for i in range(scores.shape[0]):
  nall = scores[i,:].size
  #find number of scores below this score
  nbelow = scores[i,scores[i,:] <= score].size
  #find percentile
  pct.append(100*float(nbelow)/float(nall))
 
 return pct

#Script to reproduce scoreofpercentile
def scoreofpercentile():
 return
