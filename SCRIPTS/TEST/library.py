import numpy as np

#Script to reproduce percentileofscore
def percentileofscore(scores,score):

 nall = scores.size
 #find number of scores below this score
 nbelow = scores[scores <= score].size
 #find percentile
 pct = 100*float(nbelow)/float(nall)
 
 return pct

#Script to reproduce scoreofpercentile
def scoreatpercentile(scores,pct):

 quantiles = []
 #Sort scores 
 quantiles = np.sort(scores)
 nbelow = 0
 nall = quantiles.size
 #Calculate percentiles
 percentiles = []
 for i in range(quantiles.size):
  percentiles.append(100*float(i)/float(nall-1))
 #Find where the percentile belongs
 for i in range(quantiles.size-1):
  if pct >= percentiles[i] and pct <= percentiles[i+1]:
   #Interpolate and break
   qnt = quantiles[i] + (quantiles[i+1] - quantiles[i])*(pct - percentiles[i])/(percentiles[i+1] - percentiles[i])
   break
 
 return qnt
