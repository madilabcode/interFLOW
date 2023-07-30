import numpy as np
import pandas as pd
import math
from scipy import stats
from scipy.stats import rankdata, fisher_exact
import statsmodels.stats.multitest

###Ron Sheinin implamation of  the CERNO algotitam and signature score fucntions

def signature_by_express(exp,up,down = None):
  exp_up = exp.loc[exp.index.isin(up),:]
  scores_up = exp_up.apply(np.sum,axis=0)
  scores_up = np.array(scores_up)

  scores_up = np.array(scores_up)
  
  if down is not None:
    exp_down = exp.loc[exp.index.isin(up),:]
    scores_down = exp_down.apply(np.sum,axis=1)
  
    scores_down = np.array(scores_down)
    scores = scores_up - scores_down
    
    scores = pd.Series(scores)
    scores.index = exp.columns
    return scores
  
  scores_up = pd.Series(scores_up)
  scores_up.index = exp.columns
  return scores_up
  
def make_rank_dict(exp_named,names):
  exp_named.index = names
  exp_named = exp_named.sort_values(ascending = False)
  exp_dict = {}
  for name,value in exp_named.items():
    if value not in exp_dict:
      exp_dict[value] = []
    exp_dict[value].append(name)
  
  rank_dict = {}
  counter = 1
  for key, value in exp_dict.items():
    rank_dict[counter] = value
    counter += len(value)
    
  result = {}
  for key, value in rank_dict.items():
    for item in value:
      result[item] = key
  
  return result 
  
def rank_pvalue(gs,exp,names, pval = True):
  rank_dict = make_rank_dict(exp,names)
  N = len(exp.index)
  F = 0
  na = []
  
  for gene in gs:
    try:
       F += math.log(rank_dict[gene] / N)
    except:
      na.append(gene)
      
  F *= -2
  
  if pval:
    return 1 - stats.chi2.cdf(F,2*(len(gs) - len(na)))
     
  return F


def CERNO_alg(exp,tfs_targets):
  print("stop")
  tfs = tfs_targets.keys()
  array1 = [rank_pvalue(tfs_targets[key],exp,exp.index) for key in tfs]
  M1 = len(array1)
  array1 = list(map(lambda x: x[1],map(statsmodels.stats.multitest.fdrcorrection,array1)))
  pval_dicts = {key:array1[index][0] for index,key in enumerate(tfs)}
  return {key:value for key,value in pval_dicts.items() if value <= 0.05}
