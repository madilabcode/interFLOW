import pandas as pd
import numpy as np
import random
import math
import scipy.stats as st
import zlib, json, base64
from rpy2.robjects import r
from rpy2 import robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import default_converter
import rpy2.robjects.pandas2ri as rpyp
import scipy.stats as sc
from matplotlib import pyplot as plt
import seaborn as sns
import pickle
import dill

def intersection(lst1, lst2):
    return [value for value in lst1 if value in lst2]


def reduce_table_by_key(table, key, Col_Name):
    table["is_Key"] = table[Col_Name].apply(lambda name: 1 if str(name).find(key) != -1 else 0)
    table = table.loc[table["is_Key"] == 1, :]
    table = table.drop(axis=1, labels="is_Key")
    return table


def sum_avg_ligand(ligs, exp):
    sub_exp = exp.loc[exp.index.isin(ligs), :]
    sub_exp_Avg = sub_exp.apply(np.mean, axis=1)
    return np.sum(sub_exp_Avg)


def dsa_table_update_to_lig(toExp, fromExp, DSA_Table, legRet):
    toExp = toExp.loc[toExp.index.isin(np.unique(DSA_Table["Recp"])), :]

    for rec in DSA_Table.index:
        sum_lig = sum_avg_ligand(list(legRet.loc[legRet["Receptor"] == rec, "Ligand"]), fromExp)
        DSA_Table.loc[rec, "DSA"] *= sum_lig

    for gene in toExp.index:
        toExp.loc[gene, :] = np.log1p(toExp.loc[gene, :]) * DSA_Table.loc[gene, "DSA"]

    return toExp


def dsa_with_lig(toExp, fromExp, DSA_Table, legRet):
    toExp = toExp.loc[toExp.index.isin(np.unique(DSA_Table["Recp"])), :]

    for rec in DSA_Table["Recp"]:
        sum_lig = sum_avg_ligand(list(legRet.loc[legRet["Receptor"] == rec, "Ligand"]), fromExp)
        DSA_Table.loc[DSA_Table["Recp"] == rec, "DSA"] *= sum_lig

    for gene in toExp.index:
        toExp.loc[gene, :] = np.log1p(toExp.loc[gene, :]) * float(DSA_Table.loc[DSA_Table["Recp"] == gene, "DSA"])

    return toExp.apply(np.sum, axis=0)


def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        dill.dump(obj, f, dill.HIGHEST_PROTOCOL)


def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return dill.load(f)

def sample_idents(exp,idents,sample_size):
    idents_count = idents.value_counts()
    idents_count = idents_count.apply(lambda x: max(sample_size,x))
    cts = lambda cell: idents_count[idents[cell]]
    sample_idents = pd.Series(idents.index).apply(lambda x: np.random.choice([1,0], size=1, p=[sample_size/cts(x), 1 -(sample_size/cts(x))])[0])
    sample_idents.index = idents.index
    return sample_idents[sample_idents == 1].index


def return_norm_express_per_gene(gene_col,value):
    loc, scale = sc.norm.fit(gene_col)
    return sc.norm.cdf(value,loc,scale)

def avg_dist(row,x,n):
    return sc.norm.cdf(x,row[0],row[1]/n)


def zero_inflated_cdf(x,row):
    if x == 0:
        return 0 
    n = len(row)
    non_zero_n  =len(row[row!=0])
    zero_n = len(row[row==0])
    return zero_n/n + non_zero_n/n*(sc.norm.cdf(x,*sc.norm.fit(row)))

def normalize_scale_exp(tr_path, co_path,sample_size = 200):
    with localconverter(default_converter + rpyp.converter):
        r("library(Seurat)")
        readRDS = r("readRDS")
        GetAssayData = r(f"function(obj) as.data.frame(obj[['RNA']]@data)")
        GetObjIdent = r("function(obj) as.character(obj@active.ident)")
        obj_tert = readRDS(f"InputObj/{tr_path}")
        obj_con = readRDS(f"InputObj/{co_path}")
        tr_exp = GetAssayData(obj_tert) 
        co_exp = GetAssayData(obj_con) 
        tr_idents = pd.Series(GetObjIdent(obj_tert))
        tr_idents.index = tr_exp.columns
        co_idnets = pd.Series(GetObjIdent(obj_con))
        co_idnets.index = co_exp.columns

    tr_sample = sample_idents(tr_exp,tr_idents,sample_size)
    co_sample = sample_idents(co_exp, co_idnets,sample_size)
    frames = [tr_exp[tr_sample],co_exp[co_sample]]
    unite_exp = pd.concat(frames, axis=1)
    #dist_args = unite_exp.apply(lambda row : lambda x: zero_inflated_cdf(x,row),axis=1)
    dist_args = dict(unite_exp.apply(lambda row : (row.mean(),row.std()),axis=1))
    save_obj(dist_args,r"./outputObj/dist_norm_args")
    return dist_args

def save_to_csv(table, name):
    table.to_csv(name)

def json_zip(js, ZIPJSON_KEY='base64(zip(o))'):
    js = {
        ZIPJSON_KEY: base64.b64encode(
            zlib.compress(
                json.dumps(js).encode('utf-8')
            )
        ).decode('ascii')
    }
    return json.dumps(js)


def json_unzip(zp, ZIPJSON_KEY='base64(zip(o))'):
    zp = json.loads(zp)
    assert (ZIPJSON_KEY in zp.keys())
    try:
        js = zlib.decompress(base64.b64decode(zp[ZIPJSON_KEY]))
    except:
        raise RuntimeError("Could not decode/unzip the contents")
    try:
        js = json.loads(js)
    except:
        raise RuntimeError("Could interpret the unzipped contents")
    return js
