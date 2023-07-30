import numpy as np
import pandas as pd
import os
os.environ["R_HOME"] = r"C:\Program Files\R\R-4.2.1"
os.environ['path'] += r";C:\Program Files\R\R-4.2.1"
from rpy2.robjects import r
from rpy2 import robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import default_converter
from Codes import pyDictToR as pyd
from Codes import flowGraph as tfg
#from Codes import neo_connection_db as neo
from Codes import utils
import rpy2.robjects.pandas2ri as rpyp
import concurrent.futures
import pickle

assay = "data"

def LRP(path, to_name, from_name, plot_path=None, thrshold=0.1, per_claster=False):
    with localconverter(default_converter + rpyp.converter):
        print(f"{to_name}_{from_name}")
        r["source"]("Codes/Ligand_Receptor_pipeline.R")
        createLRtable = r("createLRtable")
        createCircosPlots = r("createCircosPlots")
        DElegenedNcolor = r("DElegenedNcolor")
        DSAlegenedNcolor = r("DSAlegenedNcolor")
        DSA_PLOT_TSNE = r("DSA_PLOT_TSNE")
        readRDS = r("readRDS")
        subset = r("subset")
        FindMarkers = r("function(obj,id) FindMarkers(object = obj,ident.1=id,only.pos = TRUE,slot='data',verbose = FALSE)")


        ProtNamesInfo = pd.read_csv("files/humanProtinInfo.csv")
        ProtActions = pd.read_csv("files/humanProtinAction.csv")
        tf_dict = pyd.main_py_to_R(ProtNamesInfo, ProtActions)[0]
      
        obj = readRDS(f"{path}")
        obj = subset(obj, ident=[to_name, from_name])
        to_exp, from_exp = ExpTables(obj, to_name, from_name)
      
        lr = dict(createLRtable(obj, to_exp, from_exp, from_name, to_name, assay, thrshold=thrshold))
        if len(lr) != 3:
            return None

        lr, markerall_ligand, markerall_receptor = lr.values()
        de = dict(DElegenedNcolor(obj, from_name, to_name, lr, markerall_ligand, markerall_receptor ))
       
        if len(de) == 0:
            return None

        dsa = tfg.DSA_anaylsis(to_exp, lr["Receptor"], ProtActions, ProtNamesInfo, tfs=tf_dict)
            
        if len(dsa) == 0:
            return None

        lege_color = DSAlegenedNcolor(dsa[0])
        #sang_recp = [rec for rec,pvalue in DSA_lst[1].significant_receptors.items() if pvalue <= 0.05]
        dsa_df = dsa[0]
        dsa_Df = dsa_df.loc[dsa_df["DSA"] < np.inf, :]
        obj = dsa_score_per_cell(obj, to_exp, dsa_df)


        createCircosPlots(to_exp, lr, ro.ListVector(de), ro.ListVector(dsa[1].tfs_per_recp_in_max_flow()), ro.ListVector(lege_color),
                              from_name, to_name, 'r', obj,de_recptors=dsa[1].significant_nodes())


        DSA_PLOT_TSNE(obj, from_name, to_name)
    return lr, dsa



def dsa_score_per_cell(obj, toExpression, DSA_Table, sacle_factor=1):
    mdata = r("function(obj,col='DSA_SCORE') ifelse(col %in% (obj@meta.data %>% names()),1,0)")
    mdata_col = r("function(obj,col='DSA_SCORE') obj@meta.data[col]")
    add_to_meta = r("""function(obj,values,cells,col='DSA_SCORE') {
            names(values) = cells
           obj[[col]] = values
           return (obj)
            }""")
    with localconverter(default_converter + rpyp.converter):
        toExpression = toExpression.loc[list(filter(lambda x: x in DSA_Table.Recp, toExpression.index))]
        for row in range(toExpression.shape[0]):
            toExpression.iloc[row, :] = (np.log1p(toExpression.iloc[row, :]) * float(
                DSA_Table.loc[DSA_Table.Recp == toExpression.index[row], :].DSA))

        if not mdata(obj)[0]:
            dsa_score = pd.Series(None, toExpression.columns)
        else:
            dsa_score = pd.Series(mdata_col(obj), toExpression.columns)

        dsa_score = pd.DataFrame(dsa_score)
        dsa_score.columns = ["DSA_SCORE"]
        dsa_per_cell = toExpression.apply(sum)
        dsa_per_cell = dsa_per_cell * sacle_factor
        dsa_per_cell = pd.DataFrame(dsa_per_cell)
        dsa_per_cell.columns = ["dsa_per_cell"]
        DSA_temp = pd.merge(dsa_score, dsa_per_cell, how="left", left_index=True, right_index=True)
        DSA_temp["DSA_SCORE"] = DSA_temp.apply(
            lambda x: x["dsa_per_cell"] if np.isnan(x["DSA_SCORE"]) else x["DSA_SCORE"], axis=1)
        return add_to_meta(obj, DSA_temp["DSA_SCORE"], DSA_temp["DSA_SCORE"].index)


def _helper_LRP(x):
    if len(x) == 3:
        return LRP(x[0], x[1], x[2])
    elif len(x) == 4:
        return LRP(x[0], x[1], x[2], plot_path=x[3])

def ExpTables(obj, toName, fromName, assay=assay):
    with localconverter(default_converter + rpyp.converter):
        r("library(Seurat)")
        subset = r("subset")
        GetAssayData = r(f"function(obj) obj[['RNA']]@{assay} %>% as.data.frame()")
        tosub = subset(obj, idents=toName)
        fromsub = subset(obj, idents=fromName)
        toExpression = GetAssayData(tosub)
        fromExpression = GetAssayData(fromsub)

        return toExpression, fromExpression



def run_pipline(args, objName, max_workers=1):
    r["source"]("Codes/Ligand_Receptor_pipeline.R")
    lr_list = {}
    dsa_tables = {}
    dsa_graph = tfg.GraphsDicts()

    # with BoundedProcessPoolExecutor(max_workers=max_workers) as executor:
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        result = list(executor.map(_helper_LRP, args))
    result = list(filter(lambda x: x is not None, result))

    with localconverter(default_converter + rpyp.converter):
        obj = r(f"readRDS('{args[0][0]}')")
        for index, arg in enumerate(args):
            if arg[1] not in lr_list:
                lr_list[arg[1]] = {}
                dsa_tables[arg[1]] = {}
                dsa_graph.add_obj(arg[1], tfg.GraphsDicts())

            if result[index] is not None:
                lr_list[arg[1]][arg[2]] = result[index][0]
                dsa_tables[arg[1]][arg[2]] = result[index][1][0]
                dsa_graph.return_GD_from_key(arg[1]).add_obj(arg[2], result[index][1][1])              

        tfg.save_obj(lr_list, f"./outputObj/legRetLists_{objName}")
        tfg.save_obj(dsa_tables, f"./outputObj/DSA_Tables_{objName}")
        tfg.save_obj(dsa_graph, f"./outputObj/DSA_Graphs_{objName}")

    # print([legRetList, DSA_Tables, DSA_Graphs, DSA_Mean])
    # del result
    r("rm(list = ls())")
    r("gc()")
    return lr_list,dsa_tables, dsa_graph




def main():
    try:
        conf = pd.read_csv("config.csv")
    except:
        raise Exception("no input CSV found")

    path = list(conf.loc[conf["Var"] == "obj", "Value"])[0]



    to_vec = list(conf.loc[conf["Var"] == "Receptors Clusters", "Value"])[0].split(" ")
    from_vec = list(conf.loc[conf["Var"] == "Ligand Clusters", "Value"])[0].split(" ")
    args = [(f"InputObj/{path}", toName, fromName) for toName in to_vec for fromName in from_vec if
           fromName != toName]

    output_obj_name = list(conf.loc[conf["Var"] == "Obj Name", "Value"])[0]
    lst_Tr = run_pipline(args, output_obj_name)
   
    return lst_Tr


def format_dsa_file(x):
    annot = pd.read_csv("Annot.csv")
    try:
        type_cluster = list(annot.loc[annot['Cluster'] == int(x), 'Type'])[0]
        subtype_cluster = list(annot.loc[annot['Cluster'] == int(x), 'subType'])[0]
    except:
        raise Exception(f"cluster {x} not found")

    if subtype_cluster is np.nan:
        return f"{x}/{type_cluster}"
    return f"{x}/{type_cluster}/{subtype_cluster}"
    

if __name__ == "__main__":
   print(os.getcwd())
   main()