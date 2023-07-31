import numpy as np
import pandas as pd
import os
from rpy2.robjects import r
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import default_converter
from Codes import flowGraph as tfg
import rpy2.robjects.pandas2ri as rpyp
import concurrent.futures

def tfs_importance():
    try:
        conf = pd.read_csv("config.csv")

    except:
        raise Exception("no input CSV found")

    obj_name = list(conf.loc[conf["Var"] == "Obj Name", "Value"])[0]

    graphs = tfg.load_obj(r"./outputObj/DSA_Graphs_" + obj_name)

    to_cluster = input("To Receptor")
    from_cluster = input("Ligand Cluster")
    graph_obj = graphs.return_GD_from_key(to_cluster).return_GD_from_key(from_cluster).update_obj()
    df = graph_obj.calculate_significant_tf_Multy_sinc()
    df.to_csv(f"./files/{from_cluster}_{to_cluster}_importance.csv")


def tf_flow_graph():
    r["source"]('Codes/Ligand_Receptor_pipeline.R')
    draw_flow_graph_from_df = r("draw_flow_graph_from_df")

    try:
        conf = pd.read_csv("config.csv")

    except:
        raise Exception("no input CSV found")

    obj_name = list(conf.loc[conf["Var"] == "Obj Name", "Value"])[0]

    graphs = tfg.load_obj(r"./outputObj/DSA_Graphs_" + obj_name)

    to_cluster = input("Receptor Cluster")
    from_cluster = input("Ligand Cluster")
    tfs = input("enter all tfs")
    tfs_list = tfs.split()
    graph_obj = graphs.return_GD_from_key(to_cluster).return_GD_from_key(from_cluster)

    for tf in tfs_list:
            df = graph_obj.calculate_max_flow_for_one_tf(tf)
            df["flow"] = df["flow"].apply(lambda x: round(x, 3))
            df["flow"] = df["flow"].astype(float)
            df.to_csv(f"./files/{from_cluster}_{to_cluster}_{tf}_flow_network.csv")

            with localconverter(default_converter + rpyp.converter):
                df["flow"] = df["flow"].astype(str)

if __name__ == "__main__":
    func = input("You want to make tf flow graph [f] or calc tf importance [i]: ")
    if func == "i":
        tfs_importance()
    else:
        tf_flow_graph()
    
