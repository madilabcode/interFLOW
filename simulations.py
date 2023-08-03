import numpy as np
import pandas as pd
import os
#os.environ["R_HOME"] = r"C:\Program Files\R\R-4.2.1"
#os.environ['path'] += r";C:\Program Files\R\R-4.2.1"
from rpy2.robjects import r
from rpy2 import robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import default_converter
from Codes import pyDictToR as pyd
from Codes import flowGraph as tfg
#from Codes import neo_connection_db as neo
from Codes import utils
import rpy2.robjects.pandas2ri as rpyp
from sklearn.metrics import roc_auc_score, roc_curve, auc
import concurrent.futures
import pickle
import networkx as nx
import concurrent.futures
import matplotlib.pyplot as plt


def load_network():


    pa = pd.read_csv("files/humanProtinAction.csv")
    pi = pd.read_csv("files/humanProtinInfo.csv")
    tf_dict = pyd.main_py_to_R(pi, pa)[0]
    lr = pd.read_table(r"./files/LigandReceptorTableMouse.tsv",sep = "\t")
    receptors = lr.to.drop_duplicates().values
    tfs = list(tf_dict.keys())
    pa["Source"] = pa["Input-node Gene Symbol"]
    pa["Target"] = pa["Output-node Gene Symbol"]
    pa["capacity"] = pa["Edge direction score"]
    keeps = ["Source", "Target"]
    pa = pa[keeps]
    
    pa["Source"] = pa["Source"].apply(lambda x: x[0] + x[1:].lower()).astype(str)
    pa["Target"] = pa["Target"].apply(lambda x: x[0] + x[1:].lower()).astype(str)
    receptors = list(filter(lambda x: x in pa.Source.values or x in pa.Target.values, receptors))
    tfs = list(filter(lambda x: (x in pa.Source.values or x in pa.Target.values) and not x in receptors, tfs))

    return pa, receptors, tfs


def find_path(pa, receptors, tfs):

    g = nx.from_pandas_edgelist(pa, "Source", "Target", create_using=nx.DiGraph())    
    selcted_recep = []
    while len(selcted_recep) == 0:
        selcted_tf  = pd.Series(tfs).sample(1).values[0]
        paths = nx.single_source_shortest_path(g, selcted_tf)
        selcted_recep  = pd.Series(list(filter(lambda x: x in paths, receptors)))
        
    selcted_recep = selcted_recep.sample(5).values
    all_nodes = []
    for key in selcted_recep:
        all_nodes += paths[key]
    all_nodes = np.unique(all_nodes)
    return all_nodes
     

def build_syntatic_data(pa, path, latent_shape=100, corr_size=20, num_of_cells=1000):
    mus = np.random.normal(0,1,latent_shape)
    g = nx.from_pandas_edgelist(pa, "Source", "Target", create_using=nx.DiGraph())
    sigma =  np.random.uniform(low =0, high =0.7,size = (latent_shape, latent_shape) )
    sigma[:corr_size,:corr_size] = np.random.uniform(low =0.5, high =1,size = (corr_size, corr_size) )
    np.fill_diagonal(sigma,1)
    samples = np.random.multivariate_normal(mus, sigma, num_of_cells)
    nodes = list(g.nodes)
    indexs = [nodes.index(x) for x in path]
    A = np.random.normal(0,1, ( len(nodes), latent_shape))
    A[indexs,corr_size:] = 0
    A[indexs,:corr_size] = np.random.uniform(-0.5,5, ( len(indexs), corr_size))

    mask = np.arange(A.shape[0])
    mask = list(filter(lambda x: not x in indexs, mask))
    A[mask,:corr_size] = 0
    exp = pd.DataFrame(A@samples.T, index=nodes)
    exp = ((exp.T - exp.mean(axis=1)) / exp.std(axis=1)).T
    exp = 0.7*exp + 0.3* np.random.normal(0,1,exp.shape)
    exp *= np.random.binomial(1,0.2,exp.shape)
    return exp


def run_anaylsis(exp, pa, tfs, receptors):
    gpf = tfg.build_flowing_network_with_normlized_wights(pa, tfs, exp,wights_flag=True)
    flow_dict={recp:0 for recp in receptors}
    gd = tfg.FlowGraph(flow_dict, gpf,False, sim_flag=True)
    nodes = list(gd.multy_flow_dict.keys())
    node_flow = np.array([sum(dict(vec).values()) for vec in gd.multy_flow_dict.values()])
    flow_df = pd.DataFrame({"node":nodes, "flow":node_flow})
    flow_df = flow_df.loc[~flow_df.node.isin(["source_node","sink"])]
    return flow_df


def calculate_roc(flow_df, path):
    labels = flow_df.node.isin(path)
    fpr, tpr, _= roc_curve (labels, flow_df.flow, pos_label=1)
    auc = roc_auc_score (labels, flow_df.flow)
    return fpr, tpr, auc

def pipeline(args):
    pa, receptors, tfs = args
    path = find_path(pa, receptors, tfs)
    exp = build_syntatic_data(pa, path)
    flow_df = run_anaylsis(exp,  pd.read_csv("files/humanProtinAction.csv"), tfs, receptors)
    roc = calculate_roc(flow_df, path)
    return roc

def main(num_of_repat=10):
    pa, receptors, tfs = load_network()
    with concurrent.futures.ProcessPoolExecutor(max_workers=20) as executor:
        result = list(executor.map(pipeline,[( pa, receptors, tfs) for _ in range(num_of_repat)]))
    result = np.array(result)
    mean_fpr = np.linspace(0, 1, 1000)
    fpr_list = []
    tpr_list = []
    auc_list = []
    for fpr, tpr in zip(result[:,0], result[:,1]):
        interp_tpr = np.interp(mean_fpr, fpr, tpr)
        interp_tpr[0] = 0.0
        fpr_list.append(fpr)
        tpr_list.append(interp_tpr)
        auc_list.append(auc(fpr, tpr))
    plot_k_fold_roc(fpr_list, tpr_list, auc_list)

def plot_k_fold_roc(fpr_list, tpr_list, auc_list):

    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
            label='Chance', alpha=.8)

    mean_tpr = np.mean(tpr_list, axis = 0)
    mean_fpr = np.linspace(0, 1, 1000)
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(auc_list)

    plt.plot(mean_fpr, mean_tpr, color='b',
            label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
            lw=2, alpha=.8)

    std_tpr = np.std(tpr_list, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                    label=r'$\pm$ 1 std. dev.')

    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.xlabel('False Positive Rate',fontsize=18)
    plt.ylabel('True Positive Rate',fontsize=18)
    plt.title('Cross-Validation ROC',fontsize=18)
    plt.legend(loc="lower right", prop={'size': 15})
    plt.show()


if __name__ == "__main__":
    main()