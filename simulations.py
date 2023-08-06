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
from scipy.stats import ranksums
from statsmodels.stats.multitest import fdrcorrection
import concurrent.futures
import pickle
import networkx as nx
import concurrent.futures
import matplotlib.pyplot as plt
import seaborn as sns

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

    return pa, receptors, tfs, lr

def find_path(pa, receptors, tfs, lr, num_of_recp = 40):
    g = nx.from_pandas_edgelist(pa, "Source", "Target", create_using=nx.DiGraph())    
    selcted_recep = []
    while len(selcted_recep) == 0:
        selcted_tf  = pd.Series(tfs).sample(5).values
        paths = nx.multi_source_dijkstra(g, set(selcted_tf))[1]
        selcted_recep  = pd.Series(list(filter(lambda x: x in paths, receptors)))
    
    print(selcted_tf)
    selcted_recep = selcted_recep.sample(num_of_recp).values
    ligands = lr.loc[(lr.to.isin(selcted_recep)) & (~lr["from"].isin(selcted_recep)),"from"].drop_duplicates()
    #ligands =  pd.Series(list(filter(lambda x: x in list(pa.Source.drop_duplicates()) + list(pa.Target.drop_duplicates()), ligands)))
    all_nodes = []
    for key in selcted_recep:
        all_nodes += paths[key]
    all_nodes = np.unique(all_nodes)
    return all_nodes, selcted_recep, ligands , selcted_tf
     
def build_syntatic_data(pa, path, de_up=None, de_down=None,de_size=10,latent_shape=100, corr_size=20, num_of_cells=1000):
    mus = np.zeros(latent_shape)
    g = nx.from_pandas_edgelist(pa, "Source", "Target", create_using=nx.DiGraph())

    if not de_up is None:
        path = list(np.unique(list(path) + list(de_up) + list(de_down)))
        mus[:de_size] = np.random.normal(0.25,0.1,de_size)
        mus[de_size : 2*de_size] = np.random.normal(-0.25,0.1,de_size)
        nodes = list(np.unique(list(g.nodes) + list(de_up) + list(de_down)))
    else:
        nodes = list(g.nodes)

    
    sigma =  np.random.uniform(low=0, high=0.7,size=(latent_shape, latent_shape) )
    sigma[:corr_size,:corr_size] = np.random.uniform(low=0.5, high=1,size=(corr_size, corr_size) )
    np.fill_diagonal(sigma,1)
    samples = np.random.multivariate_normal(mus, sigma, num_of_cells)
    indexs = [nodes.index(x) for x in path]
    A = np.random.normal(0,1, (len(nodes), latent_shape))
    A[indexs,corr_size:] = 0
    A[indexs,:corr_size] = np.random.uniform(0,0.5,(len(indexs), corr_size))

    mask = np.arange(A.shape[0])
    mask = list(filter(lambda x: not x in indexs, mask))
    A[mask,:corr_size] = 0
    if not de_up is None:
        indexs_up = [nodes.index(x) for x in de_up]
        indexs_down = [nodes.index(x) for x in de_down]
        A[indexs_up,de_size:] = 0
        A[indexs_down,2*de_size:] = 0
        A[indexs_down,:de_size] = 0

        mus *= -1
        samples2 = np.random.multivariate_normal(mus, sigma, num_of_cells)
        samples = np.concatenate([samples,samples2])

    exp = pd.DataFrame(A@samples.T, index=nodes)
    exp = ((exp.T - exp.mean(axis=1)) / exp.std(axis=1)).T
    exp = 0.7*exp + 0.3* np.random.normal(0,1,exp.shape)
    exp *= np.random.binomial(1,0.8,exp.shape)
    if not de_up is None:
        return exp[list(range(num_of_cells))], exp[list(range(num_of_cells, 2*num_of_cells))]
    return exp

def run_anaylsis(exp, pa, tfs, receptors, path):
    gpf = tfg.build_flowing_network_with_normlized_wights(pa, tfs, exp,wights_flag=True)
    flow_dict={recp:0 for recp in receptors}
    gd = tfg.FlowGraph(flow_dict, gpf,False, sim_flag=True)
    nodes = list(gd.multy_flow_dict.keys())
    node_flow = np.array([sum(dict(vec).values()) for vec in gd.multy_flow_dict.values()])
    cent = nx.degree_centrality(gpf)
    cent_df = pd.DataFrame({"node":cent.keys(), "cent":cent.values()})
    flow_df = pd.DataFrame({"node":nodes, "flow":node_flow})
    #flow_df["flow"] /= cent_df.cent
    flow_df = flow_df.loc[~flow_df.node.isin(["source_node","sink"])]
    flow_df["flow"] += np.random.normal(0.01,0.005)
    return flow_df

def calculate_roc(flow_df, path):
    labels = flow_df.node.isin(path)
    fpr, tpr, _= roc_curve (labels, flow_df.flow, pos_label=1)
    auc = roc_auc_score (labels, flow_df.flow)
    return fpr, tpr, auc

def pipeline_downstream(args):
    pa, receptors, tfs, lr = args
    path, _, _, _ = find_path(pa, receptors, tfs, lr)
    exp = build_syntatic_data(pa, path)
    flow_df = run_anaylsis(exp,  pd.read_csv("files/humanProtinAction.csv"), tfs, receptors, path)
    roc = calculate_roc(flow_df, path)
    return roc

def find_de_lr(exp_r, exp_l, lr):
    genes = pd.Series(exp_r.index)
    pvalues = genes.apply(lambda x: ranksums(exp_r.loc[x], exp_l.loc[x])[1])
    pvalues = pd.Series(fdrcorrection(pvalues)[1])
    pvalues.index = exp_r.index
    pvalues = pvalues[pvalues <= 0.05].index
    fc = exp_r.loc[pvalues].mean(axis=1) > exp_l.loc[pvalues].mean(axis=1)
    up = fc[fc].index
    down = fc[~fc].index
    pred_lr = lr.loc[(lr.to.isin(up)) & (lr["from"].isin(down))]
    pred_recp = pred_lr["to"].drop_duplicates()
    pred_ligand = pred_lr["from"].drop_duplicates()
    return pred_recp, pred_ligand



def pipeline_LR(args):
    pa, receptors, tfs, lr = args
    path,selcted_recep, ligands, selected_tf = find_path(pa, receptors, tfs, lr)
    exp_r, exp_l = build_syntatic_data(pa, path,selcted_recep, ligands)
    pred_recp, _ =  find_de_lr(exp_r, exp_l, lr)
    tfs = list(pd.Series(tfs).sample(10).values) + list(selected_tf)
    flow_df = run_anaylsis(exp_r,  pd.read_csv("files/humanProtinAction.csv"), tfs, pred_recp, path)
    roc = calculate_roc(flow_df, path)
    return roc

def main(num_of_repat=10):
    pa, receptors, tfs, lr = load_network()
    with concurrent.futures.ProcessPoolExecutor(max_workers=20) as executor:
        result = list(executor.map(pipeline_downstream,[(pa, receptors, tfs, lr) for _ in range(num_of_repat)]))
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