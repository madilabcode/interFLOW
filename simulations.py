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
from scipy.stats import norm, nbinom
from statsmodels.stats.multitest import fdrcorrection
import concurrent.futures
from scipy.stats import spearmanr
import scanpy as sc
import pickle
import networkx as nx
import concurrent.futures
import matplotlib.pyplot as plt
import seaborn as sns
import time 

PA = pd.read_csv("files/humanProtinAction.csv")
PI = pd.read_csv("files/humanProtinInfo.csv")
TF_DICT = pyd.main_py_to_R(PI, PA)[0]
TF_DICT = {key : value for key, value in TF_DICT.items() if len(value) <= 100}


def load_network():
    pa = PA.copy()
    lr = pd.read_table(r"./files/LigandReceptorTableMouse.tsv",sep = "\t")
    receptors = lr.to.drop_duplicates().values
    tfs = list(TF_DICT.keys())
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


def load_tfs_target(tfs):
    targets = []
    for tf in tfs:
        targets += TF_DICT[tf]
    
    targets = np.unique(targets)
    mask = np.random.binomial(1,0.5,len(targets))
    targets = np.extract(mask, targets)
    return targets
    

def find_path(pa, receptors, tfs, lr, num_of_tfs=5, num_of_recp=15):
    g = nx.from_pandas_edgelist(pa, "Source", "Target", create_using=nx.DiGraph())    
    receptors = list(filter(lambda x: x in pa.Source.values or x in pa.Target.values, receptors))
    tfs = list(filter(lambda x: (x in pa.Source.values or x in pa.Target.values) and not x in receptors, tfs))
    flag = False

    while not flag:
        all_paths = {}
        for _ in range(num_of_tfs):
            selected_recep = []
            while len(selected_recep) == 0:
                selected_tf  = pd.Series(tfs).sample(1).values[0]
                paths = nx.single_source_shortest_path(g, selected_tf)
                paths = {key:value for key,value in paths.items() if key in receptors}
                selected_recep = paths.keys()

            all_paths[selected_tf] = paths

        selected_recep = set()
        for value in all_paths.values():
            if len(selected_recep) == 0:
                selected_recep = set(value.keys())
            else:
                selected_recep = selected_recep.intersection(set(value.keys()))
        
        if len(selected_recep) >= num_of_recp:
            flag = True
        
    selected_recep =  pd.Series(list(selected_recep)).sample(num_of_recp).values
    selected_tf = list(all_paths.keys())
    
    all_nodes = []
    for value1 in all_paths.values():
        for value2 in value1.values():
            all_nodes += list(value2)

    all_nodes = np.unique(all_nodes)
        
    print(selected_tf)
    targets = load_tfs_target(selected_tf)
    ligands = lr.loc[(lr.to.isin(selected_recep)) & (~lr["from"].isin(selected_recep)),"from"].drop_duplicates()
    #ligands =  pd.Series(list(filter(lambda x: x in list(pa.Source.drop_duplicates()) + list(pa.Target.drop_duplicates()), ligands)))
    return all_nodes, selected_recep, ligands , selected_tf, targets

def build_corraltion(nodes, path ,alpha = 0.7, num_of_cells=1000):
    exp = np.random.normal(0, 1, size = (len(nodes),num_of_cells))
    exp =  pd.DataFrame(exp, index=nodes)
    corr_exp =  np.random.normal(0, 1, size = num_of_cells)
    exp.loc[path] = alpha * corr_exp + (1 - alpha) *exp.loc[path] 
    #exp *= np.random.binomial(1,0.75,exp.shape)
    print(np.corrcoef(exp.loc[path]).mean())
    return exp

def build_syntatic_data(nodes, path, alpha = 0.7,  means = None, de_up=None , de_down=None , num_of_cells=1000):
    n, p = 5, 0.5     
    indexs = list(map(lambda x: nodes.index(x), path))
    exp = build_corraltion(nodes, path ,alpha, num_of_cells)
    df = pd.DataFrame(norm.cdf(exp.T,exp.mean(axis=1), exp.std(axis=1))).T
    if means is None:
         means = np.random.random_integers(5,100,len(nodes))
    else:
        means = pd.Series(means, index= nodes)
        means.loc[de_up] += np.random.normal(15, 5, len(de_up))
        means.loc[de_down] -= np.random.normal(15, 5, len(de_down))
        means.loc[~means.index.isin(list(de_up) + list(de_down))] += np.random.normal(0, 2, len(nodes) - len(de_up) - len(de_down))
        means = means.apply(lambda x: max(1,round(x))).values

    df["means"] = means
    df  = df.apply(lambda x: pd.Series(nbinom.ppf(x.drop("means"),x["means"],p)), axis=1 )
    df *= np.random.binomial(1,0.85,df.shape)
    df.index = nodes
    obj = sc.AnnData(df.T)
    sc.pp.normalize_total(obj, target_sum=1e4)
    sc.pp.log1p(obj)
    df = pd.DataFrame(obj.X.T, index=nodes)
    if len(path) > 0:
        correlation_matrix, p_values = spearmanr(df.loc[path], axis=1)
        return df, correlation_matrix.mean() , means
    return df, means

def run_anaylsis(exp, pa, tfs, receptors, path, path_neg=None):
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
   # flow_df["flow"] += np.random.normal(0.1,0.05)
    return flow_df

def run_anaylsis_sc(exp, pa, tfs, receptors, path, path_neg=None):
    gpf = tfg.build_flowing_network_with_normlized_wights(pa, tfs, exp,wights_flag=True)
    flow_values = []
    flow_dicts = {}

    for rec in receptors:
        try:
            flow_value, flow_dict = nx.maximum_flow(gpf, rec ,"sink", flow_func=nx.algorithms.flow.dinitz)
            flow_values.append([rec, flow_value])
            flow_dicts[rec] = flow_dict
        except:
            print("node " + rec + " is not in the graph")

    df = pd.DataFrame(list(map(lambda val: val[1], flow_values)), columns=["flow"])
    df["node"] = receptors.values
    return df

def calculate_roc(flow_df, path, path_neg=None):
    #flow_df = flow_df.loc[flow_df.node.isin(list(path)+list(path_neg))]
    labels = flow_df.node.isin(path)
    fpr, tpr, _= roc_curve (labels, flow_df.flow, pos_label=1)
    auc = roc_auc_score (labels, flow_df.flow)
    return auc

def calculate_roc_tfs(flow_df, selected_tf,  tfs):
    flow_df = flow_df.loc[flow_df.node.isin(tfs)]
    labels = flow_df.node.isin(selected_tf)
    fpr, tpr, _= roc_curve (labels, flow_df.flow, pos_label=1)
    auc = roc_auc_score (labels, flow_df.flow)
    return auc

def calculate_roc_recep(flow_df, selected_recep, receptors ):
    flow_df = flow_df.loc[flow_df.node.isin(receptors)]
    labels = flow_df.node.isin(selected_recep)
    fpr, tpr, _= roc_curve (labels, flow_df.flow, pos_label=1)
    auc = roc_auc_score (labels, flow_df.flow)
    return auc

def pipeline_downstream(args):
    np.random.seed((os.getpid() * int(time.time())) % 123456789)
    pa, receptors, tfs, lr, corr, _ = args
    path, _, _, selected_tf, _= find_path(pa, receptors, tfs, lr)
    g = nx.from_pandas_edgelist(pa, "Source", "Target", create_using=nx.DiGraph())
    nodes = list(g.nodes)
    #path_neg, _, _, _, _= find_path(pa.loc[(~pa.Source.isin(path)) & (~pa.Target.isin(path))], receptors, tfs, lr)
    path_neg, _, _, _, _= find_path(pa, receptors, tfs, lr)
    exp, corr_mean, _ = build_syntatic_data(nodes, path, corr)
    flow_df = run_anaylsis(exp, PA.copy(), tfs, receptors, path)
    #roc = calculate_roc(flow_df, path, path_neg)
    roc = calculate_roc_tfs(flow_df, selected_tf,  tfs)
    return roc, corr_mean


def cellchat_pred(exp):
    with localconverter(default_converter + rpyp.converter):
        r("library(CellChat)")
        r("library(Seurat)")
        r("library(dplyr)")

        cell_chat = r("""   
                        cell_chat_pred = function(exp){
                            n = dim(exp)[2]
                            colnames(exp) = 1:n
                            row.names(exp) = row.names(exp) %>% toupper()
                            meta = data.frame("Cells"= colnames(exp), "labels"=c(rep("rec",n/2), rep("ligand",n/2)))
                            exp = sample(exp, 200)
                            meta = meta[meta$Cells %in% colnames(exp),]
                            exp = exp[,meta$Cells]
                            obj <- createCellChat(object = as.matrix(exp), meta = meta, group.by = "labels")
                            CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
                            showDatabaseCategory(CellChatDB)
                            CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
                            obj@DB <- CellChatDB.use
                            obj <- subsetData(obj) # This step is necessary even if using the whole database
                            future::plan("multiprocess", workers = 4)
                            obj <- identifyOverExpressedGenes(obj)
                            obj <- identifyOverExpressedInteractions(obj)
                            obj <- computeCommunProb(obj)
                            #obj <- filterCommunication(obj, min.cells = 3)
                            df.net <- subsetCommunication(obj)
                            return(df.net)
                        }
                      """)
        exp.colums = list(range(0, exp.shape[1]))
        df_cc = cell_chat(exp)
        df_cc["receptor"] = df_cc["receptor"].apply(lambda x: x.split("_")[0])
        df_cc["receptor"] = df_cc["receptor"].apply(lambda x:  x[0] + x[1:].lower())

        scores = df_cc.groupby("receptor")["prob"].max()
        return scores
        



def find_de_lr(exp_r, exp_l, lr, selected_recep, ligands):
    genes = pd.Series(exp_r.index)
    pvalues = genes.apply(lambda x: ranksums(exp_r.loc[x], exp_l.loc[x])[1])
    pvalues = pd.Series(fdrcorrection(pvalues)[1])
    pvalues.index = exp_r.index
    pvalues = pvalues[pvalues <= 0.05].index
    #fc = exp_r.loc[pvalues].mean(axis=1) > exp_l.loc[pvalues].mean(axis=1)
    #up = fc[fc].index
    #down = fc[~fc].index
    up = exp_r.loc[pvalues].mean(axis=1) / exp_l.loc[pvalues].mean(axis=1) > 1.1
    up = up[up].index
    down = exp_l.loc[pvalues].mean(axis=1) / exp_r.loc[pvalues].mean(axis=1) > 1.1
    down = down[down].index
    pred_lr = lr.loc[(lr.to.isin(up)) & (lr["from"].isin(down))]
    pred_recp = pred_lr["to"].drop_duplicates()
    pred_ligand = pred_lr["from"].drop_duplicates()
    return pred_recp, pred_ligand

def roc_cellcaht(cc_scores,selected_recep, receptors ):
    recp_df = pd.DataFrame({"def":np.zeros(len(receptors)) ,"node":receptors})
    cc_scores = pd.DataFrame(cc_scores)
    cc_scores.columns = ["flow"]
    cc_scores["node"] = cc_scores.index
    cc_scores.reset_index(inplace=True,drop=True)
    cc_scores = cc_scores.merge(recp_df, on="node", how="right")
    cc_scores["flow"].fillna(0, inplace=True)
    cc_scores = cc_scores[["node","flow"]]
    return calculate_roc_recep(cc_scores, selected_recep, receptors)

def pipeline_LR(args):
    np.random.seed((os.getpid() * int(time.time())) % 123456789)
    pa, receptors, tfs, lr, corr, cc_flag = args
    path,selected_recep, ligands, selected_tf, targets = find_path(pa, receptors, tfs, lr)
    g = nx.from_pandas_edgelist(pa, "Source", "Target", create_using=nx.DiGraph())
    nodes = list(np.unique(list(g.nodes) + list(ligands)))
    exp_r, corr_mean , means  = build_syntatic_data(nodes, path, corr, means=None)
    exp_l, _ = build_syntatic_data(nodes, [], corr, means=means, de_up=ligands, de_down=selected_recep)#, targets=targets)
    exp_l.columns = list(range(exp_r.shape[1] , 2*exp_r.shape[1]))
    exp = pd.concat([exp_r, exp_l], axis=1)#.to_csv(r"./validation/sim_exp.csv")
    #exp_r.to_csv(r"./validation/scRNAseq_r.csv")
    #exp_l.to_csv(r"./validation/scRNAseq_l.csv")
    pd.Series(tfs).to_csv(r"./validation/tfs.csv")
    pred_recp, _ =  find_de_lr(exp_r, exp_l, lr, selected_recep, ligands)
   # tfs = list(pd.Series(tfs).sample(10).values) + list(selected_tf)
    flow_df = run_anaylsis_sc(exp_r,  PA.copy(), tfs, pred_recp, path)
    if cc_flag:
        cc_scores = cellchat_pred(exp)
        #recp_df = pd.DataFrame({"def":np.zeros(len(receptors)) ,"node":receptors})
        #flow_df = flow_df.merge(recp_df, on="node", how="right")
        #flow_df["flow"] = flow_df["flow"].fillna(0)
        #flow_df = flow_df[["node","flow"]]
        roc = calculate_roc_recep(flow_df, selected_recep, pred_recp)
        roc_cc = roc_cellcaht(cc_scores,selected_recep, pred_recp) 
        return roc, roc_cc, corr_mean
    
    else:
        roc = calculate_roc_recep(flow_df, selected_recep, np.unique(list(selected_recep) + list(pred_recp)))
    return roc, corr_mean


def make_box_plots(results, corr_range, num_of_repat):
    results = np.array(results)
    cors = results[:,:,-1]
    cors =  cors.mean(axis=1)
    resultsROC = np.array(results)[:,:,0].squeeze()
    df = pd.DataFrame({"ROCAUC" : resultsROC.reshape(-1), "Pathway Corr": [np.round(corr,2) for corr in cors for _ in range(num_of_repat) ]})

    if results.shape[-1] == 3:
        results_cc = np.array(results)[:,:,1].squeeze()
        df_cc = pd.DataFrame({"ROCAUC" : results_cc.reshape(-1), "Pathway Corr": [np.round(corr,2) for corr in cors for _ in range(num_of_repat) ]})
        df = pd.concat([df, df_cc])
        df["method"] = np.concatenate([np.repeat("interFLOW",num_of_repat* len(cors)), np.repeat("CellChat",num_of_repat* len(cors))])
    fig, ax = plt.subplots(figsize=[13,7])
    sns.boxenplot(ax=ax, data=df,x="Pathway Corr", y="ROCAUC", hue="method")    
    sns.set_theme(style='white',font_scale=2)
    plt.show()


def main(num_of_repat=10, corr_range=[0.2,0.25,0.3,0.5,0.7,0.9], lr_flag=True, box_plot_flag = True, cc_flag= True):
    pa, receptors, tfs, lr = load_network()
    if lr_flag:
        func = pipeline_LR
    else:
        func = pipeline_downstream

    results = []
    for corr in corr_range  :
        with concurrent.futures.ProcessPoolExecutor(max_workers=20) as executor:
            result = list(executor.map(func,[(pa, receptors, tfs, lr, corr, cc_flag) for _ in range(num_of_repat)]))
        results.append(result)
    
    tfg.save_obj(results, r"./files/test_sim")
    
    if box_plot_flag:
        make_box_plots(results, corr_range, num_of_repat)

    else:
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

def test():
    results = tfg.load_obj(r"./files/test_sim")
    make_box_plots(results, None, num_of_repat = 10)

if __name__ == "__main__":
    main(lr_flag=True, corr_range=[0.2, 0.9])   
    #main(lr_flag=False, corr_range=[0.0, 0.9])
    #test()