from cmath import isnan
from networkx.testing.test import run
import pandas as pd
import numpy as np
import math
from rpy2.robjects import NULL, r
# import matplotlib.pyplot as plt
import scipy.stats as sc
import networkx as nx
import pickle
from Codes import CERNO as ce
from Codes import utils as utils
from sklearn.metrics import mutual_info_score as mis
from Codes import pyDictToR as pyd
import concurrent.futures
import warnings
from scipy.stats import ranksums
import statsmodels.stats.multitest
import dill
warnings.filterwarnings('ignore') 

CUTOFF = 0.0
PA = pd.read_csv("./files/humanProtinAction.csv")
class GraphsDicts:
    def __init__(self):
        self.dic = {}

    def add_obj(self, key, obj):
        self.dic[key] = obj

    def return_GD_from_key(self, key):
        return self.dic[key]

    def save_obj(self, name):
        save_obj(self, name)


class FlowGraph:
    def __init__(self, flow_dics, capcity_network,do_permutation = True):
        self.flow_dics = flow_dics
        self.capcity_network = capcity_network
        self.max_flow_dict = self.calculate_max_flow_for_all_recp()
        self.max_multy_flow = None
        self.pa = self.run_multiy_source_flow()
        if do_permutation:
            self.pvalues = self.calcualte_pvalue_for_every_node(num_of_perm=10)

    @classmethod
    def make_capcity_graph(cls, network_dict):
        ProtInfo = pd.read_csv(r"files/humanProtinInfo.csv")
        pa = pd.read_csv(r"files/humanProtinAction.csv")
        pa["Target"] = pa["Output-node Gene Symbol"]
        pa["Source"] = pa["Input-node Gene Symbol"]
        pa["capacity"] = math.inf
        keeps = ["Source", "Target", "capacity"]
        pa = pa[keeps]
        genes = network_dict.keys()
        genes = list(map(lambda gene: gene.split("_")[0], genes))

        pa["Source"] = pa["Source"].apply(format_Gene_dict).astype(str)
        pa["Target"] = pa["Target"].apply(format_Gene_dict).astype(str)

        pa = pa.loc[pa["Source"].isin(genes), :]
        pa = pa.loc[pa["Target"].isin(genes), :]
        pa["Source"] += "_Source"
        pa["Target"] += "_Target"

        gp1 = nx.from_pandas_edgelist(pa, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
        gp2 = nx.from_dict_of_dicts(network_dict, create_using=nx.DiGraph())

        return nx.compose(gp1, gp2)

    @classmethod
    def make_obj_out_of_capcity_network(cls, network_dict, recptors):
        gp =FlowGraph.make_capcity_graph(network_dict)
        flow_dicts = {}

        recptors = np.unique(recptors)
        for rec in recptors:
            try:
                flow_value, flow_dict = nx.maximum_flow(gp, "s", rec + "_Source", flow_func=nx.algorithms.flow.dinitz)
                flow_dicts[rec] = flow_dict
            except:
                print("node " + rec + " is not in the graph")

        return FlowGraph(flow_dicts, gp)

    def update_obj(self):
        return FlowGraph(self.flow_dics, self.capcity_network,do_permutation=False)

    def calculate_max_flow(self, recp):
        dfFlow = self.make_df_flow(recp)
        max_flow = 0

        for row in dfFlow.iterrows():
            row = row[1]
            if row["Target"] == "sink" :
                max_flow += row["flow"]

        return max_flow

    def calculate_max_flow_for_all_recp(self):
        max_flow_dict = {}
        for recp in self.flow_dics.keys():
            max_flow_dict[recp] = self.calculate_max_flow(recp)

        return max_flow_dict
    
    def single_recptor_flow(self,recp):
        graph = self.make_df_flow(recp=recp)
        graph = graph.loc[graph.Target !="sink"]
        return graph 
    
    def make_df_flow(self, recp=None, flow=None):
        if flow is None:
            flow = self.flow_dics[recp]
        source = []
        target = []
        flows = []

        for Inode in flow.keys():
            edges = flow[Inode]
            for Onode in edges.keys():
                if edges[Onode] != 0 and Inode != Onode:
                    source.append(Inode)
                    target.append(Onode)
                    flows.append(edges[Onode])

        df = pd.DataFrame([source, target, flows])
        df = df.transpose()
        df.columns = ["Source", "Target", "flow"]
        return df

    def tfs_per_recp_in_max_flow(self):
        tf_rec = {}
        for rec in self.flow_dics.keys():
            df = self.make_df_flow(rec)
            df = df.loc[df["Target"] == "sink", :]
            tf_rec[rec] = len(df["Source"])

        return tf_rec

    @staticmethod
    def perm_for_tf(cap):
        cap_df = FlowGraph.premutate_graph(cap)
        graph =  nx.from_pandas_edgelist(cap_df, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
        flow_value,_ =  nx.maximum_flow(graph, "source_node" ,"sink", flow_func=nx.algorithms.flow.dinitz)
        return flow_value
   
    @staticmethod
    def flow_to_single_tf(args,num_of_perm=10):
        pa,tf,remove_tfs = args
        remove_tfs.remove(tf)
        pa = pa.loc[~(pa.Target.isin(remove_tfs)) & ~(pa.Source.isin(remove_tfs))]
        print(tf)
        gpf = nx.from_pandas_edgelist(pa, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
        orig_flow,_ = nx.maximum_flow(gpf, "source_node", tf, flow_func=nx.algorithms.flow.dinitz)
        flows =  np.array([FlowGraph.perm_for_tf(pa.copy()) for i in range(num_of_perm)])
        p_value = (flows >= orig_flow).sum() / len(flows)  
        return tf, p_value
    
    def flow_dict_to_single_tf(self,tf,pa=None):
        if pa is None:
            pa = self.build_multy_source_network()
            pa.columns = ["Source","Target","capacity"]
        remove_tfs = list(np.unique(list(pa.loc[pa["Target"] == "sink", "Source"])))        
        remove_tfs.remove(tf)
        pa = pa.loc[~(pa.Target.isin(remove_tfs)) & ~(pa.Source.isin(remove_tfs))]
        gpf = nx.from_pandas_edgelist(pa, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
        _,flow_dict = nx.maximum_flow(gpf, "source_node", tf, flow_func=nx.algorithms.flow.dinitz)
        df = self.make_df_flow(flow=flow_dict)
        df = df.loc[(df.Source!="source_node") & (df.Target !="sink")]
        return  df
    
    def flow_to_tfs(self,pa=None):
        if pa is None:
            pa = self.run_multiy_source_flow()
        all_tfs = np.unique(list(pa.loc[pa["Target"] == "sink", "Source"]))
        tf_flow = {}
        with concurrent.futures.ProcessPoolExecutor(max_workers=3) as executor:
            results = list(executor.map(self.flow_to_single_tf,
                                        [(pa.copy(), tf, list(all_tfs).copy()) for tf in all_tfs]))
        for element in results:
            tf_flow[element[0]] = element[1]
           
        tf_effects = pd.DataFrame({"flow":tf_flow.values()}, index =tf_flow.keys())
        tf_effects = tf_effects.sort_values(by="flow", axis=0, ascending=False)
        return tf_effects
    
    def flow_to_all_tf(self,pa=None):
        if pa is None:
            pa = self.build_multy_source_network()
        pa.columns = ["Source","Target","capacity"]
        gpf = nx.from_pandas_edgelist(pa, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
        cent = nx.degree_centrality(gpf)
        max_flow, flow_network = nx.maximum_flow(gpf, "source_node","sink", flow_func=nx.algorithms.flow.dinitz)
        flow_network = self.make_df_flow(flow=flow_network)
        
        all_tfs = flow_network.loc[flow_network["Target"] == "sink", ["Source","flow"]]
        all_tfs = all_tfs.groupby("Source")["flow"].sum()
        
        return  pa, pd.DataFrame(all_tfs.sort_values(axis=0, ascending=False))

    
    def calculate_max_flow_for_one_tf(self, tf):
        pat = self.pa.copy()
        all_tf = np.unique(pat.loc[pat["Target"] == "sink", "Source"])
        all_tf = np.delete(all_tf, np.where(all_tf == tf))
        pat["isNotTf"] = pat.apply(
                lambda row: False if row["Source"] in all_tf  or  row["Target"] in all_tf else True, axis=1)
        pat = pat.loc[pat["isNotTf"], :]
        gpf = nx.from_pandas_edgelist(pat, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
        flow_value, flow_dict  = nx.maximum_flow(gpf, "source_node", tf, flow_func=nx.algorithms.flow.dinitz)

        df = self.make_df_flow(flow=flow_dict)
        return df.loc[df.flow >= 0.1]
    
    def calculate_significant_tf_Multy_sinc(self):
        pa = self.run_multiy_source_flow()
        all_genes = np.unique(list(pa.loc[pa["Target"] == "sink", "Source"]))
        delta_dict = {}
        with concurrent.futures.ProcessPoolExecutor(max_workers=3) as executor:
            results = list(executor.map(self._calculate_flow_delta_for_df,
                                        [(gene, pa.copy(), self.max_multy_flow) for gene in all_genes]))
        for element in results:
            delta_dict[element[0]] = element[1]
           
        tf_effects = pd.DataFrame({"flow":delta_dict.values()}, index =delta_dict.keys())
        tf_effects = tf_effects.sort_values(by="flow", axis=0, ascending=False)
        return tf_effects


    @staticmethod
    def _calculate_flow_delta_for_df(args):
        gene, pat, max_flow = args
        if gene != "source_node" and gene != "sink":
            pat["isNotTf"] = pat.apply(
                lambda row: True if row["Source"] != gene  and row["Target"] != gene else False, axis=1)
            pat = pat.loc[pat["isNotTf"], :]
            gpf = nx.from_pandas_edgelist(pat, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
            flow_value = nx.maximum_flow(gpf, "source_node", "sink", flow_func=nx.algorithms.flow.dinitz)[0]
            return gene, max_flow - flow_value

    def build_multy_source_network(self):
        pa = nx.to_pandas_edgelist(self.capcity_network)
        pa.columns = ["Source", "Target", "capacity"]
    
        for recp in self.flow_dics.keys():
            pa.loc[pa.shape[0]] = ["source_node",recp,math.inf]
        
        return pa 

    def calcualte_pvalue_for_every_node(self, num_of_perm=100):
        self.pa = self.build_multy_source_network()
        self.gpf = nx.from_pandas_edgelist(self.pa, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
        _, flow_dic = nx.maximum_flow(self.gpf, "source_node","sink", flow_func=nx.algorithms.flow.dinitz)
        orig_flow = np.array([sum(dict(vec).values()) for vec in flow_dic.values()])
        flows = []
        for _ in range(num_of_perm):
            perm = self.premutate_graph(self.pa)
            gpf_perm = nx.from_pandas_edgelist(perm, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
            perm_flow_dict = nx.maximum_flow(gpf_perm, "source_node","sink", flow_func=nx.algorithms.flow.dinitz)[1]
            flow_values = [sum(dict(vec).values()) for vec in perm_flow_dict.values()]
            flows.append(flow_values)
        
        flows = np.array(flows).T
        greater_than = np.greater_equal(flows, orig_flow[:,np.newaxis])
        pvalues = np.sum(greater_than, axis=1).astype(np.float64)
        pvalues /= num_of_perm
        return pvalues 
    
    
    def significant_nodes(self):
        nodes = list(self.gpf.nodes)
        sig_pvalues = np.where(self.pvalues < 0.05)[0]
        sig_nodes = [nodes[index] for index in sig_pvalues]
        return sig_nodes

    def run_multiy_source_flow(self):
        pa = self.build_multy_source_network()
        gpf = nx.from_pandas_edgelist(pa, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
        max_flow = nx.maximum_flow(gpf, "source_node","sink", flow_func=nx.algorithms.flow.dinitz)[0]
        self.max_multy_flow = max_flow
        return pa 

    def calculate_significant_tf_Multy_sinc(self):
        pa = self.run_multiy_source_flow()
        all_genes = np.unique(list(pa.loc[pa["Target"] == "sink", "Source"]))
        delta_dict = {}
        with concurrent.futures.ProcessPoolExecutor(max_workers=20) as executor:
            results = list(executor.map(self._calculate_flow_delta_for_df,
                                        [(gene, pa.copy(), self.max_multy_flow) for gene in all_genes]))
        for element in results:
            delta_dict[element[0]] = element[1]
           
        tf_effects = pd.DataFrame({"flow":delta_dict.values()}, index =delta_dict.keys())
        tf_effects = tf_effects.sort_values(by="flow", axis=0, ascending=False)
        return tf_effects


    @staticmethod
    def premutate_graph(sub_graph,max_ineration=100):
        inf_sub_graph = sub_graph.loc[sub_graph.capacity == np.inf,:]
        sub_graph = sub_graph.loc[sub_graph.capacity < np.inf,:]

          
        sub_graph["Target"] = np.random.permutation(sub_graph["Target"])
        dup = sub_graph.drop_duplicates().index
        dup = sub_graph.loc[~sub_graph.index.isin(dup)]
       # dup2 = sub_graph.drop_duplicates().merge(sub_graph.drop_duplicates(), left_on=["Source","Target"], right_on=["Target","Source"], how="left")
       # dup2.index = sub_graph.drop_duplicates().index
       # dup2.dropna(inplace=True)
       # dup = list(dup.index) + list(dup2.index)
        inter = 0
        while len(dup)> 100 and inter < max_ineration:
            sub_graph.loc[dup,"Target"] = np.random.permutation(sub_graph.loc[dup,"Target"])
            dup = sub_graph.drop_duplicates().index
            dup = sub_graph.loc[~sub_graph.index.isin(dup)]
           # dup2 = sub_graph.drop_duplicates().merge(sub_graph.drop_duplicates(), left_on=["Source","Target"], right_on=["Target","Source"], how="left")
           # dup2.index = sub_graph.drop_duplicates().index
           # dup2.dropna(inplace=True)
           # dup = list(dup.index) + list(dup2.index)
            inter+=1
        
        return pd.concat([sub_graph,inf_sub_graph])
        
    def edge_dgree_perm(self, source_node,cap=None):
        if cap is None:
            cap = self.capcity_network
        
        cap_df = nx.to_pandas_edgelist(cap)
        cap_df.columns = ["Source", "Target", "capacity"]
        ##### 3 bins
        cap_df.sort_values("capacity",inplace=True)
        top1 = cap_df.iloc[int(cap_df.shape[0]/3)]["capacity"]
        top2 = cap_df.iloc[int(2*cap_df.shape[0]/3)]["capacity"]
        sub_cup1 = self.premutate_graph(cap_df.loc[cap_df.capacity < top1,:])
        sub_cup2 = self.premutate_graph(cap_df.loc[(cap_df.capacity >= top1) & (cap_df.capacity <= top2) ,:])
        sub_cup3 = self.premutate_graph(cap_df.loc[cap_df.capacity > top2,:])
        cap_df = pd.concat([sub_cup1,sub_cup2,sub_cup3])
        ###on bine!
        cap_df = self.premutate_graph(cap_df)
        graph =  nx.from_pandas_edgelist(cap_df, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
        flow_value, flow =  nx.maximum_flow(graph, source_node ,"sink", flow_func=nx.algorithms.flow.dinitz)
        return flow_value


    def calculate_p_value_for_tfs(self,num_of_perm=100):
        pa, flows = self.flow_to_tfs() 
        pa.columns = ["Source","Target","capacity"]
    
        for i in range(num_of_perm):
            perm_pa = self.premutate_graph(pa)
            if i == 0:
                df = self.flow_to_tfs(perm_pa)[1]
            else :
                df = df.merge(self.flow_to_all_tf(perm_pa)[1],left_index=True, right_index=True)
        
        flows = pd.DataFrame(flows)
        flows["gene"] = flows.index
        flows = flows.loc[flows.index.isin(df.index)]
        flows["p_value"] = flows.apply(lambda x: (df.loc[x[1]] > x[0]).sum() / num_of_perm ,axis=1)
        return  flows[["flow","p_value"]]
            
        
def format_Gene_dict(gene, from_dict=False):
    if from_dict:
        gene = gene[1:len(gene) - 1]

    fl = gene[0]
    sl = gene[1:]
    sl = str.lower(sl)

    return fl + sl

def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        dill.dump(obj, f, dill.HIGHEST_PROTOCOL)


def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return dill.load(f)


def mi(exp,pa, recps_for_roc=None):
    clusterex = exp.copy()
    clusterex["not_zero"] = clusterex.apply(lambda x: x.astype(bool).sum(), axis=1)
    if recps_for_roc is None:
         clusterex = clusterex.loc[clusterex.not_zero > clusterex.shape[1] * CUTOFF, :]
    else:
        clusterex = clusterex.loc[(clusterex.not_zero > clusterex.shape[1] * CUTOFF) | (clusterex.index.isin(recps_for_roc)) , :]

    clusterex.drop("not_zero",axis=1,inplace = True)

    pa = pa.loc[(pa.Source.isin(exp.index)) &  (pa.Target.isin(exp.index) )]
    pa["wights"] = pa.apply(lambda x : mis(exp.loc[x.Source],exp.loc[x.Target]),axis=1)
    pa["wights"] = pa.apply(lambda x: x.wights if x.Source in clusterex.index and x.Target in clusterex.index else 0 ,axis=1)
    pa["wights"] = (pa["wights"] - pa["wights"].min())/(pa["wights"].max() - pa["wights"].min())
    pa["capacity"] = pa["wights"]
    pa.drop("wights",axis=1,inplace = True)

    return pa

def build_flowing_network_with_normlized_wights(ProtAction, tfs_enrc, exp, recps_for_roc = None,wights_flag=True):
    pa = ProtAction
    pa["Source"] = pa["Input-node Gene Symbol"]
    pa["Target"] = pa["Output-node Gene Symbol"]
    pa["capacity"] = pa["Edge direction score"]
    keeps = ["Source", "Target", "capacity"]
    pa = pa[keeps]

    pa["Source"] = pa["Source"].apply(format_Gene_dict).astype(str)
    pa["Target"] = pa["Target"].apply(format_Gene_dict).astype(str)
    
    if wights_flag:
        pa = mi(exp,pa)
    pa.reset_index(drop=True,inplace=True)

    for tf in tfs_enrc:
        if tf in  pa["Target"].to_list():
            pa.loc[pa.shape[0]] = [tf, "sink", np.inf]
 

    gp = nx.from_pandas_edgelist(pa, "Source", "Target", edge_attr="capacity", create_using=nx.DiGraph())
    utils.save_obj(gp, "temp_grpah")

    return gp

def hypergeomtric_test(gs, markers,genes):
    N = len(genes)
    m = len(list(filter(lambda x: x in genes,gs)))
    k = len(markers)
    x = len(list(filter(lambda x: x in markers.index,gs)))
    return 1 - sc.hypergeom.cdf(x,N,m,k) 

def wilcoxon_enrcment_test(tf,gene_list,exp):
    try:
        gene_exp = exp.loc[exp.index.isin(gene_list)]
        backround_exp = exp.loc[~exp.index.isin(gene_list)]
        return ranksums(backround_exp,gene_exp,alternative="two-sided")[1]
    except:
        print(f"problem accure in {tf}")
        return 1
                                
def enrch_tfs(exp,tfs, reduce):
    expamount = exp.apply(lambda x: x.astype(bool).sum(), axis=1)
    expamount = expamount[expamount > exp.shape[1]*0.1]
    tfs = {key: value for key,value in tfs.items() if key in expamount.index}
    tfs_scores = {tf:wilcoxon_enrcment_test(tf,tfs[tf],exp.loc[expamount.index].mean(axis=1)) for tf in tfs.keys() if len(tfs[tf])>1}
    tfs_scores = {key: value for key,value in tfs_scores.items() if not np.isnan(value)}
    values = statsmodels.stats.multitest.fdrcorrection(np.array(list(tfs_scores.values())))[1]
    tfs_scores = {key:values[index] for index,key in enumerate(tfs_scores.keys())}
    if reduce:
        tfs_scores = {tf: tfs_scores[tf] for tf in tfs_scores.keys() if tfs_scores[tf]<=0.05}
    return tfs_scores


#### delete reuce and do permutation fals

def DSA_anaylsis(exp, recptors, ProtAction, ProtInfo, tfs,markers=None,recps_for_roc = None,reduce=True,do_permutation=True,wights_flag=True):
    tfs_scores = enrch_tfs(exp,tfs,reduce)
    gpf = build_flowing_network_with_normlized_wights(ProtAction, tfs_scores, exp,recps_for_roc=recps_for_roc,wights_flag=wights_flag)
    flow_values = []
    flow_dicts = {}
    recptors = np.unique(recptors)

    for rec in recptors:
        try:
            flow_value, flow_dict = nx.maximum_flow(gpf, rec ,"sink", flow_func=nx.algorithms.flow.dinitz)
            flow_values.append([rec, flow_value])
            flow_dicts[rec] = flow_dict
        except:
            print("node " + rec + " is not in the graph")

    df = pd.DataFrame(list(map(lambda val: val[1], flow_values)), columns=["DSA"])
    # df.DSA = df.DSA * scale_factor
    df.index = list(map(lambda val: val[0], flow_values))
    df["Recp"] = df.index
    gd = FlowGraph(flow_dicts, gpf,do_permutation)
    return df, gd
