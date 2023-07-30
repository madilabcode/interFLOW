import re
import numpy as np
import pandas as pd



def format_Gene_dict(gene, from_dict = False):
  if from_dict:
    gene = gene[1:len(gene) - 1]
  
  fl = gene[0]
  sl = gene[1:]
  sl = str.lower(sl)

  return fl + sl


def make_dict(data):
  tf_dict = {}
  tf_list = data.split(": [")
  first_tf = tf_list[0][1:]
  tf_list = tf_list[1:]
  new_tf_list = list()
  new_tf_list.append(format_Gene_dict(first_tf,True))

  for idex, obj in enumerate(tf_list):
      obj = obj.split("], ")
      try:
          obj[1] = format_Gene_dict(obj[1],True)
          obj[0] = [format_Gene_dict(ge,True) for ge in obj[0].split(", ")]
          new_tf_list.append(obj[0])
          new_tf_list.append(obj[1])
      except:
          print("error in " + str(idex))
          new_tf_list.append(format_Gene_dict(obj[0],True))
  
  for i in range(0, len(new_tf_list) - 1, 2):
      tf_dict[new_tf_list[i]] = new_tf_list[i + 1]
  return tf_dict

def format_graph_files(pro_name,col_names, pro_action,col_actions):
  for col_action in col_actions:
        pro_action[col_action] =  pro_action[col_action].apply(format_Gene_dict)
  for col_name in col_names:
        pro_name[col_name] =  pro_name[col_name].apply(format_Gene_dict)

def make_dist_dickt(data, tf_dict):
  tfs = tf_dict.keys()
  dist_dickt = {}

  for tf in tfs:
    try:
      meanx = np.mean(data.loc[tf,])
      sdx = np.std(data.loc[tf,])
      dist_dickt[tf] = [meanx,sdx]
    except:
      print("not valde tf" + tf)

  return dist_dickt


def main_py_to_R(pro_name, pro_action):
  data = np.load('files/TFdictBT1.npy', allow_pickle=True)
  data = str(data)
  tf_dict = make_dict(data)
  pro_action = pro_action.loc [:,["Input-node Gene Symbol","Output-node Gene Symbol","Edge direction score"]]
  format_graph_files(pro_name,["Gene Symbol"],pro_action,["Input-node Gene Symbol","Output-node Gene Symbol"])
  #dist_dickt = make_dist_dickt(r.Expression_scale,tf_dict )
  return tf_dict, pro_name, pro_action

