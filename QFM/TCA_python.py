import copy
from math import comb
import numpy as np
import pandas as pd
import random
import os
import time
from typing import Union
from multiprocessing import Pool, Process, Queue, Lock, cpu_count

import sys
sys.setrecursionlimit(10000)

class TCA_tree: # each instance is a node in the tree
    __is_leaf: bool # if this node is a leaf
    __son: list # contains all sons, every element is a TCA_tree instance; for leaves, it's a empty list
    n_cells: int # this node and its sons repersent for how many cells in total
    label: set
    __dist: float # distance to father
    name: str
    __son_shuffled: list # list[int], each number means if the son need to shuffle and the shuffled potency, 0 means don't need to shuffle


    def __init__(self, is_leaf, son,
    n_cells, label, dist, name, son_shuffled):
        self.__is_leaf = is_leaf
        self.__son = son
        self.n_cells = n_cells
        self.label = label
        self.__dist = dist
        self.name = name
        self.__son_shuffled = son_shuffled


    def __print_newick(self) -> str:
        if self.__is_leaf:
            if self.__dist == 0:
                return self.name
            return self.name + ':' + str(self.__dist)
        son_str = []
        for i in range(len(self.__son)):
            son_str.append(self.__son[i].__print_newick())
        if self.__dist == 0:
            return '(' + ','.join(son_str) + ')' + self.name
        return '(' + ','.join(son_str) + ')' + self.name + ':' + str(self.__dist)
    

    def print_newick(self) -> str:
        return self.__print_newick() + ';'



    # pseudo-overload merge
    def merge(self, arg):
        if type(arg) == set:
            self.__merge_by_label(arg)
        elif type(arg) == np.ndarray:
            self.__merge_by_ann(arg)
        elif type(arg) == list:
            self.__merge_by_ann(np.array(arg))
        else:
            raise TypeError("illegal argument, only accept set (for labels) or array (for annotations).")


    # update labels and merge nodes with same label, tip_ann defination same as update_label
    def __merge_by_ann(self, tip_ann: np.ndarray) -> None:
        if tip_ann.shape[1] != 2:
            raise ValueError("annotation has more than TWO rows.")
        if self.__is_leaf:
            rlt = tip_ann[np.argwhere(tip_ann[...,0] == self.name),1]
            if rlt.shape[0] == 0: # node not in annotation
                raise ValueError(f"node {self.name} not in annotations.")
            elif rlt.shape[0] > 1: # more than one annotation for the node
                for i in range(1, rlt.shape[0]):
                    if rlt[i] != rlt[0]:
                        raise ValueError(f"node {self.name} has more than one annotation.")
            self.label = {str(rlt[0][0])}
            return
        self.label = set()
        for i in range(len(self.__son)):
            self.__son[i].__merge_by_ann(tip_ann)
            self.label = self.label | self.__son[i].label
        if len(self.label) == 1: #merge
            self.__is_leaf = True
            self.__son = []
            self.__son_shuffled = []
    

    # merge nodes with label in merged_labels
    def __merge_by_label(self, merged_labels: set) -> None:
        if self.__is_leaf:
            return
        if self.label.issubset(merged_labels): #merge
            self.__is_leaf = True
            self.__son = []
            self.__son_shuffled = []
            return
        for i in range(len(self.__son)):
            self.__son[i].__merge_by_label(merged_labels)
    


    # pseudo-overload update_label
    def update_label(self, tip_ann: Union[np.ndarray, None] = None) -> None:
        if type(tip_ann) == np.ndarray:
            self.__update_label_by_ann(tip_ann)
        else:
            self.__update_label()


    # update labels with tip_ann
    # tip_ann is has 2 colunms, 1st col is node names for leaves, 2nd col is the label for the corresponding node
    def __update_label_by_ann(self, tip_ann: np.ndarray) -> None:
        if self.__is_leaf:
            rlt = tip_ann[np.argwhere(tip_ann[...,0] == self.name),1]
            if rlt.shape[0] == 0: # node not in annotation
                raise ValueError(f"node {self.name} not in annotations.")
            elif rlt.shape[0] > 1:# more than one annotation for the node
                for i in range(1, rlt.shape[0]):
                    if rlt[i] != rlt[0]:
                        raise ValueError(f"node {self.name} has more than one annotation")
            self.label = {str(rlt[0][0])}
            return
        label = set()
        for i in range(len(self.__son)):
            self.__son[i].update_label(tip_ann)
            label = label | self.__son[i].label
        self.label = copy.deepcopy(label)
        return
    

    # update labels with sons' labels
    def __update_label(self):
        if self.__is_leaf:
            return
        label = set()
        n_cells = 0
        for son in self.__son:
            son.__update_label()
            label = label | son.label
            n_cells = n_cells + son.n_cells
        self.label = copy.deepcopy(label)
        self.n_cells = copy.deepcopy(n_cells)
        return



    # mark all nodes with name in shuffled_names
    # return a tuple(int, list), int for father's son_shuffle, list is the nodes need to shuffle
    def mark_shuffle_by_name(self, shuffled_names: list, shuffled_nodes: list) -> tuple:
        if self.name in shuffled_names:
            shuffled_nodes.append(self)
            return (1, shuffled_nodes)
        if self.__is_leaf:
            return (0, shuffled_nodes)
        for i in range(len(self.__son)):
            flag, shuffled_nodes = self.__son[i].mark_shuffle_by_name(shuffled_names, shuffled_nodes)
            self.__son_shuffled[i] = flag
        return (0, shuffled_nodes)


    # mark all nodes with label in shuffled_labels and exact potency (length of self.label) and its label is not proper subset of target_labels
    # return same as mark_shuffle_by_name
    def mark_shuffle_by_label(self, shuffled_labels: set, target_labels: set, potency: int, shuffled_nodes: list) -> tuple:
        shuffled_labels=set(shuffled_labels)
        target_labels=set(target_labels)
        if self.label <= shuffled_labels and len(self.label) == potency: # mark if node's label in shuffled_labels and has exact potency
            shuffled_nodes.append(self)
            return (potency, shuffled_nodes)
        if self.__is_leaf or len(self.label) < potency or self.label < target_labels: 
        # len(self.label) < potency   for pruning
        # self.label < target_labels  make sure the node is not proper subset of target_labels
            return (0, shuffled_nodes)
        for i in range(len(self.__son)):
            flag, shuffled_nodes = self.__son[i].mark_shuffle_by_label(shuffled_labels, target_labels, potency, shuffled_nodes)
            if flag:
                self.__son_shuffled[i] = flag
        return (0, shuffled_nodes)



    # count shuffle marks in the tree
    # returns a dict of {potency: numbers}
    def __count_shuffle(self) -> dict:
        if self.__is_leaf:
            return {}
        rlt = []
        keypool = set()
        for i in range(len(self.__son)):
            rlt.append(self.__son[i].__count_shuffle())
            keypool = keypool | set(rlt[-1])
        rtn = {}
        for key in keypool:
            rtn[key] = 0
            for item in rlt:
                if key in item:
                    rtn[key] += item[key]
        for i in range(len(self.__son)):
            if self.__son_shuffled[i] in rtn:
                rtn[self.__son_shuffled[i]] += 1
            else:
                rtn[self.__son_shuffled[i]] = 1
        return rtn


    # pseudo-overload update_shuffle, randomly shuffle the nodes and insert them back to the tree.
    def update_shuffle(self, shuffled_nodes: list) -> None:
#         if type(shuffled_nodes) != list:
#             raise TypeError("illegal argument, only accept list of shuffled nodes.")
        random.seed(os.getpid() + int(time.time() * 100000))
        temp_node = copy.deepcopy(shuffled_nodes)
        count_shuffle = self.__count_shuffle()
        if type(temp_node[0]) == list:
            #print(count_shuffle)
            #print(temp_node)
            for i in range(len(temp_node)):
#                 if len(temp_node[i]) == 0 and not (i+1) in count_shuffle:
#                     continue
#                 if count_shuffle[i+1] != len(temp_node[i]):
#                     raise ValueError(f"number of shuffled nodes and marks in tree does not match. (Potency {i+1})")
                random.shuffle(temp_node[i])
            self.__update_shuffle_by_label(temp_node)
        elif type(temp_node[0]) == TCA_tree:
#             if count_shuffle[1] != len(temp_node):
#                 raise ValueError("number of shuffled nodes and marks in tree does not match.")
            random.shuffle(temp_node)
            self.__update_shuffle(temp_node)
        else:
            raise TypeError("illegal argument, only accept list of shuffled nodes.")
        

    # insert the nodes need to shuffle
    def __update_shuffle(self, shuffled_nodes: list) -> list:
        if self.__is_leaf:
            return shuffled_nodes
        self.n_cells = 0
        for i in range(len(self.__son)):
            if self.__son_shuffled[i]: # insert a random node and reset the son_shuffled mark
                self.__son[i] = shuffled_nodes.pop()
                self.__son_shuffled[i] = 0
            else:
                shuffled_nodes = self.__son[i].__update_shuffle(shuffled_nodes)
            self.n_cells += self.__son[i].n_cells
        return shuffled_nodes
    
    
    # insert the nodes need to shuffle with exact potency
    def __update_shuffle_by_label(self, shuffled_nodes: list) -> list:
        if self.__is_leaf:
            return shuffled_nodes
        self.n_cells = 0
        for i in range(len(self.__son)):
            if self.__son_shuffled[i]: # insert a random node with right potency and reset the son_shuffled mark
                self.__son[i] = copy.deepcopy(shuffled_nodes[self.__son_shuffled[i]-1].pop())
                self.__son_shuffled[i] = 0
            shuffled_nodes = self.__son[i].__update_shuffle_by_label(shuffled_nodes)
            self.n_cells += self.__son[i].n_cells
        return shuffled_nodes


    
    # get all n_cells and label in leaves, used in Np_Estimator
    # returns a dict of {label: [n_cells]}
    def get_leaf_n(self) -> dict:
        if self.__is_leaf:
            for item in self.label:
                rtn = {item: [self.n_cells]}
            return rtn
        rlt = []
        for i in range(len(self.__son)):
            rlt.append(self.__son[i].get_leaf_n())
        key_pool = {}
        rtn = {}
        for i in range(len(rlt)):
            key_pool.update(rlt[i])
        for key in key_pool.keys():
            rtn[key] = []
            for i in range(len(rlt)):
                if key in rlt[i]:
                    rtn[key] = rtn[key] + rlt[i][key]
        return rtn
    

    # get all n_cells and label in leaves, merge all target nodes if the labels are differnt
    # only for random shuffle
    # return same as get_leaf_n
    def get_leaf_n_no_same(self, target_labels: set) -> dict:
        if self.__is_leaf:
            if self.label <= target_labels:
                return {'target': [self.n_cells]}
            else:
                return {'others': []}
        if self.label & target_labels == set():
            return {'others': []}
        if self.label <= target_labels:
            flag_label = False
            temp_label = self.__son[0].label
            for son in self.__son:
                if son.label != temp_label:
                    flag_label = True
                    break
            if flag_label:
                return {'target': [self.n_cells]}
        rlt = []
        rtn = {'target':[], 'others':[]}
        for son in self.__son:
            rlt.append(son.get_leaf_n_no_same(target_labels))
        for item in rlt:
            if 'target' in item:
                rtn['target'] = rtn['target'] + item['target']
        return rtn

                

    # get potency for every node
    # returns a tuple of (rtn, rlt)
    # rlt is a dict of {label: [n_cells]} for this node
    # rtn is a dict of {node name: rlt} for all node
    def get_potency(self) -> tuple:
        if self.__is_leaf:
            #rlt = {self.label:[self.n_cells]}
            rlt = {list(self.label)[0]:[self.n_cells]}
            return ({self.name : rlt}, rlt)
        rlt = []
        for i in range(len(self.__son)):
            rlt.append(self.__son[i].get_potency())
        rtn = {}
        for i in range(len(rlt)):
            rtn.update(rlt[i][0])
        # merge node if all sons are monophyletic clade
        flag_len = True
        for i in range(len(rlt)):
            if len(rlt[i][1]) != 1: # merge condtion: only have one label
                flag_len = False
                break
        if flag_len:
            flag_key = True
            for i in range(len(rlt)):
                for j in range(i+1, len(rlt)):
                    if rlt[i][1].keys() != rlt[j][1].keys(): # merge condition: the labels are the same
                        flag_key = False
                        break
                if not flag_key:
                    break
        rlt_node = {}
        if flag_len and flag_key: # merge
            for key in rlt[0][1].keys():
                sum_ = 0
                for i in range(len(rlt)):
                    sum_ += rlt[i][1][key][0]
                rlt_node[key] = [sum_]
        else:
            key_pool = {}
            for i in range(len(rlt)):
                key_pool.update(rlt[i][1])
            for key in key_pool.keys():
                rlt_node[key] = []
                for i in range(len(rlt)):
                    if key in rlt[i][1]:
                        rlt_node[key] = rlt_node[key] + rlt[i][1][key]
        rtn.update({self.name : rlt_node})
        return (rtn, rlt_node)



    # get all leaves' name
    def get_leaf_names(self) -> list:
        if self.__is_leaf:
            return [self.name]
        else:
            rtn = []
            for i in range(len(self.__son)):
                rtn = rtn + self.__son[i].get_leaf_names()
            return rtn



    # get all label combinations in the tree
    def get_label_set(self, label_set: list) -> list:
        if self.label not in label_set:
            label_set.append(self.label)
        if self.__is_leaf:
            return label_set
        for i in range(len(self.__son)):
            label_set = self.__son[i].get_label_set(label_set)
        return label_set




# tranfer a string to a TCA_tree instance
def read_newick(newick_str: str) -> TCA_tree:
    def split_newick(newick_str: str) -> list:
        n_lbracket = 0
        rtn = []
        last_i = 0
        for i, char in enumerate(newick_str):
            if char == '(':
                n_lbracket += 1
            elif char == ')':
                n_lbracket -= 1
            elif char == ',' and n_lbracket == 0:
                rtn.append(newick_str[last_i:i])
                last_i = i+1
        rtn.append(newick_str[last_i:])
        return rtn
    
    newick_str = newick_str.strip(";\n")
    if newick_str[0] == '(': # not a leaf
        is_leaf = False
        node_info = newick_str[newick_str.rfind(')') + 1:]
        son_str = split_newick(newick_str[1:newick_str.rfind(')')])
        son = []
        for i in range(len(son_str)):
            son.append(read_newick(son_str[i]))
        n_cells = 0
        for i in range(len(son)):
            n_cells += son[i].n_cells
    else:
        is_leaf = True
        son = []
        node_info = newick_str
        n_cells = 1
    if ':' not in node_info: # don't have distance
        name = node_info
        dist = 0
    else:
        name = node_info.split(":")[0]
        dist = float(node_info.split(":")[1])

    return TCA_tree(
        is_leaf = is_leaf,
        son = son,
        n_cells = n_cells,
        label = set(),
        dist = dist,
        name = name,
        son_shuffled = [0] * len(son)
    )



# calculate the Np and clade size for each annotation
# returns a dataframe with 4 columns
# TipAnn is the annotation
# CladeSize is a string of clade sizes and clade numbers
# Total is the total cell number with the annotation
# EffN is the estimated Np for the annotation
def Np_Estiamtor(tree: Union[TCA_tree, str], ann: np.ndarray) -> pd.DataFrame:
    if type(tree) == TCA_tree:
        temp_tree = copy.deepcopy(tree)
    elif type(tree) == str:
        temp_tree = read_newick(tree)
    else:
        raise TypeError("illegal argument. Only accept newick string or TCA_tree class.")
    temp_tree.merge(ann)
    tree_clade = temp_tree.get_leaf_n()
    rlt = calc_eNp(tree_clade)
    #format
    clade_size = []
    total = []
    effN = []
    for ann in sorted(rlt[0]):
        effN.append(rlt[0][ann])
        total.append(rlt[2][ann])
        clade_str = []
        for n_clade in sorted(rlt[1][ann]):
            clade_str.append(f'{n_clade} ({rlt[1][ann][n_clade]})')
        clade_size.append(', '.join(clade_str))
    ann = sorted(rlt[0])
    rtn = pd.DataFrame([ann, clade_size, total, effN], dtype=object).T
    rtn.columns = ['TipAnn', 'CladeSize', 'Total', 'EffN']
    return rtn
    

# calculate the Np, clade size and total for each annotation
# returns a tuple of (Nps, mono_clades, total)
# Nps is a dict of {annotation: Np}
# mono_clades is a dict of {annotation: clade_size}
#   clade_size is a dict of {number of cells in the clade: number of clades}
# total is a dict of {annotation: total cell number with the annotation}
def calc_eNp(tree_clade: dict) -> tuple:
    rtn = {}
    rtn_total = {}
    for label in tree_clade.keys():
        temp = {}
        for item in tree_clade[label]:
            if item in temp:
                temp[item] += 1
            else:
                temp[item] = 1
        tree_clade[label] = temp
        combs = 0
        total = 0
        for clade_size in temp.keys():
            combs = combs + comb(clade_size, 2) * temp[clade_size]
            total = total + clade_size * temp[clade_size]
        if combs == 0:
            rtn[label] = np.inf
        else:
            rtn[label] = comb(total, 2) / combs
        rtn_total[label] = total
    return (rtn, tree_clade, rtn_total)



# calculate the potency for every node
# returns a dataframe with 4 colunms
# NodeName is the name of the node
# Potency is the potency of the node
# TipAnn is one of the annotations of the node
# CladeSize is the total cell number of the corresponding annotation below the node
def calc_potency(tree: TCA_tree, ann: np.ndarray) -> pd.DataFrame:
    ann = np.array(ann)
    tree.update_label(ann)
    ptc = tree.get_potency()[0]
    #format
    NodeName = []
    potency = []
    TipAnn = []
    CladeSize = []
    for node_name in ptc.keys():
        n_cells = 0
        for tip_ann in ptc[node_name].keys():
            n_cells = n_cells + sum(ptc[node_name][tip_ann])
        for tip_ann in ptc[node_name].keys():
            NodeName.append(node_name)
            TipAnn.append(tip_ann)
            CladeSize.append(sum(ptc[node_name][tip_ann]))
            if n_cells == 1:
                potency.append(0)
            else:
                potency.append(len(ptc[node_name]))
    rtn = pd.DataFrame([NodeName, potency, TipAnn, CladeSize])
    rtn = rtn.T
    rtn.columns = ['NodeName', 'Potency', 'TipAnn', 'CladeSize']
    return rtn



# shuffle all nodes with listed names
# return a shuffled tree
def shuffle_tree_by_name(tree: TCA_tree, shuffled_names: list) -> TCA_tree:
    temp_tree = copy.deepcopy(tree)
    shuffled_nodes = temp_tree.mark_shuffle(shuffled_names = shuffled_names, shuffled_nodes = [])[1]
    random.seed(os.getpid() + int(time.time()*10000))
    random.shuffle(shuffled_nodes)
    temp_tree.update_shuffle(shuffled_nodes = shuffled_nodes)
    return temp_tree


def shuffle_tree(tree: Union[TCA_tree, str],  ann: Union[np.ndarray, list, pd.DataFrame],
                 shuffled_labels: Union[set, list, None] = None, target_labels: Union[set, list, None] = None,
                 rep: int = 1, multithread: bool = False, n_thread: Union[int, None] = None) -> Union[list, pd.DataFrame]:
    if type(tree) == TCA_tree:
        temp_tree = copy.deepcopy(tree)
    elif type(tree) == str:
        temp_tree = read_newick(tree)
    else:
        raise TypeError("illegal argument. Only accept newick string or TCA_tree class.")
    if type(ann) == list:
        ann = np.array(ann)
    elif type(ann) == pd.DataFrame:
        ann = np.array(ann)
    if type(shuffled_labels) == list:
        shuffled_labels = set(shuffled_labels)
    if type(target_labels) == list:
        target_labels = set(target_labels)
    if type(rep) == float and rep % 1 == 0:
        rep = int(rep)
    if type(n_thread) == float and n_thread % 1 == 0:
        n_thread = int(n_thread)

    # traverse all label pairs
    if shuffled_labels == None:
        temp_tree.update_label(ann)
        label_set = temp_tree.get_label_set([])
        parent = []
        target = []
        shuffle_rlt = []
        for labels in label_set:
            for sub_labels in label_set:
                if len(sub_labels) == 1 or not sub_labels < labels:
                    continue
                parent.append('|'.join(sorted(labels)))
                target.append('|'.join(sorted(sub_labels)))
                rlt = shuffle_tree(temp_tree, ann, labels, sub_labels, rep, multithread, n_thread)
                shuffle_rlt.append(rlt)
        rtn = pd.DataFrame([parent, target, shuffle_rlt], dtype=object).T
        rtn.columns = ['Parent', 'Target', 'ShuffledNp']
        return rtn
    
    #mark shuffle
    for item in target_labels:
        temp_tree.merge(target_labels - {item})
    shuffled_nodes = []
    for i in range(1, len(target_labels)):
        shuffled_nodes.append(temp_tree.mark_shuffle_by_label(shuffled_labels, target_labels, i, [])[1])

    # multi-thread strategy
    if multithread:
        if n_thread == None:
            n_thread = cpu_count()
        pool = Pool(n_thread)
        rlt = []
        rtn = []
        for i in range(rep):
            #print('multi threading', len(rlt)+1)
            rlt.append(pool.apply_async(func = _run_shuffle,
                                        args = (target_labels, temp_tree, shuffled_nodes)))
        pool.close()
        pool.join()
        for item in rlt:
            rtn.append(item.get())
        return rtn
    
    # single-thread strategy
    rtn = []
    for i in range(rep):
        rtn.append(_run_shuffle(target_labels, temp_tree, shuffled_nodes))
    return rtn

def _run_shuffle(target_labels, tree: TCA_tree, shuffled_nodes):
    temp_tree = copy.deepcopy(tree)
    temp_tree.update_shuffle(shuffled_nodes)
    temp_tree.update_label()
    tree_clade = temp_tree.get_leaf_n_no_same(target_labels)
    rlt = calc_eNp(tree_clade)[0]['target']
    return rlt
