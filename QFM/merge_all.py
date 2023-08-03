import numpy as np
import rpy2.robjects as rob
import copy
import csv
import tqdm

tmpl = {"QFM":1, "sampling method":"1", "rep":1, "tip":1, "pop_size_act":1, "pop_size_source":1, "pop_size_cor":1, "bias":1, "eNp":0}

path = "D:/SYSU/helab/TCA/QFM/"

def read_act(fate_map_id, map_dict):
    path_act = path + "progenitor_population_size_act/pop_size_fate_map" + str(fate_map_id).rjust(4, '0') + ".csv"
    file_act = open(path_act)
    csv_act = csv.reader(file_act)
    for line in csv_act:
        #if line[0][0] != '-':
        #    continue
        temp_dict = {"QFM":fate_map_id, "pop_size_act":line[1], "tip":line[0]}
        temp_dict["sampling method"] = "fixed"
        for i in range(5):
            temp_dict["rep"] = i + 1
            #print(str(fate_map_id) + '/' +  str(i+1) + '/' +  str(line[0]) + 'f')
            map_dict[str(fate_map_id) + '/' +  str(i+1) + '/' +  str(line[0]) + 'f'] = copy.deepcopy(temp_dict) #QFM + rep + tip + sample_method
            #if fate_map_id == 1 and line[0] == "-1" and i == 0:
            #    print(map_dict['1/1/-1f'])
        temp_dict["sampling method"] = "proportional"
        for i in range(5):
            temp_dict["rep"] = i + 1
            map_dict[str(fate_map_id) + '/' +  str(i+1) + '/' +  str(line[0]) + 'p'] = copy.deepcopy(temp_dict)
    return map_dict

def read_source(fate_map_id, map_dict):
    path_act = path + "progenitor_population_size/tip_pop_size_fate_map" + str(fate_map_id).rjust(4, '0') + ".csv"
    file_act = open(path_act)
    csv_act = csv.reader(file_act)
    for line in csv_act:
        for i in range(5):
            map_dict[str(fate_map_id) + '/' +  str(i+1) + '/' +  str(line[0]) + 'f']["pop_size_source"] = line[1] #QFM + rep + tip + sample_method
        for i in range(5):
            map_dict[str(fate_map_id) + '/' +  str(i+1) + '/' +  str(line[0]) + 'p']["pop_size_source"] = line[1]
    return map_dict

def read_cor(fate_map_id, map_dict):
    path_act = path + "progenitor_population_size_cor/cor_pop_size_fate_map" + str(fate_map_id).rjust(4, '0') + ".csv"
    file_act = open(path_act)
    csv_act = csv.reader(file_act)
    for line in csv_act:
        for i in range(5):
            map_dict[str(fate_map_id) + '/' +  str(i+1) + '/' +  str(line[0]) + 'f']["pop_size_cor"] = line[1] #QFM + rep + tip + sample_method
        for i in range(5):
            map_dict[str(fate_map_id) + '/' +  str(i+1) + '/' +  str(line[0]) + 'p']["pop_size_cor"] = line[1]
    return map_dict

def read_bias(fate_map_id, map_dict):
    path_act = path + "progenitor_bias/tip_bias_fate_map" + str(fate_map_id).rjust(4, '0') + ".csv"
    file_act = open(path_act)
    csv_act = csv.reader(file_act)
    for line in csv_act:
        for i in range(5):
            map_dict[str(fate_map_id) + '/' +  str(i+1) + '/' +  str(line[0]) + 'f']["bias"] = line[1] #QFM + rep + tip + sample_method
        for i in range(5):
            map_dict[str(fate_map_id) + '/' +  str(i+1) + '/' +  str(line[0]) + 'p']["bias"] = line[1]
    return map_dict

def read_prop(fate_map_id, map_dict):
    path_act = path + "fraction/fraction_fate_map" + str(fate_map_id).rjust(4, '0') + ".csv"
    file_act = open(path_act)
    csv_act = csv.reader(file_act)
    for line in csv_act:
        for i in range(5):
            map_dict[str(fate_map_id) + '/' +  str(i+1) + '/' +  str(line[0]) + 'f']["pure_fraction"] = line[1] #QFM + rep + tip + sample_method
            map_dict[str(fate_map_id) + '/' +  str(i+1) + '/' +  str(line[0]) + 'f']["mixed_fraction"] = line[2] #QFM + rep + tip + sample_method
        for i in range(5):
            map_dict[str(fate_map_id) + '/' +  str(i+1) + '/' +  str(line[0]) + 'p']["pure_fraction"] = line[1]
            map_dict[str(fate_map_id) + '/' +  str(i+1) + '/' +  str(line[0]) + 'p']["mixed_fraction"] = line[2]
    return map_dict


def read_tca(tca_id, map_dict):
    path_tca = path + "TCA_population_size/" + str(tca_id).rjust(4, '0') + ".csv"
    file_tca = open(path_tca)
    csv_tca = csv.reader(file_tca)
    fate_map_id = ((tca_id - 1) % 331) + 1
    rep = ((tca_id - 1) // 662) + 1
    if ((tca_id - 1) % 662) + 1 > 331:
        flag = 'p'
    else:
        flag = 'f'
    for line in csv_tca:
        map_dict[str(fate_map_id) + '/' + str(rep) + '/' + str(line[0]) + flag]["eNp"] = line[1]
    return map_dict

def read_total(sample_id, map_dict):
    path_sample = path + "QFM2023_progenitor_npy/" + str(sample_id).rjust(4, '0') + ".npy"
    rlt = np.load(path_sample)
    fate_map_id = ((sample_id - 1) % 331) + 1
    rep = ((sample_id - 1) // 662) + 1
    if ((sample_id - 1) % 662) + 1 > 331:
        flag = 'p'
    else:
        flag = 'f'
    for item in rlt:
        map_dict[str(fate_map_id) + '/' + str(rep) + '/' + str(item[0]) + flag]["eNp"] = item[3]
        map_dict[str(fate_map_id) + '/' + str(rep) + '/' + str(item[0]) + flag]["Total"] = item[2]
        map_dict[str(fate_map_id) + '/' + str(rep) + '/' + str(item[0]) + flag]["NClade"] = calc_clade(item[1])
    return map_dict


def calc_clade(str_clade):
    clade = str_clade.split(", ")
    n_clade = 0
    for item in clade:
        temp = item.split(' (')[1].strip(')')
        n_clade = n_clade + int(temp)
    return n_clade


map_dict = {}
for i in range(331):
    map_dict = read_act(i+1, map_dict)
    map_dict = read_cor(i+1, map_dict)
    map_dict = read_bias(i+1, map_dict)
    map_dict = read_source(i+1, map_dict)
    #map_dict = read_prop(i+1, map_dict)
for i in tqdm.tqdm(range(3310)):
    #map_dict = read_tca(i+1, map_dict)
    map_dict = read_total(i+1, map_dict)

file_merge = open(path + "QFM_merge.csv", 'w')
file_merge.write("QFM,sampling method,rep,tissue,pop_size_terminal,pop_size,pop_size_corrected,bias,total,eNp,n_clade\n")
for key in map_dict.keys():
    item = map_dict[key]
    write_item = [item["QFM"], item["sampling method"], item["rep"], item["tip"], item["pop_size_act"], item["pop_size_source"], item["pop_size_cor"],
                  item["bias"], item["Total"], item["eNp"], item["NClade"]]
    write_item = [str(item) for item in write_item]
    file_merge.write(','.join(write_item) + '\n')
file_merge.close()
