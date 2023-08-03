import json
import os
import glob


def calc_double(n0, p_death, p_nd, double_round):
    n = n0
    p_double = 1 - p_death - p_nd
    for i in range(double_round):
        n_double = round(n * p_double)
        n_nd = round(n * p_nd)
        n = n_double * 2 + n_nd
    return n


def calc_node(node, n, t, tree_struct, diff_time, diff_mode_probs, life_time, tar_t, p_death, p_nd, cor_n):
    if node[0] == '-':
        t_life = tar_t - t
        double_round = int(t_life // life_time[node])
        pop_size = calc_double(n0 = n, p_death = p_death[node], p_nd = p_nd[node], double_round = double_round)
        return {node : pop_size}
    t_life = diff_time[node] - t
    double_round = int(t_life // life_time[node])
    pop_size = calc_double(n0 = n, p_death = p_death[node], p_nd = p_nd[node], double_round = double_round)
    #rtn = {node : pop_size}

    n_ii = round(pop_size * diff_mode_probs[node][0])
    n_jj = round(pop_size * diff_mode_probs[node][1])
    n_ij = round(pop_size * diff_mode_probs[node][2])
    n_i = n_ii * 2 + n_ij
    n_j = n_jj * 2 + n_ij
    if n_i + n_j < 8:
        n_i = 4
        n_j = 4
    if n_i < 4:
        n_j = n_j - 4 + n_i
        n_i = 4
    if n_j < 4:
        n_i = n_i - 4 + n_j
        n_j = 4
    cor_n_i = n_ii + n_ij
    cor_n_j = n_jj + n_ij
    if cor_n_i < 2:
        cor_n_i = 2
    if cor_n_j < 2:
        cor_n_j = 2
    pop_i = calc_node(node = tree_struct[node][0], n = n_i, t = t + life_time[node] * (double_round + 1),
                    tree_struct = tree_struct, diff_time = diff_time, cor_n = cor_n_i,
                    diff_mode_probs = diff_mode_probs, life_time = life_time, tar_t = tar_t, p_death = p_death, p_nd = p_nd)
    pop_j = calc_node(node = tree_struct[node][1], n = n_j, t = t + life_time[node] * (double_round + 1),
                    tree_struct = tree_struct, diff_time = diff_time, cor_n = cor_n_j,
                    diff_mode_probs = diff_mode_probs, life_time = life_time, tar_t = tar_t, p_death = p_death, p_nd = p_nd)
    rtn = {node : pop_size}
    rtn.update(pop_i)
    rtn.update(pop_j)
    return rtn


def calc_pop(info):
    if type(info) == str:
        info = json.loads(info)
    tree_struct = info['merge']
    diff_time = info['diff_time']
    diff_mode_probs = info['diff_mode_probs']
    life_time = info['lifetime']
    n0 = info['founder_size']
    tar_t = info['target_time']
    p_death = info['prob_loss']
    p_nd = info['prob_pm']

    node_pop = calc_node(node = info['root_id'], n = n0, t = 0, tree_struct = tree_struct, diff_time = diff_time,
                        diff_mode_probs = diff_mode_probs, life_time = life_time, tar_t = tar_t,
                        p_death = p_death, p_nd = p_nd, cor_n = 1)
    return node_pop


if __name__ == "__main__":
    read_path = './supplementary_data_1_fate_map_panel/'
    write_path = './progenitor_population_size_act/'

    for path_json in glob.glob(read_path + "fate_map*.json"):
        file_json = open(path_json)
        read_json = file_json.read()
        info = json.loads(read_json)
        pop_size = calc_pop(info)
        file_out = open(write_path + "pop_size_" + path_json.split('/')[-1].split('.')[0] + ".csv", 'w')
        for key in pop_size.keys():
            file_out.write('"' + key + '",' + str(pop_size[key]) + '\n')

