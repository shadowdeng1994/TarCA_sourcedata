import copy
import csv
import tqdm

path = "D:\\SYSU\\helab\\TCA\\QFM\\"

def read_tip(fate_map_id):
    path_map = path + "cor_population_size\\cor_pop_size_fate_map" + str(fate_map_id).rjust(4, '0') + ".csv"
    file_map = open(path_map)
    csv_map = csv.reader(file_map)
    path_bias = path + "tip_bias\\tip_bias_fate_map" + str(fate_map_id).rjust(4, '0') + ".csv"
    file_bias = open(path_bias)
    csv_bias = csv.reader(file_bias)
    rtn = {}
    map_dict = {}
    for line in csv_map:
        map_dict[line[0]] = line[1]
    for line in csv_bias:
        rtn[line[0]] = [map_dict[line[0]], line[1]]
    return rtn


def read_tca(tca_id):
    path_tca = path + "TCA_population_size\\" + str(tca_id).rjust(4, '0') + ".csv"
    file_tca = open(path_tca)
    csv_tca = csv.reader(file_tca)
    rtn = {}
    for line in csv_tca:
        rtn[line[0]] = line[1]
    return rtn


if __name__ == "__main__":
    path_fix = path + "cor_population_size_bias\\fixed\\"
    path_prop = path + "cor_population_size_bias\\proportional\\"
    fate_map_fix = {}
    for i in range(331):
        fate_map_fix[str(i + 1)] = read_tip(i + 1)
    fate_map_prop = copy.deepcopy(fate_map_fix)
    for i in tqdm.tqdm(range(3310)):
        result =  read_tca(i+1)
        map_id = str((i % 331) + 1)
        if (i % 662) + 1 > 331:
            for key in result.keys():
                fate_map_prop[map_id][key].append(result[key])
        else:
            for key in result.keys():
                #print(fate_map_fix.keys())
                #print(fate_map_fix[map_id])
                #print(map_id, key)
                #print(result)
                fate_map_fix[map_id][key].append(result[key])
    
    for map_id in fate_map_fix.keys():
        file_fix = open(path_fix + map_id.rjust(4, '0') + ".csv", 'w')
        for key in fate_map_fix[map_id].keys():
            file_fix.write(key + ',' + ','.join(fate_map_fix[map_id][key]) + '\n')
        file_fix.close()
        file_prop = open(path_prop + map_id.rjust(4, '0') + ".csv", 'w')
        for key in fate_map_prop[map_id].keys():
            file_prop.write(key + ',' + ','.join(fate_map_prop[map_id][key]) + '\n')
        file_prop.close()
