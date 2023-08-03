import glob
import re


def extract_node_info(tree_newick : str) -> dict[str:str]:
    rtn = {}
    nodes = [item.split(':')[0] for item in re.split('[(),]', tree_newick) if item != '']
    for item in nodes:
        label = item.split('_')[1]
        if label[0] == '-':
            rtn[item] = label
    return rtn


if __name__ == "__main__":
    read_path = "./supplementary_data_3_simulated_experiments/simulated_phylogeny/"
    write_path = "./node_ann/"

    for path_newick in glob.glob(read_path + "*.newick"):
        #print(path_newick)
        file_newick = open(path_newick)
        read_newick = file_newick.read()
        node_info = extract_node_info(read_newick)
        file_out = open(write_path + path_newick.split('/')[-1].split('.')[0] + ".csv", 'w')
        for key in node_info.keys():
            file_out.write(key + ',' + node_info[key] + '\n')
