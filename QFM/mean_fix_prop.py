import glob
import numpy as np 
import matplotlib.pyplot as plt
import scipy.stats

path = "D:\\SYSU\\helab\\TCA\\QFM\\"

fix_ele = {}
prop_ele = {}
for path_merge in glob.glob(path + "merge_population_size\\fixed\\*.csv"):
    read_merge = np.loadtxt(path_merge, delimiter=',')
    name_merge = path_merge.split('\\')[-1].split('.')[0]
    for row in read_merge:
        name = name_merge + str(row[0])
        mean = np.mean(row[~np.isinf(row)][2:])
        fix_ele[name] = mean
for path_merge in glob.glob(path + "merge_population_size\\proportional\\*.csv"):
    read_merge = np.loadtxt(path_merge, delimiter=',')
    name_merge = path_merge.split('\\')[-1].split('.')[0]
    for row in read_merge:
        name = name_merge + str(row[0])
        mean = np.mean(row[~np.isinf(row)][2:])
        prop_ele[name] = mean

merge_ele = np.array([[0,0]])
error_keys = []
line = []
line_keys = []
for key in fix_ele.keys():
    if np.isnan(fix_ele[key]) or np.isnan(prop_ele[key]):
        continue
    row = [[fix_ele[key], prop_ele[key]]]
    #if row[0][1] > 5.9 and row[0][1] < 6.1:
    #    continue
    merge_ele = np.concatenate((merge_ele, row), axis = 0)
    if row[0][0] / row[0][1] > 2 or row[0][0] / row[0][1] < 0.5:
        error_keys.append(key)
    if row[0][1] > 5.5 and row[0][1] < 6.5:
        line.append(row)
        line_keys.append(key)
merge_ele = merge_ele[1:,...]
merge_ele = np.log(merge_ele) / np.log(2)

#line = np.array(line)
#line = np.log(line) / np.log(2)
#plt.scatter(line[...,0], line[...,1], s = 3, alpha = 0.1, marker = '.')
plt.scatter(merge_ele[...,0], merge_ele[...,1], s = 3, alpha = 0.2, marker = '.')
plt.xlabel("log${_2}$(fixed sampling eNp)")
plt.ylabel("log${_2}$(proportional sampling eNp)")
plt.title("Estimate progenitor sizes")
plt.show()

print(scipy.stats.pearsonr(merge_ele[...,0], merge_ele[...,1]))
print(len(error_keys), len(fix_ele))
print(len(line_keys))