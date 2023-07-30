import numpy as np
import rpy2.robjects as rob
import tqdm

path = "D:/SYSU/helab/TCA/QFM/"

def read_total(sample_id):
    path_sample = path + "QFM2023/" + str(sample_id).rjust(4, '0') + ".RData"
    r_script = '''
    load("{path}")
    ann = tmp.out$EffN$TipAnn
    total = tmp.out$EffN$Total
    clade = tmp.out$EffN$Clade
    '''.format(path = path_sample)
    rob.r(r_script)
    ann = rob.r['ann']
    ann = np.array(ann, dtype=int)
    total = rob.r['total']
    total = np.array(total, dtype=int)
    clade = rob.r['clade']
    clade = np.array(clade, dtype=str)
    rlt = np.dstack((ann, total, clade))[0]
    np.save(path + "QFM2023_npy/" + str(sample_id).rjust(4, '0') + ".npy", rlt)
    return

for i in tqdm.tqdm(range(3310)):
    read_total(i+1)