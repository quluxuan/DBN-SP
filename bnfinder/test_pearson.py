import numpy as np
import xlrd
import xlwt
import math
import stats

from minepy import MINE
from pyExcelerator import *
from heapq import *

from com_mi import com_mi
if __name__ == '__main__':

    aaa=1
    f = open("ecoli100_1565.txt", "r")
    Gene = f.readlines()
    Gene1=[]
    f.close()
    for item in Gene:
        item = item.split("\t",1)
        geneno=item[0]
        gene_data=item[1].split("\t")

        gene_data = map(eval, gene_data)

        item=[geneno,gene_data]
        Gene1.append(item)
    Gene2=Gene1
    parents_num=int(Gene1.__len__()*0.3)
    ###################pearsonr###################
    parents_pearson=[]
    for item1 in Gene1:
        p_n=parents_num
        parents_potential = []
        for item2 in Gene2:
            correl = stats.pearsonr(item1[1], item2[1])
            abscorrel = abs(correl[0])
            if p_n:
                heappush(parents_potential, (abscorrel, item1[0], item2[0]))
                p_n -= 1
            elif parents_potential:
                if parents_potential[0][0] < abscorrel:
                    heapreplace(parents_potential, (abscorrel, item1[0], item2[0]))
        parents_pearson.append(parents_potential)
    ###################spearmanr###################
    # parents_spearmanr=[]
    # for item1 in Gene1:
    #     p_n = parents_num
    #     parents_potential = []
    #     for item2 in Gene2:
    #         if item1[1] == item2[1]:
    #            continue
    #         correl = stats.spearmanr(item1[1], item2[1])
    #         abscorrel = abs(correl[0])
    #         if p_n:
    #             heappush(parents_potential, (abscorrel, item1[0], item2[0]))
    #             p_n -= 1
    #         elif parents_potential:
    #             if parents_potential[0][0] < abscorrel:
    #                 heapreplace(parents_potential, (abscorrel, item1[0], item2[0]))
    #     parents_spearmanr.append(parents_potential)
    #
    #
    #
    # ###################mi###################
    # parents_mi=[]
    # for item1 in Gene1:
    #     p_n = parents_num
    #     parents_potential = []
    #     for item2 in Gene2:
    #         # aa=np.array(item1[1])
    #         # bb=np.array(item2[1])
    #         if item1[1] == item2[1]:
    #             continue
    #         aa=item1[1]
    #         bb=item2[1]
    #         correl = com_mi(item1[1], item2[1])
    #         # abscorrel = abs(correl)
    #         if p_n:
    #             heappush(parents_potential, (correl, item1[0], item2[0]))
    #             p_n -= 1
    #         elif parents_potential:
    #             if parents_potential[0][0] < correl:
    #                 heapreplace(parents_potential, (correl, item1[0], item2[0]))
    #     parents_mi.append(parents_potential)
    #
    #
    # ###################MIC###################
    # parents_mic=[]
    # for item1 in Gene1:
    #     p_n = parents_num
    #     parents_potential = []
    #     for item2 in Gene2:
    #         # aa=np.array(item1[1])
    #         # bb=np.array(item2[1])
    #         if item1[1] == item2[1]:
    #             continue
    #         aa=item1[1]
    #         bb=item2[1]
    #         mine = MINE()
    #         mine.compute_score(item1[1], item2[1])
    #         mi = mine.mic()
    #         if p_n:
    #             heappush(parents_potential, (mi, item1[0], item2[0]))
    #             p_n -= 1
    #         elif parents_potential:
    #             if parents_potential[0][0] < mi:
    #                 heapreplace(parents_potential, (mi, item1[0], item2[0]))
    #     parents_mic.append(parents_potential)
    ###################entropy###################
    # parents=[]
    # for item1 in Gene1:
    #     parents_potential = []
    #     parents_num = 15
    #     for item2 in Gene2:
    #         if item1[1] == item2[1]:
    #            continue
    #
    #         aa = np.array(item1[1])
    #         bb = np.array(item2[1])
    #
    #         # c2, p = chisquare(aa, bb)
    #         correl = scipy.stats.entropy(aa, bb)
    #         abscorrel = abs(correl)
    #         if parents_num:
    #             heappush(parents_potential, (abscorrel, item1[0], item2[0]))
    #             parents_num -= 1
    #         elif parents_potential:
    #             if parents_potential[0][0] < abscorrel:
    #                 heapreplace(parents_potential, (abscorrel, item1[0], item2[0]))
    #     parents.append(parents_potential)

    # for i_parents in parents:
    #     for i_gene in i_parents:
    #         wa=i_gene[0]
    #         wb=i_gene[1]
    #         wc=i_gene[2]
    parents=[]
    # for i in (0,Gene1.__len__()) :
    #     p = list(set(parents_pearson[i]).union(set(parents_spearmanr[i])))
    #     p = list(set(p).union(set(parents_mi[i])))
    #     p = list(set(p).union(set(parents_mic[i])))
    #     parents.append(p)

    parents = parents_pearson
    w = Workbook()
    ws = w.add_sheet('result')
    ws.write(0, 0, 'correlation')
    ws.write(0, 1, 'gene1')
    ws.write(0, 2, 'gene2')
    i = 1
    for i_parents in parents:
        for i_gene in i_parents:
            ws.write(i, 0, i_gene[1])
            ws.write(i, 1, i_gene[2])
            ws.write(i, 2, i_gene[0])
            i += 1
    w.save('eeeee.xlsx')



