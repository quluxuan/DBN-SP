#coding: utf-8
from heapq import *
import graph,stats,math,fpconst
import numpy as np
from numpy import *
from com_mi import com_mi
from minepy import MINE
from multiprocessing import Pool
import math
import numpy.linalg
###计算基因之间的关系矩阵（互信息mi_matrix，皮尔逊系数pearson_matrix，最大互信息mic_matrix）
###首先在计算某个基因与其他基因进行互信息计算时实现多线程，完成后释放/关闭？；
def compute_relation(self):
    gene_num = self.n_gen
    mi_matrix = np.zeros([gene_num,gene_num])
    mi_matrix_temp=[]
    pearson_matrix = np.zeros([gene_num,gene_num])
    pearson_matrix_temp=[]
    mic_matrix = np.zeros([gene_num,gene_num])
    mic_matrix_temp=[]
    for node in self.vertices:
        parents_select = self.regulators
        import datetime
        distrib2 = 5
        # pool = Pool(distrib2)
        ###每次的mi_temp汇成mi_matrix，mi_matrix转换成矩阵，补充对角线元素
        # rela = pool.map(compute_map,[(i,self,node) for i in range(len(parents_select))])
        rela = map(compute_map, [(i, self, node) for i in range(len(parents_select))])
        relaa=[]
        for item in rela:
            item=list(item)
            relaa.append(item)
        mi = [x[0] for x in relaa]
        pearson = [x[1] for x in relaa]
        micc = [x[2] for x in relaa]
        # pool.close()
        mi_matrix_temp.append(mi)
        pearson_matrix_temp.append(pearson)
        mic_matrix_temp.append(micc)
        nowTime = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print nowTime
        for i in range(len(parents_select)):
            if i>node.index:
                iii = mi_matrix_temp[node.index][i]
                nnn = mi_matrix[node.index][i]
                mmm = mi_matrix[i][node.index]
                mi_matrix[node.index][i] = mi_matrix_temp[node.index][i]
                mi_matrix[i][node.index] = mi_matrix_temp[node.index][i]
                pearson_matrix[node.index][i] = pearson_matrix_temp[node.index][i]
                pearson_matrix[i][node.index] = pearson_matrix_temp[node.index][i]
                mic_matrix[node.index][i] = mic_matrix_temp[node.index][i]
                mic_matrix[i][node.index] = mic_matrix_temp[node.index][i]

        ove=1
    return pearson_matrix,mi_matrix,mic_matrix
def compute_map(v):
    i, self, node = v
    parents_select = self.regulators
    mi=0
    mine = MINE()
    pearson=0
    micc=0
    if node.index < i:
        #===========compute mi
        a = np.cov(parents_select[node.index])
        b = np.cov(parents_select[i])
        AB = np.cov(parents_select[node.index], parents_select[i])
        ab = numpy.linalg.det(AB)
        mi = 0.5 * math.log(2) * ((a * b) / ab)
        # mi=com_mi(parents_select[node.index], parents_select[i])
        #===========compute pearson
        pearson = stats.pearsonr(parents_select[node.index], parents_select[i])
        pearson = abs(pearson[0])
        #===========compute MIC
        mine.compute_score(parents_select[node.index], parents_select[i])
        micc = mine.mic()
    return mi,pearson,micc
