#!/usr/bin/python
# vim: set fileencoding=utf-8 :

import networkx as nx
import scipy.io as sio
import numpy as np
import time;

#计算L_Zl

#读取.mat文件
time_start=time.time();
load_fn='/home/bubble/PycharmProjects/bnfinder/relation_gse12346_432.mat'
load_data = sio.loadmat(load_fn)

print(load_data['num'])
time_end=time.time();#time.time()为1970.1.1到当前时间的毫秒数
print time_end-time_start,
print "s"

#邻接矩阵转化为图
time_start=time.time();
G = nx.from_numpy_matrix(np.array(load_data['num']))
time_end=time.time();#time.time()为1970.1.1到当前时间的毫秒数
print time_end-time_start,
print "s"


#度中心性计算
time_start=time.time();
filename='/home/bubble/PycharmProjects/bnfinder/relation_gse12346_432degree_centrality.txt'
with open(filename,'w') as file_object:
	file_object.write(str(nx.degree_centrality(G)))
print nx.degree_centrality(G)					#度中心性计算
time_end=time.time();				#time.time()为1970.1.1到当前时间的毫秒数
print time_end-time_start,
print "s"


#计算接近中心性
time_start=time.time();
filename='/home/bubble/PycharmProjects/bnfinder/relation_gse12346_432closeness_centrality.txt'
with open(filename,'w') as file_object:
	file_object.write(str(nx.closeness_centrality(G,None,None,True,False)))
print nx.closeness_centrality(G,None,None,True,False)              #接近中心性
time_end=time.time();#time.time()为1970.1.1到当前时间的毫秒数
print time_end-time_start,
print "s"


#计算中介中心性
time_start=time.time();
filename='/home/bubble/PycharmProjects/bnfinder/relation_gse12346_432betweenness_centrality.txt'
with open(filename,'w') as file_object:
	file_object.write(str(nx.betweenness_centrality(G,None,True,None,False,None)))
print nx.betweenness_centrality(G,None,True,None,False,None)              #接近中心性
time_end=time.time();#time.time()为1970.1.1到当前时间的毫秒数
print time_end-time_start,
print "s"


#计算特征相量中心性
time_start=time.time();

filename='/home/bubble/PycharmProjects/bnfinder/relation_gse12346_432eigenvector_centrality.txt'
with open(filename,'w') as file_object:
	file_object.write(str(nx.nx.eigenvector_centrality(G)))
print nx.eigenvector_centrality(G)              #接近中心性
time_end=time.time();#time.time()为1970.1.1到当前时间的毫秒数
print time_end-time_start,
print "s"