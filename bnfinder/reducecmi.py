from conditional_mi import cmi
from numpy import mat
import numpy as np
import sys
import itertools
import math


lamda = 0.1
order0 = 2


def reduce_cmi(data_matrix1,lamda,order0,score_matrix1,gene_no):
    num=np.array(data_matrix1)
    n_gene=num.shape[0]
    G=score_matrix1
    Gval=G
    order=-1
    t=0
    gene_name=gene_no
    while t==0:
        order=order+1
        if len(sys.argv)==4:
            if order>order0:
                G=np.transpose(np.tril(G,-1))
                Gval=np.transpose(np.tril(Gval,-1))
                order=order-1
                return G,Gval,order
        #return G, Gval
        [G,Gval,t]=EdgeReduce(G,Gval,order,data_matrix1,t,lamda)
        #[G, Gval, t] = EdgeReduce(G, Gval,  data_matrix1, t, lamda)
        if t==0:

            #net=[]

            for t1 in gene_name:
                for t2 in gene_name:
                    if G[t1[0],t2[0]]!=0:
                        #print 'G'+str(t1+1) +'\t'+str(float(Gval[t1,t2]))+'\tG'+str(t2+1)
                        print t1[1] + '\t' + str(float(Gval[t1[0], t2[0]])) + '\t' + t2[1]

            print "No edge is reduce! Algorithm  finished!"
            break
        else:
            t=0
    G=np.transpose(np.tril(G,-1))
    Gval=np.transpose(np.tril(Gval,-1))
    order=order-1
    return G,Gval
def EdgeReduce(G,Gval,order,data,t,lamda):
#def EdgeReduce(G, Gval,  data, t, lamda):
    G1=np.array(G)
    Gval1=np.array(Gval)
    if order==0:
        for i in range(G1.shape[0]-1):
            for j in range(i+1,G1.shape[0]):
                if (G1[i,j]!=0)or(G1[j,i]!=0):
                    cmiv=cmi(data[i],data[j],0)
                    Gval1[i,j]=cmiv
                    Gval1[j,i]=cmiv
                    if cmiv<lamda:
                        G1[i,j]=0
                        G1[j,i]=0
        t=t+1
    else:
        for i in range(G1.shape[0]-1):
            for j in range(i+1,G1.shape[0]):
                if G1[i,j]!=0:
                    adj=[]

                    for k in range(G1.shape[0]):
                        if (G1[i,k]!=0)and(G1[j,k]!=0):
                            adj.append(k)

                    adj1=np.array(adj)
                    if adj1.shape[0]>=order:
                        combntnslist = list(itertools.combinations(adj1, order))
                        combntnslist1=np.array(combntnslist)
                        combntnsrow=combntnslist1.shape[0]
                        cmiv=0
                        v1=data[i]
                        v2=data[j]
                        for k in range(combntnsrow):
                            kk=0
                            vcs=data[combntnslist1[k,kk]]
                            kk=kk+1
                            a=cmi(v1,v2,vcs)
                            cmiv=max(cmiv,a)
                        Gval1[i,j]=cmiv
                        Gval1[j,i]=cmiv
                        if cmiv<lamda:
                            G1[i,j]=0
                            G1[j,i]=0
                        t=t+1
    return G1,Gval1,t


# def cmi(v1,v2,vcs):
#     if len(sys.argv)==2:
#         c1=np.linalg.det(np.cov(v1))
#         c2=np.linalg.det(np.cov(v2))
#         c3=np.linalg.det(np.cov(v1,v2))
#         cmiv=0.5*math.log((c1*c2)/c3)
#     elif len(sys.argv)==3:
#         c1 = np.linalg.det(np.cov(np.trans(np.concatenate((v1,vcs),axis=0))))
#         c2 = np.linalg.det(np.cov(np.trans(np.concatenate((v2, vcs), axis=0))))
#         c3 = np.linalg.det(np.cov(np.trans(vcs)))
#         c4 = np.linalg.det(np.cov(np.trans(np.concatenate((v1, v2,vcs), axis=0))))
#         cmiv=0.5*math.log((c1*c2)/(c3*c4))
#     if cmiv==float("inf"):
#         cmiv=0
#     return cmiv



f = open("E:/a/python27/bnfinder1/doc/data/dream3_10&50&100/data1001.txt", "r")
f1 = open("E:/a/python27/bnfinder1/doc/data/dream3_10&50&100/data100_20%.sif", "r")

rec=[]
gene=[]
gene_no=[]
test_p=[]
ln = f.readline()
rec = ln.strip().split()
gene.extend(rec[1:])
i=0
for item in gene:
    gene_no.append([i,item])
    i=i+1
gene_matrix=[]

score_matrix1 = [[0 for j in range(0, i)] for k in range(0, i)]
data_matrix1 = [[0 for j in range(0, 1565)] for k in range(0, i)]
data_line1=f.readlines()
ano = 0
for item in data_line1:
    bno=0
    data = item.split('\t')
    for c in data[1:]:
        data_matrix1[ano][bno]=float(c)
        bno=bno+1
    ano=ano+1


test_line1 = f1.readlines()
#######f1
aaano=0
bbbno=0
for item in test_line1:
    gene_data = item.split("\t")
    aaa=gene_data[0]
    bbb=gene_data[2].strip('\n')
    nnn = eval(gene_data[1])
    for i2 in gene_no:
        if aaa==i2[1]:
            aaano=i2[0]
            continue
        if bbb==i2[1]:
            bbbno=i2[0]
            continue
        # if aaano and bbbno:
        #     break
    score_matrix1[aaano][bbbno]=nnn
# score_matrixx1=[]
# for j in range(0,i):
#     sm=score_matrix1[j]
#     popo=sm.pop(j)
#     score_matrixx1=score_matrixx1+sm
[G,Gval]=reduce_cmi(data_matrix1,lamda,order0,score_matrix1,gene_no)