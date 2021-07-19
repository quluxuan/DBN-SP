#coding: utf-8
import sys
from scipy import interpolate

from scipy.interpolate import spline
reload(sys)
sys.setdefaultencoding('utf8')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
# from sklearn import cross_validation
from sklearn.ensemble import GradientBoostingClassifier
f = open("E:/a/python27/bnfinder2/doc/data/1565data/dream10/insilico_size10_1_timeseries.txt", "r")
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
test_matrix1 = [[0 for j in range(0, i)] for k in range(0, i)]
# score_matrix2 = [[0 for j in range(0, i)] for k in range(0, i)]
# test_matrix2 = [[0 for j in range(0, i)] for k in range(0, i)]

# score_matrix5 = [[0 for j in range(0, i)] for k in range(0, i)]
#test_matrix5 = [[0 for j in range(0, i)] for k in range(0, i)]
f1 = open("E:/a/python27/bnfinder2/doc/data/1565data/dream10/10ts_1.sif", "r")
# f2 = open("C:/Users/BUBBLE/PycharmProjects/bnfinder/doc/data/gdexp/40roc-2", "r")

#f5 = open("E:/a/python27/bnfinder1/doc/data/1565data/10roc", "r")
test_line1 = f1.readlines()
# test_line1 = f1.readlines()
# test_line2 = f2.readlines()
# test_line2 = f2.readlines()

#test_line5 = f5.readlines()
# test_line5 = f5.readlines()
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
score_matrixx1=[]
for j in range(0,i):
    sm=score_matrix1[j]
    popo=sm.pop(j)
    score_matrixx1=score_matrixx1+sm
#######f2
# aaano1=0
# bbbno1=0
# for item in test_line2:
#     gene_data = item.split("\t")
#     aaa=gene_data[0]
#     bbb=gene_data[2].strip('\n')
#     nnn = eval(gene_data[1])
#     for i2 in gene_no:
#         if aaa==i2[1]:
#             aaano=i2[0]
#             continue
#         if bbb==i2[1]:
#             bbbno=i2[0]
#             continue
#         # if aaano and bbbno:
#         #     break
#     score_matrix2[aaano][bbbno]=nnn
# score_matrixx2=[]
# for j in range(0,i):
#     sm=score_matrix2[j]
#     popo=sm.pop(j)
#     score_matrixx2=score_matrixx2+sm
#######f3

#######f5
# for item in test_line5:
#     gene_data = item.split("\t")
#     aaa=gene_data[0]
#     bbb=gene_data[2].strip('\n')
#     nnn = eval(gene_data[1])
#     for i2 in gene_no:
#         if aaa==i2[1]:
#             aaano=i2[0]
#             continue
#         if bbb==i2[1]:
#             bbbno=i2[0]
#             continue
#         # if aaano and bbbno:
#         #     break
#     score_matrix5[aaano][bbbno]=nnn
# score_matrixx5=[]
# for j in range(0,i):
#     sm=score_matrix5[j]
#     popo=sm.pop(j)
#     score_matrixx5=score_matrixx5+sm

f6 = open("E:/a/python27/bnfinder2/doc/data/1565data/dream10/10ts_1_gold", "r")
test_line1=f6.readlines()
# test_line2=test_line1

#test_line5=test_line1

aaaano=0
bbbbno=0
for item in test_line1:
    gene_data = item.split("\t")
    aaaa = gene_data[0]
    bbbb = gene_data[1]
    nnnn=gene_data[2].strip('\n')
    nnnn = eval(nnnn)
    for i2 in gene_no:
        if aaaa==i2[1]:
            aaaano=i2[0]
            continue
        if bbbb==i2[1]:
            bbbbno=i2[0]
            continue
        # if aaano and bbbno:
        #     break
    test_matrix1[aaaano][bbbbno]=nnnn
test_matrixx1=[]
for j in range(0,i):
    tm=test_matrix1[j]
    popo=tm.pop(j)
    test_matrixx1=test_matrixx1+tm

# aaaano=0
# bbbbno=0
# for item in test_line2:
#     gene_data = item.split("\t")
#     aaaa = gene_data[0]
#     bbbb = gene_data[1]
#     nnnn=gene_data[2].strip('\n')
#     nnnn = eval(nnnn)
#     for i2 in gene_no:
#         if aaaa==i2[1]:
#             aaaano=i2[0]
#             continue
#         if bbbb==i2[1]:
#             bbbbno=i2[0]
#             continue
#         # if aaano and bbbno:
#         #     break
#     test_matrix2[aaaano][bbbbno]=nnnn
# test_matrixx2=[]
# for j in range(0,i):
#     tm=test_matrix2[j]
#     popo=tm.pop(j)
#     test_matrixx2=test_matrixx2+tm



# for item in test_line5:
#     gene_data = item.split("\t")
#     aaaa = gene_data[0]
#     bbbb = gene_data[1]
#     nnnn=gene_data[2].strip('\n')
#     nnnn = eval(nnnn)
#     for i2 in gene_no:
#         if aaaa==i2[1]:
#             aaaano=i2[0]
#             continue
#         if bbbb==i2[1]:
#             bbbbno=i2[0]
#             continue
#         # if aaano and bbbno:
#         #     break
#     test_matrix5[aaaano][bbbbno]=nnnn
# test_matrixx5=[]
# for j in range(0,i):
#     tm=test_matrix5[j]
#     popo=tm.pop(j)
#     test_matrixx5=test_matrixx5+tm
#
# aa=1




y_test1=test_matrixx1
y_score1=score_matrixx1
# Compute ROC curve and ROC area for each class
fpr1,tpr1,threshold1 = roc_curve(y_test1, y_score1) ###计算真正率和假正率
roc_auc1 = auc(fpr1,tpr1)

# y_test2=test_matrixx2
# y_score2=score_matrixx2
# Compute ROC curve and ROC area for each class
# fpr2,tpr2,threshold2 = roc_curve(y_test2, y_score2) ###计算真正率和假正率
# roc_auc2 = auc(fpr2,tpr2)


# y_test5=test_matrixx5
# y_score5=score_matrixx5
# # Compute ROC curve and ROC area for each class
# fpr5,tpr5,threshold5 = roc_curve(y_test5, y_score5) ###计算真正率和假正率
# roc_auc5 = auc(fpr5,tpr5)
# plt.figure()
# xnew =np.arange(0,10,0.1)
# fpr1r=np.sort(fpr1)
# func = interpolate.interp1d(fpr1,tpr1,kind='cubic')
# ynew = func(xnew)
lw = 2.5
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
plt.figure(figsize=(10,7.5))
plt.tick_params(labelsize=19)
ax = plt.gca() # gca = 'get current axis' 获取当前坐标
ax.spines['bottom'].set_linewidth(1.5)
ax.spines['left'].set_linewidth(1.5)
ax.spines['top'].set_linewidth(1.5)
ax.spines['right'].set_linewidth(1.5)
plt.plot(fpr1, tpr1, color='#228B22',
         lw=lw, label='ROC 100ts_5(AUC = %0.4f)' % roc_auc1) ###假正率为横坐标，真正率为纵坐标做曲线
# plt.plot(fpr2, tpr2, color='maroon',
#          lw=lw, linestyle=':',label='ROC BN-SP(20%%)(AUC = %0.4f)' % roc_auc2) ###假正率为横坐标，真正率为纵坐标做曲线
# plt.plot(xnew,ynew,color='navy', lw=lw, linestyle='--')
# plt.plot(fpr5, tpr5, color='#B22222',
#          lw=lw, linestyle='--',label='ROC BN(AUC = %0.4f)' % roc_auc5) ###假正率为横坐标，真正率为纵坐标做曲线
# x_new = np.linspace(min(fpr1), max(fpr1), 50)
#
# y_smooth = spline(fpr1, tpr1, x_new)，
#
# plt.plot(x_new, y_smooth)
# plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
font2 = {'family': 'Times New Roman','weight':'normal','size': 28}
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.0])
plt.xlabel('False Positive Rate',font2)
plt.ylabel('True Positive Rate',font2)
# plt.title('Receiver operating characteristic',font2)
plt.legend(loc="lower right",prop={'family': 'Times New Roman','size': 22})
plt.show()


