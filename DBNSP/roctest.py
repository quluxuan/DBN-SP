# #coding: utf-8
# import numpy as np
#
# import matplotlib.pyplot as plt
#
# from sklearn import svm, datasets
#
# from sklearn.metrics import roc_curve, auc
#
# from sklearn.model_selection import StratifiedKFold
#
# iris = datasets.load_iris()
#
# X = iris.data
#
# y = iris.target
#
# X, y = X[y != 2], y[y != 2]  # 去掉了label为2，label只能二分，才可以。
#
# n_samples, n_features = X.shape
#
# # 增加噪声特征
#
# random_state = np.random.RandomState(0)
#
# X = np.c_[X, random_state.randn(n_samples, 200 * n_features)]
#
# cv = StratifiedKFold(n_splits=6)  # 导入该模型，后面将数据划分6份
#
# classifier = svm.SVC(kernel='linear', probability=True, random_state=random_state)  # SVC模型 可以换作AdaBoost模型试试
#
# # 画平均ROC曲线的两个参数
#
# mean_tpr = 0.0  # 用来记录画平均ROC曲线的信息
#
# mean_fpr = np.linspace(0, 1, 100)
#
# cnt = 0
#
# for i, (train, test) in enumerate(cv.split(X, y)):  # 利用模型划分数据集和目标变量 为一一对应的下标
#
#     cnt += 1
#
#     probas_ = classifier.fit(X[train], y[train]).predict_proba(X[test])  # 训练模型后预测每条样本得到两种结果的概率
#
#     fpr, tpr, thresholds = roc_curve(y[test], probas_[:, 1])  # 该函数得到伪正例、真正例、阈值，这里只使用前两个
#
#     mean_tpr += np.interp(mean_fpr, fpr, tpr)  # 插值函数 interp(x坐标,每次x增加距离,y坐标)  累计每次循环的总值后面求平均值
#
#     mean_tpr[0] = 0.0  # 将第一个真正例=0 以0为起点
#
#     roc_auc = auc(fpr, tpr)  # 求auc面积
#
#     plt.plot(fpr, tpr, lw=1, label='ROC fold {0:.2f} (area = {1:.2f})'.format(i, roc_auc))  # 画出当前分割数据的ROC曲线
#
# plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Luck')  # 画对角线
#
# mean_tpr /= cnt  # 求数组的平均值
#
# mean_tpr[-1] = 1.0  # 坐标最后一个点为（1,1）  以1为终点
#
# mean_auc = auc(mean_fpr, mean_tpr)
#
# plt.plot(mean_fpr, mean_tpr, 'k--', label='Mean ROC (area = {0:.2f})'.format(mean_auc), lw=2)
#
# plt.xlim([-0.05, 1.05])  # 设置x、y轴的上下限，设置宽一点，以免和边缘重合，可以更好的观察图像的整体
#
# plt.ylim([-0.05, 1.05])
#
# plt.xlabel('False Positive Rate')
#
# plt.ylabel('True Positive Rate')  # 可以使用中文，但需要导入一些库即字体
#
# plt.title('Receiver operating characteristic example')
#
# plt.legend(loc="lower right")
#
# plt.show()
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

#设置距离
x =np.array([0, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9,10])

#设置相似度
y =np.array([0.8579087793827057, 0.8079087793827057, 0.7679087793827057, 0.679087793827057,
    0.5579087793827057, 0.4579087793827057, 0.3079087793827057, 0.3009087793827057,
    0.2579087793827057, 0.2009087793827057, 0.1999087793827057, 0.1579087793827057,
    0.0099087793827057, 0.0079087793827057, 0.0069087793827057, 0.0019087793827057,
    0.0000087793827057])

#插值法之后的x轴值，表示从0到10间距为0.5的200个数
xnew =np.arange(0,10,0.1)

#实现函数
func = interpolate.interp1d(x,y,kind='cubic')

#利用xnew和func函数生成ynew,xnew数量等于ynew数量
ynew = func(xnew)

# 原始折线
plt.plot(x, y, "r", linewidth=1)

#平滑处理后曲线
plt.plot(xnew,ynew)
#设置x,y轴代表意思
plt.xlabel("The distance between POI  and user(km)")
plt.ylabel("probability")
#设置标题
plt.title("The content similarity of different distance")
#设置x,y轴的坐标范围
plt.xlim(0,10,8)
plt.ylim(0,1)

plt.show()
