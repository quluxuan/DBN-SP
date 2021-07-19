# DBN-SP
DBN-SP algorithm is an algorithm for learning the structure of gene regulatory network, which can be used in dynamic Bayesian and Bayesian networks, and is written in python (version > 2.4 is required, python3 is not supported yet, but we plan to finish it soon).

The main functions and functions are as follows:
BDE.py	BDE scoring function
BIC.py	BIC scoring function
bnf.py	Overall running code
com_mi.py	A complex mutual information calculation code
compute_relation.py	Three methods to calculate the code of candidate parent node set
conditional_mi.py	Calculation code of conditional mutual information
continuous.py	Processing part of continuous data
data.py	Data processing, including screening of candidate parent node sets and structure learning
network_roc.py	Draw ROC curve
sif_relation.py	Convert the output network result. sif file into adjacency matrix
reducecmi.py	Optimize the network structure and delete redundant edges

The main parameters are as follows:
score	Select the scoring function BDe or BIC (default is BDe)
limit	Size of parent node set (default 5)
sloops	Whether to allow self-circulation in dynamic Bayesian networks (false by default)
distrib	Number of processes used for multiprocessing (default is 0)
ratio	Filtering ratio when using three methods to filter candidate parent nodes (default is 10%)
lamda	Threshold of CMI calculation during structural optimization (default 0.03)
order0		Calculate the order of CMI during structural optimization (default 2)
expr	Gene expression data file
txt	Output file with suboptimal parent node set
net	Network structure output file

Demonstration
We demonstrate the DBN-SP algorithm using the data set of 10 genes in DREAM3 as an example.
First, in the bnf.py file, input the path of 10 gene expression data files, and specify the path of the output file, set each parameter, and each parameter can be modified according to the needs.
......
options.expr = ‘E:/a/python27/doc/data/DREAM3/data/10genes_data.txt’
options.txt = ‘E:/a/python27/doc/data/DREAM3/data/data10_10%.txt’ 
options.net = ‘E:/a/python27/doc/data/DREAM3/data/data10_10%.sif’ 
options.limit = 5
options.distrib = 20
options.score = BDE
options.sloops = True
......
Then, open the data.py file, find the function to find potential parent nodes, and set the screening ratio ratio when screening candidate parent nodes by three methods. The screening ratio can be modified according to needs. By default, we selected the screening ratio of 10%.
……
def get_potential_parents_matrix(self, node, parents_select,parents_matrix,sloops=True):
gene_num = self.n_gen
ratio = 0.1
pearson_matrix = parents_matrix[0]
……
Then, run bnf.py, and you can get the results of DBN structure learning.
Finally, open reducecmi.py to optimize the structure. First, set the threshold and order of CMI calculation, input gene expression data and structure learning results, and run reducecmi.py to get the experimental results of structure optimization.
......
lamda = 0.1
order0 = 2
......
f = open(“E:/a/python27/doc/data/DREAM3/data/10genes_data.txt”, “r”)
f1 = open(“E:/a/python27/doc/data/DREAM3/data/ data10_10%.sif”, “r”)
......

