#coding: utf-8
# vim: set fileencoding=utf-8 :
import graph,stats,math,fpconst
from itertools import chain
from BDE import BDE
from com_mi import com_mi
import continuous
import util
from math import log
from random import shuffle
from heapq import *
from minepy import MINE
import numpy as np
import datetime
from collections import Counter
from compute_relation import compute_relation
from ROTEST import *
from conditional_mi import cmi
class experiment:
    """Class representing an experiment, possibly containing more than one time-point
    """
    def __init__(self,name,data,perturbed,points=[],static=False):
        """Constructor

        retains the argument data as self._data, while filling the self.data with double vectors
        (duplicates for static networks, consecutive pairs for dynamic)
        """
        self._data=data
        self.perturbed=perturbed
        self.data=[]
        self.name=name
        self.points=points
        if static:
            for i,d in enumerate(self._data):
                self.data.append(self._data[i]+d)
        else:
            # in format (data from time t - 1, data from time t)
            for i,d in enumerate(self._data[1:]):
                self.data.append(self._data[i]+d)


class gene:
    """Class corresponding to a single variable

    #TODO: rename to variable or smth... gene is inappropriate
    """

    def __init__(self, name, n_disc, index):
        self.name = name
        self.n_disc = n_disc
        self.values = []
        self.index = index
        self.cpd = None

    def __str__(self):
        return self.name

    def disc_vals(self):
        """Returns discretized values of a variable

        works also for continuous variables
        """
        if self.values:
            try:
                return map(int, self.values)
            except ValueError:
                return self.values
        else:
            return [0,1]  # This should take care of continuous variables... may need a change if we introduce more complex mixtures.

    def __index__(self):
        """returns the index of a variable.
        """
        return self.index

    def __add__(self, other):
        """Implements addition on indexes, so that we can write table[var+1] instad of table[var.index+1]

        """
        return other + self.index

    def __hash__(self):
        """returns the hash value for comparisons
        """
        return self.index

    def __eq__(self, other):
        if type(other) == int:
            return self.index == other
        else:
            return self.index == other.index

    def base_weight(self):
        return max(1.5, self.n_disc)



class dataset:
    """
    A class representing a complete dataset.

    Implements input, output, data selection and also functions for learning networks given a scoring function
    """
    def __init__(self,name='no_name'):
        self.name=name
        self.n_series=0
        self.n_gen=0
        self.points=[]
        self.vertices= []
        self.vertice_names=[]
        self.static=False
        self.regulators=[]
        self.parents={}
        self.prior={}

    def fromNewFile(self, file, float_transform=continuous.BiNormMeans()):
        pert = []
        disc = {}
        cpd = {'and': [], 'or': []}
        contin = []
        self.default = None
        rec = []
        regul = []
        while not rec:
            ln = file.readline()
            rec = ln.strip().split()
        #################
        # THE PREAMBLE
        ################
        while rec[0][0] == "#":  # preamble
            if rec[0] == "#perturbed":  # experiment serie name followed by list of its perturbed vertices
                pert.append((rec[1], rec[2:]))
            elif rec[0] == "#regulators":  # list of potential regulators of all vertices (except specified in #parents command) not specified in previous or present #regulators command
                reg = filter(lambda x: True not in map(lambda r: x in r, regul), rec[1:])
                for i, r in enumerate(reg):
                    if i != reg.index(r):
                        reg[i] = None
                reg = filter(None, reg)
                regul.append(reg)
            elif rec[0] == "#parents":  # vertice name followed by list of its potential regulators
                reg = rec[2:]
                for i in range(len(reg)):
                    if i != reg.index(reg[i]):
                        reg[i] = None
                reg = filter(None, reg)
                self.parents[rec[1]] = reg
            elif rec[0]=="#discrete": #vertice name followed by list of all its possible values
                disc[rec[1]]=rec[2:]
            elif rec[0]=="#continuous": #list of vertices with continuous values
                contin.extend(rec[1:])
            elif rec[0] == "#default":  # list of default possible values or 'FLOAT' for default continuous values
                if len(rec) == 2 and rec[1] == "FLOAT":
                    self.default = []
                else:
                    self.default = rec[1:]
            elif rec[0] == "#priorvert":  # prior weight of all outgoing edges (except specified in #prioredge command) followed by list of related vertices
                val = float(rec[1])
                if val <= 0:
                    raise Exception("Non-positive prior weight in dataset")
                for name in rec[2:]:
                    self.prior[name] = val
            elif rec[0] == "#prioredge":  # prior weight of edges preceded by edge target (common for all specified edges) and followed by list of edge sources
                target_name = rec[1]
                val = float(rec[2])
                if val <= 0:
                    raise Exception("Non-positive prior weight in dataset")
                for name in rec[3:]:
                    self.prior[(target_name, name)] = val
            else:
                pass  # ignoring other stuff
            rec = []
            while not rec:
                ln = file.readline()
                rec = ln.strip().split()
        self.regulators = regul
        ################
        # condition names
        ################
        if rec[0] != 'conditions':
            self.name = rec[0]
            # todo:pbd serie not used any more?
        if len(rec[1].split(":")) == 1:
            self.static = True
        if self.static:
            conds = map(lambda x: [x], rec[1:])
            names = conds
        else:
            # todo:pbd serie not used any more?
            serie = rec[1].split(":")[0]
            # pbd: explain why conds is a singleton and not just a single standalone element
            cond = [rec[1].split(":")[1]]
            conds = [cond]
            names = [serie]
            for c in rec[2:]:
                r = c.split(":")
                if r[0] == serie:
                    cond.append(r[1])
                else:
                    serie = r[0]
                    names.append(serie)
                    cond = [r[1]]
                    conds.append(cond)
            # pbd: after that conds looks like [['1', '2', '3', ...]]
        ln = file.readline()
        ###################
        # EXPRESSION DATA #
        ###################
        self.vertices = []
        self.vertice_names = []
        pline = []  ###======================YE=================##
        # prepare the lists for expression data
        #        n_cond=sum(map(len,conds))
        exp_data = map(lambda l: [[] for el in l], conds)
        geneindex = 0
        while ln:  # expression data
            rec = ln.strip().split()
            if not rec:
                ln = file.readline()
                continue
            name = rec[0]
            vals = []
            # are the values discrete or not?
            if name in contin:  # did we say it's continuous?
                discrete = False
            elif name in disc.keys():  # did we say is's discrete?
                discrete = True
                vals = disc[name]
            elif self.default != None:  # did we specify a default
                if self.default == []:
                    discrete = False
                else:
                    discrete = True
                    vals = self.default
            else:  # We did not specify the type of values for this gene, try to guess it
                try:
                    i = int(rec[1])
                except ValueError:  # it's not an integer
                    try:
                        i = float(rec[1])
                    except ValueError:  # it's not a float -> it's an alphanumerical value
                        discrete = True
                    else:
                        discrete = False  # it's a float
                else:
                    discrete = True  # it's an int, let's map the strings to ints
            if discrete:
                line = map(int, rec[1:])
                pline.append(line)  ###======================YE=================##
                if vals == []:  # we don't know what are the possible values
                    for v in line:
                        if v not in vals:
                            vals.append(v)
                    # let's sort the values
                    vals.sort()
            else:  # not discrete
                line = map(float, rec[1:])
                pline.append(line)  ###======================YE=================##
                float_params = []
                if float_transform:  # use float_transform if specified
                    float_params = float_transform.estimate(line)
                    line = float_transform.transform(line)
            self.vertices.append(gene(name, len(vals), geneindex))
            self.vertices[-1].values = vals
            if discrete == False:
                self.vertices[-1].floatParams = float_params
            else:
                self.vertices[-1].floatParams = None
            self.vertice_names.append(name)
            # append this line to the expression data
            for el, l in zip(line, chain(*exp_data)):
                if discrete:
                    l.append(vals.index(el))
                else:
                    l.append(el)
            ln = file.readline()
            geneindex += 1
        self.regulators = pline  ###======================YE=================##
        ################
        # POSTPROCESSING#
        ################
        # process regulatory constraints
        # process the lists of perturbed genes and make the list of experiments
        for i, (c, n) in enumerate(zip(conds, names)):
            pg = []
            n="".join(n)
            for x, p in pert:
                if x == n:
                    pg.extend(
                        map(lambda n: self.vertices[self.vertice_names.index(n)], p))  # map gene names to indexes
            self.points.append(experiment(n, exp_data[i], pg, c, self.static))
        # set the number of genes
        self.n_gen = geneindex
        return self

    def learn(self, score, fpr=None,fpr_nodes=False, prior=None,distrib=False,verbose=None, n_min=1, limit=None, min_empty=None, min_optim=None):
        scr = score
        parents_matrix = compute_relation(self)
        if distrib:
            try:
                from multiprocessing import Pool#,Queue
            except ImportError:
                print "Problem invoking multiprocessing module. Running on a single CPU"
                distrib=False
                #from Queue import Queue
            else:
                pool=Pool(distrib)
        #result = map(learn_x, [(x, self, min_empty, min_optim, scr, verbose, n_min,limit) for x in self.vertices])
        if fpr:
            if fpr_nodes:
                fpr_factor = fpr
            else:
                normalizer = 0.
                n_edge = 0
                for vertex in self.vertices:
                    for parent in self.get_potential_parents1(vertex,self.regulators, scr.sloops):
                        normalizer += 1. / self.get_prior(vertex, parent)
                        n_edge += 1
                try:
                    fpr_factor = fpr * n_edge / normalizer
                    #                fpr_factor = fpr*self.n_gen/normalizer
                    #                fpr_factor = fpr/normalizer
                except ZeroDivisionError:
                    raise Exception("Potential parents set is empty for each vertex")
        else:
            fpr_factor = None
        if distrib:
            #            from itertools import izip,repeat
             result = pool.map(learn_x, [(x, self, min_empty, min_optim, scr,verbose, n_min, limit, fpr_factor,fpr_nodes,parents_matrix) for x in self.vertices])
            # result = pool.map(learn_x, [
            #     (x, self, min_empty, min_optim, scr, verbose, n_min, limit, fpr_factor, fpr_nodes) for x
            #     in self.vertices])
            # eval_func,izip([learn_x]*len(self.vertices),self.vertices))
            # pool.close() #raises weird errors...., not necessary in normal bnf operation, since the process exits after learning
        else:
             result = map(learn_x, [(x, self, min_empty, min_optim, scr, verbose, n_min,limit,fpr_factor,fpr_nodes,parents_matrix) for x in self.vertices])
            # result = map(learn_x, [
            #     (x, self, min_empty, min_optim, scr, verbose, n_min, limit, fpr_factor, fpr_nodes) for x
            #     in self.vertices])
        # print result
        self.vertices = [r[0] for r in result]
        # TODO do we need this line?^
        total_score = 0.0
        pars = {}
        subpars = {}
        for x, par, minlist, sc in result:
            pars[x] = par
            subpars[x] = minlist
            total_score += sc
        # move parents from a dictionary to a list
        par_list = []
        for v in self.vertices:
            try:
                par_list.append(pars[v.__index__()])
            except KeyError:
                print "ERRR", par_list, v.name, id(v), v.__index__(), [map(lambda x: (x.name, id(x)), l) for l in pars.values()]
        pars = par_list
        g = graph.graph()
        g.fromParents(self.vertices, pars)
        # mutuals=[]
        # alpha = 0.05  ###
        # lamda = 1
        # for regulate_gene in g.edges:
        #     targets = g.edges[regulate_gene]
        #     for target_gene in targets:
        #         ######========= mutuals用来存放具有相互关系的基因对，保证只需检测一次
        #         notmutual=0
        #         if regulate_gene.name < target_gene.name:
        #             mutual = regulate_gene.name+target_gene.name
        #         else:
        #             mutual = target_gene.name+regulate_gene.name
        #         if mutual in mutuals:
        #             continue
        #         else:
        #             mutuals.append(mutual)
        #         ######========= mutuals用来存放具有相互关系的基因对，保证只需检测一次
        #         if regulate_gene in g.edges[target_gene]:
        #             for othergene in g.edges:
        #                 if othergene == regulate_gene or othergene ==  target_gene:
        #                     continue
        #                 if regulate_gene in g.edges[othergene] or target_gene in g.edges[othergene]:
        #                     notmutual=1
        #             if len(g.edges[target_gene]) == 1 and len(g.edges[regulate_gene]) == 1 and notmutual == 0: ######=========独立的两个基因
        #                 y1 = self.regulators[target_gene.index]
        #                 y2 = self.regulators[regulate_gene.index]
        #                 X1 = self.regulators[regulate_gene.index]
        #                 X2 = self.regulators[target_gene.index]
        #                 parro1 = reoptim(y1, X1, lamda, alpha) ######========= [X1](regulate_gene)===>>y1(target_gene)
        #                 parro2 = reoptim(y2, X2, lamda, alpha)######========= [X2](target_gene)===>>y2(regulate_gene)
        #                 # cmi1 = cmi(y1, self.regulators[regulate_gene.index], XX1)
        #                 # cmi2 = cmi(y2, self.regulators[target_gene.index], XX2)
        #                 print regulate_gene,target_gene,parro1[0], parro2[0]
        #                 if abs(parro1[0]) < abs(parro2[0]):
        #                     del g.edges[regulate_gene][g.edges[regulate_gene].index(target_gene)]
        #                 else:
        #                     del g.edges[target_gene][g.edges[target_gene].index(regulate_gene)]
        #             else:######=========不独立的两个基因
        #                 X1 = []
        #                 # XX1 = []
        #                 X2 = []
        #                 # XX2 = []
        #                 i1 = 0
        #                 i2 = 0
        #                 y1 = self.regulators[target_gene.index]
        #                 y2 = self.regulators[regulate_gene.index]
        #                 for mugene1 in g.edges:
        #                     if target_gene in g.edges[mugene1]:
        #                         if mugene1 == regulate_gene:
        #                             i1i = i1
        #                         X1.append(self.regulators[mugene1.index])
        #                         i1 = i1 + 1
        #                         # if mugene1 == regulate_gene:
        #                         #     continue
        #                         # XX1.append(self.regulators[mugene1.index])
        #                 for mugene2 in g.edges:
        #                     if regulate_gene in g.edges[mugene2]:
        #                         if mugene2 == target_gene:
        #                             i2i = i2
        #                         X2.append(self.regulators[mugene2.index])
        #                         i2 = i2 + 1
        #                         # if mugene2 == target_gene:
        #                         #     continue
        #                         # XX2.append(self.regulators[mugene2.index])
        #                 cmiX1 = []
        #                 cmiX2 = []
        #                 for target1_n_gene in g.edges[target_gene]:
        #                     if target1_n_gene == regulate_gene:
        #                         continue
        #                     cmiX1.append(self.regulators[target1_n_gene.index])
        #                 for ru in g.edges:
        #                    if ru == regulate_gene:
        #                        continue
        #                    if target_gene in g.edges[ru]:
        #                        if ru in g.edges[target_gene]:
        #                            continue
        #                        cmiX1.append(self.regulators[ru.index])
        #                 for regulate1_n_gene in g.edges[regulate_gene]:
        #                     if regulate1_n_gene == target_gene:
        #                         continue
        #                     cmiX2.append(self.regulators[regulate1_n_gene.index])
        #                 for tg in g.edges:
        #                     if tg == target_gene:
        #                         continue
        #                     if regulate_gene in g.edges[tg]:
        #                         if tg in g.edges[regulate_gene]:
        #                             continue
        #                         cmiX2.append(self.regulators[tg.index])
        #                 parro1 = reoptim(y1, X1, lamda, alpha)######========= [X1](regulate_gene)===>>y1(target_gene)
        #                 parro2 = reoptim(y2, X2, lamda, alpha)######========= [X2](target_gene)===>>y2(regulate_gene)
        #                 p1 = parro1[i1i]
        #                 p2 = parro2[i2i]
        #
        #                 ########compute cmi
        #                 cmi1 = cmi(y1, self.regulators[regulate_gene.index], cmiX1)
        #                 cmi2 = cmi(y2, self.regulators[target_gene.index], cmiX2)
        #                 # if abs(parro1[i1i]) < abs(parro2[i2i]) and abs(parro1[i1i]-parro2[i2i])>0.3:
        #                 #     del g.edges[regulate_gene][g.edges[regulate_gene].index(target_gene)]
        #                 # elif abs(parro1[i1i]) > abs(parro2[i2i]) and abs(parro1[i1i]-parro2[i2i])>0.3:
        #                 #     del g.edges[target_gene][g.edges[target_gene].index(regulate_gene)]
        #                 print regulate_gene,target_gene,p1,p2,cmi1,cmi2
        #             awert=1
        #         awert = 1
        #     awert = 1

        for v, par in zip(self.vertices, pars):
            g.vertice_labelling[v] = v.name
            v_signs = self.select_1(v, par).signs()
            for p in par:
                g.edge_labelling[p, v] = v_signs[p]
        return total_score, g, subpars

    def select_1(self, vertex, parents=None, sloops=True):
        """
        given a vertice and a list of his parents,
        select parents values at time t
        and child value at time t+1,
        omitting the experiments where the child was perturbed.
        """
        if parents is None:
            parents = self.get_potential_parents(vertex, sloops)
        res = []
        for expment in self.points:
            if vertex not in expment.perturbed:
                for sample in expment.data:
                    res.append([sample[par] for par in parents] + [sample[vertex + self.n_gen]])
        return dataset_1(map(list, zip(*res)), vertex, parents)
    def select_2(self,vertex,parents_select, parents_matrix, sloops=True):
    #def select_2(self, vertex, parents_select, sloops=True):
        # s_gene=['NFKB1','NFE2L2','PIK3CA','PTEN','STAT3','IL6','TGFB1','MAPK1','TLR4','AKT1','MAPK14','RELA','TNF','MTOR']
        # if vertex.name in s_gene:
        #     parents=[]
        # else:
        #     parents = self.get_potential_parents1(vertex,parents_select, sloops)
        parents = self.get_potential_parents_matrix(vertex, parents_select, parents_matrix, sloops)
        #parents = self.get_potential_parents_single(vertex, parents_select, sloops)
        res = []
        for expment in self.points:
            if vertex not in expment.perturbed:
                for sample in expment.data:
                    res.append([sample[par] for par in parents]+[sample[vertex+self.n_gen]])
        return dataset_1(map(list,zip(*res)),vertex,parents)
    def get_potential_parents(self, node, sloops=True):
        if self.parents.has_key(node.name):  # potential parents of a node are explicitly specified
            parents = map(lambda n: self.vertices[self.vertice_names.index(n)], self.parents[node.name])
        elif self.regulators:  # regulators are specified
            parents = []
            for r in self.regulators:
                if node.name in r:
                    break
                else:
                    parents.extend(map(lambda n: self.vertices[self.vertice_names.index(n)], r))
        else:  # all vertices are potential parents
            parents = [v for v in self.vertices]
        if not sloops:
            try:
                i = parents.index(node)
                del parents[i]
            except ValueError:
                pass
        return parents

    def get_potential_parents_matrix(self, node, parents_select, parents_matrix, sloops=True):
        gene_num = self.n_gen
        ratio = 0.1
        pearson_matrix = parents_matrix[0]
        mi_matrix = parents_matrix[1]
        mic_matrix = parents_matrix[2]
        if self.parents.has_key(node.name): # potential parents of a node are explicitly specified
            parents=map(lambda n: self.vertices[self.vertice_names.index(n)],self.parents[node.name])
        else:
            parents = []
            parents_potential = []
            ############==========================pearsonr=====================###################
            parents_num = round(gene_num * ratio)
            parents_pearson = []
            for i in range(len(parents_select)):
                abscorrel=pearson_matrix[node.index][i]
                if parents_num > 0:
                    heappush(parents_potential, (abscorrel, i))
                    parents_num -= 1
                elif parents_potential:
                    if parents_potential[0][0]<abscorrel:
                        heapreplace(parents_potential, (abscorrel, i))
            for (coe,no) in parents_potential:
                parents_pearson.append(self.vertices[no])
            # print "gene:", node.name, "parents_pearson===",
            # for gene123 in parents_pearson:
            #     print gene123.name,
            print ""
            ############==========================pearsonr=====================###################
            # ############==========================mi=====================###################
            parents_potential = []
            parents_mi = []
            parents_num = round(gene_num * ratio)
            for i in range(len(parents_select)):
                mi = mi_matrix[node.index][i]
                if parents_num > 0:
                    heappush(parents_potential, (mi, i))
                    parents_num -= 1
                elif parents_potential:
                    if parents_potential[0][0] < mi:
                        heapreplace(parents_potential, (mi, i))
            for (coe, no) in parents_potential:
                parents_mi.append(self.vertices[no])
            # print "gene:", node.name, "parents_mi==="
            # for gene123 in parents_mi:
            #     print gene123.name,
            # ############==========================mi=====================###################
            ############==========================MIC=====================###################
            parents_potential = []
            parents_MIC = []
            parents_num = round(gene_num * ratio)
            for i in range(len(parents_select)):
                mic = mic_matrix[node.index][i]
                if parents_num > 0:
                    heappush(parents_potential, (mic, i))
                    parents_num -= 1
                elif parents_potential:
                    if parents_potential[0][0] < mic:
                        heapreplace(parents_potential, (mic, i))
            for (coe, no) in parents_potential:
                parents_MIC.append(self.vertices[no])
            # print "gene:", node.name, "parents_MIC==="
            # for gene123 in parents_MIC:
            #     print gene123.name,
        pars = parents_pearson + parents_mi + parents_MIC
        parents = list(set(pars))
        if not sloops:
            try:
                i = parents.index(node)
                del parents[i]
            except ValueError:
                pass
        print node.name, parents.__len__()
        return parents
            ############==========================MIC=====================###################

    def get_potential_parents3(self, node, parents_select, sloops=True):
        """selects potential parents for node"""
        gene_num=self.n_gen
        ratio=0.4
        # from multiprocessing import Pool
        # distrib2 = 4
        # pool2 = Pool(distrib2)
        if self.parents.has_key(node.name): # potential parents of a node are explicitly specified
            parents=map(lambda n: self.vertices[self.vertice_names.index(n)],self.parents[node.name])
        else: # all vertices are potential parents
            parents = []
            parents_potential = []
            ############==========================pearsonr=====================###################
            parents_num=round(gene_num*ratio)
            parents_pearson = []
            # parents_potential = compute_pearson(node, parents_num, parents_select)

            for i in range(len(parents_select)):
                if parents_select[node.index] == parents_select[i]:
                    continue
                correl=stats.pearsonr(parents_select[node.index],parents_select[i])
                abscorrel=abs(correl[0])
                if parents_num > 0:
                    heappush(parents_potential, (abscorrel, i))
                    parents_num -= 1
                elif parents_potential:
                    if parents_potential[0][0]<abscorrel:
                        heapreplace(parents_potential, (abscorrel, i))
            for (coe,no) in parents_potential:
                parents_pearson.append(self.vertices[no])
            print "gene:", node.name, "parents_pearson===",
            for gene123 in parents_pearson:
                print gene123.name,
            print ""
            ############==========================pearsonr=====================###################
            # ############==========================mi=====================###################
            parents_potential = []
            parents_mi = []
            parents_num = round(gene_num * ratio)
            for i in range(len(parents_select)):
                if parents_select[node.index] == parents_select[i]:
                    continue
                mi = com_mi(parents_select[node.index], parents_select[i])
                if parents_num > 0:
                    heappush(parents_potential, (mi, i))
                    parents_num -= 1
                elif parents_potential:
                    if parents_potential[0][0] < mi:
                        heapreplace(parents_potential, (mi, i))
            for (coe, no) in parents_potential:
                parents_mi.append(self.vertices[no])
            print "gene:", node.name, "parents_mi==="
            for gene123 in parents_mi:
                print gene123.name,
            # ############==========================mi=====================###################
            ############==========================MIC=====================###################
            parents_potential = []
            parents_MIC = []
            parents_num = round(gene_num * ratio)
            for i in range(len(parents_select)):
                if parents_select[node.index] == parents_select[i]:
                    continue
                mine = MINE()
                mine.compute_score(parents_select[node.index], parents_select[i])
                mi = mine.mic()
                if parents_num > 0:
                    heappush(parents_potential, (mi, i))
                    parents_num -= 1
                elif parents_potential:
                    if parents_potential[0][0] < mi:
                        heapreplace(parents_potential, (mi, i))
            for (coe, no) in parents_potential:
                parents_MIC.append(self.vertices[no])
            print "gene:", node.name, "parents_MIC==="
            for gene123 in parents_MIC:
                print gene123.name,

            ############==========================MIC=====================###################
            ############==========================spearman=====================###################
            # parents_potential = []
            # parents_spearman = []
            # parents_num = round(gene_num * ratio)
            # for i in range(len(parents_select)):
            #     if parents_select[node.index] == parents_select[i]:
            #         continue
            #     correl=stats.spearmanr(parents_select[node.index], parents_select[i])
            #     abscorrel = abs(correl[0])
            #     if parents_num:
            #         heappush(parents_potential, (abscorrel, i))
            #         parents_num -= 1
            #     elif parents_potential > 0:
            #         if parents_potential[0][0] < abscorrel:
            #             heapreplace(parents_potential, (abscorrel, i))
            # for (coe, no) in parents_potential:
            #     parents_spearman.append(self.vertices[no])
            ############==========================spearman=====================###################
        # if not sloops:
        #     try:
        #         i=parents.index(node)
        #         del parents[i]
        #     except ValueError:
        #         pass
        pars=parents_pearson+parents_mi+parents_MIC
        parents = list(set(pars))

        www=90
        print node.name,parents.__len__()
        return parents
    def get_potential_parents1(self, node, parents_select, sloops=True):
        """selects potential parents for node"""
        gene_num=self.n_gen
        ratio=0.4
        if self.parents.has_key(node.name): # potential parents of a node are explicitly specified
            parents=map(lambda n: self.vertices[self.vertice_names.index(n)],self.parents[node.name])
        else: # all vertices are potential parents
            parents = []
            parents_potential = []
            ############==========================pearsonr=====================###################
            parents_num=round(gene_num*ratio)
            parents_pearson = []
            parents_pearson_aft = []
            for i in range(len(parents_select)):
                if parents_select[node.index] == parents_select[i]:
                    continue
                correl=stats.pearsonr(parents_select[node.index],parents_select[i])
                abscorrel=abs(correl[0])
                if parents_num > 0:
                    heappush(parents_potential, (abscorrel, i))
                    parents_num -= 1
                elif parents_potential:
                    if parents_potential[0][0]<abscorrel:
                        heapreplace(parents_potential, (abscorrel, i))
            parents_potential.sort(reverse=True)
            mid_ratio=round(0.5*gene_num*ratio)
            for (coe,no) in parents_potential:
                if mid_ratio:
                    parents_pearson.append(self.vertices[no])
                    mid_ratio=mid_ratio-1
                else:
                    parents_pearson_aft.append(self.vertices[no])
            print "gene:", node.name, "parents_pearson==="
            for gene123 in parents_pearson:
                print gene123.name
            print "gene:", node.name, "parents_pearson_aft==="
            for gene123 in parents_pearson_aft:
                print gene123.name
            ############==========================pearsonr=====================###################
            # ############==========================mi=====================###################
            parents_potential = []
            parents_mi = []
            parents_mi_aft = []
            parents_num = round(gene_num * ratio)
            for i in range(len(parents_select)):
                if parents_select[node.index] == parents_select[i]:
                    continue
                mi = com_mi(parents_select[node.index], parents_select[i])
                if parents_num > 0:
                    heappush(parents_potential, (mi, i))
                    parents_num -= 1
                elif parents_potential:
                    if parents_potential[0][0] < mi:
                        heapreplace(parents_potential, (mi, i))
            parents_potential.sort(reverse=True)
            mid_ratio = round(0.5 * gene_num * ratio)
            for (coe, no) in parents_potential:
                if mid_ratio:
                    parents_mi.append(self.vertices[no])
                    mid_ratio = mid_ratio - 1
                else:
                    parents_mi_aft.append(self.vertices[no])
            print "gene:", node.name, "parents_mi==="
            for gene123 in parents_mi:
                print gene123.name
            print "gene:", node.name, "parents_mi_aft==="
            for gene123 in parents_mi_aft:
                print gene123.name
            # ############==========================mi=====================###################
            ############==========================MIC=====================###################
            parents_potential = []
            parents_MIC = []
            parents_MIC_aft = []
            parents_num = round(gene_num * ratio)
            for i in range(len(parents_select)):
                if parents_select[node.index] == parents_select[i]:
                    continue
                mine = MINE()
                mine.compute_score(parents_select[node.index], parents_select[i])
                mi = mine.mic()
                if parents_num > 0:
                    heappush(parents_potential, (mi, i))
                    parents_num -= 1
                elif parents_potential:
                    if parents_potential[0][0] < mi:
                        heapreplace(parents_potential, (mi, i))
            parents_MIC_aft.sort(reverse=True)
            mid_ratio = round(0.5 * gene_num * ratio)
            for (coe, no) in parents_potential:
                if mid_ratio:
                    parents_MIC.append(self.vertices[no])
                    mid_ratio = mid_ratio - 1
                else:
                    parents_MIC_aft.append(self.vertices[no])
            print "gene:", node.name, "parents_MIC==="
            for gene123 in parents_MIC:
                print gene123.name
            print "gene:", node.name, "parents_MIC_aft==="
            for gene123 in parents_MIC_aft:
                print gene123.name
            ############==========================MIC=====================###################
            ############==========================spearman=====================###################
            parents_potential = []
            parents_spearman = []
            parents_spearman_aft = []
            parents_num = round(gene_num * ratio)
            for i in range(len(parents_select)):
                if parents_select[node.index] == parents_select[i]:
                    continue
                correl=stats.spearmanr(parents_select[node.index], parents_select[i])
                abscorrel = abs(correl[0])
                if parents_num:
                    heappush(parents_potential, (abscorrel, i))
                    parents_num -= 1
                elif parents_potential > 0:
                    if parents_potential[0][0] < abscorrel:
                        heapreplace(parents_potential, (abscorrel, i))
            parents_spearman_aft.sort(reverse=True)
            mid_ratio = round(0.5 * gene_num * ratio)
            for (coe, no) in parents_potential:
                if mid_ratio:
                    parents_spearman.append(self.vertices[no])
                    mid_ratio = mid_ratio - 1
                else:
                    parents_spearman_aft.append(self.vertices[no])
            print "gene:", node.name, "parents_spearman==="
            for gene123 in parents_spearman:
                print gene123.name
            print "gene:", node.name, "parents_spearman_aft==="
            for gene123 in parents_spearman_aft:
                print gene123.name
            ############==========================spearman=====================###################
        # if not sloops:
        #     try:
        #         i=parents.index(node)
        #         del parents[i]
        #     except ValueError:
        #         pass
        pars=parents_pearson+parents_mi+parents_MIC+parents_spearman
        parents = list(set(pars))
        parents_aft=parents_spearman_aft+parents_MIC_aft+parents_pearson_aft
        c = dict(Counter(parents_aft))
        pars_aft = [key for key, value in c.items() if value > 1]
        parents=parents+pars_aft
        parents = list(set(parents))
        www=90
        print node.name,parents.__len__()
        return parents

    def get_potential_parents_single(self, node, parents_select, sloops=True):
        """selects potential parents for node"""
        gene_num=self.n_gen
        ratio=0.2
        if self.parents.has_key(node.name): # potential parents of a node are explicitly specified
            parents=map(lambda n: self.vertices[self.vertice_names.index(n)],self.parents[node.name])
        else: # all vertices are potential parents
            parents = []
            ############==========================pearsonr=====================###################
            parents_potential = []
            parents_mi = []
            parents_num = round(gene_num * ratio)
            for i in range(len(parents_select)):
                if parents_select[node.index] == parents_select[i]:
                    continue
                correl=stats.pearsonr(parents_select[node.index],parents_select[i])
                abscorrel=abs(correl[0])
                if parents_num > 0:
                    heappush(parents_potential, (abscorrel, i))
                    parents_num -= 1
                elif parents_potential:
                    if parents_potential[0][0]<abscorrel:
                        heapreplace(parents_potential, (abscorrel, i))
            for (coe,no) in parents_potential:
                parents_mi.append(self.vertices[no])
            ############==========================pearsonr=====================###################
        if not sloops:
            try:
                i=parents.index(node)
                del parents[i]
            except ValueError:
                pass
        parents = parents_mi
        www=90
        print 'hey',node.name,parents.__len__()
        return parents




    def get_prior(self, vertex, parent):
        return self.prior.get((vertex.name, parent.name), self.prior.get(parent.name, 1))
    def write_txt(self,subpars,file_name):
        """Outputs suboptimal parents sets into a text file
        """

        f=open(file_name,"w")
        for v,minlist in subpars.items():
            f.write('\n%s'% v.name)
            for prob, pars in minlist:
                prob_s=util.safe_exponent(prob)
                f.write('\n %s '% prob_s)
                for p in pars:
                    f.write(' %s'% p.name)
        f.close()


def learn_x(v):
    #x, self, min_empty, min_optim, scr, verbose, n_min, limit, fpr_factor, fpr_nodes = v
    x,self,min_empty,min_optim,scr,verbose,n_min,limit,fpr_factor,fpr_nodes,parents_matrix=v
    selected_data = self.select_2(x, self.regulators, parents_matrix, sloops=scr.sloops)
    #selected_data = self.select_2(x, self.regulators, sloops=scr.sloops)
    x_priors = {p : self.get_prior(x,p) for p in selected_data.parents}
    selected_data.weight_parents(x_priors,fpr_factor,fpr_nodes,scr)
    if min_empty:
        score_empty=scr.data_score(selected_data.subset([]))+scr.graph_score(len(selected_data.parents),x,[],len(selected_data))
        score_max=score_empty-math.log(min_empty,2)
    else:
        score_max=fpconst.PosInf
    if min_optim:
        score_delta=-math.log(min_optim,2)
    else:
        score_delta=fpconst.PosInf
    (par,sc),minlist = scr.learn_1(selected_data,verbose,n_min,limit,score_max,score_delta)
    ###plus the delete edges
    # y = self.regulators[x.index]
    # X=[]
    # for item in par:
    #     X.append(self.regulators[item.index])
    # alpha = 0.05  ###
    # lamda = 1
    # parro = reoptim(y, X, lamda, alpha)
    # pars = []
    # for no in parro:
    #     pars.append(par[no])
    # par=pars
    return x,par,minlist,sc



class dataset_1:
    """Class representing data selected for learning parents of one vertex
    """

    def __init__(self, data, vertex, parents):
        self.data = data
        self.vertex = vertex
        self.parents = parents
        self.weights = {}

    def weight_parents(self, priors, fpr_factor=None,fpr_nodes=False, score=None):
        """Assigns weights to potential parents
        Required before learning
        """
        if self.parents:
            if fpr_factor:  # calculate weights with fpr procedure
                if fpr_nodes:
                    normalizer = 0.
                    for parent in self.parents:
                        normalizer += 1. / priors[parent]
                    norm_factor = fpr_factor / normalizer
                else:
                    norm_factor = fpr_factor
                n_parents = len(self.parents)
                n_data = len(self)
                no_tries = max(100, int(round(log(2 + n_parents) * (max(priors.values()) / norm_factor))))
                # if max_tries:
                #     no_tries = min(no_tries, max_tries)

                self_empty = self.subset()
                empty_score = score.data_score(self_empty) + score.graph_score(n_parents, self.vertex, [], n_data)
                sample_scores = [[score.lower_bound_for_data_score(self_empty)] * n_parents]

                # no_tries times shuffle vertex data and add scores to the distributions
                v_data = [d for d in self.data[-1]]  # create a copy before shuffling
                for i in range(no_tries):
                    shuffle(self.data[-1])
                    sample_scores.append([score.data_score(self.subset([par])) for par in self.parents])
                distributions = [sorted(scores) for scores in zip(*sample_scores)]
                self.data[-1] = v_data

                for i, par in enumerate(self.parents):
                    # compute data_score thresholds
                    thr_ind = no_tries * norm_factor / priors[par]
                    ind1 = min(no_tries, int(thr_ind))
                    ind2 = min(no_tries, ind1 + 1)
                    p = thr_ind - ind1
                    ds_thr = (1 - p) * distributions[i][ind1] + p * distributions[i][ind2]

                    # compute required parents' graph_score
                    g_score = max(0, empty_score - ds_thr)

                    # compute parent's weight
                    r = 2.001
                    r_score = score.graph_score(n_parents, self.vertex, [r], n_data)
                    while r_score < g_score:
                        r *= 2
                        r_score = score.graph_score(n_parents, self.vertex, [r], n_data)
                    l = r / 2
                    l_score = score.graph_score(n_parents, self.vertex, [l], n_data)
                    while r - l > 1e-9 and r_score - l_score > 1e-5:
                        w = (l + r) / 2
                        w_score = score.graph_score(n_parents, self.vertex, [w], n_data)
                        if w_score < g_score:
                            l = w
                            l_score = w_score
                        else:
                            r = w
                            r_score = w_score
                    self.weights[par] = w

            else:  # calculate weights based on priors only
                for par in self.parents:
                    self.weights[par] = par.base_weight() ** priors[par]
            # sort parents and data according to weights
            order = [i for i, pp in sorted(enumerate(self.parents), key=lambda (i, par): self.weights[par])]
            self.parents = [self.parents[i] for i in order]
            self.data = [self.data[i] for i in order] + [self.data[-1]]

    def __len__(self):
        return len(self.data[-1])

    def signs(self):
        """Calculates signs of correlation coefficients between a child and its parents
        """
        res = {}
        for i, p in enumerate(self.parents):
            try:
                sign = stats.pearsonr(self.data[i], self.data[-1])[0]
            except ValueError:
                print "problem finding out sign of interaction"  # self.data[i],self.data[-1]
                res[p] = "*"
            else:
                if sign >= 0:
                    res[p] = "+"
                else:
                    res[p] = "-"
        return res

    def subset(self, sub=[]):
        """Creates a new dataset_1 object with a subset of parents
        """
        ssub = set(sub)
        sub_ind = [i for i, p in enumerate(self.parents) if p in ssub]
        parents = [self.parents[i] for i in sub_ind]
        data = [self.data[i] for i in sub_ind]
        data.append(self.data[-1])
        d1 = dataset_1(data, self.vertex, parents)
        if self.weights:
            d1.weights = {par: self.weights[par] for par in parents}
        return d1

    def stats(self):
        """Counts frequencies of value vectors of
        parents and child (stats_all) and
        parents only (stats_parents)
        """

        def stats_disc(data):
            stats_all = {}
            stats_par = {}
            for d in zip(*data):
                d_par = tuple(d[:-1])
                d_all = tuple(d)
                stats_par[d_par] = stats_par.get(d_par, 0) + 1
                stats_all[d_all] = stats_all.get(d_all, 0) + 1
            return stats_all, stats_par

        def strec(key, prob, d):
            first, first_disc = d[0]
            rest = d[1:]
            if rest == []:
                stats_par[key] = stats_par.get(key, 0) + prob
                if first_disc:
                    stats_all[key + (first,)] = stats_all.get(key + (first,), 0) + prob
                else:
                    stats_all[key + (0,)] = stats_all.get(key + (0,), 0) + prob * (1 - first)
                    stats_all[key + (1,)] = stats_all.get(key + (1,), 0) + prob * first
            else:
                if first_disc:
                    strec(key + (first,), prob, rest)
                else:
                    strec(key + (0,), prob * (1 - first), rest)
                    strec(key + (1,), prob * first, rest)

        d_disc = map(lambda p: p.n_disc, self.parents + [self.vertex])
        if 0 not in d_disc:  # are all variables discrete?
            return stats_disc(self.data)  # -yes
        # -no
        stats_all = {}
        stats_par = {}
        for d in zip(*self.data):
            strec((), 1, zip(d, d_disc))
        return stats_all, stats_par