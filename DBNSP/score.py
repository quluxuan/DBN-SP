import graph,fpconst
from heapq import heappush, heappop, heapreplace


class score:
    """
    Abstract class implementing a scoring function.

    To obtain a working implementation, one has to define methods graph_score and data_score
    """

    def __init__(self, data_factor=1.0, prior=None, sloops=False, **kwds):
        self.data_factor = data_factor
        self.sloops = sloops  # do we allow self-loops in the networks ?

    def graph_score(self, number_of_potential_parents, vertices, weights_of_parents, number_of_data_points):
        """
        Method for computing graph score factor - ABSTRACT method

        it computes the score given :
        number_of_parents - number of potential parents of a given gene
        vertices - given gene
        par_disc - list of weights_of_parents of parents
        number_of_data_points - number of data points
        """
        pass

    def lower_bound_for_data_score(self, selected_data_empty):
        """
        method for computing the lower bound of the data score factor

        parameters:
        vertices - number of the vertice
        data - list of datapoints :

        each datapoint is a vector containing 2*n_vert values (expression levels for all genes for time points t and t+1
        """
        return 0.0

    def data_score(self, selected_data):
        """
        abstract method for computing the data score factor

        parameters:
        gene_vertex - number of the vertice
        numbers_of_parents - list of numbers of parents
        data - list of datapoints :

        each datapoint is a vector containing 2*n_vert values (expression levels for all genes for time points t and t+1
        """

        pass

    def subsets(self, list, k):
        """
        generate all the k-long subsets of a list
        """

        s = graph.stack()
        s.put((list, k, []))
        try:
            while True:  # exception used to get out of the loop
                list, k, acc = s.get()
                if k == 0:
                    yield acc
                elif list == []:
                    pass
                else:
                    s.put((list[1:], k, acc))
                    s.put((list[1:], k - 1, acc + [list[0]]))
        except IndexError:
            pass  # we've exhausted all options

    def learn_1(self, selected_data, verbose=None, n_min=1, limit=None, score_max=fpconst.PosInf,
                score_delta=fpconst.PosInf):

        if verbose:
            print 'Learning parents of', selected_data.vertex.name, '...',

        #        if not self.sloops:
        #            selected_data.rm_sloops()
        v = selected_data.vertex
        nd = len(selected_data)
        parents = selected_data.parents
        p_weights = selected_data.weights
        n = len(parents)
        try:
            lim = int(limit)  ##Limit of the size of parents subsets
        except TypeError:  # limit was None
            lim = 3

        selected_data_empty = selected_data.subset([])
        mindata = self.lower_bound_for_data_score(selected_data_empty)
        min_set = minset(n_min, score_max, score_delta,
                         self.data_score(selected_data_empty) + self.graph_score(n, v, [], nd))  # empty parents set
        if n:  # are there any potential parents?
            w_min = p_weights[parents[0]]
            w_max = p_weights[parents[-1]]
            if w_min == w_max:  # we can use algorithm 2
                weight = w_min
                size = 1
                mg = self.graph_score(n, v, [weight], nd)
                while min_set.accepts(mg + mindata) and (size <= lim):  # we can possibly add (sub-)optimal scores
                    for sub in self.subsets(parents, size):
                        # print "sub.size ", len(sub)
                        selected_data_sub = selected_data.subset(sub)
                        min_set.add(mg + self.data_score(selected_data_sub), sub)
                    size += 1

                    mg = self.graph_score(n, v, [weight] * size, nd)
            else:  # we have to use algorithm 1
                subsets = []  # successors of considered yet potential parents sets
                for parent in parents:  # one-element parents sets
                    # print "one parent"
                    heappush(subsets, (self.graph_score(n, v, [p_weights[parent]], nd), [p_weights[parent]], [parent]))
                while subsets:
                    # print subsets
                    mg, weights, sub = heappop(subsets)
                    # print sub
                    if not min_set.accepts(mg + mindata):  # we cannot improve the score
                        break
                    selected_data_sub = selected_data.subset(sub)
                    min_set.add(mg + self.data_score(selected_data_sub), sub)
                    # insert sub's successors
                    if len(sub) < lim:
                        last_parent = parents.index(sub[-1])
                        for parent in parents[last_parent + 1:]:
                            sub_succ = sub + [parent]
                            weights_succ = weights + [p_weights[parent]]
                            mg_succ = self.graph_score(n, v, weights_succ, nd)
                            heappush(subsets, (mg_succ, weights_succ, sub_succ))
        if verbose:
            print 'done', min_set
        return min_set.optimal, min_set.tolist()

    #    def learn_all(self,vertices,data,n_points):
    #        par = {}
    #        for v in vertices:
    #            par[v] = self.learn_1(v,vertices,data,n_points)
    #        return par

    def score_graph(self, g, data):
        s = 0.0
        n_vert = len(g.vertices)
        for i, v in enumerate(g.vertices):
            p = g.parents(v)
            selected_data = data.select_1(i, p)
            sg = self.graph_score(n_vert, v, map(lambda par: par.n_disc, p), len(selected_data))
            sd = self.data_score(selected_data)
            s += sg + sd
        return s


class minset:
    def __init__(self, size, score_max, score_delta, emptyscore):
        self.optimal = ([], emptyscore)
        self.escore = emptyscore
        self.free = size
        self.mscore = score_max
        self.dscore = score_delta
        self.mset = []
        self.add(emptyscore, [])

    def add(self, score, parents):
        if score < min(self.optimal[1] + self.dscore, self.mscore):
            if self.free:
                heappush(self.mset, (-score, parents))
                self.free -= 1
                if not self.free:
                    self.mscore = -self.mset[0][0]
            elif self.mset:
                heapreplace(self.mset, (-score, parents))
                self.mscore = -self.mset[0][0]
        if score < self.optimal[1]:
            self.optimal = (parents, score)
            while self.mset and -self.mset[0][0] > self.optimal[1] + self.dscore:
                tmp = heappop(self.mset)
                self.free += 1

    def __str__(self):
        self.mset.sort()
        minstr = ''
        for (score, parents) in reversed(self.mset):
            minstr += '\n Score ' + str(-score) + ' for parents set:'
            for p in parents:
                minstr += ' ' + p.name
        return minstr

    def accepts(self, score):
        return score < min(self.optimal[1] + self.dscore, self.mscore)

    #    def optimal(self):
    #        negscore, parset= max(self.mset)
    #        return parset, -negscore

    def tolist(self):
        self.mset.sort()
        minlist = []
        for (score, parents) in self.mset:
            minlist.append((score + self.escore, parents))
        minlist.reverse()
        return minlist