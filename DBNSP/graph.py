class stack:
    def __init__(self):
        self.stack = []
    def put(self,item):
        self.stack.append(item)
        return self
    def get(self):
        res = self.stack[-1]
        self.stack=self.stack[:-1]
        return res
    def empty(self):
        return len(self.stack)==0
class graph:
    def __init__(self):
        self.vertices = []
        self.edges = {}
        self.edge_labelling = {}
        self.vertice_labelling = {}

    def add_edge(self, st, end, label=" "):
        try:
            self.edges[st].append(end)
        except KeyError:
            if st not in self.vertices:
                raise Exception('Wrong starting vertice')
            self.edges[st] = [end]

        self.edge_labelling[st, end] = label
    def add_vert(self, vert, label=" "):
        self.vertices.append(vert)
        self.vertice_labelling[vert] = label
        self.edges[vert] = []
    def __repr__(self):
        rep = "Graph: \n"
        for v in self.vertices:
            rep += "\t%s" % v
            try:
                rep += "(%s) => " % str(self.vertice_labelling[v])
            except KeyError:
                rep += " => "
            for v2 in self.edges[v]:
                rep += "%s" % v2
                try:
                    rep += "(%s), " % str(self.edge_labelling[(v, v2)])
                except KeyError:
                    rep += ", "
            rep += "\n"
        # print self.edge_labelling, self.vertice_labelling
        return rep
    def parents(self, vert):
        par = []
        for v in self.vertices:
            if vert in self.edges[v]:
                par.append(v)
        return par
    def fromParents(self, verts, pars):
        self.__init__()
        self.vertices = []
        for v in verts:
            self.add_vert(v)
            self.edges[v] = []
        for v, par in zip(verts, pars):
            for p in par:
                self.add_edge(p, v)
        return self
    def num_edges(self):
        return sum(map(len, self.edges.values()))
    def weighted_edges(self, weights, fraction=0.0):
        # print "XX",fraction
        res = graph()
        res.vertices = self.vertices
        res.vertice_labelling = self.vertice_labelling
        for u in self.vertices:
            res.edges[u] = []
        for u in self.vertices:
            if weights:  # we have weights for all edges
                s = 0.0
                w_dict = {}
                if weights[u]:
                    w0, l = weights[u][0]
                    for w, l in weights[u]:  # calculate the weights for edges
                        # print w,l
                        w2 = 2 ** (w - w0)
                        s += w2
                        for node in l:
                            try:
                                w_dict[node] += w2
                            except KeyError:
                                w_dict[node] = w2
                    for v, w in w_dict.items():
                        # print w/s,fraction,v,u
                        if w / s > fraction:
                            res.add_edge(v, u, str(w / s))
        return res
    def to_SIF(self, weights=None):
        result = ""
        if weights:
            g1 = self.weighted_edges(weights)
        else:
            g1 = self
        for u in g1.vertices:
            for v in g1.edges[u]:
                label = g1.edge_labelling[(u, v)]
                result += "%s\t%s\t%s\n" % (g1.vertice_labelling[u], label, g1.vertice_labelling[v])
        # print "XX",result,g1.edges
        return result