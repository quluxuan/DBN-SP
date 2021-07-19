from BDE import BDE
from BIC import BIC

from data import dataset
import util
import datetime
start=datetime.datetime.now()
if __name__=="__main__":
    nowTime = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print nowTime
    import sys
    #from dispatch import *
    #test(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]))
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-s", "--score", dest="score",default="BDE",help="Scoring function: BDE (default) or BIC or MDL or MIT")
    parser.add_option("-e", "--expr", dest="expr",help="Expression data file")
    parser.add_option("-t", "--txt", dest="txt",help="Output file with the suboptimal parents sets")
    parser.add_option("-n", "--net", dest="net",help="Output file with the network structure")
    parser.add_option("-c", "--cpd", dest="cpd",help="Output file with the conditional probability distributions")
    parser.add_option("-b", "--bif", dest="bif",help="Output file with the optimal network in BIF format")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False, help="Print comments")
    parser.add_option("-p", "--prior-pseudocount", dest="prior",help="Total pseudocount number in Dirichlet priors of BDE score (default: 1 pseudocount for each possible value vector)")
    parser.add_option("-i", "--suboptimal", dest="n_min", default=1, type=int,help="Number of computed suboptimal parents sets")
    parser.add_option("-m", "--min-empty", dest="min_empty", default=None, type=float,help="Reported suboptimal parents sets' minimal probability ratio to the empty set")
    parser.add_option("-o", "--min-optimal", dest="min_optim", default=None, type=float,help="Reported suboptimal parents sets' minimal probability ratio to the optimal set")
    parser.add_option("-l", "--limit", dest="limit",help="Limit of the size of parents subsets")
    parser.add_option("-f", "--fraction", dest="frac", type=float,default=None,help="Minimum weight of a parent to be considered in a parent set.")
    parser.add_option("-r", "--fpr", dest="fpr", type=float, default=None, help="Expected proportion of false positives (i.e. FP edges unless -u is specified; switched off by default)")
    parser.add_option("-u", "--fpr-nodes", dest="fpr_nodes", action="store_true", default=False, help="Interpret -r parameter as the expected proportion of regulons having false positive regulators")
    parser.add_option("-g", "--sloops", dest="sloops", action='store_true', default=False, help="Allow self-loops in Dynamic Bayesian Networks (no self-loops by default)")
    parser.add_option("-k", "--cpu", dest="distrib", type=int, default=0, help="Number of processes to use for multiprocessing, defaults to 0 (no multiprocessing)")
    (options, args) = parser.parse_args()
    # options.expr = '/home/omnisky/bnfinder/doc/data/1565data/bd100_1565'
    # options.txt = "/home/omnisky/bnfinder/doc/data/1565data/outbd100_1565_10%.txt"
    # options.net = "/home/omnisky/bnfinder/doc/data/1565data/outbd100_1565_10%.sif"

    options.expr = 'E:/a/python27/bnfinder/doc/data/DREAM3/data/10genes_data.txt'
    options.txt = "E:/a/python27/bnfinder/doc/data/DREAM3/result/10genes/data10_10%.txt"
    options.net = "E:/a/python27/bnfinder/doc/data/DREAM3/result/10genes/data10_10%.sif"
    options.limit = 5
    options.distrib = 20
    # options.score = BIC
    # options.sloops = True


    if options.expr == None:
        print "You must provide an expression file name. Run bnf --help for more options."
    else:
        if options.verbose:
            print 'Loading data from:', options.expr, '...',
            d = dataset(options.expr).fromNewFile(open(options.expr))
        if options.verbose:
            print 'done\nLearning network ', d.name, ' with', options.score, 'scoring function'
        score = eval(options.score)(prior=options.prior, sloops=options.sloops)
        score, g, subpars = d.learn(score=score,fpr=options.fpr,fpr_nodes=options.fpr_nodes,prior=options.prior,verbose=options.verbose,n_min=options.n_min,limit=options.limit,min_empty=options.min_empty,min_optim=options.min_optim,distrib=options.distrib,)
        if options.verbose:
            print 'Total score of optimal network:', score
            # print g
            # print g.to_SIF(subpars)

        if options.txt:
            d.write_txt(subpars, options.txt)
            if options.verbose:
                print "Suboptimal parents sets written to:", options.txt

        if options.net:
            if options.n_min == 1:
                net_str = g.to_SIF()
            else:
                net_str = g.to_SIF(subpars)
            f = open(options.net, "w")
            f.write(net_str)
            f.close()
            # g.toDot().write_dot(options.net)
            # if options.verbose:
                # print "Network graph written to:", options.net
                # print net_str

        if options.cpd:
            f_cpd = open(options.cpd, "w")
            # d.write_cpd(g,f_cpd,options.prior)
            # f_cpd.write("XXXXXXXXXXXXXXXXX\n")
            # print "XXXXX",repr(options.frac)
            if options.frac:
                g1 = g.weighted_edges(subpars, float(options.frac))
                util.pprint(f_cpd, d.to_cpd(g1, options.prior))
            else:  # outpu the best network
                util.pprint(f_cpd, d.to_cpd(g, options.prior))
            f_cpd.close()
            if options.verbose:
                print "Conditional probability distributions written to:", options.cpd
        if options.bif:
            d.write_bif(g, options.bif, options.prior, ['Network graph learned with %s scoring function from expression data file \'%s\'' % (options.score, options.expr)])
            if options.verbose:
                print "Bayesian network written to:", options.bif
        end = datetime.datetime.now()
        print end-start