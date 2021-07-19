#coding: utf-8
import numpy as np
from scipy.sparse import csr_matrix
from scipy import sparse
from scipy import optimize
def lp_tgn(Y, X, lamda):

    lamda = 1
    n, m = Y.shape
    p, q = X.shape
    c = n*m
    h = n*p
    f1 = np.ones([2*c, 1])
    f2 = lamda*np.ones([2*h, 1])
    f = np.vstack((f1, f2))
    # f=f.transpose()
    Y_1 = Y.transpose()
    b_eq = np.array(Y_1)###缺少一步把Y_1每列合并为一个列向量
    row = []
    for i in range(0, c):
        row.append(i)
    col = []
    for i in range(0, c):
        col.append(i)
    data = np.ones(c)
    A1 = csr_matrix((data, (np.array(row), np.array(col))), shape=(c, c)).toarray()
    data1 = -np.ones(c)
    A2 = csr_matrix((data1, (np.array(row), np.array(col))), shape=(c, c)).toarray()
    Z = np.transpose(X)
    Z_1 = Z
    for i in range(0, n-1):
        Z_1=np.blkdiag(Z_1, Z)###缺少对Z_1的对角线操作
    Z_2 = -Z_1
    Aeq = np.hstack((A1, A2))
    Aeq = np.hstack((Aeq, Z_1))
    Aeq = np.hstack((Aeq, Z_2))
    A_eq = np.array(Aeq)

    lb = np.zeros((2*(c+h), 1))
    bound = []
    num = 2*(c+h)
    for i in range(0, num):
        bound.append((0, None))
    # print bound
    bound = tuple(bound)
    x = optimize.linprog(f, A_eq=A_eq, b_eq=b_eq, bounds=bound)
    # print x

    x = x.x
    if isinstance(x,float):
        return 100
    # print type(x)
    x.transpose()
    s = np.zeros((n, p))
    t = np.zeros((n, p))
    J = np.zeros((n, p))
    for i in range(1, n+1):
        for k in range(0, p):
            s[i-1][k] = x[2*n*m+(i-1)*p+k]
            t[i-1][k] = x[2*n*m+n*p+(i-1)*p+k]
    J = s-t
    return J
def reoptim(y, X, lamda, alpha):
# if __name__=="__main__":
    par=[]
    if X:
        aaa=1
    else:
        return par
    alpha = 0.01###
    lamda = 1
    # X = [[7.27115,7.34202,7.47145,7.40772,7.34687,7.2781,6.90899,7.24638,7.22903,7.44396,7.48013,7.46915,7.25278,7.25721,7.23643,7.40374,7.40313,7.35249,7.60202,7.3977,7.35311,7.92033,7.34215,7.50883,7.59041,7.51859,7.76208,7.25042,7.40896,7.71962,7.60947,6.61947,7.1655,7.11698,6.79979,7.5137,7.51889,6.70122,7.09738,7.06708,6.83741,6.71678,6.78401,6.63479,7.4503,7.05799,6.74564,6.59339,7.56551,7.49477,7.48145,7.75082,7.50698,7.45544,7.35131,7.3174,7.43963,7.31067,7.60301,7.37945,7.95273,7.43532,7.71917,7.5359,7.60282,7.46779,7.11646,7.24586,6.82614,6.82578,6.79895,7.07746,6.61795,6.41606,7.19077,7.09761,7.07811,7.18269,6.69542,6.66914,7.2801,7.06089,7.90371,7.36667,7.4576,7.45522,7.67412,7.72422,7.66936,7.50939,7.39191,7.22094,7.64488,7.61062,7.55139,7.20251,7.77469,7.77282,7.90681,7.36765],
    #      [8.52764,8.37024,8.35357,8.45709,8.35158,8.49377,8.61388,8.53042,8.71845,8.31984,8.43194,8.36732,8.47082,8.53232,8.62372,8.3487,8.38408,8.45807,8.42419,8.61958,8.7405,8.62279,8.67056,8.77746,7.83429,8.24566,7.71816,8.01624,7.85634,7.95388,8.27466,8.62334,8.44718,8.30232,8.6473,8.49359,8.57868,8.60725,8.46676,8.46954,8.46625,8.60997,8.56635,8.66132,8.50955,8.56761,8.75657,8.73129,8.05512,7.96947,7.76601,8.15491,8.11892,7.59101,8.60536,8.55159,8.61417,8.22201,7.96002,7.8583,8.82929,8.71717,8.29886,8.18211,8.22845,7.71834,8.49893,8.57698,8.37918,8.46206,8.83276,8.78178,8.96293,8.90819,8.484,8.46072,8.60546,8.34651,9.38133,9.03446,8.11361,8.36858,8.14509,8.45093,8.39702,8.37998,8.17759,8.3616,8.34502,7.91751,7.60726,7.93757,7.92077,8.22992,7.63906,7.94753,8.39679,7.85365,8.73734,8.86213],
    #      [7.7489,7.61698,8.15653,6.30678,6.24277,7.73988,7.65934,8.44874,7.01438,6.29407,6.04706,6.09265,6.44004,6.20284,6.46298,6.05785,6.1301,6.14035,6.33775,6.56551,6.3039,6.1779,6.30212,6.03285,8.78342,8.28687,8.20916,8.54521,8.37137,8.56552,6.42871,9.68673,6.96584,8.33473,9.13122,7.11395,7.06914,6.5429,7.12989,6.60463,7.30119,6.9582,7.12134,6.94934,8.86325,8.32163,9.12657,9.05212,8.21608,8.38398,8.15911,8.09951,8.40803,8.38745,6.49829,6.32687,6.46488,8.4857,8.33361,8.4397,6.35109,6.28565,6.11301,8.3308,8.26803,8.67195,6.03099,5.82696,6.84671,6.93233,7.46455,7.40057,11.1585,11.3168,5.969,6.00334,6.26399,6.51524,6.5572,6.61811,8.11618,7.73068,8.49789,7.24661,7.33591,7.00443,7.0042,6.81912,6.87528,7.7599,7.88375,8.29217,7.59873,7.92289,8.33539,7.49384,8.01658,8.27273,6.09911,6.17532],
    #      [9.83726,10.0249,10.1733,9.6876,9.58425,9.62205,9.53728,9.606,9.45068,9.78349,9.87048,9.78274,9.46495,9.4187,9.54709,9.47923,9.44928,9.51089,9.16018,9.08763,9.57599,9.286,9.19433,9.72894,10.1048,10.0478,9.90265,9.96097,9.90369,9.97809,9.61355,9.8538,9.5496,9.73097,9.71611,9.16198,9.16284,9.35784,9.20732,9.07899,9.14724,9.1199,9.16069,9.2872,9.10307,8.86483,9.04735,9.03689,9.95859,9.96618,9.65973,9.82436,10.1547,9.81076,9.53948,9.70638,9.95062,9.84873,9.848,9.8271,9.40693,9.45029,9.56888,10.088,10.0958,9.99864,9.22991,9.43759,9.16024,9.19017,9.50339,9.33109,9.98849,9.93914,9.20235,9.35565,9.15356,9.02365,8.59928,8.50121,9.83984,9.9888,10.0447,9.96643,9.89377,9.80635,9.84757,9.77705,9.73088,9.96373,9.94257,9.92507,9.80318,9.93152,9.85732,10.0941,10.1053,9.85397,9.62556,9.55723],
    #      [7.88838,7.85292,7.85123,8.66684,8.73663,9.92755,9.83971,9.18583,10.474,8.60158,8.58624,8.40628,8.23063,8.24371,8.62643,8.08057,8.5439,8.29213,7.684,7.76227,7.91382,7.82883,7.8367,7.79397,8.39875,8.26595,8.22852,8.36906,8.1956,8.29014,7.59201,8.80183,7.70004,7.36154,7.95836,7.46791,7.55572,7.95843,7.53832,7.64967,7.63222,8.02716,7.87287,7.73284,7.58018,7.68444,8.05662,8.68884,8.13995,7.6639,8.23044,8.12156,8.25615,8.30963,8.38346,8.53554,7.79223,8.19334,8.16685,8.08302,7.74668,7.76914,7.42555,8.18487,8.26772,8.22851,7.60281,7.99036,8.08931,8.03154,8.015,7.92402,12.666,12.3254,7.63517,7.78476,7.6346,7.69691,9.39457,9.30185,8.10381,8.12273,8.33254,8.03944,7.91476,7.80593,7.94367,7.47716,7.73631,8.10512,7.8244,8.13823,7.96777,7.90357,7.9838,7.95051,8.10066,8.09111,7.54098,7.91157],
    #      [11.7597,11.8792,11.8129,11.4809,11.4977,11.5549,11.3079,11.2384,11.2143,11.5325,11.6479,11.6041,11.6118,11.6685,11.4439,11.3062,11.6782,11.4532,11.2696,11.3135,11.231,11.0243,11.439,11.2659,11.1473,11.3365,11.3597,11.1585,11.3182,11.2067,11.6924,11.1865,11.3643,11.2452,11.1256,11.3081,11.1506,10.7177,11.0252,11.1918,10.9154,11.0308,11.1248,11.0152,11.2025,10.5963,10.539,10.2648,11.215,10.8223,11.3835,11.2417,11.037,11.4664,11.8353,11.3831,11.9039,11.2868,11.4501,11.3493,11.3864,11.0547,11.6679,11.4006,11.5386,11.3433,11.7424,11.9388,11.6598,11.7523,11.8063,11.7359,10.4556,10.9168,11.5838,11.722,11.4984,11.7015,9.6447,9.53304,11.3828,11.4326,11.4327,11.3139,11.3575,11.2897,11.2056,11.3429,11.2161,11.3393,10.2715,11.2806,11.2915,11.4183,11.1688,11.435,10.9278,11.3045,11.6427,11.481],
    #      [9.49539,9.61674,9.84158,9.7527,9.75633,9.83815,9.6036,9.57055,9.81554,9.73316,9.66057,9.72355,9.4983,9.44919,9.62046,9.42263,9.78235,9.32088,8.62775,8.55286,9.01239,9.11415,8.59154,8.88111,8.92372,8.92983,9.38281,8.9647,9.04214,9.05968,8.92219,9.2311,8.76079,8.72174,9.24771,8.50813,8.46138,8.19075,8.44123,8.48332,8.39237,8.22683,8.38778,8.54887,8.96368,8.65682,8.71379,8.79047,9.10808,9.05741,8.73819,9.05375,9.52234,8.98785,9.3243,9.41744,8.86225,9.11652,9.25051,8.67627,9.09499,8.96118,9.41186,9.394,9.39861,8.64858,9.49995,9.42737,9.13577,9.15815,9.41142,9.48629,9.67883,9.66697,9.12537,9.10919,9.20373,9.31547,9.68337,9.42766,9.10681,9.29891,9.22878,9.27456,9.09827,9.23715,9.22206,9.07714,9.09604,9.02465,9.26149,8.8956,8.96637,9.25486,8.95956,9.16029,9.36774,9.14015,9.28995,9.37912],
    #      [7.56622,7.44385,7.71532,7.78181,7.58769,7.6816,7.43434,7.67846,7.80509,7.56678,7.38434,7.58599,7.46619,7.45173,7.54859,7.55168,7.56983,7.55162,7.70735,7.63315,7.96859,7.92736,7.81097,7.83316,7.69116,7.52406,7.54426,7.74107,7.67009,7.7177,7.48534,7.53013,7.55401,7.67232,7.78262,7.17593,7.36834,7.67026,7.27854,7.55526,7.44469,7.47424,7.71729,7.50123,7.50731,7.25726,7.69704,7.33412,7.76436,7.83865,7.62938,7.82471,7.70933,7.70363,7.67226,7.60105,7.60489,7.36672,7.63202,7.82698,7.93497,7.90024,7.64105,7.5709,7.6269,7.61778,7.53225,7.39839,7.37251,7.48399,7.55381,7.66426,7.9718,7.75408,7.61275,7.49794,7.51168,7.94033,8.2577,8.11646,7.51446,7.38216,7.49068,7.68992,7.50147,7.74236,7.69905,7.79297,7.62681,7.611,7.87143,7.80644,7.52815,7.41674,7.71738,7.4544,7.41394,7.49773,7.78231,7.87418],
    #      [4.33721,4.28538,4.21563,4.31541,4.37142,4.30921,5.74099,5.56373,5.48113,4.94044,4.90572,5.17374,5.40787,4.86138,5.37891,4.95992,5.0748,4.88242,4.27578,4.58095,4.38204,4.79897,4.7428,4.79778,4.73438,4.55062,4.44763,4.47446,4.4717,4.64452,4.22001,4.89995,4.46194,5.75855,6.95577,3.99353,3.90766,4.38751,4.4032,4.32033,5.17854,4.73839,4.79533,4.53558,5.59744,7.00824,7.68237,8.2886,4.65915,4.41773,4.48247,5.01584,4.66077,4.50582,6.17221,5.29862,5.28864,4.7356,4.55541,4.48671,5.49401,4.84895,4.5174,4.98046,4.70298,4.45166,4.21544,4.5405,4.6063,4.83322,5.8035,6.19621,4.95001,5.3222,4.27609,4.3377,4.31941,4.45307,5.35159,4.76967,4.69503,4.32817,5.067,4.60243,4.50294,4.5129,4.60917,4.47737,4.7543,4.4438,4.39856,4.48087,4.43127,4.72673,4.50182,4.34487,4.94673,4.41748,4.8362,4.60717],
    #      [7.66698,7.66822,7.72285,7.73853,7.75253,7.89045,7.80683,8.18775,7.91553,7.6836,7.63611,7.77916,7.84375,7.87908,7.60833,7.47311,7.80765,7.65053,13.0931,12.5181,12.1647,12.4149,12.5968,12.0008,9.13184,9.55152,8.95199,8.93042,9.11013,8.30692,8.29314,10.233,10.6742,10.2265,10.0733,8.30078,8.09201,8.6887,9.4783,9.53319,9.62957,9.90406,9.38692,9.35165,8.34586,8.99162,10.4698,9.77803,9.08072,9.44773,8.86079,8.822,9.38566,9.08814,8.05157,8.00368,7.66898,9.05751,9.49254,8.65946,8.43091,8.19398,8.14368,9.60246,9.3362,8.64487,6.70604,6.7804,6.61744,6.4466,7.79807,7.63144,7.20038,7.55481,6.8531,6.62767,7.36033,7.38374,7.86565,8.07585,9.10713,9.37329,9.05362,8.77216,8.84181,8.85673,8.79537,9.01628,8.99129,8.96988,8.97657,8.5774,9.19824,9.02568,8.7471,9.48147,9.31569,8.78291,8.30442,7.82618],
    #      [10.2011,10.1187,10.2137,10.084,10.1735,9.98081,7.94797,7.97867,7.98438,7.88327,8.32524,8.09839,8.119,8.189,7.9668,8.1514,8.16326,8.17208,10.1719,10.002,10.3602,10.4448,10.3221,10.4848,10.2968,10.0882,10.6885,10.3498,10.2082,10.7655,10.4397,10.3029,10.5093,10.326,10.1598,10.5337,10.1974,9.77016,10.5038,10.5035,10.4764,10.4394,10.0604,10.0259,10.4634,10.4528,10.3198,10.3717,10.126,9.92693,10.6045,10.2744,9.76125,10.6691,11.3168,11.3234,10.8717,10.1389,10.2283,10.557,10.1427,10.2311,10.446,10.3244,10.4194,10.4603,10.33,10.6457,10.4179,10.3075,10.3951,10.3189,11.8881,11.7655,10.1985,10.4303,10.0307,10.0189,11.1711,10.9617,10.1981,10.2959,10.4742,10.1557,10.1258,10.2306,10.1046,10.1875,10.2926,10.2142,9.89872,10.3544,10.4521,10.5221,10.5233,10.4024,10.1653,10.3834,10.6213,10.7198],
    #      [7.67362,7.99922,8.04277,7.73782,7.78105,7.85629,7.58657,7.63438,7.76606,7.76635,7.63799,7.79022,7.81681,7.76652,7.72399,7.68172,7.66182,7.51844,7.43628,7.34464,7.39224,7.38544,7.22989,7.38907,7.46618,7.42858,7.58486,7.15312,7.35585,7.4197,7.88651,7.29866,7.58963,7.50043,7.3598,6.96044,6.80492,6.9059,6.93083,6.9756,7.16011,6.81462,6.8499,7.00176,7.25668,6.94071,6.95953,7.0922,7.45326,7.42676,7.22827,7.42672,7.17748,7.2551,7.37981,7.29451,7.54922,7.24403,7.35311,7.31642,7.41503,7.32518,7.34694,7.27876,7.22803,7.34266,7.60072,7.86046,7.53088,7.59193,7.52513,7.71467,7.00648,7.56228,7.59927,7.71729,7.63655,7.60404,7.14653,7.30809,7.114,6.509,7.15973,7.58611,7.59022,7.38481,7.59683,7.73215,7.47628,7.27207,7.25768,7.2892,7.18595,7.22831,7.22507,7.02321,7.27453,7.30944,7.35253,7.46777],
    #      [9.50863,9.57927,9.42919,9.26574,9.13601,8.72151,9.44167,9.67699,9.37589,9.31234,9.38271,9.30456,9.59538,9.50413,9.65508,9.35271,9.36942,9.37,8.71084,8.84085,8.89887,8.91355,8.75436,8.80081,9.10982,9.20639,9.39406,8.7909,9.11028,8.82842,8.58875,8.27029,8.62423,8.59333,8.23819,9.10325,8.98207,8.54689,9.14496,8.88773,8.76911,9.02655,8.47045,8.57391,8.77483,8.98201,8.51186,8.4672,8.81826,9.14872,8.81691,8.69311,9.10157,8.80786,8.71462,8.67865,8.90373,8.96262,9.13469,8.84692,8.85885,8.39516,8.83335,8.93336,9.21791,9.09931,8.78472,8.60907,8.2672,8.39071,9.2947,9.048,8.43742,8.5487,8.641,8.85297,9.01029,9.14361,7.96335,8.24758,9.16688,9.18495,9.10607,8.64607,8.54834,8.55555,8.85675,8.63855,8.66984,9.11361,9.56387,8.84756,9.34894,9.15894,9.28882,9.34437,9.33545,9.0186,8.81499,8.8062],
    #      [8.09553,8.23207,8.30046,8.53182,8.39968,8.42921,8.13926,8.14979,8.27252,8.27224,8.2775,8.35988,8.16878,8.19471,8.4538,8.44271,8.18776,8.33696,8.69515,8.5555,9.13925,8.94184,8.56477,9.07736,8.74843,8.65278,8.75065,8.77346,8.86862,8.77464,8.50701,8.77417,8.79944,8.85061,8.72961,8.6065,8.56608,9.15338,8.62432,8.75647,8.58302,8.55933,8.70815,8.77726,8.48179,8.37101,8.43983,8.7675,8.72908,8.90816,8.83733,8.78764,8.68166,8.88163,8.79,8.67532,8.42913,8.52476,8.9737,8.73775,8.8093,9.15553,8.89979,8.58313,8.64811,8.8179,7.91407,8.20684,8.3976,8.373,8.13021,8.0424,8.86022,8.57197,8.44664,8.66903,8.49073,8.66208,9.04492,9.1667,8.67595,8.42024,8.67141,8.73628,8.74156,8.67308,8.6996,8.63768,8.61845,8.70995,8.53023,8.79862,8.9252,8.61981,8.9669,8.65724,8.26731,8.5409,8.83901,9.01167]]
    # y = [11.6948000000000,11.6217000000000,11.4917000000000,11.3815000000000,11.3900000000000,11.0106000000000,11.1130000000000,10.7701000000000,10.9724000000000,11.4267000000000,11.5404000000000,11.4645000000000,11.2708000000000,11.2179000000000,11.0796000000000,11.2799000000000,11.5255000000000,11.3076000000000,11.2463000000000,11.4694000000000,10.8372000000000,11.3573000000000,11.5740000000000,10.8398000000000,11.3164000000000,11.3641000000000,11.2718000000000,11.3360000000000,11.4813000000000,11.2767000000000,11.7839000000000,10.7376000000000,11.3629000000000,11.1976000000000,10.5982000000000,11.9811000000000,11.9583000000000,9.86958000000000,11.6155000000000,11.9504000000000,11.1996000000000,11.3951000000000,10.4489000000000,10.4473000000000,11.6601000000000,11.2729000000000,10.9692000000000,10.4439000000000,11.4315000000000,11.4351000000000,11.5194000000000,11.3462000000000,11.2329000000000,11.3233000000000,11.0291000000000,10.5557000000000,11.2373000000000,11.4438000000000,11.5153000000000,11.4231000000000,11.6503000000000,10.9384000000000,11.0440000000000,11.5522000000000,11.6176000000000,11.4395000000000,11.4243000000000,11.4930000000000,11.4366000000000,11.3836000000000,11.4332000000000,11.5109000000000,9.86137000000000,10.3323000000000,11.6038000000000,11.5329000000000,11.8851000000000,11.9697000000000,8.54514000000000,8.74382000000000,11.6530000000000,11.4590000000000,11.6205000000000,11.7587000000000,11.9109000000000,11.8187000000000,11.8401000000000,11.8901000000000,11.8927000000000,11.5929000000000,11.4200000000000,11.4298000000000,11.7924000000000,11.4847000000000,11.5061000000000,11.7240000000000,11.3013000000000,11.4260000000000,11.1340000000000,10.9064000000000]
    X = np.matrix(X)
    y = np.matrix(y)
    X = X[:, 0:100]
    y = y[:, 0:100]
    at=X[-1,3]
    bt = X[-1, :5]
    J = lp_tgn(y, X, lamda)
    p, q = X.shape
    if isinstance(J, int) and J == 100:
        e = np.zeros([q,1])
        return e
    J_sparse = np.zeros(J.shape)
    J_value = np.zeros(J.shape)
    J_index = [j for j in J]
    J_index = J_index[0].tolist()
    index = [i for i, x in enumerate(J_index) if abs(x) >= alpha]
    index_c = [i for i, x in enumerate(J_index) if abs(x) < alpha]
    for ind in index:
        J_sparse[0][ind] = J[0][ind]
    # for inc in index_c:
    #     J_value[0][inc] = J[0][inc]
    # XX = None
    # while index_c:
    #     for i in index:
    #         if XX is None:
    #             XX = X[i]
    #         else:
    #             XX = np.vstack((XX, X[i]))
    #     J1 = lp_tgn(y, XX, lamda)
    #     J_index = [j for j in J1]
    #     J_index = J_index[0].tolist()
    #     index1 = [i for i, x in enumerate(J_index) if abs(x) >= alpha]
    #     index1_c = [i for i, x in enumerate(J_index) if abs(x) < alpha]
    #     index_cc = []
    #     for i in index1_c:
    #         index_cc.append(index[i])
    #     index_c = index_cc
    #     index_x = []
    #     for i in index1:
    #         index_x.append(index[i])
    #     index = index_x
    #     i = 1
    #     j = 1
    #
    #     for ind1, ind in zip(index1,index):
    #         J_sparse[0][ind] = J1[0][ind1]
    #     for indc in index_c:
    #         J_sparse[0][indc] = 0
    #     for ind1_c, ind_c in zip(index1_c,index_c):
    #         J_value[0][ind_c] = J1[0][ind1_c]
    #     YYY = 90
    J_value = J_sparse+J_value
    i=0

    # for value in J_value[0]:
    #     if abs(value) > 0.2:
    #         par.append(i)
    #     i=i+1

    opl = 0
    # print J_value
    return J_value[0]

