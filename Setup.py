import sys
import time
from hashlib import sha3_256, shake_128
from sage.all import *
# from utilities import generate_ldt_lists
from collections import defaultdict
from griffin import griffin_parameters

# Constants
P = 2 ** 127 - 1
L = 32768
B = 128
S = 128
ETA = 2
R = 4
RHO_STAR = 1 / 16
# polynomial degree parameter m decide SIZE_U
SIZE_U_OPTIONS = {32: (2 ** 12), 64: (2 ** 13)}
# (LOQUAT_STAR, LDT_SEC) decide LDT query complexity kappa
KAPPA_OPTIONS = {(1, 100): 50, (1, 128): 64, (0, 100): 25, (0, 128): 32}


def validate_params(loquat_star, algebraic_hash, ldt_sec_choice, m):
    if loquat_star not in [0, 1]:
        print(f"Invalid choice of loquat_star, require 0 or 1, given {loquat_star}")
        sys.exit()
    if algebraic_hash not in [0, 1]:
        print(f"Invalid choice of algebraic_hash, require 0 or 1, given {algebraic_hash}")
        sys.exit()
    if (loquat_star, ldt_sec_choice) not in KAPPA_OPTIONS.keys():
        print(f"Invalid choice of ldt secure parameters, require 100 or 128, given {ldt_sec_choice}")
        sys.exit()
    if m not in SIZE_U_OPTIONS.keys():
        print(f"Invalid choice of m, require 32 or 64, given {m}")
        sys.exit()


def setup_griffin_parameters(t, c):
    (p, t, capacity, security_level, d, dinv, N, mat, alphas, betas,
     round_constants) = griffin_parameters(P ** 2, t, c, S)
    parameters_hash = (p, t, capacity, security_level, d, dinv, N, mat,
                       alphas, betas, round_constants)
    return parameters_hash


def generate_ldt_lists(Fp2, size_H, size_U):
    ldt_lists = []
    g = Fp2.multiplicative_generator()
    g1 = g ** ((P ** 2 - 1) // size_H)
    listH = [(g1 ** i) * g for i in range(size_H)]
    g2 = g ** ((P ** 2 - 1) // size_U)
    dictU = {g2 ** i: () for i in range(size_U)}
    ldt_lists.append(dictU)

    # compute dictU1,...dictUr
    for _ in range(R):
        dict_i_plus_1 = defaultdict(list)
        for item in ldt_lists[-1]:
            item_new = item ** (2 ** ETA)
            dict_i_plus_1[item_new].append(item)
        ldt_lists.append(dict(dict_i_plus_1))
    # sorting ldt_lists for easier LDT
    sorted_ldt_lists = [list(ldt_lists[R].keys())]
    for i in range(R):
        temp_list = sorted_ldt_lists[0]
        temp_dic = ldt_lists[R - i]
        U_list = []
        for item in temp_list:
            if isinstance(item, list):
                U_list.extend([temp_dic[item[ind]] for ind in range(2**ETA)])
            else:
                U_list.extend([temp_dic[item]])
        sorted_ldt_lists.insert(0, U_list)
    return g, listH, sorted_ldt_lists


def loquat_setup(loquat_star=0, algebraic_hash=0, ldt_sec_choice=100, m=32):
    st = time.time()
    validate_params(loquat_star, algebraic_hash, ldt_sec_choice, m)

    kappa = KAPPA_OPTIONS[(loquat_star, ldt_sec_choice)]
    size_U = SIZE_U_OPTIONS[m]
    Fp = GF(P)
    Fp2 = GF(P ** 2, 'a')
    listI = [Fp.random_element() for _ in range(L)]
    n = B // m
    size_H = 2 * m
    # set parameters for univariate sumcheck and ldt
    g, listH, ldt_lists = generate_ldt_lists(Fp2, size_H, size_U)
    tree_cap = floor(log(kappa, 2) - 1)

    if algebraic_hash:
        # Griffin capacity
        c = ceil(2 * 128 / log(P ** 2, 2))
        loquat_hash = setup_griffin_parameters(4, c)
        loquat_expand = setup_griffin_parameters(3, c)
    else:
        loquat_hash = 0
        loquat_expand = 0

    pp = {
        "Fp": Fp,
        "Fp2": Fp2,
        "p": P,
        # "L": L,   # L can be find by len(listI)
        # "B": B,   # B can be find by m * n
        "listI": listI,
        "m": m,
        "n": n,
        "eta": ETA,
        "kappa": kappa,
        "rho_star": RHO_STAR,
        "r": R,
        "tree_cap": tree_cap,
        "g": g,
        "listH": listH,
        "ldt_lists": ldt_lists,
        "algebraic_hash": algebraic_hash,
        "loquat_hash": loquat_hash,
        "loquat_expand": loquat_expand
    }
    print("Setup running time: {} seconds".format(time.time() - st))
    print("-" * 50)
    return pp
