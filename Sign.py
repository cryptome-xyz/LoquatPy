import time

from sage.all import *
import numpy as np
import array as arr
from Utils import *


def loquat_sign_phase_1(sk, n, m, p, Fp, Fp2, listH, g, eta, kappa, ldt_lists, algebraic_hash, loquat_hash, tree_cap,
                        salt, poly_def, poly_Z_H):
    st = time.time()
    size_H = len(listH)
    ldt_eta = 2 ** eta
    # initialize vectors
    T, list_r, poly_c, poly_c_prime = ([] for _ in range(4))
    # compute witness polynomial
    for j in range(n):
        vec_c_j = []
        rs = [Fp.random_element() for _ in range(m)]
        list_r.extend(rs)
        T.extend([floor(1 / 2 * (1 - kronecker(r, p))) for r in rs])
        vec_c_j.extend([Fp2(sk * r), Fp2(r)] for r in rs)
        poly_c.append((poly_def.lagrange_polynomial(zip(listH, vec_c_j))))

    # convert list_r to a vector (for easier algebraic computation later)
    vec_r = vector(Fp, list_r)

    # generate degree kappa * 2^ldt_eta random polynomial
    poly_r = poly_def([Fp2.random_element() for _ in range(kappa * ldt_eta + 1)])
    # multiply the vanish polynomial with random polynomial (this is for masking the witness polynomial)
    poly_Z_H_r = poly_Z_H * poly_r
    # compute poly_c_prime that is the masked witness polynomial
    poly_c_prime = [poly + poly_Z_H_r for poly in poly_c]
    # evaluate the masked witness polynomial over elements of listU
    list_U = ldt_lists[0]
    eval_c_prime = [eval_poly_over_listU(poly, list_U) for poly in poly_c_prime]
    # print(eval_c_prime)
    # print(salt)
    leaf_c_prime = get_leaf_hash(True, eval_c_prime, algebraic_hash, loquat_hash, salt, Fp2)

    root_c, tree_c = MT_commit(algebraic_hash, loquat_hash, leaf_c_prime, tree_cap, salt, Fp2)
    sigma_1 = [root_c, T]

    print("Signing phase 1 running time: {} seconds".format(time.time() - st))
    print("-" * 50)

    return eval_c_prime, tree_c, poly_c_prime, vec_r, sigma_1


def loquat_sign_phase_2(Fp2, Fp, loquat_hash, loquat_expand, algebraic_hash, sigma_1, salt, msg, B, vec_r, listI, sk):
    st = time.time()
    h_1, phase_1_challenge = get_phase_1_challenge(Fp2, loquat_hash, loquat_expand, algebraic_hash, sigma_1, salt, msg,
                                                   B)

    # challenges contain B random selected elements from listI
    vec_I = vector(Fp, [listI[ch] for ch in phase_1_challenge])
    vec_o = vec_r.pairwise_product(vector([element + sk for element in vec_I]))

    sigma_2 = vec_o

    print("Signing phase 2 running time: {} seconds".format(time.time() - st))
    print("-" * 50)

    return h_1, vec_I, sigma_2


def loquat_sign_phase_3(Fp, Fp2, B, n, m, algebraic_hash, loquat_hash, loquat_expand, sigma_2, h_1, salt, listH,
                        ldt_lists, vec_I, poly_def, poly_c_prime, tree_cap):
    st = time.time()
    h_2, phase_2_challenge = get_phase_2_challenge(Fp, Fp2, B, n, algebraic_hash, loquat_hash, loquat_expand, sigma_2,
                                                   h_1, salt)

    # split phase 2 challenge first B elements are lambdas while the last loquat_n elements are epsilon
    vec_lambda, vec_epsilon = vector(phase_2_challenge[:B]), vector(phase_2_challenge[B:])
    vec_lambda_I = vec_lambda.pairwise_product(vec_I)

    poly_f_j, poly_q = ([] for _ in range(2))

    for j in range(n):
        vec_q_j = []
        for i in range(m):
            vec_q_j.append(Fp2(vec_lambda[i + j * m]))
            vec_q_j.append(Fp2(vec_lambda_I[i + j * m]))
        # interpolate the polynomial
        poly_q_j = poly_def.lagrange_polynomial(zip(listH, vec_q_j))
        poly_q.append(poly_q_j)
        poly_f_j.append(poly_q_j * poly_c_prime[j])

    poly_f = sum(Fp2(vec_epsilon[ind]) * poly_f_j[ind] for ind in range(n))

    # randomly samples 4m + kappa * 2**eta coefficients for polynomial poly_s
    poly_s = poly_def([Fp2.random_element() for _ in range(poly_f.degree())])

    # compute sum over listH of S
    value_S = sum(poly_s(item) for item in listH)
    eval_s = eval_poly_over_listU(poly_s, ldt_lists[0])
    leaf_s = get_leaf_hash(False, eval_s, algebraic_hash, loquat_hash, salt, Fp2)
    # get Merkle tree commitment over eval_list_U_s
    root_s, tree_s = MT_commit(algebraic_hash, loquat_hash, leaf_s, tree_cap, salt, Fp2)
    sigma_3 = [root_s, value_S]

    print("Signing phase 3 running time: {} seconds".format(time.time() - st))
    print("-" * 50)

    return sigma_3, tree_s, eval_s, h_2, poly_f, poly_s, vec_lambda, vec_epsilon


def loquat_sign_phase_4(sigma_3, algebraic_hash, loquat_hash, h_2, poly_f, poly_s, poly_Z_H, ldt_lists, tree_cap, salt,
                        Fp2):
    st = time.time()
    h_3 = loquat_hash_func(algebraic_hash, loquat_hash, sigma_3 + [h_2], salt, Fp2)
    # print("h_3", h_3)
    # z = h_3
    poly_f_prime = h_3[0] * poly_f + poly_s
    # compute poly_g and poly_h such that poly_f_prime = poly_g + poly_Z_H * poly_h
    poly_h, poly_g = poly_f_prime.quo_rem(poly_Z_H)
    # evaluate poly_h over listU
    eval_h = eval_poly_over_listU(poly_h, ldt_lists[0])
    leaf_h = get_leaf_hash(False, eval_h, algebraic_hash, loquat_hash, salt, Fp2)
    root_h, tree_h = MT_commit(algebraic_hash, loquat_hash, leaf_h, tree_cap, salt, Fp2)
    sigma_4 = [root_h]
    print("Signing phase 4 running time: {} seconds".format(time.time() - st))
    print("-" * 50)
    return sigma_4, tree_h, eval_h, poly_f_prime, h_3, poly_h


def loquat_sign_phase_5(Fp2, algebraic_hash, loquat_hash, loquat_expand, sigma_4, h_3, n, m, salt, sigma_2, vec_lambda,
                        vec_epsilon, x, poly_Z_H, poly_f_prime, poly_h, value_S, poly_def, ldt_lists, eval_c_prime, eta,
                        eval_s, eval_h, poly_c_prime, poly_s, rho_star):
    st = time.time()
    h_4, vec_e = get_phase_4_challenges(algebraic_hash, loquat_hash, loquat_expand, sigma_4 + [h_3], salt, Fp2)
    vec_lambda_o = vec_lambda.pairwise_product(sigma_2)
    mu = sum(vec_epsilon[j] * sum(vec_lambda_o[ind + j * m] for ind in range(m)) for j in range(n))

    # rational constraint poly_p
    size_H = poly_Z_H.degree()
    denominator_poly = poly_def(size_H * x)
    numerator_poly = size_H * poly_f_prime - size_H * poly_Z_H * poly_h - (h_3[0] * Fp2(mu) + value_S)
    poly_p, poly_p_rem = numerator_poly.quo_rem(denominator_poly)
    eval_p = eval_poly_over_listU(poly_p, ldt_lists[0])
    ldt_eta = 2 ** eta
    size_U = len(ldt_lists[0]) * ldt_eta
    # stack all codewords together as a (n + 3) * (size_U//2**ldt_eta) matrix
    matrix_pi = Matrix(Fp2, 0, size_U)

    for item in eval_c_prime:
        matrix_pi = matrix_pi.stack(Matrix(flatten(item)))
    matrix_pi = matrix_pi.stack(Matrix(flatten(eval_s)))
    matrix_pi = matrix_pi.stack(Matrix(flatten(eval_h)))
    matrix_pi = matrix_pi.stack(Matrix(flatten(eval_p)))
    flattened_list_U = flatten(ldt_lists[0])
    exp_list = [Fp2(int((rho_star - (poly_c_prime[j].degree() + 1) / size_U) * size_U)) for j in range(n)]
    exp_list.extend([Fp2(int((rho_star - (poly_s.degree() + 1) / size_U) * size_U)),
                     Fp2(int((rho_star - (poly_h.degree() + 1) / size_U) * size_U)),
                     Fp2(int((rho_star - (poly_p.degree() + 1) / size_U) * size_U))])

    # all_rows = []
    for ind in range(matrix_pi.nrows()):
        row = matrix_pi.row(ind)
        new_row = [flattened_list_U[ind_prime] ** exp_list[ind] * matrix_pi.row(ind)[ind_prime] for ind_prime in
                   range(len(row))]
        matrix_pi = matrix_pi.stack(Matrix(new_row))

    temp = (Matrix(vec_e) * matrix_pi).row(0)
    num_elements = len(ldt_lists[0])
    vec_f_0 = [list(temp[ind * ldt_eta: (ind + 1) * ldt_eta]) for ind in range(num_elements)]

    print("Signing phase 5 running time: {} seconds".format(time.time() - st))
    print("-" * 50)

    return vec_f_0, h_4


def loquat_sign_phase_6(vec_f_0, h_4, loquat_hash, algebraic_hash, Fp2, salt, tree_cap, ldt_lists, eta, poly_def,
                        rho_star, r):
    st = time.time()
    vec_f_i = vec_f_0.copy()
    h_i = h_4
    tree_f = []
    root_f = []
    coe_f_r = []
    vec_ldt = []
    ldt_eta = 2 ** eta
    for ind in range(r + 1):
        leaf_f_i = get_leaf_hash(False, vec_f_i, algebraic_hash, loquat_hash, salt, Fp2)
        # MT commit to vec_f_i
        root_f_i, tree_f_i = MT_commit(algebraic_hash, loquat_hash, leaf_f_i, tree_cap, salt, Fp2)
        sigma_i = root_f_i
        h_i_plus_1 = loquat_hash_func(algebraic_hash, loquat_hash, sigma_i, salt, Fp2)
        list_U_i = ldt_lists[ind]
        vec_f_i_plus_1 = []
        tree_f.append(tree_f_i)
        root_f.append(root_f_i)

        if ind == r:
            y_values = []
            for s in range(len(vec_f_i)):
                for my_eta in range(ldt_eta):
                    y_values.append(vec_f_i[s][my_eta])
            poly_f_r = poly_def.lagrange_polynomial(zip(ldt_lists[r], y_values))
            # print(poly_f_r)
            degree_f_r = int(rho_star * len(ldt_lists[r]) - 1)
            coe_f_r = poly_f_r.coefficients()[:degree_f_r + 1]
            # need to hash them with previous challenge h_i_plus_1
            sigma_i = []
            for item in coe_f_r:
                sigma_i.append(item)
            sigma_i.append(h_i_plus_1)
            h_i_plus_1 = loquat_hash_func(algebraic_hash, loquat_hash, sigma_i, salt, Fp2)
            h_i = h_i_plus_1
            break
        else:
            for y in range(0, len(list_U_i), 4):
                temp = []
                for ind_prime in range(ldt_eta):
                    poly_P_y = poly_def.lagrange_polynomial(zip(list_U_i[y + ind_prime], vec_f_i[y + ind_prime]))
                    temp.append(poly_P_y(h_i_plus_1))
                vec_f_i_plus_1.append(temp)

        vec_ldt.append([item for item in vec_f_i_plus_1])
        h_i = h_i_plus_1
        vec_f_i.clear()
        vec_f_i = vec_f_i_plus_1
    print("Signing phase 6 running time: {} seconds".format(time.time() - st))
    print("-" * 50)
    return tree_f, root_f, coe_f_r, vec_ldt, h_i


def loquat_sign_phase_7(algebraic_hash, loquat_expand, h_i, kappa, Fp2, n, eval_c_prime, eval_s, eval_h, tree_c,
                        tree_s, tree_h, ldt_lists, tree_cap, tree_f, r, vec_ldt):
    st = time.time()
    # expand h_i to get query sets
    # print(h_i)
    queries = get_ldt_query(algebraic_hash, loquat_expand, h_i[0], kappa, Fp2)
    print(queries)
    # get authentication paths for tree c, s and h

    query_leaf_c = [[eval_c_prime[j][item] for j in range(n)] for item in queries]
    query_leaf_s = [eval_s[item] for item in queries]
    query_leaf_h = [eval_h[item] for item in queries]
    auth_c = [MT_open(item, tree_c) for item in queries]
    auth_s = [MT_open(item, tree_s) for item in queries]
    auth_h = [MT_open(item, tree_h) for item in queries]

    # append additional final tree nodes for tree_c, s, h if it cannot be recomputed from MT
    final_tree_indices = set(range(len(tree_c[-1])))
    used_indices = {floor(item / 2 ** (log(len(ldt_lists[0]), 2) - tree_cap)) for item in queries}
    unused_indices = final_tree_indices - used_indices
    additional_node_us = {index: [tree_c[-1][index], tree_s[-1][index], tree_h[-1][index]] for index in unused_indices}

    auth_r, query_r, additional_node_ldt = [], [], {i: {} for i in range(r)}
    powers_of_4 = [4 ** i for i in range(r + 1)]

    # get authentication paths for LDT trees
    for i in range(r + 1):
        indices = set(floor(queries[ind] / powers_of_4[i]) for ind in range(len(queries)))
        auth_r.append([MT_open(index, tree_f[i]) for index in indices])
        if i > 0:
            query_r.append({index: vec_ldt[i - 1][index] for index in indices})

    # append additional final tree nodes for ldt trees if it cannot be recomputed from MT
    for i in range(r):
        final_indices = set(range(len(tree_f[i][-1])))
        used_indices = set(
            floor(floor(item / powers_of_4[i]) / 2 ** (log(len(ldt_lists[i]), 2) - tree_cap)) for item in queries)
        unused_indices = final_indices - used_indices
        if unused_indices:
            additional_node_ldt[i] = {index: tree_f[i][-1][index] for index in unused_indices}

    print("Signing phase 7 running time: {} seconds".format(time.time() - st))
    print("-" * 50)

    return query_leaf_c, query_leaf_s, query_leaf_h, auth_c, auth_s, auth_h, additional_node_us, auth_r, query_r, additional_node_ldt


def cal_sig_size(sig):
    total_bits = 0
    (sigma_1, sigma_2, sigma_3, sigma_4, root_f, auth_c, auth_s, auth_h, query_leaf_c, query_leaf_s,
     query_leaf_h, additional_node_us, coe_f_r, auth_r, query_r,
     additional_node_ldt, salt) = (sig[i] for i in ("sigma_1", "sigma_2", "sigma_3", "sigma_4", "root_f", "auth_c", "auth_s",
                                              "auth_h", "query_leaf_c", "query_leaf_s", "query_leaf_h",
                                              "additional_node_us", "coe_f_r", "auth_r", "query_r",
                                              "additional_node_ldt", "salt"))
    # convert signature elements to bits
    # bits for sigma_1
    total_bits += sigma_1[0][0].to_integer().bit_length()
    total_bits += sum(bin(element)[2:].count('1') + bin(element)[2:].count('0') for element in sigma_1[1])
    # bits for sigma_2
    for item in sigma_2:
        total_bits += item.lift().bit_length()
    # bits for sigma_3
    total_bits += sigma_3[0][0].to_integer().bit_length()
    total_bits += sigma_3[1].to_integer().bit_length()
    # bits for sigma_4
    total_bits += sigma_4[0][0].to_integer().bit_length()
    # bits for root_f
    for item in root_f:
        total_bits += item[0].to_integer().bit_length()
    # bits for authentication paths of us trees
    for ind in range(len(auth_c)):
        for ind_prime in range(len(auth_c[ind])):
            total_bits += auth_c[ind][ind_prime][0].to_integer().bit_length()
            total_bits += auth_s[ind][ind_prime][0].to_integer().bit_length()
            total_bits += auth_h[ind][ind_prime][0].to_integer().bit_length()
    # bits for authentication paths of ldt trees
    for item in auth_r:
        for ind in range(len(item)):
            # print(item[ind])
            total_bits += item[ind][0][0].to_integer().bit_length()
    # bits for coe_f_r
    for item in coe_f_r:
        total_bits += item.to_integer().bit_length()
    # bits for leaf nodes
    for item in query_leaf_h:
        for items in item:
            total_bits += items.to_integer().bit_length()
    for item in query_leaf_s:
        for items in item:
            total_bits += items.to_integer().bit_length()
    for item in query_leaf_c:
        for items in item:
            for item_prime in items:
                total_bits += item_prime.to_integer().bit_length()
    for items in query_r:
        for item in items.values():
            for i in item:
                total_bits += i.to_integer().bit_length()
    if len(additional_node_us.keys()) > 0:
        for key, value in additional_node_us.items():
            total_bits += key.bit_length()
            for item in value:
                total_bits += item[0].to_integer().bit_length()
        for key, value in additional_node_ldt.items():
            total_bits += key.bit_length()
            for key1, value1 in value.items():
                total_bits += key1.bit_length()
                for item1 in value1:
                    total_bits += item1.to_integer().bit_length()
    total_bits += salt.to_integer().bit_length()
    print("Signature size is (bits)", total_bits)
    print('-' * 50)


def loquat_sign(pp, sk, msg):
    st = time.time()
    (Fp, Fp2, p, listI, m, n, eta, kappa, rho_star, r, tree_cap, g, listH, ldt_lists,
     algebraic_hash, loquat_hash, loquat_expand) = (pp[key] for key in ("Fp", "Fp2", "p", "listI", "m",
                                                                        "n", "eta", "kappa", "rho_star", "r",
                                                                        "tree_cap", "g", "listH",
                                                                        "ldt_lists", "algebraic_hash", "loquat_hash",
                                                                        "loquat_expand"))
    poly_def = PolynomialRing(Fp2, 'x')
    x = var('x')
    size_H = len(listH)
    salt = Fp2.random_element()
    # print(salt)
    B = m * n
    # compute vanish polynomial for multiplicative coset H
    poly_Z_H = poly_def(x) ** size_H - g ** size_H

    eval_c_prime, tree_c, poly_c_prime, vec_r, sigma_1 = loquat_sign_phase_1(sk, n, m, p, Fp, Fp2, listH, g,
                                                                             eta,
                                                                             kappa, ldt_lists, algebraic_hash,
                                                                             loquat_hash, tree_cap, salt,
                                                                             poly_def, poly_Z_H)

    h_1, vec_I, sigma_2 = loquat_sign_phase_2(Fp2, Fp, loquat_hash, loquat_expand, algebraic_hash, sigma_1, salt, msg,
                                              B, vec_r,
                                              listI, sk)
    sigma_3, tree_s, eval_s, h_2, poly_f, poly_s, vec_lambda, vec_epsilon = loquat_sign_phase_3(Fp, Fp2, B, n, m,
                                                                                                algebraic_hash,
                                                                                                loquat_hash,
                                                                                                loquat_expand, sigma_2,
                                                                                                h_1, salt, listH,
                                                                                                ldt_lists, vec_I,
                                                                                                poly_def, poly_c_prime,
                                                                                                tree_cap)

    sigma_4, tree_h, eval_h, poly_f_prime, h_3, poly_h = loquat_sign_phase_4(sigma_3, algebraic_hash, loquat_hash, h_2,
                                                                             poly_f,
                                                                             poly_s, poly_Z_H, ldt_lists, tree_cap,
                                                                             salt, Fp2)

    vec_f_0, h_4 = loquat_sign_phase_5(Fp2, algebraic_hash, loquat_hash, loquat_expand, sigma_4, h_3, n, m, salt,
                                       sigma_2, vec_lambda,
                                       vec_epsilon, x, poly_Z_H, poly_f_prime, poly_h, sigma_3[1], poly_def, ldt_lists,
                                       eval_c_prime, eta, eval_s, eval_h, poly_c_prime, poly_s, rho_star)
    tree_f, root_f, coe_f_r, vec_ldt, h_i = loquat_sign_phase_6(vec_f_0, h_4, loquat_hash, algebraic_hash, Fp2, salt,
                                                                tree_cap, ldt_lists, eta, poly_def,
                                                                rho_star, r)

    query_leaf_c, query_leaf_s, query_leaf_h, auth_c, auth_s, auth_h, additional_node_us, auth_r, query_r, additional_node_ldt = loquat_sign_phase_7(
        algebraic_hash, loquat_expand, h_i, kappa, Fp2, n, eval_c_prime, eval_s, eval_h, tree_c,
        tree_s, tree_h, ldt_lists, tree_cap, tree_f, r, vec_ldt)

    sig = {
        "sigma_1": sigma_1,
        "sigma_2": sigma_2,
        "sigma_3": sigma_3,
        "sigma_4": sigma_4,
        "root_f": root_f,
        "auth_c": auth_c,
        "auth_s": auth_s,
        "auth_h": auth_h,
        # "T": sigma_1[1],
        "query_leaf_c": query_leaf_c,
        "query_leaf_s": query_leaf_s,
        "query_leaf_h": query_leaf_h,
        "additional_node_us": additional_node_us,
        "coe_f_r": coe_f_r,
        "auth_r": auth_r,
        "query_r": query_r,
        "additional_node_ldt": additional_node_ldt,
        "salt": salt
    }

    print("Sign running time: {} seconds".format(time.time() - st))
    print("-" * 50)

    cal_sig_size(sig)
    return sig
