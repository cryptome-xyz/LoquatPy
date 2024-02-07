from Utils import *
import time


def loquat_verify_step_1(Fp, Fp2, loquat_hash, loquat_expand, algebraic_hash, sigma_1, salt, msg, B, n, sigma_2,
                         sigma_3, root_f, coe_f_r, kappa, sigma_4, r):
    h_1, phase_1_challenge = get_phase_1_challenge(Fp2, loquat_hash, loquat_expand, algebraic_hash, sigma_1, salt, msg,
                                                   B)
    h_2, phase_2_challenge = get_phase_2_challenge(Fp, Fp2, B, n, algebraic_hash, loquat_hash, loquat_expand, sigma_2,
                                                   h_1, salt)

    h_3 = loquat_hash_func(algebraic_hash, loquat_hash, sigma_3 + [h_2], salt, Fp2)
    h_4, vec_e = get_phase_4_challenges(algebraic_hash, loquat_hash, loquat_expand, sigma_4 + [h_3], salt, Fp2)
    h_i = h_4
    h_i_list = []
    for ind in range(r + 1):
        sigma_i = [root_f[ind]]
        h_i_plus_1 = loquat_hash_func(algebraic_hash, loquat_hash, sigma_i, salt, Fp2)
        h_i_list.append(h_i_plus_1)
        # print("h_i", h_i)
        if ind == r:
            sigma_i = []
            for item in coe_f_r:
                sigma_i.append(item)
            sigma_i.append(h_i_plus_1)
            h_i_plus_1 = loquat_hash_func(algebraic_hash, loquat_hash, sigma_i, salt, Fp2)
            h_i = h_i_plus_1
            break
        h_i = h_i_plus_1
    queries = get_ldt_query(algebraic_hash, loquat_expand, [h_i[0]], kappa, Fp2)

    return phase_1_challenge, phase_2_challenge, vec_e, queries, h_3, h_i_list


def loquat_verify_step_2(Fp2, phase_2_challenge, phase_1_challenge, listI, m, n, listH, query_leaf_c, kappa, eta,
                         ldt_lists, queries, poly_Z_H, h_3, sigma_2, sigma_3, B, poly_def, query_leaf_s, query_leaf_h,
                         rho_star, vec_e):
    vec_lambda = vector(phase_2_challenge[:B])
    vec_epsilon = vector(phase_2_challenge[B:])
    ldt_eta = 2 ** eta
    size_U = len(ldt_lists[0]) * ldt_eta
    size_H = len(listH)
    vec_I = vector(
        [listI[ind] for ind in phase_1_challenge])  # challenges contain B random selected elements from listI
    # poly_q = []
    vec_lambda_I = vec_lambda.pairwise_product(vec_I)
    listU = ldt_lists[0]
    # query list
    query_list = [listU[queries[ind]] for ind in range(len(queries))]
    value_f = []
    for j in range(n):
        vec_q_j = []
        value_f_j = []
        for i in range(m):
            vec_q_j.append(Fp2(vec_lambda[i + j * m]))
            vec_q_j.append(Fp2(vec_lambda_I[i + j * m]))
        poly_q_j = poly_def.lagrange_polynomial(zip(listH, vec_q_j))
        for ind in range(kappa):
            value_f_j_ind = []
            for ind_prime in range(ldt_eta):
                value_f_j_ind.append(poly_q_j(query_list[ind][ind_prime]) * query_leaf_c[j][ind][ind_prime])
            value_f_j.append(value_f_j_ind)
        value_f.append(value_f_j)

    eval_f = []
    for i in range(kappa):
        temp_list = []
        for j in range(ldt_eta):
            temp_list.append(sum(Fp2(vec_epsilon[ind]) * value_f[ind][i][j] for ind in range(n)))
        eval_f.append(temp_list)
    vec_lambda_o = vec_lambda.pairwise_product(sigma_2)
    mu = sum(vec_epsilon[j] * sum(vec_lambda_o[ind + j * m] for ind in range(m)) for j in range(n))
    S = sigma_3[1]
    z = h_3[0]
    # val_Z_H = [[poly_Z_H(item) for item in val] for val in query_list]
    vec_p = []
    for ind in range(kappa):
        temp = []
        for ind_prime in range(ldt_eta):
            val_Z_H = poly_Z_H(query_list[ind][ind_prime])
            denominator_value = size_H * query_list[ind][ind_prime]
            numerator_value = size_H * (z * eval_f[ind][ind_prime] + query_leaf_s[ind][ind_prime]) - size_H * val_Z_H * \
                              query_leaf_h[ind][ind_prime] - (z * Fp2(mu) + S)
            # print(numerator_value)
            # print(denominator_value)
            temp.append(numerator_value / denominator_value)
        vec_p.append(temp)
    # reconstruct interleaved code matrix
    matrix_pi = Matrix(Fp2, 0, kappa * ldt_eta)
    for item in query_leaf_c:
        matrix_pi = matrix_pi.stack(Matrix(flatten(item)))
    matrix_pi = matrix_pi.stack(Matrix(flatten(query_leaf_s)))
    matrix_pi = matrix_pi.stack(Matrix(flatten(query_leaf_h)))
    matrix_pi = matrix_pi.stack(Matrix(flatten(vec_p)))
    flattened_query = flatten(query_list)
    exp_list = []
    for j in range(n):
        exp_list.append(Fp2(int((rho_star - (2 * m + kappa * ldt_eta + 1) / size_U) * size_U)))
    exp_list.append(Fp2(int((rho_star - (4 * m + kappa * ldt_eta - 1) / size_U) * size_U)))
    exp_list.append(Fp2(int((rho_star - (2 * m + kappa * ldt_eta) / size_U) * size_U)))
    exp_list.append(Fp2(int((rho_star - (2 * m - 1) / size_U) * size_U)))

    # print(exp_list)

    for ind in range(matrix_pi.nrows()):
        row = matrix_pi.row(ind)
        new_row = [flattened_query[ind_prime] ** exp_list[ind] * matrix_pi.row(ind)[ind_prime] for ind_prime in
                   range(len(row))]
        matrix_pi = matrix_pi.stack(Matrix(new_row))

    '''for ind in range(num_rows):
        temp = []
        for ind_prime in range(len(matrix_pi.row(ind))):
            temp.append(flattened_query[ind_prime] ** exp_list[ind] * matrix_pi.row(ind)[ind_prime])
        matrix_pi = matrix_pi.stack(Matrix(temp))'''

    matrix_e = Matrix(vec_e)
    temp = (matrix_e * matrix_pi).row(0)
    # print(temp)
    vec_f_0 = []
    for ind in range(kappa):
        temp_vec = []
        for ind_prime in range(ldt_eta):
            temp_vec.append(temp[ind_prime + ind * ldt_eta])
        vec_f_0.append(temp_vec)
    return vec_f_0


def loquat_verify_step_3(phase_1_challenge, sigma_2, sigma_1, p, pk, sigma_3, sigma_4, ldt_lists,
                         query_leaf_h, query_leaf_s, tree_cap, auth_h, auth_s, auth_c, queries, additional_node_us, eta,
                         query_leaf_c, algebraic_hash, rho_star, poly_def, coe_f_r, loquat_hash, salt, Fp2, r, auth_r,
                         vec_f_0, additional_node_ldt, root_f, query_r, h_i_list):
    # check the response of Legendre symbols
    for ind in range(len(phase_1_challenge)):
        if floor(1 / 2 * (1 - kronecker(sigma_2[ind], p))) != ((pk[phase_1_challenge[ind]] + sigma_1[1][ind]) % 2):
            return 0
    # check Merkle tree
    ldt_eta = 2 ** eta
    leaf_c = get_leaf_hash(True, query_leaf_c, algebraic_hash, loquat_hash, salt, Fp2)
    leaf_s = get_leaf_hash(False, query_leaf_s, algebraic_hash, loquat_hash, salt, Fp2)
    leaf_h = get_leaf_hash(False, query_leaf_h, algebraic_hash, loquat_hash, salt, Fp2)
    additional_node_c = {}
    additional_node_s = {}
    additional_node_h = {}
    if len(additional_node_us.keys()) > 0:
        for key, value in additional_node_us.items():
            additional_node_c[key] = additional_node_us[key][0]
            additional_node_s[key] = additional_node_us[key][1]
            additional_node_h[key] = additional_node_us[key][2]
    # print(additional_node_c)
    root_c_prime = MT_recompute_root(algebraic_hash, loquat_hash, queries, auth_c, leaf_c, tree_cap,
                                     additional_node_c, ldt_lists, 0, salt, Fp2)

    root_s_prime = MT_recompute_root(algebraic_hash, loquat_hash, queries, auth_s, leaf_s, tree_cap,
                                     additional_node_s, ldt_lists, 0, salt, Fp2)
    root_h_prime = MT_recompute_root(algebraic_hash, loquat_hash, queries, auth_h, leaf_h, tree_cap,
                                     additional_node_h, ldt_lists, 0, salt, Fp2)

    if root_c_prime == 0 or root_c_prime != sigma_1[0] or root_s_prime != sigma_3[
        0] or root_s_prime == 0 or root_h_prime != sigma_4[0] or root_h_prime == 0:
        return 0
    # print(vec_f_0[ind] for ind in range(len(queries)))
    for i in range(r + 1):
        # check LDT Merkle tree
        index = []
        if i == 0:
            leaf_f_0 = get_leaf_hash(False, vec_f_0, algebraic_hash, loquat_hash, salt, Fp2)
            root_f_i_prime = MT_recompute_root(algebraic_hash, loquat_hash, queries, auth_r[0], leaf_f_0, tree_cap,
                                               additional_node_ldt[0], ldt_lists, 0, salt, Fp2)
            if root_f_i_prime != root_f[i] or root_f_i_prime == 0:
                print("MT f0 error")
                return 0
        elif 0 < i < r:
            for ind in range(len(queries)):
                if floor(queries[ind] / 4 ** i) not in index:
                    index.append(floor(queries[ind] / 4 ** i))
            vec_f_i = query_r[i - 1].values()
            leaf_f_i = get_leaf_hash(False, vec_f_i, algebraic_hash, loquat_hash, salt, Fp2)
            root_f_i_prime = MT_recompute_root(algebraic_hash, loquat_hash, index, auth_r[i], leaf_f_i,
                                               tree_cap, additional_node_ldt[i], ldt_lists, i, salt, Fp2)
            if root_f_i_prime != root_f[i] or root_f_i_prime == 0:
                print("MT sub-vector error")
                return 0
        else:
            vec_f_i = list(query_r[r - 1].values())
            index = []
            degree_f_r = int(rho_star * len(ldt_lists[r]) - 1)
            for ind in range(len(queries)):
                if floor(queries[ind] / 4 ** i) not in index:
                    index.append(floor(queries[ind] / 4 ** i))
            for i in range(len(index)):
                x_values = ldt_lists[r][index[i]]
                y_values = vec_f_i[i]
                poly_f_r = poly_def.lagrange_polynomial(zip(x_values, y_values))
                if poly_f_r.degree() > degree_f_r:
                    print("Reject: Invalid LDT answer with degree of polynomial f_r greater than", degree_f_r)
                    return 0
                coe_f_r_prime = poly_f_r.coefficients()[:degree_f_r + 1]
                for i in range(len(coe_f_r_prime)):
                    if coe_f_r_prime[i] != coe_f_r[i]:
                        print("Reject: Invalid LDT coefficieints at round r =", r)
                        return 0
            break
        # check ldt consistency
        for ind in range(len(index)):
            if i == 0:
                y_values = vec_f_0[ind]
            elif i < r and i != 0:
                y_values = list(query_r[i - 1].values())[ind]
            x_values = ldt_lists[i][index[ind]]
            poly_P_y = poly_def.lagrange_polynomial(zip(x_values, y_values))
            value = poly_P_y(h_i_list[i])
            next_query = query_r[i]
            if floor(index[ind] / ldt_eta) in next_query.keys():
                if value not in next_query.get(floor(index[ind] / ldt_eta)):
                    print("ldt consistency error")
                    return 0
            else:
                print("ldt consistency error")
                return 0

    return 1


def loquat_verify(pp, pk, msg, sig):
    st = time.time()
    (Fp, Fp2, p, listI, m, n, eta, kappa, rho_star, r, tree_cap, g, listH, ldt_lists,
     algebraic_hash, loquat_hash, loquat_expand) = (pp[key] for key in ("Fp", "Fp2", "p", "listI", "m",
                                                                        "n", "eta", "kappa", "rho_star", "r",
                                                                        "tree_cap", "g", "listH",
                                                                        "ldt_lists", "algebraic_hash", "loquat_hash",
                                                                        "loquat_expand"))
    (sigma_1, sigma_2, sigma_3, sigma_4, root_f, auth_c, auth_s, auth_h, query_leaf_c, query_leaf_s,
     query_leaf_h, additional_node_us, coe_f_r, auth_r, query_r,
     additional_node_ldt, salt) = (sig[i] for i in
                                   ("sigma_1", "sigma_2", "sigma_3", "sigma_4", "root_f", "auth_c", "auth_s",
                                    "auth_h", "query_leaf_c", "query_leaf_s", "query_leaf_h",
                                    "additional_node_us", "coe_f_r", "auth_r", "query_r",
                                    "additional_node_ldt", "salt"))
    B = n * m
    size_H = len(listH)
    poly_def = PolynomialRing(Fp2, 'x')
    x = var('x')
    poly_Z_H = poly_def(x) ** size_H - g ** size_H

    phase_1_challenge, phase_2_challenge, vec_e, queries, h_3, h_i_list = loquat_verify_step_1(Fp, Fp2, loquat_hash,
                                                                                               loquat_expand,
                                                                                               algebraic_hash,
                                                                                               sigma_1, salt, msg, B, n,
                                                                                               sigma_2,
                                                                                               sigma_3, root_f, coe_f_r,
                                                                                               kappa,
                                                                                               sigma_4, r)

    vec_f_0 = loquat_verify_step_2(Fp2, phase_2_challenge, phase_1_challenge, listI, m, n, listH, query_leaf_c, kappa,
                                   eta,
                                   ldt_lists, queries, poly_Z_H, h_3, sigma_2, sigma_3, B, poly_def, query_leaf_s,
                                   query_leaf_h,
                                   rho_star, vec_e)

    result = loquat_verify_step_3(phase_1_challenge, sigma_2, sigma_1, p, pk, sigma_3, sigma_4, ldt_lists,
                                  query_leaf_h, query_leaf_s, tree_cap, auth_h, auth_s, auth_c, queries,
                                  additional_node_us, eta,
                                  query_leaf_c, algebraic_hash, rho_star, poly_def, coe_f_r, loquat_hash, salt, Fp2, r,
                                  auth_r,
                                  vec_f_0, additional_node_ldt, root_f, query_r, h_i_list)

    print("Verify running time: {} seconds".format(time.time() - st))
    print("-" * 50)

    return result
