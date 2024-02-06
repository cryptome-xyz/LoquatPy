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
    queries = get_ldt_query(algebraic_hash, loquat_expand, h_i[0], kappa, Fp2)

    return phase_1_challenge, phase_2_challenge, vec_e, queries, h_3


def loquat_verify_step_2(Fp2, phase_2_challenge, phase_1_challenge, listI, m, n, listH, query_leaf_c, kappa, eta,
                         ldt_lists, queries, poly_Z_H, h_3, sigma_2, sigma_3):
    vec_lambda = phase_2_challenge[:B]
    vec_epsilon = phase_2_challenge[B:]
    vec_I = [listI[ind] for ind in phase_1_challenge]  # challenges contain B random selected elements from listI
    poly_q = []
    for j in range(n):
        vec_q_j = []
        vec_q_j.extend(
            [Fp_2(vec_lambda[ind + j * m]), Fp_2(vec_lambda[ind + j * m] * vec_I[ind + j * m])] for ind in range(m))
        # interpolate the polynomial
        poly_q_j = poly_def.lagrange_polynomial(zip(listH, vec_q_j))
        poly_q.append(poly_q_j)
    vec_f_j = []
    ldt_eta = 2 ** eta
    size_H = len(listH)
    for ind in range(kappa):
        temp = []
        for j in range(ldt_eta):
            temp_list = []
            for j_prime in range(n):
                temp_list.append(poly_q[j_prime](ldt_lists[0][queries[ind]][j]) * query_leaf_c[ind][j_prime][j])
            temp.append(temp_list)
        vec_f_j.append(temp)
    vec_f = []
    for item in vec_f_j:
        temp_list = []
        for j_prime in range(ldt_eta):
            temp = Fp2(0)
            for j in range(n):
                temp = temp + Fp_2(vec_epsilon[j]) * item[j_prime][j]
            temp_list.append(temp)
        vec_f.append(temp_list)
    z = h_3
    vec_p = []
    mu = 0
    # recompute the claimed sumcheck value mu
    for j in range(n):
        value = vec_epsilon[j]
        temp = 0
        for ind in range(m):
            temp = temp + (vec_lambda[ind + j * m] * sigma_2[ind + j * m])
        mu = mu + (temp * value)
    S = sigma_3[1]

    for ind in range(kappa):
        temp = []
        for ind_prime in range(ldt_eta):
            val_Z_H = poly_Z_H(ldt_lists[0][queries[ind]][ind_prime])
            denominator_value = size_H * ldt_lists[0][queries[ind]][ind_prime]
            numerator_value = size_H * (z * vec_f[ind][ind_prime] + query_leaf_s[ind][ind_prime]) - size_H * val_Z_H * \
                              query_leaf_h[ind][ind_prime] - (z * mu + S)
            temp.append(numerator_value / denominator_value)
        vec_p.append(temp)

    # reconstruct interleaved code matrix
    matrix_pi = Matrix(Fp_2, 0, kappa * ldt_eta)
    for j in range(n):
        temp = []
        for ind in range(kappa):
            for ind_prime in range(ldt_eta):
                temp.append(query_leaf_c[ind][j][ind_prime])
        matrix_pi = matrix_pi.stack(Matrix(temp))

    big_list_query = []
    for ind in range(kappa):
        for item in ldt_lists[0][queries[ind]]:
            big_list_query.append(item)

    for ind in range(3):
        temp = []
        for ind_prime in range(kappa):
            for ind_prime_prime in range(2 ** ldt_eta):
                if ind == 0:
                    temp.append(query_leaf_s[ind_prime][ind_prime_prime])
                elif ind == 1:
                    temp.append(query_leaf_h[ind_prime][ind_prime_prime])
                else:
                    temp.append(vec_p[ind_prime][ind_prime_prime])
        matrix_pi = matrix_pi.stack(Matrix(temp))

    exp_list = []
    for j in range(loquat_n):
        exp_list.append(Fp_2(int((rho_star - (2 * m + kappa * 2 ** ldt_eta + 1) / size_U) * size_U)))
    exp_list.append(Fp_2(int((rho_star - (4 * m + kappa * 2 ** ldt_eta - 1) / size_U) * size_U)))
    exp_list.append(Fp_2(int((rho_star - (2 * m + kappa * 2 ** ldt_eta) / size_U) * size_U)))
    exp_list.append(Fp_2(int((rho_star - (2 * m - 1) / size_U) * size_U)))

    num_rows = matrix_pi.nrows()
    for ind in range(num_rows):
        temp = []
        for ind_prime in range(len(matrix_pi.row(ind))):
            temp.append(big_list_query[ind_prime] ** exp_list[ind] * matrix_pi.row(ind)[ind_prime])
        matrix_pi = matrix_pi.stack(Matrix(temp))

    matrix_e = Matrix(vec_e)
    temp = (matrix_e * matrix_pi).row(0)
    vec_f_0 = []
    for ind in range(kappa):
        temp_vec = []
        for ind_prime in range(2 ** ldt_eta):
            temp_vec.append(temp[ind_prime + ind * (2 ** ldt_eta)])
        vec_f_0.append(temp_vec)


def loquat_verify(pp, pk, msg, sig):
    st = time.time()
    result = 1
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
    poly_Z_H = poly_def(x) ** size_H - g ** size_H

    phase_1_challenge, phase_2_challenge, vec_e, queries, h_3 = loquat_verify_step_1(Fp, Fp2, loquat_hash,
                                                                                     loquat_expand, algebraic_hash,
                                                                                     sigma_1, salt, msg, B, n, sigma_2,
                                                                                     sigma_3, root_f, coe_f_r, kappa,
                                                                                     sigma_4, r)

    return result
KeyGen.py