import pickle
import sys

from griffin import griffin_sponge
from sage.all import *
from hashlib import sha3_256
from hashlib import shake_256


def eval_poly_over_listU(poly, ldt_list):
    return [[poly(item) for item in val] for val in ldt_list]


def sha3_256_func(Fp2, inputs):
    sha256 = sha3_256()
    if isinstance(inputs, bytes):
        sha256.update(inputs)
    elif isinstance(inputs, int):
        sha256.update(repr(inputs).encode('utf-8'))
    else:
        sha256.update(inputs.encode('utf-8'))
    byte_output = sha256.digest()
    return convert_bytes_to_Fp2(Fp2, byte_output)


def shake_256_func(input_value, output_bytes):
    shake = shake_256()
    if isinstance(input_value, bytes):
        shake.update(input_value)
    elif isinstance(input_value, int):
        shake.update(repr(input_value).encode('utf-8'))
    else:
        shake.update(input_value.encode('utf-8'))
    return shake.digest(output_bytes)


def loquat_hash_func(algebraic_hash, loquat_hash, inputs, salt, Fp2):
    hash_input = sum(convert_to_Fp2_list(inputs))
    if algebraic_hash:  # griffin hash
        output = griffin_sponge(loquat_hash, [hash_input] + [salt], 1)
        return output
    else:  # standard hash
        # convert input to byte string
        byte_input = convert_to_bytes(hash_input)
        byte_salt = convert_to_bytes(salt)
        return [sha3_256_func(Fp2, byte_input + byte_salt)]


def loquat_expand_func(byte_output, algebraic_hash, loquat_expand, input_seed, output_length, Fp2):
    if algebraic_hash:
        output = griffin_sponge(loquat_expand, input_seed, output_length)
        if byte_output:
            return convert_to_bytes(output)
        else:
            return output
    else:
        byte_seed = convert_to_bytes(input_seed)
        output_byte = shake_256_func(byte_seed, output_length * 32)
        if byte_output:
            return output_byte
        else:
            # convert output bytes to Fp2 elements
            return [convert_bytes_to_Fp2(Fp2, output_byte[i * 32: (i + 1) * 32]) for i in range(output_length)]


def convert_to_Fp2_list(inputs):
    if isinstance(inputs, list):
        return [element for sublist in inputs for element in convert_to_Fp2_list(sublist)]
    else:
        return [inputs]


def convert_to_bytes(inputs):
    if isinstance(inputs, list):
        return b"".join(convert_to_bytes(item) for item in inputs)
    else:
        return pickle.dumps(inputs)


# The verifier should be able to recompute the leaf hash outside of the circuit and
# append the resulting leaf node to SNARK witness
def get_leaf_hash(witness_poly, eval_poly, algebraic_hash, loquat_hash, salt, Fp2):
    hashed_leaf = []
    if witness_poly:  # the case of witness polynomial
        list_length = len(eval_poly[0])
        for i in range(list_length):
            leaf_vec = [poly[i] for poly in eval_poly]
            hashed_leaf.append(loquat_hash_func(algebraic_hash, loquat_hash, leaf_vec, salt, Fp2))
    else:
        hashed_leaf = [loquat_hash_func(algebraic_hash, loquat_hash, poly, salt, Fp2) for poly in eval_poly]
    return hashed_leaf


def MT_commit(algebraic_hash, loquat_hash, leaf, tree_cap, salt, Fp2):
    tree = [leaf]
    # tree_size = len(leaf)
    num_layers = int(log(len(leaf), 2)) - tree_cap
    for layer in range(num_layers):
        new_layer = [
            loquat_hash_func(algebraic_hash, loquat_hash, leaf[ind] + leaf[ind + 1], salt, Fp2) for ind in
            range(0, len(leaf), 2)
        ]
        leaf = new_layer
        tree.append(leaf)
    root = loquat_hash_func(algebraic_hash, loquat_hash, [h[0] for h in leaf], salt, Fp2)
    return root, tree


def convert_bytes_to_Fp2(Fp2, data):
    if len(data) != 32:
        print("Length error when converting bytes to Fp2 element")
        sys.exit()
    else:
        byte_part1 = int.from_bytes(data[:16], byteorder='big')
        byte_part2 = int.from_bytes(data[16:], byteorder='big')
        # print(Fp2)
        return Fp2([byte_part2, byte_part1])


def get_phase_1_challenge(Fp2, loquat_hash, loquat_expand, algebraic_hash, sigma_1, salt, msg, B):
    # lift bit string T to Fp2 element
    field_T = Fp2(int(''.join(map(str, sigma_1[1])), 2))
    # hash the message
    hashed_message = sha3_256_func(Fp2, msg)
    h_1 = loquat_hash_func(algebraic_hash, loquat_hash, sigma_1[0] + [field_T] + [hashed_message], salt, Fp2)
    output = loquat_expand_func(True, algebraic_hash, loquat_expand, h_1, 8, Fp2)
    # convert bytes to bits
    bits = ''.join(format(byte, '08b') for byte in output)
    # get integer with 15 bits
    result = [int(bits[i:i + 15], 2) for i in range(B)]
    return h_1, result


def get_phase_2_challenge(Fp, Fp2, B, n, algebraic_hash, loquat_hash, loquat_expand, sigma_2, h_1, salt):
    a = Fp2.gen()
    # converts sigma_2 elements to Fp_2
    input_list = [Fp2(sigma_2[i] + sigma_2[i + 1] * a) for i in range(0, len(sigma_2), 2)]
    input_list.append(h_1)

    # hash h_1 || sigma_2
    h_2 = loquat_hash_func(algebraic_hash, loquat_hash, input_list, salt, Fp2)

    # expand h_2 with one padding to 0
    output = loquat_expand_func(False, algebraic_hash, loquat_expand, h_2, (B + n) // 2, Fp2)
    result = [Fp(coeff) for item in output for coeff in item.polynomial().coefficients(sparse=False)]

    return h_2, result


def get_phase_4_challenges(algebraic_hash, loquat_hash, loquat_expand, inputs, salt, Fp2):
    h_4 = loquat_hash_func(algebraic_hash, loquat_hash, inputs, salt, Fp2)
    result = loquat_expand_func(False, algebraic_hash, loquat_expand, h_4, 14, Fp2)
    return h_4, result


def get_ldt_query(algebraich_hash, loquat_expand, h_i, kappa, Fp2):
    output_seq = ceil(10 * kappa / 254)
    output = loquat_expand_func(True, algebraich_hash, loquat_expand, h_i, output_seq, Fp2)
    # convert the output to bit string
    bit_string = ''.join(format(byte, '08b') for byte in output)
    bit_length = 10
    result = [int(bit_string[i:i + bit_length], 2) for i in range(0, kappa * bit_length, bit_length)]

    return result


def MT_open(index, merkle_tree):
    path = []
    if index % 2 == 0:
        path.append(merkle_tree[0][index + 1])
    else:
        path.append(merkle_tree[0][index - 1])

    for i in range(1, len(merkle_tree) - 1):
        index = floor(index / 2)
        if index % 2 == 0:
            path.append(merkle_tree[i][index + 1])
        else:
            path.append(merkle_tree[i][index - 1])
    return path