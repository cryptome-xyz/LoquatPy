import sys

from Setup import loquat_setup
from KeyGen import loquat_keygen
from Sign import loquat_sign
from Verify import loquat_verify
def main():
    # loquat_star = 0/1; 1 for stronger security
    # algebraic_hash = 0/1; 1 use griffin, 0 use SHA
    # ldt_choice: 100/128 100 for 100-bit security ldt
    pp = loquat_setup(1, 1, 100)
    sk, pk = loquat_keygen(pp)
    msg = "Hello World"
    test_num = 1
    for _ in range(test_num):
        sig = loquat_sign(pp, sk, msg)
        result = loquat_verify(pp, pk, msg, sig)
        if result == 1:
            print("Accept: Valid signature")
        else:
            print("Reject: Invalid Signature")
            sys.exit()


if __name__ == '__main__':
    main()