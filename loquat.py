from Setup import loquat_setup
from KeyGen import loquat_keygen
from Sign import loquat_sign
from Verify import loquat_verify
def main():
    pp = loquat_setup(0, 0, 100, 32)
    sk, pk = loquat_keygen(pp)
    msg = "Hello World"
    sig = loquat_sign(pp, sk, msg)
    # result = loquat_verify(pp, pk, msg, sig)

if __name__ == '__main__':
    main()