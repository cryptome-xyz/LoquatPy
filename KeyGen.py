import time
from sage.all import kronecker
from sage.all import floor

def loquat_keygen(pp):
    st = time.time()
    Fp = pp.get("Fp")
    listI = pp.get("listI")
    p = pp.get("p")

    sk = Fp.random_element()

    pk = [floor(1 / 2 * (1 - kronecker(sk + item, p))) for item in listI]
    
    print("Keygen running time: {} seconds".format(time.time() - st))
    print("-" * 50)
    return sk, pk
