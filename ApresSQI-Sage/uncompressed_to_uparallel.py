from sage.all import *
from ApresSQI import SQIsign
from xonly import sqrtNonConstTime, customTwoIsogeny
from sys import argv

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print("Usage: sage uncompressed_to_smart.sage uncompressed smart")
        sys.exit(1)
    param = sys.argv[1]
    assert param == 'p248'
    in_variant = sys.argv[2]
    assert in_variant in ['uncompressed']
    out_variant = sys.argv[3]
    assert in_variant != out_variant
    assert out_variant in ['uparallel']

    in_filename = f'KATs/PQCsignKAT_{param}_lvl1_{in_variant}.rsp'
    out_filename = f'KATs/PQCsignKAT_{param}_lvl1_{out_variant}.rsp'

    in_file = open(in_filename, "r")
    out_file = open(out_filename, "w")
    print(f'# lvl1_{param}\n', file=out_file)

    topstr = f'############# Using {param} prime #############'
    print(topstr)

    ZIP_CHAIN_LEN = 4
    ZIP_CHAIN_HALFBYTES = 2*31
    FP2_ENCODED_BYTES = 64
    TORSION_2POWER_SECPAR_HALFBYTES = 2*16
    SIGNATURE_LEN = 320

    def pk_encode(pk):
        E_A = pk.E_A
        pk = E_A.a_invariants()[1]
        pk_bytes = f'{int(pk[0]):0{FP2_ENCODED_BYTES}X}'
        pk_bytes = f'{int(pk[1]):0{FP2_ENCODED_BYTES}X}' + pk_bytes
        pk_bytes = bytearray.fromhex(pk_bytes)[::-1]

        return pk_bytes.hex().upper()

    def pk_decode(signer, pk_bytes):
        pk = [pk_bytes[i:i+FP2_ENCODED_BYTES] for i in range(0, len(pk_bytes), FP2_ENCODED_BYTES)]
        pk_re = bytearray.fromhex(pk[0])[::-1]
        pk_re = Integer('0x' + pk_re.hex())
        pk_im = bytearray.fromhex(pk[1])[::-1]
        pk_im = Integer('0x' + pk_im.hex())
        A = signer.F([pk_re, pk_im])
        E_A = EllipticCurve(signer.F, [0,A,0,1,0])
        return E_A

    complete = lambda input: '0'*(len(input) % 2) + input
    def sig_decode_uncompressed(signer, sig):
        siglen = 2*(ZIP_CHAIN_LEN+1)*FP2_ENCODED_BYTES
        assert(len(sig) == siglen)
        gens_E1 = [sig[i:i + 2*FP2_ENCODED_BYTES] for i in range(0, len(sig), 2*FP2_ENCODED_BYTES)]
        GENS_E1 = []
        for elt_bytes in gens_E1:
            elt = [elt_bytes[i:i + FP2_ENCODED_BYTES] for i in range(0, len(elt_bytes), FP2_ENCODED_BYTES)]
            elt_re = bytearray.fromhex(elt[0])[::-1]
            elt_re = Integer('0x' + elt_re.hex())
            elt_im = bytearray.fromhex(elt[1])[::-1]
            elt_im = Integer('0x' + elt_im.hex())
            elt = signer.F([elt_re, elt_im])
            GENS_E1.append(elt)

        return GENS_E1

    def sig_msg_encode_uparallel(sigma, msg):
        sigma_bytes = ''
        for i in range(len(sigma)):
            gens_i_re = f'{int(sigma[i][0]):0{FP2_ENCODED_BYTES}X}'
            gens_i_im = f'{int(sigma[i][1]):0{FP2_ENCODED_BYTES}X}'
            sigma_bytes = gens_i_im + gens_i_re + sigma_bytes
        sigma_msg_bytes = bytearray.fromhex(sigma_bytes)[::-1] + bytearray.fromhex(msg)

        return sigma_msg_bytes.hex().upper()

    signer = SQIsign(param)
    signer.load_privkey(param=param)

    count = 0
    for line in in_file:
        if line.startswith("pk = "):
            pk_bytes = line.split("pk = ")[1].strip().replace("\n","")
            E_A = pk_decode(signer, pk_bytes)
            assert(signer.pk.E_A == E_A)
        elif line.startswith("smlen = "):
            smlen = int(line.split("smlen = ")[1].strip().replace("\n",""))
        elif line.startswith("sm = "):
            sm = line.split("sm = ")[1].strip().replace("\n","")
            sig = sm[:2*SIGNATURE_LEN]
            GENS_E1 = sig_decode_uncompressed(signer, sig)
            msg = sm[2*SIGNATURE_LEN:]
            mlen = (smlen - SIGNATURE_LEN)
            assert(len(msg) == 2*mlen)

            print(f'count = {count} -- converting {in_variant} to {out_variant} for message: \n{msg}')
            sigma = signer.getUparallel_from_uncompressed(GENS_E1)

            print('#'*len(topstr) + '\n')
            print(f'Done! signature was {sigma}' + '\n')

            print("Verifying...")
            verified = signer.verify(msg, sigma, signer.pk)
            print('Signature was ' + ('CORRECT' if verified else 'WRONG'))
            print('\n\n')

            pk_bytes = pk_encode(signer.pk)
            sigma_msg_bytes = sig_msg_encode_uparallel(sigma, msg)
            print(f'count = {count}', file=out_file)
            print('seed =', file=out_file)
            print(f'mlen = {mlen}', file=out_file)
            print(f'msg = {msg}', file=out_file)
            print(f'pk = {pk_bytes}', file=out_file)
            print('sk =', file=out_file)
            print(f'smlen = {len(sigma_msg_bytes)//2}', file=out_file)
            print(f'sm = {sigma_msg_bytes}\n', file=out_file)

            count += 1
