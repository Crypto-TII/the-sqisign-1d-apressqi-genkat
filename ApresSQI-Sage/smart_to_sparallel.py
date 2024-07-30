from sage.all import *
from ApresSQI import SQIsign
from xonly import sqrtNonConstTime
from sys import argv

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print("Usage: sage uncompressed_to_smart.sage uncompressed smart")
        sys.exit(1)
    param = sys.argv[1]
    assert param == 'p248'
    in_variant = sys.argv[2]
    assert in_variant in ['smart']
    out_variant = sys.argv[3]
    assert in_variant != out_variant
    assert out_variant in ['sparallel']

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
    SIGNATURE_LEN = 157

    def pk_encode(pk):
        E_A = pk.E_A
        A = E_A.a_invariants()[1]
        delta = A**2 - 4
        sdelta = sqrtNonConstTime(delta)
        alpha = (-A-sdelta)/2
        pk = alpha
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
        pk = signer.F([pk_re, pk_im])
        A = -(pk + 1/pk)
        E_A = EllipticCurve(signer.F, [0,A,0,1,0])
        return E_A

    complete = lambda input: '0'*(len(input) % 2) + input

    def sig_decode_smart(sig):
        zip_chain_bytes = ZIP_CHAIN_LEN*ZIP_CHAIN_HALFBYTES
        siglen = zip_chain_bytes + 2 + 2*TORSION_2POWER_SECPAR_HALFBYTES
        assert(len(sig) == siglen)
        zip_chain = [sig[i:i + ZIP_CHAIN_HALFBYTES] for i in range(0, zip_chain_bytes, ZIP_CHAIN_HALFBYTES)]
        ZIP_CHAIN = []
        for scalar_bytes in zip_chain:
            scalar = bytearray.fromhex(scalar_bytes)[::-1]
            scalar = Integer('0x' + scalar.hex())
            ZIP_CHAIN.append(scalar)
        b_byte = sig[zip_chain_bytes:zip_chain_bytes + 2]
        b = Integer('0x' + b_byte)
        zip = [b, ZIP_CHAIN]
        r_bytes = sig[zip_chain_bytes+2:zip_chain_bytes+2 + TORSION_2POWER_SECPAR_HALFBYTES]
        r = bytearray.fromhex(r_bytes)[::-1]
        r = Integer('0x' + r.hex())
        s_bytes = sig[-TORSION_2POWER_SECPAR_HALFBYTES:]
        s = bytearray.fromhex(s_bytes)[::-1]
        s = Integer('0x' + s.hex())

        return zip, r, s

    def sig_msg_encode_sparallel(sigma, msg):
        swap, zip_chain, r, s, alpha2, alpha_chal = sigma
        sigma_bytes = ''
        for i in range(ZIP_CHAIN_LEN):
            zip_chain_i = f'{int(zip_chain[i]):0{ZIP_CHAIN_HALFBYTES}X}'
            sigma_bytes = zip_chain_i + sigma_bytes
        sigma_bytes = complete(hex(swap)[2:]) + sigma_bytes
        r_bytes = f'{int(r):0{TORSION_2POWER_SECPAR_HALFBYTES}X}'
        s_bytes = f'{int(s):0{TORSION_2POWER_SECPAR_HALFBYTES}X}'
        sigma_bytes = s_bytes + r_bytes + sigma_bytes
        for alpha in [alpha2, alpha_chal]:
            alpha_re = f'{int(alpha[0]):0{FP2_ENCODED_BYTES}X}'
            alpha_im = f'{int(alpha[1]):0{FP2_ENCODED_BYTES}X}'
            sigma_bytes = alpha_im + alpha_re + sigma_bytes
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
            zip, r, s = sig_decode_smart(sig)
            msg = sm[2*SIGNATURE_LEN:]
            mlen = (smlen - SIGNATURE_LEN)
            assert(len(msg) == 2*mlen)

            print(f'count = {count} -- converting {in_variant} to {out_variant} for message: \n{msg}')
            sigma = signer.getSparallel_from_smart(zip, s, msg)

            print('#'*len(topstr) + '\n')
            print(f'Done! signature was {sigma}' + '\n')

            print("Verifying...")
            verified = signer.verify(msg, sigma, signer.pk)
            print('Signature was ' + ('CORRECT' if verified else 'WRONG'))
            print('\n\n')

            pk_bytes = pk_encode(signer.pk)
            sigma_msg_bytes = sig_msg_encode_sparallel(sigma, msg)
            print(f'count = {count}', file=out_file)
            print('seed =', file=out_file)
            print(f'mlen = {mlen}', file=out_file)
            print(f'msg = {msg}', file=out_file)
            print(f'pk = {pk_bytes}', file=out_file)
            print('sk =', file=out_file)
            print(f'smlen = {len(sigma_msg_bytes)//2}', file=out_file)
            print(f'sm = {sigma_msg_bytes}\n', file=out_file)

            count += 1
