from ApresSQI import SQIsign
from xonly import sqrtNonConstTime
import time
import os

from sys import argv

if __name__=="__main__":
    param = 'p248'
    if len(argv) > 1:
        param = argv[1]
    seeded = True
    if len(argv) > 2:
        seeded = not argv[2] == 'False'
    compressed = True
    if len(argv) > 3:
        compressed = not argv[3] == 'False'

    assert param in ['NIST', '7-block', '4-block', 'p248']
    assert  not seeded

    topstr = f'############# Using {param} prime #############'
    print(topstr)

    if param == 'NIST':
        ZIP_CHAIN_LEN = 14
        ID2ISO_COMPRESSED_LONG_TWO_ISOG_ZIP_CHAIN_BYTES = 2*10
        TORSION_2POWER_BYTES = 2*10
        TORSION_3POWER_BYTES = 2*8
        TORSION_23POWER_BYTES = 2*17
        kat_name = '782'
        variant_name = 'p1913'
    elif param == '4-block':
        ZIP_CHAIN_LEN = 4
        ID2ISO_COMPRESSED_LONG_TWO_ISOG_ZIP_CHAIN_BYTES = 2*31
        TORSION_2POWER_BYTES = 2*31
        kat_name = 'p4'
        variant_name = 'p4'
    elif param == '7-block':
        ZIP_CHAIN_LEN = 7
        ID2ISO_COMPRESSED_LONG_TWO_ISOG_ZIP_CHAIN_BYTES = 2*19
        TORSION_2POWER_BYTES = 2*19
        kat_name = 'p7'
        variant_name = 'p7'
    elif param == 'p248':
        ZIP_CHAIN_LEN = 4
        ID2ISO_COMPRESSED_LONG_TWO_ISOG_ZIP_CHAIN_BYTES = 2*31
        TORSION_2POWER_BYTES = 2*31
        kat_name = 'p248'
        variant_name = 'p248'
    FP2_ENCODED_BYTES = 64
    TORSION_2POWER_SECPAR_BYTES = 2*16

    def pk_encode(pk):
        E_A = pk.E_A
        A = E_A.a_invariants()[1]
        if compressed:
            delta = A**2 - 4
            sdelta = sqrtNonConstTime(delta)
            alpha = (-A-sdelta)/2
            pk = alpha
        else:
            pk = A
        pk_bytes = f'{int(pk[0]):0{FP2_ENCODED_BYTES}X}'
        pk_bytes = f'{int(pk[1]):0{FP2_ENCODED_BYTES}X}' + pk_bytes
        pk_bytes = bytearray.fromhex(pk_bytes)[::-1]

        return pk_bytes.hex().upper()

    complete = lambda input: '0'*(len(input) % 2) + input
    def sig_msg_encode(sigma, msg):
        zip, r, s = sigma
        zip_chain = zip[1]
        sigma_bytes = ''
        for i in range(ZIP_CHAIN_LEN):
            zip_chain_i = f'{int(zip_chain[i]):0{ID2ISO_COMPRESSED_LONG_TWO_ISOG_ZIP_CHAIN_BYTES}X}'
            sigma_bytes = zip_chain_i + sigma_bytes
        sigma_bytes = complete(hex(zip[0]).upper()[2:]) + sigma_bytes
        if param == 'NIST':
            r_bytes = f'{int(r):0{TORSION_23POWER_BYTES}X}'
            select23 = hex(2*s[2]+s[0]).upper()[2:]
            s2_bytes = f'{int(s[1]):0{TORSION_2POWER_BYTES}X}'
            s3_bytes = f'{int(s[3]):0{TORSION_3POWER_BYTES}X}'
            sigma_msg_bytes = s3_bytes + s2_bytes + complete(select23) + r_bytes + sigma_bytes
        # elif param == '4-block' or param == '7-block' or param == 'p248':
        else:
            r_bytes = f'{int(r):0{TORSION_2POWER_SECPAR_BYTES}X}'
            s2_bytes = f'{int(s):0{TORSION_2POWER_SECPAR_BYTES}X}'
            sigma_msg_bytes = s2_bytes + r_bytes + sigma_bytes
        sigma_msg_bytes = bytearray.fromhex(sigma_msg_bytes)[::-1] + bytearray.fromhex(msg)

        return sigma_msg_bytes.hex().upper()

    def sig_msg_encode_uncompressed(sigma, msg):
        gens, E1_A, _, _ = sigma
        sigma_bytes = ''
        for i in range(ZIP_CHAIN_LEN):
            gens_i_re = f'{int(gens[i][0]):0{FP2_ENCODED_BYTES}X}'
            gens_i_im = f'{int(gens[i][1]):0{FP2_ENCODED_BYTES}X}'
            sigma_bytes = gens_i_im + gens_i_re + sigma_bytes
        E1_A_bytes_re = f'{int(E1_A[0]):0{FP2_ENCODED_BYTES}X}'
        E1_A_bytes_im = f'{int(E1_A[1]):0{FP2_ENCODED_BYTES}X}'
        sigma_bytes = E1_A_bytes_im + E1_A_bytes_re + sigma_bytes
        sigma_msg_bytes = bytearray.fromhex(sigma_bytes)[::-1] + bytearray.fromhex(msg)

        return sigma_msg_bytes.hex().upper()

    signer = SQIsign(param)
    signer.load_privkey(param=param)

    kat_file = open(f'KATs/PQCsignKAT_{kat_name}_lvl1.rsp', 'w')
    print(f'# lvl1_{variant_name}\n', file=kat_file)

    mlen = 0
    for count in range(100):
        mlen = mlen + 33
        msg_bytes = os.urandom(mlen)
        msg = msg_bytes.hex().upper()
        # msg = "Ska'ru værra med i bakken?"

        print(f'count = {count} -- Signing message: \n{msg}')

        tstart = time.time()

        sigma = signer.Sign(msg, seeded = seeded, compressed = compressed)

        total_time = time.time() - tstart

        print('#'*len(topstr) + '\n')
        print(f'Done! signature was {sigma}' + '\n')

        print("Verifying...")
        if compressed:
            verified = signer.verify(msg, sigma, signer.pk)
        else:
            verified = signer.verify_uncompressed(msg, sigma, signer.pk)
        print('Signature was ' + ('CORRECT' if verified else 'WRONG'))
        print('\n\n')

        pk_bytes = pk_encode(signer.pk)
        if compressed:
            sigma_msg_bytes = sig_msg_encode(sigma, msg)
        else:
            sigma_msg_bytes = sig_msg_encode_uncompressed(sigma, msg)
        print(f'count = {count}', file=kat_file)
        print('seed =', file=kat_file)
        print(f'mlen = {mlen}', file=kat_file)
        print(f'msg = {msg}', file=kat_file)
        print(f'pk = {pk_bytes}', file=kat_file)
        print('sk =', file=kat_file)
        print(f'smlen = {len(sigma_msg_bytes)//2}', file=kat_file)
        print(f'sm = {sigma_msg_bytes}\n', file=kat_file)
