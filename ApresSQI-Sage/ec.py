
####################################################################################################################
##################################  This file contains elliptic curve functions   ##################################
####################################################################################################################

from sage.all import *
from xonly import isMontgomery, xPoint, sqrtDeterministic, sqrtNonConstTime
import hashlib

#############################
#                           #
#     Basis generation      #
#                           #
#############################
    
def CompleteBasis(R, D, x = 1):
    E = R.curve()
    i = E.base_field().gens()[0]
    cof = sqrt(E.order())/D 
    facD = factor(D)
    Drad = radical(D)
    Rsmalls = []
    Rsmall = R*(D/Drad)
    for (l, e) in facD:
        Rsmalli = Rsmall*(Drad/l)
        assert Rsmalli
        Rsmalls.append(Rsmalli)

    for _ in range(1000):
        x += i
        try:
            S = E.lift_x(x)*cof
        except:
            continue
        Ssmall = S*(D/Drad)
        basis = True
        for index, (l, e) in enumerate(facD):
            Ssmalli = Ssmall*(Drad/l)
            if (not Ssmalli) or (Rsmalls[index].weil_pairing(Ssmalli, l) == 1):
                basis = False
                break
        if basis:
            RmS = point_difference(R, S)
            if RmS != (R - S).xy()[0]:
                S = -S
                assert RmS == (R - S).xy()[0]
            S.set_order(D)
            return S
    assert False, "Something went wrong in Complete Basis..."

def TorsionBasis(E, D, xOnly = 0, seeds = None, small_ns = None, small_s = None):
    i = E.base_field().gens()[0]
    x = Integer(1)
    cof = sqrt(E.order())/D
    facD = factor(D)
    Drad = radical(D)    
    ## case 1: With seeds ##
    if seeds:
        p, q = seeds
        # Multiply by cofactor after, because of verification shenanigans
        P = E.lift_x(p)
        Q = E.lift_x(q)
        PmQ = point_difference_new(P, Q)
        if PmQ != (P - Q).xy()[0]:
            Q = -Q
            assert PmQ == (P - Q).xy()[0]
        P, Q = P*cof, Q*cof
        return P, Q
    ## case 2: Generate basis + seeds ##
    if small_ns and small_s:
        #P point
        for n, x in enumerate(small_ns):
            try:
                P = E.lift_x(x)
            except:
                continue
            Pexp = P*cof
            Psmall = Pexp*(D/Drad)
            fullorder = True
            for (l, e) in facD:
                if not Psmall*(Drad/l):
                    fullorder = False
                    break
            if fullorder:
                Pexp.set_order(D)
                break
        Psmalls = []
        for (l, e) in facD:
            Psmalli = Psmall*(Drad/l)
            assert Psmalli
            Psmalls.append(Psmalli)
        # Q point
        for m, x in enumerate(small_s):
            try:
                Q = E.lift_x(x)
            except:
                continue
            Qexp = Q*cof
            Qsmall = Qexp*(D/Drad)
            basis = True
            for index, (l, e) in enumerate(facD):
                Qsmalli = Qsmall*(Drad/l)
                if (not Qsmalli) or (Psmalls[index].weil_pairing(Qsmalli, l) == 1):
                    basis = False
                    break
            if basis:
                Qexp.set_order(D)
                break
        PmQ = point_difference_new(P, Q)
        if PmQ != (P - Q).xy()[0]:
            Qexp = -Qexp
        return Pexp, Qexp, (n, m)
    
    ## Case 3: No seeds involved
    while True:
        x += i
        try:
            P = E.lift_x(x)*cof
        except:
            continue
        Psmall = P*(D/Drad)
        fullorder = True
        for (l, e) in facD:
            if not Psmall*(Drad/l):
                fullorder = False
                break
        if fullorder:
            P.set_order(D)
            break
    Q = CompleteBasis(P, D)
    if xOnly == 1:
        return P, Q
    PmQ = P-Q
    xP, xQ, xPmQ = xPoint(P.xy()[0], E), xPoint(Q.xy()[0], E), xPoint(PmQ.xy()[0], E)
    if xOnly == 2:
        return P, Q, xP, xQ, xPmQ
    return xP, xQ, xPmQ

def TorsionBasis_2tof(E, alpha, f, xOnly = 0):
    F = alpha.parent()
    i = F.gens()[0]
    p = F.characteristic()
    cof = (p+1)/2**f

    x = F(1)
    while True:
        x += i
        if x.is_square():
            continue
        try:
            P = E.lift_x(x)
            break
        except:
            continue
    P = cof * P

    z = F(2)
    while True:
        z += i
        # print(f'z = {z} -- is not square z: {not z.is_square()} -- is square z-1: {(z-1).is_square()}')
        if not z.is_square() or (z-1).is_square():
            continue
        try:
            Q = E.lift_x(z*alpha)
            break
        except:
            continue
    Q = cof * Q

    PmQ = point_difference_new(P, Q)
    if PmQ != (P - Q).xy()[0]:
        Q = -Q
        assert PmQ == (P - Q).xy()[0]
    P.set_order(2**f)
    Q.set_order(2**f)

    if xOnly == 1:
        return P, Q
    PmQ = P-Q
    xP, xQ, xPmQ = xPoint(P.xy()[0], E), xPoint(Q.xy()[0], E), xPoint(PmQ.xy()[0], E)
    if xOnly == 2:
        return P, Q, xP, xQ, xPmQ
    return xP, xQ, xPmQ

def point_difference(P, Q, x_only = False):
	#follows code from SQIsign NIST version 1.0
	#input: affine points P,Q
	#		affine curve constant ProjA = [A, [1,0]]
	#output: x-coordinate of P-Q = PmQ

	#check if all inputs are affine
    if x_only:
        E = P.curve
        assert isMontgomery(E)
        A = E.a2()
        Px, Qx = P.X, Q.X
    else:
        E = P.curve()
        assert isMontgomery(E)
        A = E.a2()
        Px, Qx = P.xy()[0], Q.xy()[0]
    PmQZ = Px - Qx
    t2 = Px*Qx
    t3 = t2 - 1
    t0 = PmQZ*t3
    PmQZ = PmQZ**2
    t0 = t0**2
    t1 = t2 + 1
    t3 = Px + Qx
    t1 = t1*t3
    t2 = t2*A
    t2 = 2*t2
    t1 = t1 + t2
    t2 = t1**2
    t0 = t2-t0
    assert t0.is_square()
    t0 = sqrtDeterministic(t0)
    PmQX = t0 + t1
    return PmQX/PmQZ

def point_difference_new(P, Q, x_only=False):
    #check if all inputs are affine
    if x_only:
        E = P.curve
        assert isMontgomery(E)
        A = E.a2()
        Px, Qx = P.X, Q.X
    else:
        E = P.curve()
        assert isMontgomery(E)
        A = E.a2()
        Px, Qx = P.xy()[0], Q.xy()[0]

    t0 = Px * Qx
    t1 = 1 * 1
    Bxx = t0 - t1
    Bxx = Bxx**2
    Bxz = t0 + t1
    t0 = Px * 1
    t1 = 1 * Qx
    Bzz = t0 + t1
    Bxz = Bxz * Bzz
    Bzz = t0 - t1
    Bzz = Bzz**2
    t0 = t0 * t1
    t0 = t0 * A
    t0 = t0 + t0
    Bxz = Bxz + t0

    t0 = Bxz**2
    t1 = Bxx * Bzz
    t0 = t0 - t1
    t0 = sqrtNonConstTime(t0)
    PmQx = Bxz + t0
    PmQz = Bzz
    return PmQx/PmQz

#############################
#                           #
#     Hashing to chall      #
#                           #
#############################

FP2_ENCODED_BYTES = 64
def hashToPoint(D, msg, E, seeds=None, small_ns=None, small_s=None):
    J = E.j_invariant()
    E_bytes = f'{int(J[0]):0{FP2_ENCODED_BYTES}x}'
    E_bytes = f'{int(J[1]):0{FP2_ENCODED_BYTES}x}' + E_bytes
    E_bytes = bytearray.fromhex(E_bytes)[::-1]

    E_msg_bytes = E_bytes + bytearray.fromhex(msg)

    H = hashlib.shake_256()
    H.update(E_msg_bytes)
    s = int.from_bytes(H.digest(32), 'little')

    if small_ns and small_s:
        P, Q, seeds = TorsionBasis(E, D, small_ns=small_ns, small_s=small_s)
        return P + s*Q, seeds
    else:
        if seeds:
            P, Q = TorsionBasis(E, D, seeds=seeds, small_ns=small_ns, small_s=small_s)
        else:
            P, Q = TorsionBasis(E, D, xOnly = 1)
        return P + s*Q

RADIX = 64
NWORDS_FIELD = 4
def hashToPointNewFast(E, msg):
    i = E.base_field().gens()[0]
    p = E.base_field().characteristic()
    cof = (p+1)//(2**128)

    J = E.j_invariant()
    E_bytes = f'{int(J[0]):0{FP2_ENCODED_BYTES}x}'
    E_bytes = f'{int(J[1]):0{FP2_ENCODED_BYTES}x}' + E_bytes
    # print('E_bytes =', E_bytes)

    E_bytes = bytearray.fromhex(E_bytes)[::-1]

    E_msg_bytes = E_bytes + bytearray.fromhex(msg)

    z = 1
    while True:
        z += i
        if not z.is_square():
            break
    for inc in range(256):
        inc_bytes = f'{inc:0{2}x}'
        input_bytes = E_msg_bytes + bytearray.fromhex(inc_bytes)[::-1]

        H = hashlib.shake_256()
        H.update(input_bytes)
        s = int.from_bytes(H.digest(NWORDS_FIELD*RADIX//8), 'little')

        x = z*s #This step fails when not using p4

        try:
            K = E.lift_x(x)
        except Exception as e:
            continue
        K = cof * K
        return K
    # Failed to hash after 2**8 attempts
    assert False
