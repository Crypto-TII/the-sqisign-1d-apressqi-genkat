
####################################################################################################################
################################## This file contains the functions needed to do  ##################################
##################################  compute the isogenies from the torsion basis  ##################################
####################################################################################################################



import cProfile, pstats
from io import StringIO
from copy import *
from math import floor
from tors_basis import *


profile = False

if profile:
	pr = cProfile.Profile()
	pr.enable()


################################################
#Isogenies and isogeny chains
################################################

def two_iso_curve(P):
	"""
	Function that computes the 2-isogeneous curve (kernel generated by P, where [2]P != (0,0))
	Follows code from SQIsign NIST version 1.0
 
	Input: point P of order 2, P != (0,0) -- kernel generator
	Output: curve parameter APlusProj = [A24plus, C24]
			constants K0, K1 reused in two_iso_eval
	"""

	assert(P[0] != [0, 0])

	A24plus = fp2_sqr(P[0])
	C24 = fp2_sqr(P[1])
	A24plus = fp2_sub(C24, A24plus)
	K0 = fp2_add(P[0], P[1])
	K1 = fp2_sub(P[0], P[1])

	return [A24plus, C24], [K0, K1]


def two_iso_eval(Q, K):
	"""
	Function that evaluates a point Q at a 2-isogeny 
	Follows code from SQIsign NIST version 1.0
 
	Input: point Q to evaluate, constants K0, K1
	Output: point QEval = phi(Q) on the image curve
	"""
	t0 = fp2_add(Q[0], Q[1])
	t1 = fp2_sub(Q[0], Q[1])
	t2 = fp2_mul(K[0], t1)
	t1 = fp2_mul(K[1], t0)
	t0 = fp2_add (t2, t1)
	t1 = fp2_sub(t2, t1)
	QEvalX = fp2_mul(Q[0], t0)
	QEvalZ = fp2_mul(Q[1], t1)

	return [QEvalX, QEvalZ]

def four_iso_curve(P):
	"""
	Function that computes the 4-isogeneous curve (kernel generated by P, where [2]P != (0,0))
	Follows code from SQIsign NIST version 1.0
 
	Input: point P of order 4, where [2]P != (0,0)
	Output: curve param APlusProj = [A24plus, C24]
	 		constants K1, K2, K3, reused in four_iso_eval
	"""

	K0 = fp2_sqr(P[0])
	K1 = fp2_sqr(P[1])
	K2 = fp2_add(K1, K0)
	K3 = fp2_sub(K1, K0)
	A24plus = fp2_mul(K2, K3)
	C24 = fp2_sqr(K1)

	K4 = fp2_add(P[0], P[1])
	K2 = fp2_sub(P[0], P[1])
	K0 = fp2_add(K1, K1)
	K0 = fp2_add(K0, K0)

	return [A24plus, C24], [K0, K1, K2, K3, K4]


def four_iso_eval(Q, K):
	"""
	Function that evaluates a point Q at a 4-isogeny (generated by P with 2[P] != (0,0))
	Follows code from SQIsign NIST version 1.0
 
	Input: point Q to evaluate, constants K=[K0,K1,K2,K3,K4] from four_isog_curve
	Output: point QEval = phi(Q) on the image curve
	"""

	t0 = fp2_add(Q[0], Q[1])
	t1 = fp2_sub(Q[0], Q[1])
	QEvalX = fp2_mul(t0, K[2])
	QEvalZ = fp2_mul(t1, K[4])
	t0 = fp2_mul(t0, t1)
	t0 = fp2_mul(t0, K[0])
	t1 = fp2_add(QEvalX, QEvalZ)
	QEvalZ = fp2_sub(QEvalX, QEvalZ)
	t1 = fp2_sqr(t1)
	QEvalZ = fp2_sqr(QEvalZ)
	QEvalX = fp2_add(t0, t1)
	t0 = fp2_sub(t0, QEvalZ)
	QEvalX = fp2_mul(QEvalX, t1)
	QEvalZ = fp2_mul(QEvalZ, t0)

	return [QEvalX, QEvalZ]		#6M + 2S


def four_iso_curve_singular(P, ProjA24plus):
	"""
	Function that computes the 4-isogeneous curve (kernel generated by P where [2]P = (0,0))
	Follows code from SQIsign NIST version 1.0
 
	Input: point P of order 4, where [2]P = (0,0)
	Output: curve param APlusProj = [A24plus, C24]
	 		constants K1, K2, K3, reused in four_iso_eval
	"""

	K1 = fp2_copy(ProjA24plus[1])

	assert(P[0] == P[1] or P[0] == fp2_neg(P[1]))

	if P[0] == P[1]:
		K0 = fp2_copy(ProjA24plus[0])
		K2 = fp2_sub(ProjA24plus[0], ProjA24plus[1])
		C24 = fp2_neg(K2)

	elif P[0] == fp2_neg(P[1]):
		K2 = fp2_copy(ProjA24plus[0])
		K0 = fp2_sub(ProjA24plus[0], ProjA24plus[1])
		C24 = fp2_neg(K0)
		C24 = fp2_copy(K2)

	A24plus = fp2_copy(K1)

	return [A24plus, C24], [K0, K1, K2]


def four_iso_eval_singular(P, Q, K):
	"""
	Function that evaluates a point Q at a 4-isogeny (generated by point with [2]P = (0,0))
	Follows code from SQIsign NIST version 1.0
 
	Input: point P of order 4, where [2]P = (0,0), point Q to evaluate, constants K=[K0,K1,K2] from four_isog_curve_singular
	Output:  point QEval = phi(Q) on the image curve
	"""

	t0 = fp2_add(Q[0], Q[1])
	t2 = fp2_sub(Q[0], Q[1])
	t0 = fp2_sqr(t0)
	t2 = fp2_sqr(t2)
	QEvalZ = fp2_sub(t0, t2)

	if P[0] == P[1]:
		t1 = fp2_copy(t2)
	else:
		t1 = fp2_copy(t0)
		t0 = fp2_copy(t2)

	QEvalX = fp2_mul(QEvalZ, K[0])
	QEvalZ = fp2_mul(QEvalZ, K[2])
	QEvalZ = fp2_mul(QEvalZ, t1)
	t1 = fp2_mul(t1, K[1])
	QEvalX = fp2_add(QEvalX, t1)
	QEvalX = fp2_mul(QEvalX, t0)

	return [QEvalX, QEvalZ]


def four_iso_chain_naive(P, eval_pts, ProjA, e):
	"""
	Computes a isogeny of degree 2^e using a chain of 4-isogenies 
 
 	Input: P -- kernel point of order 2^e
			eval_pts -- list of points to be pushed through
			ProjA -- curve constant A = ProjA[0]/ProjA[1]
			e --isogeny degree 2^e
	Output: image curve coefficients
	# 		optional: image of points in eval_pts
	"""


	for i in range(e-2, -2, -2):
		T = xDBLe(P, ProjA, i)
		A24plus, K = four_iso_curve(T)

		if i != 0:
			P = four_iso_eval(T, P, K)

		for j in range(len(eval_pts)):
			eval_pts[j] = four_iso_eval(T, eval_pts[j], K)

		ProjA = ProjACplus_to_ProjA(A24plus)

	return ProjA, eval_pts


def four_iso_chain_strategy(P, eval_pts, A24plus, e, strategy):
	"""
	Function that computes chain of 4-isogenies using strategy
	Algorithm follows the latest SIKE C code
	Input: P -- kernel point of order 2^e
			eval_pts -- list of points to be pushed through
			A24 -- curve constant A = (A+2C : 4C)
			e -- isogeny degree 2^e
			strategy
	Output: image curve coefficients
	 		optional: image of points in eval_pts
	"""

	# if e is odd, start with a 2-isogeny
	if e%2 == 1:
		T = xDBLeA24(P, A24plus, e-1)
		A24plus, K = two_iso_curve(T)
		P = two_iso_eval(P, K)
		for j in range(len(eval_pts)):
			eval_pts[j] = two_iso_eval(eval_pts[j], K)
		e -= 1

	MAX = len(strategy) + 1
	pts = [0] * 50
	pts_index = [0] * 50
	npts = 0
	ii = 0

	#go through tree according to strategy
	index = 0
	for row in range(1, MAX):
		while index < (MAX - row):
			pts[npts] = [P[0], P[1]]
			pts_index[npts] = index
			npts += 1
			m = strategy[ii]
			ii += 1
			P = xDBLeA24(P, A24plus, 2*m)
			index += m

		A24plus, K = four_iso_curve(P)

		for i in range(npts):
			pts[i] = four_iso_eval(pts[i], K)
		for i in range(len(eval_pts)):
			eval_pts[i] = four_iso_eval(eval_pts[i], K)

		P = [pts[npts-1][0],  pts[npts-1][1]]
		index = pts_index[npts-1]
		npts -= 1

	A24plus, K = four_iso_curve(P)
	ProjA = ProjACplus_to_ProjA(A24plus)

	for i in range(len(eval_pts)):
		eval_pts[i] = four_iso_eval(eval_pts[i], K)

	return ProjA, eval_pts


def four_iso_chain_opt(P, eval_pts, ProjA, e, strategy):
	"""
	Function used to call correct optimal strategy depending on e being even or odd
	"""
	if e%2 == 0:
		return four_iso_chain_strategy_opt_even(P, eval_pts, ProjA, e, strategy)
	
	else:
		return four_iso_chain_strategy_opt_odd(P, eval_pts, ProjA, e, strategy)

def four_iso_chain_strategy_opt_even(P, eval_pts, ProjA, e, strategy):
	"""
	Function that computes a 2^e degree isogeny using a chain of 4-isogenies with optimal strategy (e even)
	Algorithm follows the latest SIKE C code
	Input: P -- kernel point of order 2^e with e even
			eval_pts -- list of points to be pushed through
			A24 -- curve constant A = (A+2C : 4C)
			e -- isogeny degree 2^e
			strategy -- defines strtegy and isog degree
	Output: image curve coefficients
	 		optional: image of points in eval_pts
	"""

	assert ProjA[1] == [1,0]

	#convert to curve constant A24plus = (A+2C : 4C)
	A24plus = ProjA_to_ProjACplus(ProjA)

	MAX = len(strategy) + 1
	pts = [0] * 50
	pts_index = [0] * 50
	npts = 0
	ii = 0
	isog_count = 0

	#go through tree according to strategy
	index = 0
	for row in range(1, MAX):
		while index < (MAX - row):
			pts[npts] = [P[0], P[1]]
			pts_index[npts] = index
			npts += 1
			m = strategy[ii]
			ii += 1
			#affine mul before first isogeny
			if isog_count == 0:
				P = xDBLe_aff(P, ProjA, 2*m)
			else:
				P = xDBLeA24(P, A24plus, 2*m)
			index += m

		# in the first step the kernel may contain (0,0)
		if isog_count == 0 and xDBLaff(P, ProjA)[0] == [0, 0]:

			#convert to curve constant A24plus = (A+2C : 4C)
			A24plus = ProjA_to_ProjACplus(ProjA)

			A24plus, K = four_iso_curve_singular(P, A24plus)

			for i in range(npts):
				pts[i] = four_iso_eval_singular(P, pts[i], K)
			for i in range(len(eval_pts)):
				eval_pts[i] = four_iso_eval_singular(P, eval_pts[i], K)

		else:
			A24plus, K = four_iso_curve(P)

			for i in range(npts):
				pts[i] = four_iso_eval(pts[i], K)
			for i in range(len(eval_pts)):
				eval_pts[i] = four_iso_eval(eval_pts[i], K)

		isog_count += 1

		P = [pts[npts-1][0],  pts[npts-1][1]]
		index = pts_index[npts-1]
		npts -= 1

	A24plus, K = four_iso_curve(P)
	ProjA = ProjACplus_to_ProjA(A24plus)

	for i in range(len(eval_pts)):
		eval_pts[i] = four_iso_eval(eval_pts[i], K)

	return ProjA, eval_pts


def four_iso_chain_strategy_opt_odd(P, eval_pts, ProjA, e, strategy):
	"""
	Function that computes a 2^e degree isogeny using a chain of 4-isogenies with optimal strategy (e odd)
	Algorithm follows the latest SIKE C code
	Input: P -- kernel point of order 2^e with e even
			eval_pts -- list of points to be pushed through
			A24 -- curve constant A = (A+2C : 4C)
			e -- isogeny degree 2^e
			strategy -- defines strtegy and isog degree
	Output: image curve coefficients
	 		optional: image of points in eval_pts
	"""

	assert ProjA[1] == [1,0]

	#convert to curve constant A24plus = (A+2C : 4C)
	A24plus = ProjA_to_ProjACplus(ProjA)

	MAX = len(strategy) + 1
	pts = [0] * 50
	pts_index = [0] * 50
	npts = 0
	ii = 0
	isog_count = 0

	#Special case for 2-isogeny and 2^3-isogeny. Doesnt work in the singular case
	if MAX == 1:
		assert e == 1 or e == 3
		if e == 3:
			P4 = xDBLaff(P, ProjA)				
			A24plus, K = four_iso_curve(P4)
			#point for computing 2-isogeny
			P = four_iso_eval(P, K)
			for i in range(npts):
				pts[i] = four_iso_eval(pts[i], K)
			for i in range(len(eval_pts)):
				eval_pts[i] = four_iso_eval(eval_pts[i], K)


		A24plus, K = two_iso_curve(P)

		for i in range(npts):
			pts[i]= two_iso_eval(pts[i], K)
		for i in range(len(eval_pts)):
			eval_pts[i] = two_iso_eval(eval_pts[i], K)
		
		ProjA = ProjACplus_to_ProjA(A24plus)
		return ProjA, eval_pts

	#go through tree according to strategy
	index = 0

	for row in range(1, MAX):
		while index < (MAX - row):
			pts[npts] = [P[0], P[1]]
			pts_index[npts] = index
			npts += 1
			m = strategy[ii]
			ii += 1
			#affine mul before first isogeny
			if isog_count == 0:
				P = xDBLe_aff(P, ProjA, 2*m)
			else:
				P = xDBLeA24(P, A24plus, 2*m)
			index += m

		# in the first step the kernel may contain (0,0)
		if isog_count == 0:
			
			#eliminate extra factor 2 for odd e
			P4 = xDBLaff(P, ProjA)

			# in the first step the kernel may contain (0,0)
			if xDBLaff(P4, ProjA)[0] == [0, 0]:

				#convert to curve constant A24plus = (A+2C : 4C)
				A24plus = ProjA_to_ProjACplus(ProjA)

				A24plus, K = four_iso_curve_singular(P4, A24plus)

				#point for computing 2-isogeny
				P = four_iso_eval_singular(P4, P, K)

				for i in range(npts):
					pts[i] = four_iso_eval_singular(P4, pts[i], K)
				for i in range(len(eval_pts)):
					eval_pts[i] = four_iso_eval_singular(P4, eval_pts[i], K)

			#otherwise, compute standard 4-isogeny
			else:
				A24plus, K = four_iso_curve(P4)
				
				#point for computing 2-isogeny
				P = four_iso_eval(P, K)

				for i in range(npts):
					pts[i] = four_iso_eval(pts[i], K)
				for i in range(len(eval_pts)):
					eval_pts[i] = four_iso_eval(eval_pts[i], K)

			#compute 2-isogeny as the second step
			#P is point of order 2

			A24plus, K = two_iso_curve(P)

			for i in range(npts):
				pts[i]= two_iso_eval(pts[i], K)
			for i in range(len(eval_pts)):
				eval_pts[i] = two_iso_eval(eval_pts[i], K)

		else:
			A24plus, K = four_iso_curve(P)

			for i in range(npts):
				pts[i] = four_iso_eval(pts[i], K)
			for i in range(len(eval_pts)):
				eval_pts[i] = four_iso_eval(eval_pts[i], K)

		isog_count += 1

		P = [pts[npts-1][0],  pts[npts-1][1]]
		index = pts_index[npts-1]
		npts -= 1

	A24plus, K = four_iso_curve(P)
	ProjA = ProjACplus_to_ProjA(A24plus)

	for i in range(len(eval_pts)):
		eval_pts[i] = four_iso_eval(eval_pts[i], K)

	return ProjA, eval_pts


def four_iso_chain_nist(P, eval_pts, ProjA, e, strategy):
	"""
	Function that computes a 2^e degree isogeny using a chain of 4-isogenies with optimal strategy (e odd)
	Algorithm follows code from SQIsign NIST version 1.0

	We take the following approach.
		- first simple 4-isogeny, may be singular ([2]P=(0,0) for kernel generator P)
		- then 2-isog if e is odd
		- then 4-isog chain with strategies for the remaining e//2-1 4-isogs
 
	Input: P -- kernel point of order 2^e with e even
			eval_pts -- list of points to be pushed through
			ProjA --  curve constant A = ProjA[0]/ProjA[1]
			e -- isogeny degree 2^e
			strategy -- strategy for 4-iso chain of e//2-1 steps
	Output: image curve coefficients (A:C)
	 		optional: image of points in eval_pts
	"""
 

	#convert to curve constant A24plus = (A+2C : 4C)
	A24plus = ProjA_to_ProjACplus(ProjA)

	#first 4-isog
	T = xDBLeA24(P, A24plus, e-2)	#point of order 4

	#check if [2]T = (0,0)
	if xDBLA24(T, A24plus)[0] != [0, 0]:
		A24plus, K = four_iso_curve(T) 
		P = four_iso_eval(P, K)
		for j in range(len(eval_pts)):
			eval_pts[j] = four_iso_eval(eval_pts[j], K)

	else:
		A24plus, K = four_iso_curve_singular(T, A24plus)
		P = four_iso_eval_singular(T, P, K)
		for j in range(len(eval_pts)):
			eval_pts[j] = four_iso_eval_singular(T, eval_pts[j], K)

	e -= 2

	#call 4-isog chain, includes fist 2-isogeny if e is odd
	ProjA, eval_pts = four_iso_chain_strategy(P, eval_pts, A24plus, e, strategy)

	return ProjA, eval_pts


def four_iso_chain_nist_aff(P, eval_pts, A, e, strategy):
	"""
	Function that computes a 2^e degree isogeny using a chain of 4-isogenies with optimal strategy (e odd)
	Algorithm follows code from SQIsign NIST version 1.0

	We take the following approach.
		- first simple 4-isogeny, may be singular ([2]P=(0,0) for kernel generator P)
		- then 2-isog if e is odd
		- then 4-isog chain with strategies for the remaining e//2-1 4-isogs
 
	Input: P -- kernel point of order 2^e with e even
			eval_pts -- list of points to be pushed through
			ProjA -- affine curve constant A = ProjA[0]/1
			e -- isogeny degree 2^e
			strategy -- strategy for 4-iso chain of e//2-1 steps
	Output: image curve coefficients (A:C)
	 		optional: image of points in eval_pts
	"""

	#convert to curve constant A24plus = (A+2C : 4C)
	A24plus = ProjA_to_ProjACplus(A)

	assert A[1] == [1,0]
	#first 4-isog
	#T = xDBLe(P, A24plus, e-2)	#point of order 4
	T = xDBLe_aff(P, A, e - 2)

	#check if [2]T = (0,0)
	if xDBLaff(T, A)[0] != [0, 0]:
		A24plus, K = four_iso_curve(T)
		P = four_iso_eval(P, K)
		for j in range(len(eval_pts)):
			eval_pts[j] = four_iso_eval(eval_pts[j], K)

	else:
		A24plus, K = four_iso_curve_singular(T, A24plus)
		P = four_iso_eval_singular(T, P, K)
		for j in range(len(eval_pts)):
			eval_pts[j] = four_iso_eval_singular(T, eval_pts[j], K)

	e -= 2

	#call 4-isog chain, includes fist 2-isogeny if e is odd
	ProjA, eval_pts = four_iso_chain_strategy(P, eval_pts, A24plus, e, strategy)

	return ProjA, eval_pts


def two_iso_chain_naive(P, eval_pts, ProjA, e):
	"""
	Naive computation of a chain of two isogenies 

	Input: kernel point P of order 2^e
			eval_pts: list of points to be pushed through
			Aproj: curve constant A = ProjA[0]/ProjA[1]
			e: isogeny degree 2^e
	Output: image curve coefficients
	 		optional: image of points in eval_pts
	"""

	for i in range(e-1, -1, -1):
		T = xDBLe(P, ProjA, i)
		A24plus = two_iso_curve(T)

		if i != 0:
			P = two_iso_eval(T, P)

		for j in range(len(eval_pts)):
			eval_pts[j] = two_iso_eval(T, eval_pts[j])

		ProjA = ProjACplus_to_ProjA(A24plus)

	return ProjA, eval_pts


def three_iso_curve(P):
	"""
	Function that computes the 3-isogenous curve (kernel generated by P, of order 3) 
	Follows SIKE spec

	Input: point P of order 3
	Output: curve param APlusMinusProj = [A24plus, A24minus]
			constants K1, K2 reused in three_iso_eval
	"""

	K1 = fp2_sub(P[0], P[1])
	t0 = fp2_sqr(K1)
	K2 = fp2_add(P[0], P[1])
	t1 = fp2_sqr(K2)
	t2 = fp2_add(t0, t1)
	t3 = fp2_add(K1, K2)
	t3 = fp2_sqr(t3)
	t3 = fp2_sub(t3, t2)
	t2 = fp2_add(t1, t3)
	t3 = fp2_add(t3, t0)
	t4 = fp2_add(t3, t0)
	t4 = fp2_add(t4, t4)
	t4 = fp2_add(t1, t4)
	A24minus = fp2_mul(t2, t4)
	t4 = fp2_add(t1, t2)
	t4 = fp2_add(t4, t4)
	t4 = fp2_add(t0, t4)
	A24plus = fp2_mul(t3, t4)

	return [A24plus, A24minus], [K1, K2]


def three_iso_eval(Q, K):
	"""
	Function that evaluates a point Q at the 3-isogeny with kernel generated by P of order 3
	Follows SIKE spec

	Input: point Q to evaluate, constants K=[K1,K2] from three_isog_curve
	Output: point QEval = phi(Q) on the image curve
	"""

	t0 = fp2_add(Q[0], Q[1])
	t1 = fp2_sub(Q[0], Q[1])
	t0 = fp2_mul(K[0], t0)
	t1 = fp2_mul(K[1], t1)
	t2 = fp2_add(t0, t1)
	t0 = fp2_sub(t1, t0)
	t2 = fp2_sqr(t2)
	t0 = fp2_sqr(t0)
	QEvalX = fp2_mul(Q[0], t2)
	QEvalZ = fp2_mul(Q[1], t0)

	return [QEvalX, QEvalZ]


def three_iso_chain_strategy(P, eval_pts, ProjA, e, strategy):
	"""
	Function that computes a chain of 3-isogenies using a strategy 
	Algorithm follows the latest SIKE C code
 
	Input: P -- kernel point of order 3^e
			eval_pts -- list of points to be pushed through
			ProjA -- curve constant A = ProjA[0]/ProjA[1]
			e -- isogeny degree 3^e
			strategy
	Output: image curve coefficients
			optional: image of points in eval_pts
	"""
	
	MAX = len(strategy) + 1
	pts = [0] * 50
	pts_index = [0] * 50
	npts = 0
	ii = 0

	A24pm = ProjA_to_ProjACplusminus(ProjA)

	#go through tree according to strategy
	index = 0
	for row in range(1, MAX):
		while index < (MAX - row):
			pts[npts] = [P[0], P[1]]
			pts_index[npts] = index
			npts += 1
			m = strategy[ii]
			ii += 1
			P = xTPLeA24pm(P, A24pm, m)
			index += m

		A24pm, K = three_iso_curve(P)

		for i in range(npts):
			pts[i] = three_iso_eval(pts[i], K)
		for i in range(len(eval_pts)):
			eval_pts[i] = three_iso_eval(eval_pts[i], K)

		P = [pts[npts-1][0],  pts[npts-1][1]]
		index = pts_index[npts-1]
		npts -= 1

	A24pm, K = three_iso_curve(P)
	ProjA = ProjACplusminus_to_ProjA(A24pm)

	for i in range(len(eval_pts)):
		eval_pts[i] = three_iso_eval(eval_pts[i], K)

	return ProjA, eval_pts


def three_iso_chain_naive(P, eval_pts, ProjA, e):
	"""
	Function that computes a chain of 3-isogenies without strategy 
	Algorithm follows the latest SIKE C code
 
	Input: P -- kernel point of order 3^e
			eval_pts -- list of points to be pushed through
			ProjA -- curve constant A = ProjA[0]/ProjA[1]
			e -- isogeny degree 3^e
	Output: image curve coefficients
			optional: image of points in eval_pts
	"""

	A24pm = ProjA_to_ProjACplusminus(ProjA)

	for i in range(e-1, -1, -1):
		T = xTPLeA24pm(P, A24pm, i)
		A24pm, K = three_iso_curve(T)

		if i != 0:
			P = three_iso_eval(P, K)

		for j in range(len(eval_pts)):
			eval_pts[j] = three_iso_eval(eval_pts[j], K)

		ProjA = ProjACplusminus_to_ProjA(A24pm)

	return ProjA, eval_pts



################################################
# Misc functions
################################################

def cost(counter, metric):
	"""
	Function that takes in the counters and calculates the cost (according to an input metric)
 
	Input: array of fp_costs in [mul, sqr, add, sub, inv, Legendre, sqrt]
	Output: cost in fp_ops with regards to the metric
	"""
 
	assert(len(counter) == len(metric))
	res = 0

	for i in range(len(counter)):
		res += counter[i]*metric[i]

	return floor(res)


def print_error(message):

	print("ERROR:")
	print(message)
	print('''\
	  +---\
\           _o/
 \     +---/ \
  \           \__ |
   \         \ \ \|
	\         \/  |
	 \         \  |
	  \         \
	   `.

	''')

""" 

The below is for testing -- comment out when done

"""


if __name__ == "__main__":
	################################################
	#Parameters
	################################################
	# p = 2^248 * 5 - 1
	p = 2261564242916331941866620800950935700259179388000792266395655937654553313279
	update_p(p)
	ProjA = [[6,0],[1,0]]
	#point of full order 2^248 * 5
	P = [[1773546402379325112226663529801404653321449577126288419713404067877613711532, 1433550135797697868049241062915569883468141391825533592943161127007286690774],[1,0]]
	Q=[[1773546402379325112226663529801404653321449577126288419713404067877613711532, 1433550135797697868049241062915569883468141391825533592943161127007286690774],[1,0]]
	P = xMUL(P, 5, ProjA)

	test = xMUL(P, 2**6*3*67*503, ProjA)
	real = xMUL_noDACs(P, 2**6*3*67*503, ProjA)

	if proj_point_normalize(test) == proj_point_normalize(real):
		print("New mul passed")
	else:
		print("New mul DID NOT passed")

 
	strategy = [59, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 27, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1]
	strategy2 = [48, 28, 19, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 7, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 20, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1]
	strategy3 = [109, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 96, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 84, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 73, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 63, 9, 8, 7, 6, 5, 4, 3, 2, 1, 54, 8, 7, 6, 5, 4, 3, 2, 1, 46, 7, 6, 5, 4, 3, 2, 1, 39, 6, 5, 4, 3, 2, 1, 33, 5, 4, 3, 2, 1, 28, 4, 3, 2, 1, 24, 3, 2, 1, 22, 1, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
	strategy248 = [60, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 28, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 12, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1] 
	A, eval_pts = four_iso_chain_strategy_opt_even(P, [Q], ProjA, 248, strategy248)
	#A2, eval_pts2 = four_iso_chain_nist(P, [Q], ProjA, 248, strategy2)
	#A3, eval_pts2 = four_iso_chain_nist(P, [Q], ProjA, 248, strategy3)
	#A, eval_pts = four_iso_chain_naive(P, [], ProjA, 248)
	#print(proj_point_normalize(A) == proj_point_normalize(A2))
	#print(j_invariant(A) == j_invariant(A2))
	#print(proj_point_normalize(A) == proj_point_normalize(A3))
	#print(j_invariant(A) == j_invariant(A3))
	#print(proj_point_normalize(A2) == proj_point_normalize(A3))
	#print(j_invariant(A2) == j_invariant(A3))


	if j_invariant(A) != [1005398237125227574239782691732783210009464266310061134537048278699034927901, 1370507504720702308022614524592851585299189237225148676571172873730749479201]:
		print_error("wrong output curve")
	else:
		print("even 4-isogeny chain test passed")

	#print(proj_point_normalize(A))
	#print(j_invariant(A))
	#Anew, isom = normalize_curve(A)
	#Q = ec_isomorphism_eval(eval_pts[0], isom)
	#print(Anew)
	#print(j_invariant(Anew))
	#print(Q)
	#print(xMUL(Q, 5, Anew))

	#P, Q, PmQ = basis_two_torsion(ProjA, 248)


	P=[[1141137256681557842965575901453345171670607752348919894450361273321971809511,595242365908776507002867825567342700584228436229705025443455832192758761782], [1,0]]
	Q=[[1490097333908319625384015451201834628520151192627517060347973771868071750688,64046556013838021704173044896656232854927427822603392072816728017896082922], [1,0]]
	PmQ=[[1665114252890774976774923640609911541753831812199311166262411535823770966582,1401161614135734264496166816599744494792343057106115563420914145388297201420], [1,0]]

	R = ladder_3pt(P,Q,PmQ,1234567,ProjA)
	R = proj_point_normalize(R)
 
	if R != [[2224230809325394377324749729923901119842877587012945901825662125704538991483, 1284913618136205709971398293573979761263299741630739978464054009149579865088], [1,0]]:
		print_error("wrong output point")
	else:
		print("3pt ladder test passed")

	#print(counter_fp)
	#print(counter_fp2)

	#3-isogeny chain test
	# p = 2^3 * 3^152 * 5 * 7 * 13 - 1
	p = 12120822769750759633013197633117582063548858768767324622185441117685272637239
	update_p(p)
	ProjA = [[6,0],[1,0]]

	#points of order 3^152
	P = [[6339434100391704367219122257368377957419474873595142657584846995035893263157,2464918024073825684093989098490259707650641022375467513054280103206693467727],[1,0]]
	Q = [[4059275945121653528930292969970815053938232637863183768162618361191697955968,3639743966419002013258087048543560288009325574740710081807281905261048086748],[1,0]]

	strategy = [66, 38, 21, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 17, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32, 16, 8, 4, 2, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1]
	ProjA, eval_pts = three_iso_chain_strategy(P, [Q], ProjA, 152, strategy)
	# ProjA, eval_pts = three_iso_chain_naive(P, [Q], ProjA, 152)

	if j_invariant(ProjA) != [1792045146026259896633838605571577017029052098617428379722180755779948922631, 7969154354288370266682464356092764029397112129948386704883684222958524751135]:
		print_error("wrong output curve")
	else:
		print("3-isogeny chain test passed")

	assert(xTPLe(eval_pts[0],ProjA,151)[1] != [0,0])
	assert(xTPLe(eval_pts[0],ProjA,152)[1] == [0,0])

	#midim mul test
	ProjA = [[6,0],[1,0]]
	P = [[526272772743314460161259066546239743803819788904242655092284234378541516556, 2483025186817906569986084049509515687629377703967006388428599553213354504238], [1,0]]
	Q = [[12077669291719591556821436394299819867826444136316972427738258953806746517703, 6355363792140223122313185126760147755617510995955973646150609649439533819536], [1,0]]
	PmQ = [[4579411811330125902492361540009507071606785799463667378208339286938891481412, 5992181484094561776965758706216585395432963484050629171702043490800095530028], [1,0]]
	R = [[1527890450134353601163154659524161610397744488783156200725572479848804949796, 8571405644709434338526515621250308010428647319966025686545586633498058337506], [1,0]]

	S = xMULbidim(P, Q, PmQ, 123456789012345678901234567890, 987654321098765432109876543210, ProjA)

	
	if proj_point_normalize(S) != R:
		print_error("wrong output in bidim mul")
	else:
		print("bidim multiplication test passed")

	p = 158309497004143235930663456066565499018142557160055458647695915635818731929599
	update_p(p)
	ProjA = [[6,0],[1,0]]
	strategy246 = [59, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 27, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1]
	strategy250 = [61, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 29, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 13, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 5, 4, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1]
	P = [[149821601174543526846007568325654721171402155692424641875815905523459940724605, 102137030952563629965960638781577900045082929991982384843452457269573691092224], [1,0]]
	A, eval_pts = four_iso_chain_strategy_opt_odd(P, [], ProjA, 249, strategy248)
	#A, eval_pts = four_iso_chain_nist(P, [], ProjA, 249, strategy246)
	if j_invariant(A) != [117709873713877940346095209274172838550315445518050768191102653827479627824039, 84200665789192361271695362957459910980696937866096852495318131758557958176303]:
		print_error("wrong output curve")
	else:
		print("odd 4-isogeny chain test passed")

	if profile:
		pr.disable()
		s = StringIO()
		sortby = 'cumulative'
		ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
		ps.print_stats()
		print(s.getvalue())
