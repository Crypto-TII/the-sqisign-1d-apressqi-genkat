
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

	if e == 0:
		return ProjACplus_to_ProjA(A24plus), eval_pts

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
	if e > 0:
		ProjA, eval_pts = four_iso_chain_strategy(P, eval_pts, A24plus, e, strategy)
	else:
		ProjA = ProjACplus_to_ProjA(A24plus)

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
