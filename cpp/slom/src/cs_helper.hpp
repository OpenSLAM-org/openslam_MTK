/*
 *  Copyright (c) 2008--2011, Universitaet Bremen
 *  All rights reserved.
 *
 *  Author: Christoph Hertzberg <chtz@informatik.uni-bremen.de>
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the Universitaet Bremen nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */
/**
 * @file slom/src/cs_helper.hpp
 * @brief CXSparse extensions
 */

#ifndef CS_HELPER_HH_
#define CS_HELPER_HH_

#include "cs.h"

namespace SLOM {

namespace internal {



/* x(0:j) = x(0:j) + beta * A(0:j,j), where x is a dense vector and A(:,j) is sparse */
CS_INT cs_scatter_upper (const cs *A, CS_INT j, CS_ENTRY beta, CS_INT *w, CS_ENTRY *x, CS_INT mark,
		cs *C, CS_INT nz)
{
	CS_INT i, p, *Ap, *Ai, *Ci ;
	CS_ENTRY *Ax ;
	if (!CS_CSC (A) || !w || !CS_CSC (C)) return (-1) ;     /* check inputs */
	Ap = A->p ; Ai = A->i ; Ax = A->x ; Ci = C->i ;
	for (p = Ap [j] ; p < Ap [j+1] ; p++)
	{
		i = Ai [p] ;                            /* A(i,j) is nonzero */
		if(i>=mark) break; // FIXME: this assumes columns of A are sorted, so we can break here
		if (w [i] < mark)
		{
			w [i] = mark ;                      /* i is new entry in column j */
			Ci [nz++] = i ;                     /* add i to pattern of C(:,j) */
			x [i] = beta * Ax [p] ;             /* x(i) = beta*A(i,j) */
		}
		else
		{
			x [i] += beta * Ax [p] ;    /* i exists in C(:,j) already */
		}
	}
	return (nz) ;
}


/* 
 * calculate upper half of A * B, (assuming columns of A are sorted)
 * first entry in result is the diagonal entry
 */
cs *cs_JtJ (const cs *A, const cs *B)
{
	CS_INT p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values, *Bi ;
	CS_ENTRY *x, *Bx, *Cx ;
	cs *C ;
	if (!CS_CSC (A) || !CS_CSC (B)) return (NULL) ;      /* check inputs */
	if (A->n != B->m) return (NULL) ;
	m = A->m ; anz = A->p [A->n] ;
	n = B->n ; Bp = B->p ; Bi = B->i ; Bx = B->x ; bnz = Bp [n] ;
	w = (CS_INT*)cs_calloc (m, sizeof (CS_INT)) ;                    /* get workspace */
	values = (A->x != NULL) && (Bx != NULL) ;
	x = values ? (CS_ENTRY*) cs_malloc (m, sizeof (CS_ENTRY)) : NULL ; /* get workspace */
	C = cs_spalloc (m, n, anz + bnz, values, 0) ;        /* allocate result */
	if (!C || !w || (values && !x)) return (cs_done (C, w, x, 0)) ;
	Cp = C->p ;
	for (j = 0 ; j < n ; j++)
	{
		if (nz + m > C->nzmax && !cs_sprealloc (C, 2*(C->nzmax)+m))
		{
			return (cs_done (C, w, x, 0)) ;             /* out of memory */
		} 
		Ci = C->i ; Cx = C->x ;         /* C->i and C->x may be reallocated */
		Cp [j] = nz ;                   /* column j of C starts here */
		
		w[j] = j+1;    // mark entry C(j,j)
		Ci [nz++] = j; // add index for entry j
		x[j] = 0;      // reset entry
		for (p = Bp [j] ; p < Bp [j+1] ; p++)
		{
			nz = cs_scatter_upper(A, Bi [p], Bx ? Bx [p] : 1, w, x, j+1, C, nz) ;
		}
		if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
	}
	Cp [n] = nz ;                       /* finalize the last column of C */
//	cs_sprealloc (C, 0) ;               /* remove extra space from C */
	return (cs_done (C, w, x, 1)) ;     /* success; free workspace, return C */
}

}  // namespace internal

}  // namespace SLOM

#endif /*CS_HELPER_HH_*/
