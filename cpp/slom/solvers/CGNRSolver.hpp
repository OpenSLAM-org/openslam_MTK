/**
 * @file /mtk-trunk/slom/solvers/CGNRSolver.hpp
 * @brief Brief description
 * 
 */

#ifndef CGNRSOLVER_HPP_
#define CGNRSOLVER_HPP_

#include "IterativeSolver.hpp"


namespace SLOM {

/*
 *
 */
class CGNRSolver : public IterativeSolver {
	
	VectorType p;   // current direction of descent
	VectorType z;   // current preconditioned gradient
	double n2s;
	void init(VectorType &z, const VectorType &delta);
	bool recycle;
public:
	CGNRSolver(int max_it, bool recycle=true) : IterativeSolver(max_it), n2s(-1), recycle(recycle) {}
	
	
	bool solve(VectorType& delta, const double & lambda);
	size_t memory() const {return sizeof(double) * (getM() + 3*getN());}

};

} /* namespace SLOM */
#endif /* CGNRSOLVER_HPP_ */
