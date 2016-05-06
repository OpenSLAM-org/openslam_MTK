/**
 * @file /mtk/slom/solvers/IterativeSolver.hpp
 * @brief Brief description
 * 
 */

#ifndef ITERATIVESOLVER_HPP_
#define ITERATIVESOLVER_HPP_

#include "../Solver.hpp"

#include <cassert>

namespace SLOM {

/*
 *
 */
class IterativeSolver : public Solver {
protected:
	int max_it;
public:
	IterativeSolver(int max_it = 1000) : max_it(max_it) {}
	
	int get_max_it() const { return max_it;}
	void set_max_it(int maxIt) { 
		assert(maxIt > 0 && "Need at least one iteration!");
		max_it = maxIt;
	}
};

} /* namespace SLOM */
#endif /* ITERATIVESOLVER_HPP_ */
