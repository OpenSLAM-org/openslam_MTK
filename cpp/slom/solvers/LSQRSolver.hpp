/**
 * @file LSQRSolverMatlab.hpp
 * @brief Brief description
 * 
 */

#ifndef LSQRSOLVER_HPP_
#define LSQRSOLVER_HPP_


#include "GolubKahanLanczos.hpp"


/*
 *
 */
namespace SLOM {

class LSQRSolver : public GolubKahanLanczos {
	
	bool recycle;
	VectorType w;
	double nextGamma;
public:
	LSQRSolver(int max_it = 1000, bool recycle = false) : GolubKahanLanczos(max_it), recycle(recycle) {}
	
	bool solve(VectorType& x, const double & lambda);
	size_t memory() const {return sizeof(double) * (2*getM() + 3*getN());}
};

}

#endif /* LSQRSOLVER_HPP_ */
