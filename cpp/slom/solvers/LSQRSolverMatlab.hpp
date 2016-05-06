/**
 * @file LSQRSolverMatlab.hpp
 * @brief Brief description
 * 
 */

#ifndef LSQRSOLVERMATLAB_HPP_
#define LSQRSOLVERMATLAB_HPP_


#include "GolubKahanLanczos.hpp"


/*
 *
 */
namespace SLOM {

class LSQRSolverMatlab : public GolubKahanLanczos {
	
public:
	LSQRSolverMatlab(int max_it = 1000) : GolubKahanLanczos(max_it) {}
	
	bool solve(VectorType& x, const double & lambda);
	size_t memory() const {return sizeof(double) * (2*getM() + 3*getN());}
};

}

#endif /* LSQRSOLVERMATLAB_HPP_ */
