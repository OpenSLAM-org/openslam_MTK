/**
 * @file LSMRSolver.hpp
 * @brief Brief description
 * 
 */

#ifndef LSMRSOLVER_HPP_
#define LSMRSOLVER_HPP_

#include "IterativeSolver.hpp"

#include "GolubKahanLanczos.hpp"


namespace SLOM {

/*
 *
 */
class LSMRSolver : public GolubKahanLanczos {
	// TODO get more parameters as well
	
public:
	LSMRSolver(int max_it = 1000) : GolubKahanLanczos(max_it) {};
	bool solve(Eigen::VectorXd& x, const double & lambda = 0.0);
	size_t memory() const {return sizeof(double) * (2*getM() + 4*getN());}
};

} /* namespace SLOM */
#endif /* LSMRSOLVER_HPP_ */
