/**
 * @file /mtk-trunk/slom/solvers/GolubKahanLanczos.hpp
 * @brief Brief description
 * 
 */

#ifndef GOLUBKAHANLANCZOS_HPP_
#define GOLUBKAHANLANCZOS_HPP_


#include "IterativeSolver.hpp"


namespace SLOM {

/*
 *
 */
class GolubKahanLanczos : public IterativeSolver {

	VectorType u, v, vt, z;
	double alpha, beta;
	
	
public:
	GolubKahanLanczos(int max_it) : IterativeSolver(max_it) { };
	
	
	/**
	 * Initializes the GKL process
	 */
	void init();
	
	
	/**
	 * Runs the next bidiagonalization step, updates z, alpha and beta
	 */
	void step();
	
	const VectorType& getZ() const { return z;}
	const double& getAlpha() const {return alpha;}
	const double& getBeta()  const {return beta; }
	
	/**
	 * Finds a rotation [c s] which rotates [a b] to [r 1], with r=norm([a b]).
	 * @param c
	 * @param s
	 * @param a
	 * @param b
	 * @return radius r
	 */
	static inline double rotation(double &c, double &s, const double &a, const double &b) {
		double r = std::sqrt(a*a + b*b);
		c = a/r;
		s = b/r;
		return r;
	}

	
};

} /* namespace SLOM */
#endif /* GOLUBKAHANLANCZOS_HPP_ */
