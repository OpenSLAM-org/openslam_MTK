/**
 * @file /mtk-trunk/slom/solvers/GolubKahanLanczos.cpp
 * @brief Brief description
 * 
 */

#include "GolubKahanLanczos.hpp"

namespace SLOM {


void GolubKahanLanczos::init() {
	const SparseType& A = getJ();
	const VectorType& b = getResiduum();
	beta = std::sqrt(getRSS()); // == b.norm();
	u = (1.0/beta) * b;
	v = A.transpose() * u;
	applyPreconditioner(v, true);
	alpha = v.norm();
	z = v = (1.0/alpha) * v;
	applyPreconditioner(z, false);
}

void GolubKahanLanczos::step() {
	const SparseType& A = getJ();
	u = A * z - alpha * u;
	beta = u.norm();
	u = (1.0/beta) * u;
	vt = A.transpose() * u;
	applyPreconditioner(vt, true);
	v = vt - beta * v;
	alpha = v.norm();
	z = v = (1.0/alpha) * v;
	applyPreconditioner(z, false);
}

} /* namespace SLOM */
