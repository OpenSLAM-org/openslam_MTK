/**
 * @file /mtk-trunk/slom/src/Optimizer.cpp
 * @brief Brief description
 * 
 */



#include "../Optimizer.hpp"

#include <slom/Estimator.hpp>
#include <slom/CallBack.hpp>

#include <slom/solvers/LSMRSolver.hpp>
#include <slom/solvers/LSQRSolverMatlab.hpp>
#include <slom/solvers/LSQRSolver.hpp>
#include <slom/solvers/CGNRSolver.hpp>

#ifdef SLOM_CHOLMOD_SUPPORT
#include <Eigen/CholmodSupport>
#endif

#include <slom/algorithms/PowellsDogLeg.hpp>



#include <mtk/build_manifold.hpp>



#define SLOM_OPTIONS_ADDOPTION(type, name, abbr, deflt, description) (#name "," #abbr, po::value(&name)->default_value(name), description \
	)
#define SLOM_OPTIONS_CONSTRUCT(type, name, abbr, deflt, description) name=deflt;
	


// syntax:
// SLOM_OPTIMIZER_OPTIONS( ((type, name, abbr, deflt, description)) ((...)) )




#define SLOM_OPTIMIZER_OPTIONS(name, options_list) \
	void OptimizerOptions::make_options() { \
		options.add_options() \
		MTK_TRANSFORM(SLOM_OPTIONS_ADDOPTION, options_list); \
	} \
	void OptimizerOptions::setDefaultOptions() { \
		MTK_TRANSFORM(SLOM_OPTIONS_CONSTRUCT, options_list) \
	}



namespace SLOM {



SLOM_OPTIMIZER_OPTIONS("Optimizer Options",
		((char,solver, S, 'C', "Use 'C'holesky, LS'Q'/'q'R, LS'M'R or C'g/G'LS solver"))
		((char,preconditioner, K, 'J', 
				"Preconditioner for iterative Solvers: n'0'ne, 'j'acobi, Block-'J'acobi, incomplete 'c'holesky1, or incomplete 'C'holesky2"))
		((int, inner_iteration, i, 100, "Number of iterations for LSQR/LSMR, negative value means -I/dim(measurements)"))
		((int, optimization_steps, s, 10, "number of optimization steps"))
		((double, lambda, l, -1, "Starting lambda for Levenberg optimization. Negative values imply Gauss-Newton optimization"))
		((int, verbosity, v, 2, "Verbosity of Estimator"))
)





Optimizer::Optimizer(const OptimizerOptions& opt, std::ostream& out) : opt(opt), e(new Estimator), cb(new CallBack(opt.verbosity, out)) {
	if(opt.lambda > 0) {
		algo = new LevenbergMarquardt(opt.lambda);
		e->changeAlgorithm(algo);
	} else algo = 0;
	
	int inner_its = opt.inner_iteration;
	switch(opt.solver) {
	case 'c': 
#ifdef SLOM_CHOLMOD_SUPPORT
		iterativeSolver = 0;
		e->changeSolver(new SLOM::CholeskySolver<Eigen::CholmodDecomposition<SparseInterface::SparseType, Eigen::Upper> >(Eigen::CholmodSupernodalLLt));
		break;
#else
		// no break
#endif
	default:
		std::cerr << "Unknown solver option '" << opt.solver << "', using Cholesky instead.\n";
		// no break. Fall through to CholeskySolver
	case 'C': 
		iterativeSolver = 0;
		e->changeSolver(new SLOM::CholeskySolver<>); 
		break;
	case 'Q': 
		iterativeSolver = new SLOM::LSQRSolver(inner_its, true); 
		e->changeSolver(iterativeSolver); 
		break;
	case 'q':
		iterativeSolver = new SLOM::LSQRSolver(inner_its, false);
		e->changeSolver(iterativeSolver);
		break;
	case 'M': 
		iterativeSolver = new SLOM::LSMRSolver(inner_its); 
		e->changeSolver(iterativeSolver); 
		break;
	case 'G': // CGLS with recycling
		iterativeSolver = new SLOM::CGNRSolver(inner_its, true);
		e->changeSolver(iterativeSolver);
		break;
	case 'g': // CGLS w/o recycling
		iterativeSolver = new SLOM::CGNRSolver(inner_its, false);
		e->changeSolver(iterativeSolver);
		break;
	}
	switch(opt.preconditioner){
	case 'J':
		e->changePreconditioner(new SLOM::BlockJacobiPreconditioner());
		break;
	case 'j':
		e->changePreconditioner(new SLOM::JacobiPreconditioner());
		break;
	case 'C':
		e->changePreconditioner(new SLOM::CholIncPreconditioner());
		break;
	default:
		e->changePreconditioner(0);
		break;
	}
	if(opt.verbosity > 0) {
		e->setCallBack(*cb);
	}
	
}

Optimizer::~Optimizer() {
	
}

void Optimizer::optimize() {
	for(int i=0; i<opt.optimization_steps; ++i){
		e->optimizeStep();
	}
}

void Optimizer::adaptInnnerIteration(int iterations) { 
	if(!iterativeSolver) return;
	int it = iterations < 0 ? -iterations / (e->getN()+1) + 1 : iterations;
//		std::cerr << iterations << " " << e.getN() << " " << it << "\n";
	iterativeSolver->set_max_it( it);
}

double Optimizer::getLambda() const {
	if(algo) return algo->getLambda();
	return 0.0/0.0;
}
void Optimizer::setLambda(const double &lambda) {
	if(algo) algo->setLambda(lambda);
}



#ifdef SLOM_CHOLMOD_SUPPORT
template<>
size_t CholeskySolver<Eigen::CholmodDecomposition<SparseInterface::SparseType, Eigen::Upper> >::memory() const {
	// FIXME very rough estimate
	return 3*(get_JtJ().nonZeros()) * (sizeof(SparseType::Scalar) + sizeof(SparseType::Index));
}
#endif


}  // namespace SLOM




