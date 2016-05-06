/**
 * @file /mtk/slom/Optimizer.hpp
 * @brief Brief description
 * 
 */

#ifndef OPTIMIZER_HPP_
#define OPTIMIZER_HPP_

#include <iostream>
#include <boost/program_options.hpp>
#include <boost/scoped_ptr.hpp>


namespace SLOM {





namespace po = boost::program_options;


class OptimizerOptions {
	
	po::options_description options;
	bool initialized;

public:
	char solver;
	char preconditioner;
	int inner_iteration;
	int optimization_steps;
	double lambda;
	int verbosity;
	
	
	
private:
	void setDefaultOptions();
	void make_options();
	
	
public:
	OptimizerOptions() : options("OptimizerOptions"), initialized(false) {
		setDefaultOptions();
	}
	po::options_description & getOptions() {
		if(!initialized) {
			make_options();
			initialized = true;
		}
		return options;
	}
	
	operator po::options_description& () {
		return getOptions();
	}
};


class Estimator;
class CallBack;
class LevenbergMarquardt;
class IterativeSolver;

/*
 *
 */
class Optimizer {
	const OptimizerOptions &opt;
	
	boost::scoped_ptr<Estimator> e;
	
	LevenbergMarquardt *algo;
	IterativeSolver* iterativeSolver;
	
	boost::scoped_ptr<CallBack> cb;

	
public:
	Optimizer(const OptimizerOptions& opt = OptimizerOptions(), std::ostream& out = std::cout);
	
	~Optimizer();
	
	
	
	void optimize();
	
	void adaptInnnerIteration(int iterations);

	double getLambda() const;
	void setLambda(const double &lambda);
	
	Estimator & getEstimator() {return *e;}
};

} /* namespace SLOM */
#endif /* OPTIMIZER_HPP_ */
