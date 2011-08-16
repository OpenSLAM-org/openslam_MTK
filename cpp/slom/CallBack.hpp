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
 * @file slom/CallBack.hpp
 * @brief CallBack class for user interaction.
 * Currently this class is only used to get timing informations for subroutines
 * (with controllable verbosity level).
 */

#ifndef CALLBACK_H_
#define CALLBACK_H_


#include <iostream>
#include <cassert>

#include <boost/preprocessor/seq.hpp>
#include <boost/preprocessor/stringize.hpp>

#include "TicToc.hpp"



// These are the possible events:
#define CALLBACK_EVENTS (createSparse) (optimizeStep) (evaluate) (calculateJacobian) (choleskyDecompose) (choleskySolve) (qrSolve) (calculate_delta)



// Helper macros for BOOST_SEQ
#define CALLBACK_COMMA_SEPARATE(unused, data, elem) elem, 
#define CALLBACK_CASE(unused,data,elem) case elem: return BOOST_PP_STRINGIZE(elem);



namespace SLOM {

class Estimator;

/**
 * CallBack class to get runtime information from the algorithm
 * Use Estimator::setCallBack to register a CallBack object.
 */
class CallBack {
	
public:	
	enum Event { 
		BOOST_PP_SEQ_FOR_EACH(CALLBACK_COMMA_SEPARATE,~,CALLBACK_EVENTS)
		numberOfEvents
	};
	struct nullstream : std::ostream {
		struct nullbuf: std::streambuf {
			int overflow(int c) { return traits_type::not_eof(c); }
		} m_sbuf;
		nullstream(): std::ios(&m_sbuf), std::ostream(&m_sbuf) {}
	} nul;
	
protected:
	// calling Estimator object:
	const Estimator &parent;
	
	int verbosity_level;
	
	// output stream:
	std::ostream &out;
	
	
	int depth;
	
	
	inline const char * eventString(Event event){
		switch(event){
		BOOST_PP_SEQ_FOR_EACH(CALLBACK_CASE,~,CALLBACK_EVENTS)
		default: assert(false && "Unknown Event");
		}
		return 0;
	}
	
	// A timer for each event:
	TicToc timer[numberOfEvents];
	
	
public:
	CallBack(const Estimator &parent, int verbosity = 1, std::ostream &out= std::cout) : parent(parent), verbosity_level(verbosity), out(verbosity > 0 ? out : nul), depth(0) {};
	virtual ~CallBack() {};
	
	virtual void start(Event event){
		assert(event < numberOfEvents);
		if(depth<verbosity_level) {
			for(int i=0; i<depth; i++) out << "\t";
			out << "starting " << eventString(event) << std::endl;
			timer[event].tic();
		}
		++depth;
	}
	
	virtual void finish(Event event){
		assert(event < numberOfEvents);
		--depth;
		if(depth < verbosity_level) {
			for(int i=0; i<depth; i++) out << "\t";
			out << "finished " << eventString(event) << " in "  
					<< timer[event].toc() << "s" << std::endl;
		}
	}
	
	std::ostream &output(int verbosity){
		if(verbosity <= verbosity_level)
			return out;
		else
			return nul;
	}
	
	
};


class scoped_callback{
	CallBack *cb;
	CallBack::Event e;
	CallBack::nullstream nul;
	friend class Estimator;
	scoped_callback(CallBack *cb, CallBack::Event e) : cb(cb), e(e) {
		if(cb)
			cb->start(e);
	}
	~scoped_callback() {
		if(cb)
			cb->finish(e);
	}
	std::ostream &output(int verbosity){
		if(cb)
			return cb->output(verbosity);
		else 
			return nul;
	}
	
};

}

// undefine internal macros:
#undef CALLBACK_EVENTS
#undef CALLBACK_TO_STR
#undef CALLBACK_COMMA_SEPARATE 


// user access macros
#define CALLBACK_START(event) scoped_callback cb(callback, CallBack:: event);
#define CALLBACK_FINISH(event) 


#endif /* CALLBACK_H_ */
