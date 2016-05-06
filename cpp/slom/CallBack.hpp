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

#include <stack>

#include "TicToc.hpp"



namespace SLOM {

class Estimator;

/**
 * CallBack class to get runtime information from the algorithm
 * Use Estimator::setCallBack to register a CallBack object.
 */
class CallBack {
	
public:	
	struct nullstream : std::ostream {
		struct nullbuf: std::streambuf {
			int overflow(int c) { return traits_type::not_eof(c); }
		} m_sbuf;
		nullstream(): std::ios(&m_sbuf), std::ostream(&m_sbuf) {}
	} nul;
	
protected:
	std::stack<TicToc> stck;
	
	int verbosity_level;
	
	// output stream:
	std::ostream &out;
	
	
	int depth;
	
	
public:
	__attribute__((deprecated))
	CallBack(const Estimator & /*parent*/, int verbosity = 1, std::ostream &out= std::cout) : verbosity_level(verbosity), out(out), depth(0) {};
	CallBack(int verbosity = 1, std::ostream &out= std::cout) : verbosity_level(verbosity), out(out), depth(0) {};
	virtual ~CallBack() {}
	
	virtual void start(const std::string& event) {
		if(depth < verbosity_level) {
			for(int i=0; i<depth; i++) out << "\t";
			out << "starting " << event << std::endl;
			stck.push(TicToc());
		}
		++depth;
	}
	
	virtual void finish(const std::string& event) {
		--depth;
		if(depth < verbosity_level) {
			for(int i=0; i<depth; i++) out << "\t";
			out << "finished " << event << " in "  
					<< stck.top().toc() << "s" << std::endl;
			stck.pop();
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
	const std::string e;
	CallBack::nullstream nul;
public:
	scoped_callback(CallBack *cb, const std::string& e_) : cb(cb), e(e_) {
		if(cb)
			cb->start(e);
	}
	~scoped_callback() {
		if(cb)
			cb->finish(e);
	}
	std::ostream &output(int verbosity) {
		if(cb)
			return cb->output(verbosity);
		else 
			return nul;
	}
	
	std::ostream & operator()(int verbosity=1) {
		return output(verbosity);
	}
	
};

}

// user access macros
#define SLOM_STRINGIFY(x) SLOM_STRINGIFY_I(x)
#define SLOM_STRINGIFY_I(x) #x

#define CALLBACK_START(event) scoped_callback cb(getCallback(), SLOM_STRINGIFY(event));
#define CALLBACK_START_M() scoped_callback cb(getCallback(), __PRETTY_FUNCTION__);
#define CALLBACK_FINISH(event)


#endif /* CALLBACK_H_ */
