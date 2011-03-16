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
 * @file slom/src/indexed_list.hpp
 * @brief Intrusive double linked list.
 * An intrusive double linked list, (based on boost::intrusive::list) which destructs its members after erasing
 * and equips its members with an index
 */

#ifndef INDEXED_LIST_HH_
#define INDEXED_LIST_HH_

#include <boost/intrusive/list.hpp>

namespace SLOM {

namespace internal {

/**
 * \c indexed_list_hook implements an intrusive list hook which also provides an index.
 */
struct indexed_list_hook : public boost::intrusive::list_base_hook<> {
	int idx;
	indexed_list_hook() : idx(-1) {}
};

/**
 * indexed_list is an intrusive double link list, which also provides an index, 
 * and takes care of deallocation.
 */
template<class type>
struct indexed_list : boost::intrusive::list<type, boost::intrusive::constant_time_size<false> > {
	
	typedef boost::intrusive::list<type, boost::intrusive::constant_time_size<false> > Container;
	
	
	static void delete_A(type* a) { delete a;}
	
	//Container list; // TODO inheritance or aggregation?
	int last_index;
	
	
	indexed_list() : last_index(0) {}
	
	~indexed_list() {
		Container::clear_and_dispose(delete_A);
	}
	
	
	
	template<class X>
	struct id {
		friend class indexed_list;
		typedef typename type::template holder<X>* ptr_type;
		id(ptr_type ptr) : ptr(ptr) {}
		id() : ptr(0) {}
	protected:
		ptr_type ptr;
	};
	
	template<class X>
	id<X> insert(typename type::template holder<X> *x){
		Container::push_back(*x);
		assert(x->idx == -1);
		x->idx = last_index;
		last_index += x->DIM;
		return id<X>(x);
	}
	
	template<class X>
	typename id<X>::ptr_type unhook(id<X> &x) {
		typename Container::iterator it = Container::s_iterator_to(*x.ptr);
		int old_idx = it->idx;
		it = Container::erase(it);
		x.ptr->idx = -1;
		//x.ptr = 0;
		int idx_diff = it->idx - old_idx;
		for( ; it != Container::end(); ++it){
			it->idx -= idx_diff;
		}
		last_index -= idx_diff;
		return x.ptr;
	}
	
	template<class X> 
	void remove(id<X> &x){
		typename Container::iterator it = Container::s_iterator_to(*x.ptr);
		int old_idx = it->idx;
		it = Container::erase_and_dispose(it, delete_A);
		x.ptr = 0;
		int idx_diff = it->idx - old_idx;
		for( ; it != Container::end(); ++it){
			it->idx -= idx_diff;
		}
		last_index -= idx_diff;
	}
	
	int getDim() const {
		return last_index;
	}
	
};

}  // namespace internal

} // namespace SLOM

#endif /* INDEXED_LIST_HH_ */
