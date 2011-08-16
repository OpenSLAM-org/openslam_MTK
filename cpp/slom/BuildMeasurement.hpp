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
 * @file slom/BuildMeasurement.hpp
 * @brief Macros to construct random variables and measurements.
 */
#ifndef SLOM_AUTOCONSTRUCT_H_
#define SLOM_AUTOCONSTRUCT_H_


#include <boost/type_traits.hpp>

#include "../mtk/build_manifold.hpp"


/** define a variable id named @c name_id for the variable @c name */
#define SLOM_RANDOMVAR_ID(name) \
	typedef SLOM_VAR_ID_TYPE(name) name ## _id;


/**
 * Construct a manifold @c name with given entries. 
 * Syntax is like #MTK_BUILD_MANIFOLD. Additionally a variable id named @c name_id is defined.
 */
#define SLOM_BUILD_RANDOMVAR(name, entries) \
	MTK_BUILD_MANIFOLD(name, entries) \
	SLOM_RANDOMVAR_ID(name)

/**
 * Macro to generate Measurements.
 * @param name is the class-name of the Measurement, 
 * @param dim is its dimension.
 * @param variables are the RVs the measurement depends on and 
 * @param data can be some arbitrary data variables.
 * 
 * Both variables and data must be given in a list like this:
 * @code
 * SLOM_BUILD_MEASUREMENT(OdoMeas, 3,
 *    ((Pose, start)) ((Pose, end)),
 *    ((double, linear)) ((double, theta))
 * )
 * @endcode
 * Whitespace is optional, but the double parentheses are necessary.
 * Construction is done entirely in preprocessor.
 * After declaration the function void @c name::eval(MTK::vectview<double, dim> ret) const
 * has to be implemented, which can be done using the macro
 * @c SLOM_IMPLEMENT_MEASUREMENT.
 */
#define SLOM_BUILD_MEASUREMENT(name, dim, variables, data) \
SLOM_BUILD_MEASUREMENT_EXTENDIBLE(name, dim, variables, data) \
{ eval(__ret);} \
void eval(MTK::vectview<double, DIM> __ret) const \
SLOM_CLOSE_MEASUREMENT_EXTENDIBLE(name)

/**
 * 
 */
#define SLOM_IMPLEMENT_MEASUREMENT(name, ret_name) \
	void name::eval(MTK::vectview<double, name::DIM> ret_name) const


/**
 * @warning Experimental feature
 */
#define SLOM_BUILD_MEASUREMENT_JACOBIAN(name, dim, variables, data) \
SLOM_BUILD_MEASUREMENT_EXTENDIBLE(name, dim, variables, data) \
SLOM_CLOSE_MEASUREMENT_EXTENDIBLE(name)

#define SLOM_IMPLEMENT_MEASUREMENT_JACOBIAN(name, ret_name) \
	void name::eval(MTK::vectview<double, name::DIM> ret_name) const


/**
 * SLOM_BUILD_MEASUREMENT_EXTENDIBLE(name, dim, variables, data)
 * generates extendible measurements.
 * Syntax is equal to SLOM_BUILD_MEASUREMENT, but can be followed by further
 * enums, functions and variables (also static) and further constructors
 * It must be finished by SLOM_CLOSE_MEASUREMENT_EXTENDIBLE()
 * @warning Experimental feature
 */
#define SLOM_BUILD_MEASUREMENT_EXTENDIBLE(name, dim, variables, data) \
struct name { \
	enum {DIM = dim}; \
	SLOM_VARREFLIST(variables) \
	MTK_TRANSFORM(SLOM_MSRMNT_GENERATE_DATALIST, data)     \
	name( \
		MTK_TRANSFORM_COMMA(SLOM_MSRMNT_CONSTRUCTOR_ARG, variables) \
		MTK_TRANSFORM(SLOM_MSRMNT_CONSTRUCTOR_ARG_D, data) \
		) : \
		MTK_TRANSFORM_COMMA(SLOM_MSRMNT_GENERATE_CONSTRUCTOR, variables) \
		MTK_TRANSFORM(SLOM_MSRMNT_GENERATE_CONSTRUCTOR_D, data) {}\
	template<class Func> \
	void traverse_variables(Func __function) const { \
		MTK_TRANSFORM(SLOM_MSRMNT_GENERATE_TRAVERSE, variables)} \
	void eval(MTK::vectview<double, DIM> __ret, bool numerical_jacobian) const

#ifndef PARSED_BY_DOXYGEN
//////// internals //////

#define SLOM_VARREFLIST(seq) \
BOOST_PP_FOR_1( \
		( \
				BOOST_PP_SEQ_SIZE(seq), \
				BOOST_PP_SEQ_HEAD(seq), \
				BOOST_PP_SEQ_TAIL(seq) (~), \
				0 ),\
		MTK_ENTRIES_TEST, MTK_ENTRIES_NEXT, SLOM_VARREF_OUTPUT)

#define SLOM_VARREF_OUTPUT(r, state) SLOM_VARREF_OUTPUT_I state
#define SLOM_VARREF_OUTPUT_I(s, head, seq, idx) \
	MTK_APPLY_MACRO_ON_TUPLE(~, \
		BOOST_PP_IF(BOOST_PP_DEC(s), SLOM_OUPUT_VARREF, SLOM_OUPUT_VARREF_AND_ENUM), \
		( BOOST_PP_TUPLE_REM_2 head, idx)) 

#define SLOM_OUPUT_VARREF(type, id, idx) \
	SLOM_VAR_REF_TYPE(type, idx) id; 
#define SLOM_OUPUT_VARREF_AND_ENUM(type, id, idx) \
	SLOM_OUPUT_VARREF(type, id, idx) \
	enum {DEPEND = type::DOF + idx}; \


namespace SLOM {
namespace internal {

template<class T>
struct add_const_ref{
	typedef typename boost::add_reference<typename boost::add_const<T>::type>::type type;
};

}  // namespace internal
}  // namespace SLOM

#define SLOM_CLOSE_MEASUREMENT_EXTENDIBLE(name); \
	typedef SLOM::MeasID<name> id; \
};

#define SLOM_VAR_REF_TYPE(type, idx) SLOM::internal::IMeasurement_Holder::VarRef<type, idx>
#define SLOM_VAR_ID_TYPE(type) SLOM::VarID<type>

#define SLOM_MSRMNT_GENERATE_VARLIST(type, id) SLOM_VAR_REF_TYPE(type, 0) id;
#define SLOM_MSRMNT_GENERATE_DATALIST(type, id) type id;
#define SLOM_MSRMNT_CONSTRUCTOR_ARG(type, id) SLOM_VAR_ID_TYPE(type) id
#define SLOM_MSRMNT_CONSTRUCTOR_ARG_D(type_, id) , SLOM::internal::add_const_ref<type_ >::type id
#define SLOM_MSRMNT_GENERATE_CONSTRUCTOR(type, id) id(id)
#define SLOM_MSRMNT_GENERATE_CONSTRUCTOR_D(type, id) , id(id)
#define SLOM_MSRMNT_GENERATE_DEPEND(type, id) + type::DOF
#define SLOM_MSRMNT_GENERATE_TRAVERSE(type, id) __function(id);

#endif /* not PARSED_BY_DOXYGEN */


#endif /*SLOM_AUTOCONSTRUCT_H_*/
