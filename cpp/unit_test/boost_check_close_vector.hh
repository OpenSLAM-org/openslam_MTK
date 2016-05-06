/*
 * boost_check_close_vector.hh
 *
 *  Created on: 11.05.2010
 *      Author: chtz
 */

#ifndef BOOST_CHECK_CLOSE_VECTOR_HH_
#define BOOST_CHECK_CLOSE_VECTOR_HH_


#include <boost/test/floating_point_comparison.hpp>

#include <Eigen/Core>

#define MTK_CHECK_CLOSE(left, right, tol) BOOST_CHECK_CLOSE(left, right, 100*tol)


namespace boost {

namespace test_tools {

using unit_test::readonly_property;


template<class ftype, int dim1, int dim2>
class close_at_tolerance<Eigen::Matrix<ftype, dim1, dim2> >{
public:
	// Public typedefs
	typedef bool result_type;
	typedef ftype FPT;
	
	typedef Eigen::Matrix<ftype, dim1, dim2> comp_type;
	
	// Constructor
	template<typename ToleranceBaseType>
	explicit    close_at_tolerance( percent_tolerance_t<ToleranceBaseType>  tolerance, 
			floating_point_comparison_type          fpc_type = FPC_STRONG ) 
	: p_fraction_tolerance( tt_detail::fpt_abs( static_cast<FPT>(0.01)*tolerance.m_value ) )
	, p_strong_or_weak( fpc_type ==  FPC_STRONG )
	{}
	template<typename ToleranceBaseType>
	explicit    close_at_tolerance( fraction_tolerance_t<ToleranceBaseType> tolerance, 
			floating_point_comparison_type          fpc_type = FPC_STRONG ) 
	: p_fraction_tolerance( tt_detail::fpt_abs( tolerance.m_value ) )
	, p_strong_or_weak( fpc_type ==  FPC_STRONG )
	{}
	
	bool        operator()( comp_type left, comp_type right ) const
	{
		FPT diff = ( left - right ).norm();
		FPT d1   = tt_detail::safe_fpt_division( diff, 1 + right.norm() );  // modification: 1 + ... to check for absolute closeness
		FPT d2   = tt_detail::safe_fpt_division( diff, 1 + left.norm() );
		
		return p_strong_or_weak 
		? (d1 <= p_fraction_tolerance && d2 <= p_fraction_tolerance) 
				: (d1 <= p_fraction_tolerance || d2 <= p_fraction_tolerance);
	}
	
	// Public properties
	readonly_property<FPT>  p_fraction_tolerance;
	readonly_property<bool> p_strong_or_weak;
};

}
}

#endif /* BOOST_CHECK_CLOSE_VECTOR_HH_ */
