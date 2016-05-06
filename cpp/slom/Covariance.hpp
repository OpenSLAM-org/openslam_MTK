/**
 * @file /mtk-trunk/slom/Covariance.hpp
 * @brief Brief description
 * 
 */


#ifndef COVARIANCE_HPP_
#define COVARIANCE_HPP_

#include "CholeskyCovariance.hpp"

#include <boost/type_traits/is_arithmetic.hpp>

namespace SLOM {

namespace internal {

template<class scalar>
struct InvDeviation {
	scalar scale;
	InvDeviation(const scalar& scale) : scale(scale) {}
	template<int dim>
	void apply(MTK::vectview<scalar, dim> arr) const {
		arr *= scale;
	}
};

template<>
struct InvDeviation<void>{
	template<int dim, class scalar>
	void apply(MTK::vectview<scalar, dim> /* arr */) const {
		// do nothing
	}
};

template<int dim, class scalar>
struct InvDeviation<Eigen::Matrix<scalar, dim, 1> >{
	typedef Eigen::Matrix<scalar, dim, 1> vect;
	vect dev;
	InvDeviation(const vect& dev) : dev(dev) {}
	void apply(MTK::vectview<scalar, dim> arr) const {
		arr.array() *= dev.array();
	}
};

template<int dim, class scalar>
struct InvDeviation<SLOM::CholeskyCovariance<dim, scalar> > {
	typedef SLOM::CholeskyCovariance<dim, scalar> DevType;
	
	DevType dev;
	
	InvDeviation(const DevType& dev) : dev(dev) {}
	
	void apply(MTK::vectview<scalar, dim> arr) const {
		dev.apply(arr);
	}
};

}  // namespace internal


#if 0
template<class scalar>
static internal::InvDeviation<scalar> Covariance(const scalar& scale){
	return internal::InvDeviation<scalar>(1.0/std::sqrt(scale));
}
template<class scalar>
static internal::InvDeviation<scalar> InvCovariance(const scalar& scale){
	return internal::InvDeviation<scalar>(std::sqrt(scale));
}
template<class scalar>
static internal::InvDeviation<scalar> Deviation(const scalar& scale){
	return internal::InvDeviation<scalar>(1.0/scale);
}
template<class scalar>
static internal::InvDeviation<scalar> InvDeviation(const scalar& scale){
	return internal::InvDeviation<scalar>(scale);
}

#endif

namespace internal {

template<class Type, bool isScalar = boost::is_arithmetic<Type>::value>
struct DeviationChooser;


// By default it is assumed that Type is a scalar
template<class Type>
struct DeviationChooser<Type, true>
{
	typedef InvDeviation<Type> type;
	
	template<bool sqrt_, bool inv_>
	static type compute(const Type& cov) {
		// default: assume Type is a scalar
		using std::sqrt;
		Type val = cov;
		if(sqrt_) val = sqrt(val);
		if(inv_)  val = 1.0/val;
		return type(val);
	}
};


/**
 * Default MatrixVectorChooser is supposed to fail to compile
 */
template<class scalar, int rows, int cols>
struct MatrixVectorChooser{
	typedef typename scalar::ONLY_VECTORS_OR_SQUARE_MATRICES_ARE_ALLOWED type;
	template<bool sqrt_, bool inv_, class storage>
	static type compute(const storage& ){
		return 0;
	}
	
};
//FIXME: static assertion ONLY_VECTORS_OR_SQUARE_MATRICES_ARE_ALLOWED;

template<class scalar, int rows>
struct MatrixVectorChooser<scalar, rows, 1>{
	typedef Eigen::Matrix<scalar, rows, 1> storage;
	typedef InvDeviation<storage> type;
	
	template<bool sqrt_, bool inv_>
	static type compute(const storage& cov){
		type ret(cov);
		if(sqrt_) ret.dev = ret.dev.cwiseSqrt();
		if(inv_)  ret.dev = ret.dev.cwiseInverse();
		return ret;
	}
};

template<class scalar, int rows>
struct MatrixVectorChooser<scalar, rows, rows>{
	typedef SLOM::CholeskyCovariance<rows, scalar> storage;
	typedef InvDeviation<storage> type;
	
	template<bool sqrt_, bool invert_>
	static type compute(const Eigen::Matrix<scalar, rows, rows>& cov){
		typedef Eigen::Matrix<scalar, rows, rows> Matrix;
		Matrix temp;
		const scalar *data = cov.data();
		if(!invert_ && sqrt_) {
			// HACK which is necessary, because the inverse of a Cholesky factor is not the Cholesky factor of the inverse.
			temp = cov.llt().solve(Matrix::Identity());
			data = temp.data();
		}
		storage ret(data, sqrt_ ? CholeskyMode::CHOLESKY_FULL : CholeskyMode::COPY_UPPER_FULL);
		if(invert_ || sqrt_)  ret.invert();
		return ret;
	}
};




template<class Derived>
struct DeviationChooser<Eigen::MatrixBase<Derived>, false>{
	typedef Eigen::MatrixBase<Derived> Matrix;
	typedef MatrixVectorChooser<typename Matrix::Scalar, Matrix::RowsAtCompileTime, Matrix::ColsAtCompileTime> Chooser;
	typedef typename Chooser::type type;
	
	template<bool sqrt_, bool inv_>
	static type compute(const Matrix& cov){
		return Chooser::template compute<sqrt_, inv_>(cov);
	}
};


}  // namespace internal


#if 0
template<int dim, class scalar>
static internal::InvDeviation<scalar> Covariance(const scalar& scale){
	return internal::InvDeviation<scalar>(1.0/std::sqrt(scale));
}
template<int dim, class scalar>
static internal::InvDeviation<scalar> InvCovariance(const scalar& scale){
	return internal::InvDeviation<scalar>(std::sqrt(scale));
}
template<int dim, class scalar>
static internal::InvDeviation<scalar> Deviation(const scalar& scale){
	return internal::InvDeviation<scalar>(1.0/scale);
}
template<int dim, class scalar>
static internal::InvDeviation<scalar> InvDeviation(const scalar& scale){
	return internal::InvDeviation<scalar>(scale);
}
#endif


#define SLOM_MAKE_DEVIATIONCHOOSER(templ, arg_type, what, sqrt, inv) \
template<templ> \
static typename internal::DeviationChooser<arg_type >::type what(const arg_type& cov) { \
	return internal::DeviationChooser<arg_type >::template compute<sqrt, inv>(cov); \
}

#define SLOM_MAKE_ALL_DEVIATIONCHOOSER(templ, arg_type) \
		SLOM_MAKE_DEVIATIONCHOOSER(templ, arg_type, Covariance, true, true) \
		SLOM_MAKE_DEVIATIONCHOOSER(templ, arg_type, InvCovariance, true, false) \
		SLOM_MAKE_DEVIATIONCHOOSER(templ, arg_type, StandardDeviation, false, true) \
		SLOM_MAKE_DEVIATIONCHOOSER(templ, arg_type, InvStandardDeviation, false, false)


SLOM_MAKE_ALL_DEVIATIONCHOOSER(class Type, Type)
SLOM_MAKE_ALL_DEVIATIONCHOOSER(class Derived, Eigen::MatrixBase<Derived>)

// TODO does this make sense?
SLOM_MAKE_ALL_DEVIATIONCHOOSER(int dim, SLOM::CholeskyCovariance<dim>)




}  // namespace SLOM



#endif /* COVARIANCE_HPP_ */
