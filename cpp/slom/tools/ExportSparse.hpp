/**
 * @file /mtk-trunk/slom/tools/ExportSparse.hpp
 * @brief Brief description
 * 
 */

#ifndef EXPORTSPARSE_HPP_
#define EXPORTSPARSE_HPP_

#include "../src/Sparse.hpp"

namespace SLOM {

/*
 *
 */
class ExportSparse {
public:
	typedef internal::SparseType SparseType;
	struct Color {
		char r,g,b;
		Color(char r, char g, char b) : r(r), g(g), b(b) {}
		Color(char k = 0) : r(k), g(k), b(k) {}
		operator char*() { return &r; }
	};
	
	typedef std::pair<int,std::vector<int> > ColumnSet;
	typedef std::vector<std::pair<Color,ColumnSet> > ColoringSet;
	
	static void toPGM(const std::string& fname, const SparseType& A, const ColoringSet& cs=ColoringSet(), 
			bool sorted = false, bool sym = false, const Color& bg=Color(255), 
			const Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic,SparseType::Index>* p = 0
			);
	
};

} /* namespace SLOM */
#endif /* EXPORTSPARSE_HPP_ */
