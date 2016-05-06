/**
 * @file /mtk-trunk/slom/tools/ExportSparse.cpp
 * @brief Brief description
 * 
 */

#include "ExportSparse.hpp"

#include <fstream>
#include <vector>

namespace SLOM {


void ExportSparse::toPGM(const std::string& fname, const SparseType& A, const ColoringSet& cs, 
		bool sorted, bool sym, const Color& bg,
		const Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic,SparseType::Index>* p) {
	
	int rows = A.rows(), cols = A.cols();
	std::vector<Color> contents(rows*cols, bg); // all white image)
	
	
	if(cs.empty()){
		for(int j=0; j<A.outerSize(); ++j) {
			for(SparseType::InnerIterator iit(A, j); iit; ++iit){
				contents[iit.row()*cols + iit.col()] = Color(0);
			}
		}
	} else {
		int curCol = 0;
		std::vector<std::pair<int, Color> > p_inv;
		p_inv.resize(cols);
		for(ColoringSet::const_iterator cit = cs.begin(); cit != cs.end(); ++cit) {
			const Color& col = cit->first;
			const ColumnSet& cms = cit->second;
			int dof = cms.first;
			for(ColumnSet::second_type::const_iterator it = cms.second.begin(); it!= cms.second.end(); ++it) {
				for(int j=0; j<dof; ++j) {
					int idx = *it+j;
					if(p) idx = p->indices()[idx];
					p_inv[idx] = std::make_pair(curCol++, col);
				}
			}
		}
		assert(cols == curCol);
		
		
		
//		std::vector<Color> contents(rows*cols, bg); // all white image
	 	//	char *data = contents.data();
		for(int j=0; j<A.outerSize(); ++j) {
			for(SparseType::InnerIterator iit(A, j); iit; ++iit){
				int rw = sorted && sym ? p_inv[iit.row()].first : iit.row();
				int cl = sorted        ? p_inv[iit.col()].first : iit.col();
				Color col = p_inv[iit.col()].second;
				if(sym)
					reinterpret_cast<int&>(col) |= reinterpret_cast<const int &>(p_inv[iit.row()].second);
				contents[rw*cols + cl] = col;
			}
		}
#if 0
		for(ColoringSet::const_iterator cit = cs.begin(); cit != cs.end(); ++cit) {
			const Color& col = cit->first;
			const ColumnSet& cms = cit->second;
			int dof = cms.first;
			for(ColumnSet::second_type::const_iterator it = cms.second.begin(); it!= cms.second.end(); ++it) {
				for(int j=0; j<dof; ++j) {
					for(SparseType::InnerIterator iit(A, *it + j); iit; ++iit){
						if(sorted) {
							//							contents[iit.row() * cols + curCol] = col;
						} else {
							contents[iit.row()*cols + iit.col()] = col;
						}
					}
					//					++curCol;
				}
			}
			
		}
#endif
	}
	std::ofstream file(fname.c_str(), std::ios_base::binary);
	file << "P6\n" << cols << " " << rows << " 255\n";
	
	file.write(contents.front(), contents.size()*sizeof(Color));
	
}


} /* namespace SLOM */
