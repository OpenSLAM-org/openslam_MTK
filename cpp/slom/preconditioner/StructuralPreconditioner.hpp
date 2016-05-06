/**
 * @file /mtk/slom/preconditioner/StructuralPreconditioner.hpp
 * @brief Brief description
 * 
 */

#ifndef STRUCTURALPRECONDITIONER_HPP_
#define STRUCTURALPRECONDITIONER_HPP_

#include <deque>

#include "Preconditioner.hpp"
#include "../../mtk/types/vect.hpp"
#include <slom/CallBack.hpp>

namespace SLOM {


// FIXME rename this and make it globally usable?
template<class LmMeas>
struct StructuralPreconditionerTraits;
/*
 *
 */
template<class LandmarkMeasurement>
class StructuralPreconditioner : public Preconditioner {
	typedef StructuralPreconditionerTraits<LandmarkMeasurement> traits;
	typedef typename traits::PoseType PoseType;
	typedef typename traits::LandmarkType LandmarkType;
//	const static typename traits::PoseRefType poseRef = traits::poseRef;
//	const static typename traits::LmRefType lmRef  = traits::lmRef;
	typedef LandmarkMeasurement LmMeas;
	typedef internal::Measurement_Holder<LmMeas> LmMeasHolder;
	
	typedef typename PoseType::scalar Scalar;
	
	typedef SLOM::VarID<PoseType> PoseID;
	typedef SLOM::VarID<LandmarkType> LandmarkID;
	
	
	typedef MTK::vect<PoseType::DOF, Scalar> VectType;
	
	typedef Eigen::Matrix<Scalar, PoseType::DOF, PoseType::DOF> PoseBlock;
	typedef Eigen::Matrix<Scalar, LandmarkType::DOF, LandmarkType::DOF> LmBlock;
	
	// TODO if (Pose::DOF+LM::DOF)*Meas::DIM < Pose::DOF*LM::DOF, it *might* be more efficient to store original Jacobian block
	typedef Eigen::Matrix<Scalar, LandmarkType::DOF, PoseType::DOF> PoseLmBlock;
	
	
	template<class Matrix> 
	struct InplaceLLT {
		typedef typename Matrix::Scalar Scalar;
		enum {DIM = Matrix::RowsAtCompileTime};
		// Theoretically, mat and chol could be stored in a union
		Matrix mat;
		Eigen::LLT<Matrix, Eigen::Upper> chol;
		
		void reset() {
			mat.setZero();
		}
		template<class Derived>
		void rankUpdate(const Eigen::MatrixBase<Derived>& m, const Scalar& scale = Scalar(1.0)){
			assert(m.cols() == DIM);
//			std::cout << mat << "\n" << m << std::endl;
			mat.template selfadjointView<Eigen::Upper>().rankUpdate(m.transpose(), scale);
//			std::cout << mat << "\n-------\n";
		}
		
		void compute() {
			bool success = (chol.compute(mat).info() == Eigen::Success);
			if(!success) {
				std::cerr << "Cholesky decomposition failed:\n" << mat << std::endl;
			}
			assert(success && "Cholesky decomposition failed!");
		}
		template<bool transpose,class Mtrx>
		void apply(const Mtrx& ref) const { // FIXME replace constref by Ref-class
			// TODO find out when to use which:
			if(transpose)
				chol.matrixL().solveInPlace(ref);
			else
				chol.matrixU().solveInPlace(ref);
		}
	};
	
	struct PoseLink{
		PoseID pose; // required?
		InplaceLLT<PoseBlock> llt;
		PoseBlock link_t; // link from current to previous pose, i.e. A_(t-1, t) = A(t, t-1)^T; or L(t,t-1)^T = U(t-1,t)
//		PoseLink* prev; // TODO could be replaced by parent to allow for a tree-like back-bone
		
		PoseLink(const PoseID& pose) : pose(pose) {}
		
		void reset() {
			llt.reset();
			link_t.setZero();
		}
	};
	
	typedef std::deque<PoseLink, Eigen::aligned_allocator<PoseLink> > PoseContainer;
	
	struct PoseLmLink {
		PoseLmBlock mat;
		const LmMeasHolder* meas;
		PoseLink* pose;
		int poseNumber;  // for quick check of neighboring poses
		
		PoseLmLink(const LmMeasHolder *meas, PoseLink* pose, int poseNumber) : meas(meas), pose(pose), poseNumber(poseNumber) {}
		
	};
	
	typedef std::deque<PoseLmLink, Eigen::aligned_allocator<PoseLmLink> > PoseLmLinks;
	
	struct LandMarkEntry {
		LandmarkID landmark;
		InplaceLLT<LmBlock> llt;
		
		PoseLmLinks links;
		
		LandMarkEntry(const LandmarkID& id) : landmark(id) {}
		
		bool operator<(const LandmarkID &id) { return landmark.index() < id.index(); }
		
	};
	
	typedef std::deque<LandMarkEntry, Eigen::aligned_allocator<LandMarkEntry> > LandmarkContainer;
	
	LandmarkContainer landmarks;
	
	
	PoseContainer poses;
	
	
	void update(){
#ifdef SLOM_JACOBI_BLOCKS
		CALLBACK_START_M();
		// FIXME FATAL Check for getsOptimized() everywhere!
		getJ();
		// reset poses
		typename PoseContainer::iterator ps = poses.begin(), ps_end=poses.end();
		for(; ps != ps_end; ++ps) {
			ps->reset();
			// TODO additionally add odometry links if available!
		}
		// TODO this loop could be run in parallel, if PoseLinks are properly guarded
		for(typename LandmarkContainer::iterator lm = landmarks.begin(); lm!=landmarks.end(); ++lm) {
			InplaceLLT<LmBlock> &llt = lm->llt; 
			llt.reset();
			typename PoseLmLinks::iterator ms, last_ms, ms_end=lm->links.end();
			for(ms = lm->links.begin(); ms != ms_end; ++ms) {
				const LmMeasHolder& meas = *ms->meas;
				// FIXME replace (&LandmarkMeasurement::...) by traits members
				llt.rankUpdate(meas.jacSubblock(&LandmarkMeasurement::lm));
				if(ms->pose->pose.getsOptimized()){
					ms->pose->llt.rankUpdate(meas.jacSubblock(&LandmarkMeasurement::p));
					ms->mat = meas.jacSubblock(&LandmarkMeasurement::lm).transpose() * meas.jacSubblock(&LandmarkMeasurement::p);
				} else {
					cb(10) << "Hit fixed pose\n";
				}
			}
			llt.compute();
			if(!0){
				// FIXME dirty pointer hack for testing
				internal::RVHolder<LandmarkType> *lmh = *reinterpret_cast<internal::RVHolder<LandmarkType>**>(&lm->landmark);
				lmh->calcBlockJacobi();
				LmBlock R1 = llt.chol.matrixU(), R2 = lmh->jacobiBlock.matrixU();
				// NOTE This is identical:
				cb(10) << R1 << ". Diff: " << (R1 - R2).array().abs().maxCoeff() << std::endl;
			}
			
			for(last_ms = ms = lm->links.begin(); ms != ms_end; last_ms = ms, ++ms) {
				if(!ms->pose->pose.getsOptimized()) continue;
				llt.template apply<true>(ms->mat);
				if(ms->poseNumber == last_ms->poseNumber+1) {
					cb(10) << "LM " << lm - landmarks.begin() << " links pose " << ms->poseNumber << " to pose " << last_ms->poseNumber;
					if(last_ms->pose->pose.getsOptimized()){
						ms->pose->link_t -= last_ms->mat.transpose() * ms->mat;
						cb(10) << "\t link_t=" << ms->pose->link_t  << std::endl;
					} else {
						cb(10) << "\t pose " << last_ms->poseNumber << " is fixed\n";
					}
				}
				ms->pose->llt.rankUpdate(ms->mat, -1);
			}
		}
		
		for(ps = poses.begin(); ps!=ps_end && !ps->pose.getsOptimized();++ps);
		
		if(ps == ps_end){
			cb(10) << "All poses are fixed!\n";
			return;
		}
		ps->llt.compute();
//		cb(10) << "Optimized? " << ps->pose.getsOptimized() << "\n"; 
		
		for(++ps; ps != poses.end(); ++ps) {
			assert((ps-1)->pose.getsOptimized());
			int idx = ps - poses.begin();
			cb(10) << idx - 1 << ": " << ps->link_t;
			(ps-1)->llt.template apply<false>(ps->link_t);
			ps->llt.rankUpdate(ps->link_t, -1);
			cb(10) << " --> " << ps->link_t << "\tllt= " << ps->llt.mat << std::endl; 
			ps->llt.compute();
		}
#else
	assert(false && "You must enable SLOM_JACOBI_BLOCKS for this to work");
#endif
	}
	
	template<bool transpose>
	void applyJL(Eigen::VectorXd& res) const {
		CALLBACK_START_M();
		// TODO can be run parallel
		for(typename LandmarkContainer::const_iterator lm = landmarks.begin(); lm!=landmarks.end(); ++lm) {
			const InplaceLLT<LmBlock> &llt = lm->llt;
			MTK::vectview<double, LandmarkType::DOF> lmRes(res.data() + lm->landmark.index());
			//= res.template segment<LandmarkType::DOF>(lm->landmark.index());
			cb(15) << lmRes;
			llt.template apply<transpose>(lmRes);
			cb(15) << " " << lmRes << "\n";
		}
	}
	template<bool transpose>
	void applyJC(Eigen::VectorXd& res) const {
		CALLBACK_START_M();
		for(typename LandmarkContainer::const_iterator lm = landmarks.begin(); lm!=landmarks.end(); ++lm) {
			cb(15) << "Landmark " << lm - landmarks.begin() << std::endl;
			for(typename PoseLmLinks::const_iterator ms = lm->links.begin(); ms != lm->links.end(); ++ms) {
				if(!ms->pose->pose.getsOptimized()) {cb(15) << "\tSkipping fixed pose\n"; continue; }
				MTK::vectview<Scalar, PoseType::DOF> posRes(res.data()+ms->pose->pose.index());
				MTK::vectview<Scalar, LandmarkType::DOF> lmRes(res.data()+lm->landmark.index());
				cb(15) << "ms->mat " << ms->mat << ", posRes " << posRes << ", lmRes " << lmRes;
				if(transpose) {
					posRes -= ms->mat.transpose() * lmRes;
					cb(15) << ", posRes " << posRes << std::endl;
				} else {
					lmRes -= ms->mat * posRes;
					cb(15) << ", lmRes " << lmRes << std::endl;
				}
			}
		}
	}
	template<bool transpose>
	void applyJP(Eigen::VectorXd& res) const {
		CALLBACK_START_M()
		if(poses.empty()) return;
		if(transpose){
			typename PoseContainer::const_iterator ps=poses.begin(), ps_prev = ps, ps_end= poses.end();
			for( ; ps != ps_end; ps_prev = ps, ++ps) {
				if(!ps->pose.getsOptimized()) { cb(15) << "Skipping pose " << ps - poses.begin() << std::endl; continue;}
				MTK::vectview<Scalar, PoseType::DOF> curRes(res.data() + ps->pose.index());
				cb(15) << "curRes=" << curRes;
				if(ps!=ps_prev){
					MTK::vectview<Scalar, PoseType::DOF> prevRes(res.data() + ps_prev->pose.index());
					cb(15) << " - " << ps->link_t.transpose() << " * " << prevRes;
					curRes -= ps->link_t.transpose() * prevRes;
					cb(15) << " = " << curRes;
				}
				ps->llt.template apply<true>(curRes);
				cb(15) << "; --llt--> " << curRes << std::endl;
			}
		} else {
			typename PoseContainer::const_reverse_iterator ps=poses.rbegin(), ps_end= poses.rend();
			for( ; ps != ps_end; ++ps) {
				if(!ps->pose.getsOptimized()) { cb(15) << "Skipping pose " << poses.rend() - ps << std::endl; continue;}
				MTK::vectview<Scalar, PoseType::DOF> curRes(res.data() + ps->pose.index());
				cb(15) << "curRes = " << curRes;
				ps->llt.template apply<false>(curRes);
				cb(15) << " --llt--> " << curRes;
				if(ps+1==ps_end || !(ps+1)->pose.getsOptimized()) { cb(15) << "Stop\n"; break;}
				MTK::vectview<Scalar, PoseType::DOF> nextRes(res.data() + (ps+1)->pose.index());
				cb(15) << "nextRes = " << nextRes << " - " << ps->link_t << " * " << curRes;
				nextRes -= ps->link_t*curRes;
				cb(15) << " = " << nextRes << std::endl;
			}
			
		}
	}
	
public:
	
	
	
	StructuralPreconditioner() {};
	
	void compute() {
		if(haveValuesChanged()){
			update();
		}
	}
	
	void addPose(SLOM::VarID<PoseType> id) {
		assert(id.getsOptimized() && "Don't add fixed variables");
		if(!poses.empty() && (!(poses.back().pose < id) || !true)) {
			CALLBACK_START_M();
			cb(15) << poses.back().pose.index() << " >= "<< id.index() << std::endl;
//			assert((poses.empty() || poses.back().pose < id ) && "Poses must be inserted in order!");
		}
		poses.push_back(id);
	}
	
	void addLandmark(SLOM::VarID<LandmarkType> id) {
		if(!landmarks.empty() && !(landmarks.back().landmark < id)){
			CALLBACK_START_M();
			cb(15) << "StrPre> addLandmark(id= " << id.index() << " )";
			cb(15) << "\t" << landmarks.back().landmark.index() << "\n";
		}
//		assert((landmarks.empty() || landmarks.back().landmark < id ) && "Landmarks must be inserted in order!");
		landmarks.push_back(id);
	}
	
	void addLandmarkMeasurement(SLOM::MeasID<LandmarkMeasurement> id) {
		const LandmarkMeasurement& meas = id.getHolder().m_; // FIXME add getter into holder or directly into MeasID
		assert(meas.p == poses.back().pose);
//		assert(meas.lm < landmarks.back().landmark);
		typename LandmarkContainer::iterator lm = std::lower_bound(landmarks.begin(), landmarks.end(), meas.lm);
		assert(lm != landmarks.end() && "Did not find landmark");
		if(lm->landmark != meas.lm){
			CALLBACK_START_M();
			cb(15) << "StrPre> " << lm->landmark.index() << " != " << meas.lm.index() << std::endl;
		}
		lm->links.push_back(PoseLmLink(&id.getHolder(), &poses.back(), poses.size()));
	}
	
	void apply(Eigen::VectorXd& vect, bool transpose) const {
		CALLBACK_START_M();
		// FIXME ***RETHINK*** which order these should be applied
		if(transpose){
			applyJL<true>(vect);
			applyJC<true>(vect);
			applyJP<true>(vect);
		} else {
			applyJP<false>(vect);
			applyJC<false>(vect);
			applyJL<false>(vect);
		}
	}
};

} /* namespace SLOM */
#endif /* STRUCTURALPRECONDITIONER_HPP_ */
