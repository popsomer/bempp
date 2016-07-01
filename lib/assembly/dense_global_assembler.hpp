#ifndef bempp_dense_global_assembler_hpp
#define bempp_dense_global_assembler_hpp

#include "../common/common.hpp"
#include "../common/armadillo_fwd.hpp"
#include "../common/scalar_traits.hpp"

#include <memory>


#include "../fiber/explicit_instantiation.hpp"
#include "assembly_options.hpp"
#include "discrete_dense_boundary_operator.hpp"
#include "context.hpp"

#include "../common/auto_timer.hpp"
#include "../common/multidimensional_arrays.hpp"
#include "../common/not_implemented_error.hpp"
#include "../fiber/explicit_instantiation.hpp"
#include "../fiber/serial_blas_region.hpp"
#include "../fiber/local_assembler_for_integral_operators.hpp"
#include "../grid/entity.hpp"
#include "../grid/entity_iterator.hpp"
#include "../grid/grid.hpp"
#include "../grid/grid_view.hpp"
#include "../grid/mapper.hpp"
#include "../space/space.hpp"

#include "../common/armadillo_fwd.hpp"
#include "../common/complex_aux.hpp"
#include <stdexcept>
#include <iostream>

#include <tbb/parallel_for.h>
#include <tbb/spin_mutex.h>
#include <tbb/task_scheduler_init.h>

#include "../fiber/default_local_assembler_for_integral_operators_on_surfaces.hpp"
#include "../fiber/quadrature_strategy.hpp"

#include "../fiber/geometrical_data.hpp"
#include "../common/common.hpp"


namespace Fiber
{

/** \cond FORWARD_DECL */
template <typename ResultType> class LocalAssemblerForIntegralOperators;
template <typename ResultType> class LocalAssemblerForPotentialOperators;
/** \endcond */

} // namespace Fiber

namespace Bempp
{

namespace
{

// Body of parallel loop for matrix (without rois)
template <typename BasisFunctionType, typename ResultType>
class DWFALBPeter
{
public:
    typedef tbb::spin_mutex MutexType;
	typedef typename ScalarTraits<ResultType>::RealType CoordinateType;
    DWFALBPeter(std::string strIn, const std::vector<int>& testIndices, const std::vector<std::vector<GlobalDofIndex> >& testGlobalDofs, const std::vector<std::vector<GlobalDofIndex> >& trialGlobalDofs, const std::vector<std::vector<BasisFunctionType> >& testLocalDofWeights, const std::vector<std::vector<BasisFunctionType> >& trialLocalDofWeights, Fiber::LocalAssemblerForIntegralOperators<ResultType>& assembler, arma::Mat<ResultType>& result, MutexType& mutex, std::vector< Point3D<CoordinateType> > testPos, std::vector< Point3D<CoordinateType> > trialPos, arma::Mat<ResultType>& rois, ResultType globMax, arma::Mat<ResultType>& rowMax) : 
        str(strIn), m_testIndices(testIndices), m_testGlobalDofs(testGlobalDofs), m_trialGlobalDofs(trialGlobalDofs), m_testLocalDofWeights(testLocalDofWeights), m_trialLocalDofWeights(trialLocalDofWeights), m_assembler(assembler), m_result(result), m_mutex(mutex), m_testPos(testPos), m_trialPos(trialPos), m_rois(rois), m_globMax(globMax), m_rowMax(rowMax) { }  


    void operator() (const tbb::blocked_range<size_t>& r) const {
	const int elementCount = m_testIndices.size();
	std::vector<arma::Mat<ResultType> > localResult;
	CoordinateType percDecay = 0.65;
	CoordinateType thrp = 0.1;
	std::string::size_type sz; // alias of size_t

	for (size_t trialIndex = r.begin(); trialIndex != r.end(); ++trialIndex) {
	    // Evaluate integrals over pairs of the current trial element and all the test elements	
            m_assembler.evaluateLocalWeakFormsPeter(str,TEST_TRIAL, m_testIndices, trialIndex, ALL_DOFS, localResult);	
            const int trialDofCount = m_trialGlobalDofs[trialIndex].size();
	    CoordinateType maxNow = std::abs(thrp*m_globMax);
	    if(str.at(2) == 'l') {
		maxNow = std::abs(thrp*m_rowMax(trialIndex,0));
	    }
            // Global assembly
	    {
                MutexType::scoped_lock lock(m_mutex);
		for (int testIndex = 0; testIndex < elementCount; ++testIndex) {
		    CoordinateType wind = 1.0;
		    
		    if( (str.at(0) == 'n') | (str.at(0) == 't') | (str.at(0) == 'k')  | (str.at(0) == 'j')  | (str.at(0) == 'b') ) {
			// Do nothing if str.at(0) == 'n' | 't' | 'k' | 'j' | 'b' : standard BEM or multiply kernel. Include to avoid 'else'-case later for wrong str
		    } else if ((str.at(0) == 'f')  | (str.at(0) == 'a') ) {
			CoordinateType b = std::stof(str.substr(4),&sz);		
			CoordinateType a = (1-percDecay)*b;
			CoordinateType dist = std::sqrt( std::pow(m_testPos[testIndex].y-m_trialPos[trialIndex].y, 2.0) + std::pow(m_testPos[testIndex].z-m_trialPos[trialIndex].z, 2.0) ); // Distance without x to include stationary points (but also shadow when coll in illuminated...)
			if (dist >= b) {
			   continue; // Window is zero
			}
			else if (dist > a) {
			    wind = exp(2*exp(-(b-a)/(dist-a) )/((dist-a)/(b-a) -1));
			}
		    } else if (str.at(0) == 'i') { // only compression on the illuminated side, 'j' is the version with only modification in the kernel
			CoordinateType b = std::stof(str.substr(4),&sz);		
			CoordinateType a = (1-percDecay)*b;
			CoordinateType dist = std::sqrt( std::pow(m_testPos[testIndex].x-m_trialPos[trialIndex].x, 2.0) + std::pow(m_testPos[testIndex].y-m_trialPos[trialIndex].y, 2.0) + std::pow(m_testPos[testIndex].z-m_trialPos[trialIndex].z, 2.0) ); // Distance without x to include stationary points (but also shadow when coll in illuminated...)
			if (dist > b) {
			   continue; 
			}
			else if ((dist > a) && (m_testPos[testIndex].x <= -0.15) ) {
			    // collocation point x negative is illuminated side
			    wind = exp(2*exp(-(b-a)/(dist-a) )/((dist-a)/(b-a) -1));
			}
			if(wind < std::abs(m_testPos[testIndex].x + 0.15)) {
				wind = std::abs(m_testPos[testIndex].x + 0.15);
			}
		    } else if (str.at(0) == 'd') { // Correlations through physical distance
			CoordinateType b = std::stof(str.substr(4),&sz);
			CoordinateType minDist = 1e9; // Find nearest point
			for (int nei = 0; nei < elementCount; ++nei) {
			    CoordinateType dist = std::sqrt(std::pow(m_testPos[testIndex].x-m_trialPos[nei].x, 2.0) + std::pow(m_testPos[testIndex].y-m_trialPos[nei].y, 2.0) + std::pow(m_testPos[testIndex].z-m_trialPos[nei].z, 2.0) );
			    if ((std::abs(m_rois(nei,trialIndex)) >= std::abs(maxNow) ) && (dist < minDist) ) {
				minDist = dist;
				if (dist < 10.0^(-13) ) {
					wind = 1; // testIndex is itself above threshold
				}
				else if (dist <= b) {
					wind = exp(2.0*exp(-b/dist)/(dist/b -1.0));
				}
			    }
			}
			if(abs(wind) < 10.0^(-13) ) {
			    continue;
			}
			CoordinateType dist = std::sqrt(std::pow(m_testPos[testIndex].x-m_trialPos[trialIndex].x, 2.0) + std::pow(m_testPos[testIndex].y-m_trialPos[trialIndex].y, 2.0) + std::pow(m_testPos[testIndex].z-m_trialPos[trialIndex].z, 2.0) );
			if (dist < b) { // Enforce the Green singularity. Formula below is at most 1: when dist == 0
			    wind = std::max<CoordinateType>(wind,exp(2.0*exp(-b/dist)/(dist/b -1.0)));
			}			
		     } else if( (str.at(0) == 'c') | (str.at(0) == 'e') ) { // Correlations through distance to correlation threshold: might give patches with window not identically one inside but faster
			CoordinateType dist = std::abs(maxNow/m_rois(testIndex,trialIndex));
			CoordinateType b = std::stof(str.substr(4),&sz)*thrp;
			if(dist >= b) {
			   continue;
			} else if (dist > thrp) {
		           wind = exp(2*exp(-(b-thrp)/(dist-thrp) )/((dist-thrp)/(b-thrp) -1));
			}
		    } else {
			throw 20;
		    }
		    const int testDofCount = m_testGlobalDofs[testIndex].size();
		    for (int trialDof = 0; trialDof < trialDofCount; ++trialDof) {
                        int trialGlobalDof = m_trialGlobalDofs[trialIndex][trialDof];
                        if (trialGlobalDof < 0)
                            continue;
                        for (int testDof = 0; testDof < testDofCount; ++testDof) {
                            int testGlobalDof = m_testGlobalDofs[testIndex][testDof];
                            if (testGlobalDof < 0)
                                continue;
                            assert(std::abs(m_testLocalDofWeights[testIndex][testDof]) > 0.);
                            assert(std::abs(m_trialLocalDofWeights[trialIndex][trialDof]) > 0.);
                            m_result(testGlobalDof, trialGlobalDof) += conj(m_testLocalDofWeights[testIndex][testDof]) * m_trialLocalDofWeights[trialIndex][trialDof] *localResult[testIndex](testDof, trialDof)*wind;
			}
		    }
		} // end loop testIndex
	    }
	}
    } // end operator()
private:
    std::string str;
    const std::vector<int>& m_testIndices;
    const std::vector<std::vector<GlobalDofIndex> >& m_testGlobalDofs;
    const std::vector<std::vector<GlobalDofIndex> >& m_trialGlobalDofs;
    const std::vector<std::vector<BasisFunctionType> >& m_testLocalDofWeights;
    const std::vector<std::vector<BasisFunctionType> >& m_trialLocalDofWeights;
    // mutable OK because Assembler is thread-safe. (Alternative to "mutable" here:
    // make assembler's internal integrator map mutable)
    typename Fiber::LocalAssemblerForIntegralOperators<ResultType>& m_assembler;
    // mutable OK because write access to this matrix is protected by a mutex
    arma::Mat<ResultType>& m_result;
    arma::Mat<ResultType>& m_rois;
    // mutex must be mutable because we need to lock and unlock it
    MutexType& m_mutex;
    std::vector< Point3D<CoordinateType> > m_testPos;
    std::vector< Point3D<CoordinateType> > m_trialPos;
    ResultType m_globMax;
    arma::Mat<ResultType> m_rowMax;
};



// Body of parallel loop for computing rois
template <typename BasisFunctionType, typename ResultType>
class DWFALBRois
{
public:
    typedef tbb::spin_mutex MutexType;
	typedef typename ScalarTraits<ResultType>::RealType CoordinateType;
    DWFALBRois(std::string strIn, const std::vector<int>& testIndices, const std::vector<std::vector<GlobalDofIndex> >& testGlobalDofs, const std::vector<std::vector<GlobalDofIndex> >& trialGlobalDofs, const std::vector<std::vector<BasisFunctionType> >& testLocalDofWeights, const std::vector<std::vector<BasisFunctionType> >& trialLocalDofWeights, Fiber::LocalAssemblerForIntegralOperators<ResultType>& assembler, MutexType& mutex, arma::Col<ResultType> * solV, arma::Mat<ResultType> * wm, std::vector< Point3D<CoordinateType> > testPos, std::vector< Point3D<CoordinateType> > trialPos, arma::Mat<ResultType>& rois) : 
        str(strIn), m_testIndices(testIndices), m_testGlobalDofs(testGlobalDofs), m_trialGlobalDofs(trialGlobalDofs), m_testLocalDofWeights(testLocalDofWeights), m_trialLocalDofWeights(trialLocalDofWeights), m_assembler(assembler), m_mutex(mutex), m_solV(solV), m_wm(wm), m_testPos(testPos), m_trialPos(trialPos), m_rois(rois) {} 

    void operator() (const tbb::blocked_range<size_t>& r) const {
        const int elementCount = m_testIndices.size();
        std::vector<arma::Mat<ResultType> > localResult;
	arma::Col<ResultType> solV = *m_solV;
	arma::Mat<ResultType> wm = *m_wm;

	CoordinateType a = 0.1;
	CoordinateType b = 0.3;
        for (size_t testIndex = r.begin(); testIndex != r.end(); ++testIndex) {
            {
                MutexType::scoped_lock lock(m_mutex);
                for(int windIdx=0; windIdx < elementCount; ++windIdx) {
		    CoordinateType dist = std::sqrt(std::pow(m_testPos[testIndex].x-m_trialPos[windIdx].x, 2.0) + std::pow(m_testPos[testIndex].y-m_trialPos[windIdx].y, 2.0) + std::pow(m_testPos[testIndex].z-m_trialPos[windIdx].z, 2.0) );
		    if(dist > b) continue;
		    CoordinateType wind = 1;
		    if(dist > a) wind = exp(2*exp(-(b-a)/(dist-a) )/((dist-a)/(b-a) -1));
		    for (int trialIndex = 0; trialIndex < elementCount; ++trialIndex) {
			m_rois(testIndex, trialIndex) += wm(windIdx, trialIndex)*solV(windIdx)*wind;
		    }
		}
            }
        }
    }
private:
    std::string str;
    const std::vector<int>& m_testIndices;
    const std::vector<std::vector<GlobalDofIndex> >& m_testGlobalDofs;
    const std::vector<std::vector<GlobalDofIndex> >& m_trialGlobalDofs;
    const std::vector<std::vector<BasisFunctionType> >& m_testLocalDofWeights;
    const std::vector<std::vector<BasisFunctionType> >& m_trialLocalDofWeights;
    typename Fiber::LocalAssemblerForIntegralOperators<ResultType>& m_assembler;
    arma::Mat<ResultType>& m_rois;
    MutexType& m_mutex;
    arma::Col<ResultType> * m_solV;
    arma::Mat<ResultType> * m_wm;
    std::vector< Point3D<CoordinateType> > m_testPos;
    std::vector< Point3D<CoordinateType> > m_trialPos;
};


// Build a list of lists of global DOF indices corresponding to the local DOFs on each element of space.grid().
template <typename BasisFunctionType> void ggdsPeter( const Space<BasisFunctionType>& space, std::vector<std::vector<GlobalDofIndex> >& globalDofs, std::vector<std::vector<BasisFunctionType> >& localDofWeights)
{
    // Get the grid's view so that we can iterate over elements
    const GridView& view = space.gridView();
    const int elementCount = view.entityCount(0);
    // Global DOF indices corresponding to local DOFs on elements
    globalDofs.clear();
    globalDofs.resize(elementCount);
    // Weights of the local DOFs on elements
    localDofWeights.clear();
    localDofWeights.resize(elementCount);
    // Gather global DOF lists
    const Mapper& mapper = view.elementMapper();
    std::unique_ptr<EntityIterator<0> > it = view.entityIterator<0>();
    while (!it->finished()) {
        const Entity<0>& element = it->entity();
        const int elementIndex = mapper.entityIndex(element);
        space.getGlobalDofs(element, globalDofs[elementIndex],localDofWeights[elementIndex]);
        it->next();
    }
}
} // namespace Peter */


/** \cond FORWARD_DECL */
class AssemblyOptions;
class EvaluationOptions;
template <typename ValueType> class DiscreteBoundaryOperator;
template <typename BasisFunctionType> class Space;
template <typename BasisFunctionType, typename ResultType> class Context;
/** \endcond */

template <typename BasisFunctionType, typename ResultType>
class DenseGlobalAssembler
{
public:
	typedef typename ScalarTraits<BasisFunctionType>::RealType CoordinateType;
	typedef Fiber::LocalAssemblerForIntegralOperators<ResultType> LocalAssemblerForIntegralOperators;
	typedef Fiber::LocalAssemblerForPotentialOperators<ResultType> LocalAssemblerForPotentialOperators;
	static std::unique_ptr<DiscreteBoundaryOperator<ResultType> > assembleDetachedWeakForm(
                            const Space<BasisFunctionType>& testSpace,
                            const Space<BasisFunctionType>& trialSpace,
                            LocalAssemblerForIntegralOperators& assembler,
                            const Context<BasisFunctionType, ResultType>& context);
	static std::unique_ptr<DiscreteBoundaryOperator<ResultType> > assemblePotentialOperator(
                            const arma::Mat<CoordinateType>& points,
                            const Space<BasisFunctionType>& trialSpace,
                            LocalAssemblerForPotentialOperators& assembler,
                            const EvaluationOptions& options);

static std::unique_ptr<DiscreteBoundaryOperator<ResultType> >
assembleDetachedWeakFormPeter(std::string str, const Space<BasisFunctionType>& testSpace, const Space<BasisFunctionType>& trialSpace, LocalAssemblerForIntegralOperators& assembler, const Context<BasisFunctionType, ResultType>& context, arma::Col<ResultType> * solV, std::vector<ResultType> * rhsV, arma::Mat<ResultType> * wm)
{
	std::vector< Point3D<CoordinateType> > testPos;
	testSpace.getGlobalDofPositions(testPos);
	std::vector< Point3D<CoordinateType> > trialPos;
	trialSpace.getGlobalDofPositions(trialPos);

	const AssemblyOptions& options = context.assemblyOptions();
	// Global DOF indices corresponding to local DOFs on elements
	std::vector<std::vector<GlobalDofIndex> > testGlobalDofs, trialGlobalDofs;
	std::vector<std::vector<BasisFunctionType> > testLocalDofWeights, trialLocalDofWeights;
	ggdsPeter(testSpace, testGlobalDofs, testLocalDofWeights);
	if (&testSpace == &trialSpace) {
	    trialGlobalDofs = testGlobalDofs;
            trialLocalDofWeights = testLocalDofWeights;
	} else {ggdsPeter(trialSpace, trialGlobalDofs, trialLocalDofWeights);}
	const int testElementCount = testGlobalDofs.size();
	const int trialElementCount = trialGlobalDofs.size();

	// Enumerate the test elements that contribute to at least one global DOF
	std::vector<int> testIndices;
	testIndices.reserve(testElementCount);
	for (int testIndex = 0; testIndex < testElementCount; ++testIndex) {
            const int testDofCount = testGlobalDofs[testIndex].size();
	    if (testIndex % (testElementCount/10) == 0) {
	    }
            for (int testDof = 0; testDof < testDofCount; ++testDof) {
		int testGlobalDof = testGlobalDofs[testIndex][testDof];
		if (testGlobalDof >= 0) {
		    testIndices.push_back(testIndex);
		    break;
		}
	    }
	}
	int roisSiz = 2;
	if ((str.at(0) == 'd') | (str.at(0) == 'c') ) {
	    roisSiz = testSpace.globalDofCount();
	} 
	arma::Mat<ResultType> rois(roisSiz,roisSiz);
	rois.fill(0.);

	typename DWFALBPeter<BasisFunctionType, ResultType>::MutexType mutex;
	typename DWFALBRois<BasisFunctionType, ResultType>::MutexType mutexR;
	const ParallelizationOptions& parallelOptions = options.parallelizationOptions();
	int maxThreadCount = 1;
	if (!parallelOptions.isOpenClEnabled()) {
	   if (parallelOptions.maxThreadCount() == ParallelizationOptions::AUTO)
		maxThreadCount = tbb::task_scheduler_init::automatic;
	   else maxThreadCount = parallelOptions.maxThreadCount();
	}
	ResultType globMax = 0;

	arma::Mat<ResultType> rowMax(roisSiz,1);
	rowMax.fill(0.);
	
	if ((str.at(0) == 'd') || (str.at(0) == 'c') ) {
	    tbb::task_scheduler_init scheduler(maxThreadCount);
	    std::cout << "Starting computing rois\n";
	    {
        	Fiber::SerialBlasRegion region;
        	tbb::parallel_for(tbb::blocked_range<size_t>(0, trialElementCount), DWFALBRois<BasisFunctionType, ResultType>(str, testIndices, testGlobalDofs, trialGlobalDofs,testLocalDofWeights, trialLocalDofWeights, assembler, mutex, solV, wm, testPos, trialPos, rois));
	    }
	    std::cout << "Computing maxima\n"; 
	    for(int i = 0; i < roisSiz; i ++) { // Better to make this parallel as well
		for(int j = 0; j < roisSiz; j ++) {
		    if(abs(rois(i,j)) > abs(rowMax(i)) ) {
			rowMax(i) = rois(i,j);
			if(abs(rowMax(i)) > abs(globMax) ) {
			    globMax = rowMax(i);
			}
		    }
		}
	    }
	}
	if ((str.at(0) != 't') & (str.at(0) != 'a') & (str.at(0) != 'b') & (str.at(0) != 'e') ) { // Full solution
	    arma::Mat<ResultType> result(testSpace.globalDofCount(), trialSpace.globalDofCount());
	    result.fill(0.);
	    tbb::task_scheduler_init scheduler(maxThreadCount);
	    std::cout << "Starting assembling system matrix\n";
	    {
        	Fiber::SerialBlasRegion region;
        	tbb::parallel_for(tbb::blocked_range<size_t>(0, trialElementCount), DWFALBPeter<BasisFunctionType, ResultType>(str, testIndices, testGlobalDofs, trialGlobalDofs,testLocalDofWeights, trialLocalDofWeights, assembler, result, mutex, testPos, trialPos, rois, globMax, rowMax));
	    }
	    std::cout << "Ended assembling system matrix\n";
	    return std::unique_ptr<DiscreteBoundaryOperator<ResultType> >(new DiscreteDenseBoundaryOperator<ResultType>(result));
	}
	arma::Mat<ResultType> oneRow(1, trialSpace.globalDofCount() );
	oneRow.fill(0.);
	std::vector<arma::Mat<ResultType> > localResult;

	std::string::size_type sz;
	size_t testIndex = std::stof(str.substr(9),&sz);
	testIndices.resize(1);
	testIndices[0] = testIndex;

	CoordinateType percDecay = 0.8;

	for(size_t trialIndex = 0; trialIndex != trialElementCount; ++trialIndex) {
		assembler.evaluateLocalWeakFormsPeter(str,TEST_TRIAL, testIndices, trialIndex, ALL_DOFS, localResult);
		const int trialDofCount = trialGlobalDofs[trialIndex].size();
		const int testDofCount = testGlobalDofs[testIndex].size();
		CoordinateType wind = 1.0;
		if (str.at(0) == 'a') {
			std::string::size_type szz; // alias of size_t
			CoordinateType b = std::stof(str.substr(4),&szz);		
			CoordinateType a = (1-percDecay)*b;
			CoordinateType dist = std::sqrt( std::pow(testPos[testIndex].y-trialPos[trialIndex].y, 2.0) + std::pow(testPos[testIndex].z-trialPos[trialIndex].z, 2.0) ); // Distance without x to include stationary points (but also shadow when coll in illuminated...)
			if (dist > b) {
			   wind = 0.0;
			}
			else if (dist > a) {
			    wind = exp(2*exp(-(b-a)/(dist-a) )/((dist-a)/(b-a) -1));
			}
		} else if(str.at(0) == 'e') { // Correlations through distance to correlation threshold: might give patches with window not identically one inside but faster
			CoordinateType dist = std::abs(arma::max(arma::max(rois))/rois(0,trialIndex));
			CoordinateType thrp = std::stof(str.substr(4),&sz);
			CoordinateType b = 1.5*thrp;
			if(dist < thrp) {
			   wind = 1;
			} else if (dist < b) {
		           wind = exp(2*exp(-(b-thrp)/(dist-thrp) )/((dist-thrp)/(b-thrp) -1));
			}
		    }
	    	for (int trialDof = 0; trialDof < trialDofCount; ++trialDof) {
			int trialGlobalDof = trialGlobalDofs[trialIndex][trialDof];
			if (trialGlobalDof < 0)
			    continue;
			for (int testDof = 0; testDof < testDofCount; ++testDof) {
			    int testGlobalDof = testGlobalDofs[testIndex][testDof];
			    if (testGlobalDof < 0)
				continue;
			    oneRow(0,trialGlobalDof) += conj(testLocalDofWeights[testIndex][testDof]) * trialLocalDofWeights[trialIndex][trialDof] *localResult[0](testDof, trialDof)*wind; // localResult only has one meaningful element (in the std::vector) which corresponds to the chosen row
			}
		}
	}
	return std::unique_ptr<DiscreteBoundaryOperator<ResultType> >(new DiscreteDenseBoundaryOperator<ResultType>(oneRow));
}


};

} // namespace Bempp

#endif
