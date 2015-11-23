#ifndef bempp_dense_global_assembler_hpp
#define bempp_dense_global_assembler_hpp

#include "../common/common.hpp"
#include "../common/armadillo_fwd.hpp"
#include "../common/scalar_traits.hpp"

#include <memory>


//Peter old:
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

//Peter:
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

// Body of parallel loop
template <typename BasisFunctionType, typename ResultType>
class DWFALBPeter
{
public:
    typedef tbb::spin_mutex MutexType;
	typedef typename ScalarTraits<ResultType>::RealType CoordinateType;
    DWFALBPeter(std::string strIn, const std::vector<int>& testIndices, const std::vector<std::vector<GlobalDofIndex> >& testGlobalDofs, const std::vector<std::vector<GlobalDofIndex> >& trialGlobalDofs, const std::vector<std::vector<BasisFunctionType> >& testLocalDofWeights, const std::vector<std::vector<BasisFunctionType> >& trialLocalDofWeights, Fiber::LocalAssemblerForIntegralOperators<ResultType>& assembler, arma::Mat<ResultType>& result, MutexType& mutex, arma::Col<ResultType> * solV, std::vector<ResultType> * rhsV, arma::Mat<ResultType> * wm, std::vector< Point3D<CoordinateType> > testPos, std::vector< Point3D<CoordinateType> > trialPos, arma::Mat<ResultType>& rois ) : 
        str(strIn), m_testIndices(testIndices), m_testGlobalDofs(testGlobalDofs), m_trialGlobalDofs(trialGlobalDofs), m_testLocalDofWeights(testLocalDofWeights), m_trialLocalDofWeights(trialLocalDofWeights), m_assembler(assembler), m_result(result), m_mutex(mutex), m_solV(solV), m_rhsV(rhsV), m_wm(wm), m_testPos(testPos), m_trialPos(trialPos), m_rois(rois) { }  


    void operator() (const tbb::blocked_range<size_t>& r) const {
        const int elementCount = m_testIndices.size();
        std::vector<arma::Mat<ResultType> > localResult;
	arma::Col<ResultType> solV = *m_solV;
	arma::Mat<ResultType> wm = *m_wm;
std::stringstream toPrint;
toPrint << str << " aposdjf " << r.begin() << " asf " << r.end() << "\n";
//std::cout << toPrint.str();
//std::cout << str << " aposdjf " << r.begin() << " asf " << r.end() << "\n";
        for (size_t trialIndex = r.begin(); trialIndex != r.end(); ++trialIndex) {
            // Evaluate integrals over pairs of the current trial element and all the test elements
//            m_assembler.evaluateLocalWeakForms(TEST_TRIAL, m_testIndices, trialIndex, ALL_DOFS, localResult);	
            m_assembler.evaluateLocalWeakFormsPeter(str,TEST_TRIAL, m_testIndices, trialIndex, ALL_DOFS, localResult);	
            const int trialDofCount = m_trialGlobalDofs[trialIndex].size();
            // Global assembly
            {
                MutexType::scoped_lock lock(m_mutex);
		CoordinateType a = 0.1;
		CoordinateType b = 0.3;
		CoordinateType percDecay = 0.8;
		ResultType maxi = 0.0;
		int idxMax = -1;	
//	int method = 1; // 1=fixedWindows, 4 = correlations, ...
	if (str.at(0) == 'c') {
                // Loop over test indices
                for (int testIndex = 0; testIndex < elementCount; ++testIndex) {
// Correlations:
			for(int windIdx=0; windIdx < elementCount; ++windIdx) {
				CoordinateType dist = std::sqrt(std::pow(m_testPos[testIndex].x-m_trialPos[windIdx].x, 2.0) + std::pow(m_testPos[testIndex].y-m_trialPos[windIdx].y, 2.0) + std::pow(m_testPos[testIndex].z-m_trialPos[windIdx].z, 2.0) );
				ResultType wind = 0;
				if(dist < a) {
					wind = 1;
				}
				else if (dist < b) {
					wind = exp(2*exp(-(b-a)/(dist-a) )/((dist-a)/(b-a) -1));
				}
				m_rois(testIndex, trialIndex) += wm(windIdx, trialIndex)*solV(windIdx)*wind;

if ( ( (trialIndex == 0) & (testIndex == 1) ) | ( (trialIndex == 1) & (testIndex == 1) )| ( (trialIndex == 1) & (testIndex == 0) ) ) {
	std::stringstream s;
	s << trialIndex+testIndex << " " << wind << " " << windIdx << " " << dist << " / " << m_result(testIndex, trialIndex) << std::endl;
}
			} // end comput corr
			if (std::abs(m_rois(testIndex,trialIndex)) > std::abs(maxi) ) {
				maxi = m_rois(testIndex,trialIndex);
				idxMax = testIndex;
			}
                } // All correlations computed in rois, make result
	}
	ResultType thrp = 0.04;
	maxi = maxi*thrp;
	for (int testIndex = 0; testIndex < elementCount; ++testIndex) {
		// Find nearest
		CoordinateType minDist = 1e9;
		ResultType wind = 0.0;
	if (str.at(0) == 'n') {
		wind = 1.0;
	}
	else if (str.at(0) == 't') {
		wind = 1.0;
	}
	else if (str.at(0) == 'f') {
//		a = 0.6;
//		b = 1.1;
//		a = 0.3;
//		b = 0.5;
		percDecay = 0.8;
		std::string::size_type sz;     // alias of size_t
//		b = std::stof(str,&sz);		
//std::cout << str.substr(2) << std::endl;
		b = std::stof(str.substr(2),&sz);		
		a = (1-percDecay)*b;
//		if ((trialIndex == 0) && (testIndex == 0)) {
//			std::cout << b << "=b, a=" << a << "sadf" << std::stof("-1.9",&sz) << std::endl;
//		}
		CoordinateType dist = std::sqrt( std::pow(m_testPos[testIndex].y-m_trialPos[trialIndex].y, 2.0) + std::pow(m_testPos[testIndex].z-m_trialPos[trialIndex].z, 2.0) );
		// Distance without x to include stationary points (but also shadow when coll in illuminated...)
//	if ((m_testPos[testIndex].x < 0) || (m_trialPos[trialIndex].x > 0)) {
		if (dist <= a) {
			wind = 1;
//			std::cout << "Window 1 for fixedWindows with x=" << m_testPos[testIndex].x << std::endl;
		}
		else if (dist <= b) {
			wind = exp(2*exp(-(b-a)/(dist-a) )/((dist-a)/(b-a) -1));
		}
		wind = 1.0; // Multiply by window in kernel
//	}
	}
	else if (str.at(0) == 'c') {
		for (int nei = 0; nei < elementCount; ++nei) {
			CoordinateType dist = std::sqrt(std::pow(m_testPos[testIndex].x-m_trialPos[nei].x, 2.0) + std::pow(m_testPos[testIndex].y-m_trialPos[nei].y, 2.0) + std::pow(m_testPos[testIndex].z-m_trialPos[nei].z, 2.0) );
			if ((std::abs(m_rois(nei,trialIndex)) >= std::abs(maxi) ) && (dist < minDist) ) {
				minDist = dist;
				if (dist <= a) {
					wind = 1;
				}
				else if (dist <= b) {
					wind = exp(2*exp(-(b-a)/(dist-a) )/((dist-a)/(b-a) -1));
				}
			}
		}
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
    // mutable OK because Assembler is thread-safe. (Alternative to "mutable" here:
    // make assembler's internal integrator map mutable)
    typename Fiber::LocalAssemblerForIntegralOperators<ResultType>& m_assembler;
    // mutable OK because write access to this matrix is protected by a mutex
    arma::Mat<ResultType>& m_result;
    arma::Mat<ResultType>& m_rois;
    // mutex must be mutable because we need to lock and unlock it
    MutexType& m_mutex;
	arma::Col<ResultType> * m_solV;
	std::vector<ResultType> * m_rhsV;
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

std::cout << "entered wakFormPeter.\n";

	std::vector< Point3D<CoordinateType> > testPos;
	testSpace.getGlobalDofPositions(testPos);
	std::vector< Point3D<CoordinateType> > trialPos;
	trialSpace.getGlobalDofPositions(trialPos);
std::cout << "got pos.\n";
	std::stringstream s;
	for(int i =0; i < 2; ++i) { 
//		std::cerr << i  << " : " << testPos[i].x << ", " << testPos[i].y << ", " << testPos[i].z << std::endl;
	}
	const AssemblyOptions& options = context.assemblyOptions();
    // Global DOF indices corresponding to local DOFs on elements
    std::vector<std::vector<GlobalDofIndex> > testGlobalDofs, trialGlobalDofs;
    std::vector<std::vector<BasisFunctionType> > testLocalDofWeights,
        trialLocalDofWeights;
std::cout << "starting ggs.\n";
    ggdsPeter(testSpace, testGlobalDofs, testLocalDofWeights);
std::cout << "ended ggs.\n";
    if (&testSpace == &trialSpace) {
        trialGlobalDofs = testGlobalDofs;
        trialLocalDofWeights = testLocalDofWeights;
    } else
        ggdsPeter(trialSpace, trialGlobalDofs, trialLocalDofWeights);
    const int testElementCount = testGlobalDofs.size();
    const int trialElementCount = trialGlobalDofs.size();

    // Enumerate the test elements that contribute to at least one global DOF
    std::vector<int> testIndices;
    testIndices.reserve(testElementCount);
    for (int testIndex = 0; testIndex < testElementCount; ++testIndex) {
        const int testDofCount = testGlobalDofs[testIndex].size();
        for (int testDof = 0; testDof < testDofCount; ++testDof) {
            int testGlobalDof = testGlobalDofs[testIndex][testDof];
            if (testGlobalDof >= 0) {
                testIndices.push_back(testIndex);
                break;
            }
        }
    }

//std::cout << "making operator matrix.\n";
    arma::Mat<ResultType> result;//(testSpace.globalDofCount(), trialSpace.globalDofCount());
//    result.fill(0.); // Create and fill the operator's matrix
//std::cout << "created  operator matrix.\n";

//arma::Mat<ResultType> rois(testSpace.globalDofCount(), trialSpace.globalDofCount());
//    arma::Mat<ResultType> rois;
    int roisSiz = 2;
    if (str.at(0) == 'c') {
	roisSiz = testSpace.globalDofCount();
//	rois(testSpace.globalDofCount(), trialSpace.globalDofCount());
    } else {
//	rois(2,2);
   }
   arma::Mat<ResultType> rois(roisSiz,roisSiz);
   rois.fill(0.);

std::cout << "created rois.\n";
    typedef DWFALBPeter<BasisFunctionType, ResultType> Body;
    typename Body::MutexType mutex;
    const ParallelizationOptions& parallelOptions = options.parallelizationOptions();
    int maxThreadCount = 1;
    if (!parallelOptions.isOpenClEnabled()) {
        if (parallelOptions.maxThreadCount() == ParallelizationOptions::AUTO)
            maxThreadCount = tbb::task_scheduler_init::automatic;
        else
            maxThreadCount = parallelOptions.maxThreadCount();
    }
	
if (str.at(0) == 't') {
arma::Mat<ResultType> oneRow(testSpace.globalDofCount(), 1);
oneRow.fill(0.);
std::vector<arma::Mat<ResultType> > localResult;

std::string::size_type sz;
size_t trialIndex = std::stof(str.substr(2),&sz);	
const int elementCount = testIndices.size();

assembler.evaluateLocalWeakFormsPeter(str,TEST_TRIAL, testIndices, trialIndex, ALL_DOFS, localResult);
const int trialDofCount = trialGlobalDofs[trialIndex].size();

for (int testIndex = 0; testIndex < elementCount; ++testIndex) {
	const int testDofCount = testGlobalDofs[testIndex].size();
	for (int trialDof = 0; trialDof < trialDofCount; ++trialDof) {
		int trialGlobalDof = trialGlobalDofs[trialIndex][trialDof];
		if (trialGlobalDof < 0)
			continue;
		for (int testDof = 0; testDof < testDofCount; ++testDof) {
			int testGlobalDof = testGlobalDofs[testIndex][testDof];
			if (testGlobalDof < 0)
				continue;
//		oneRow(testGlobalDof, trialGlobalDof) += conj(testLocalDofWeights[testIndex][testDof]) * trialLocalDofWeights[trialIndex][trialDof] *localResult[testIndex](testDof, trialDof);
			oneRow(testGlobalDof, 0) += conj(testLocalDofWeights[testIndex][testDof]) * trialLocalDofWeights[trialIndex][trialDof] *localResult[testIndex](testDof, trialDof);
		}
	}
}

return std::unique_ptr<DiscreteBoundaryOperator<ResultType> >(new DiscreteDenseBoundaryOperator<ResultType>(oneRow));
/*
arma::Mat<ResultType> result(testSpace.globalDofCount(), 1);
result.fill(0.);

	tbb::task_scheduler_init scheduler(maxThreadCount);
   {
        Fiber::SerialBlasRegion region;
        tbb::parallel_for(tbb::blocked_range<size_t>(testRow,testRow), Body("n ", testIndices, testGlobalDofs, trialGlobalDofs,testLocalDofWeights, trialLocalDofWeights, assembler, result, mutex, solV, rhsV, wm, testPos, trialPos, rois));
    }
//std::cout << "resRow" << oneRow << " asdf\n";
//   return std::unique_ptr<DiscreteBoundaryOperator<ResultType> >(new DiscreteDenseBoundaryOperator<ResultType>(oneRow));
   return std::unique_ptr<DiscreteBoundaryOperator<ResultType> >(new DiscreteDenseBoundaryOperator<ResultType>(result));
*/
}
else {

    arma::Mat<ResultType> result(testSpace.globalDofCount(), trialSpace.globalDofCount());
    result.fill(0.);

    tbb::task_scheduler_init scheduler(maxThreadCount);
    {
        Fiber::SerialBlasRegion region;
        tbb::parallel_for(tbb::blocked_range<size_t>(0, trialElementCount), Body(str, testIndices, testGlobalDofs, trialGlobalDofs,testLocalDofWeights, trialLocalDofWeights, assembler, result, mutex, solV, rhsV, wm, testPos, trialPos, rois));
   }
   return std::unique_ptr<DiscreteBoundaryOperator<ResultType> >(new DiscreteDenseBoundaryOperator<ResultType>(result));
}

	std::vector<ResultType> rhsVe = *rhsV;
	ResultType tmpErr, re, rh;
	ResultType globalNor = 0.0;
	ResultType globalNorRers = 0.0;
	ResultType globalErr = 0.0;
	ResultType zeroR = 0.0;
	int locmax = -1;
	arma::Col<float> locmaxs(testSpace.globalDofCount());
	locmaxs.fill(-1.0);
	std::stringstream m;
/*
std::cout << testSpace.globalDofCount() << " " << result.size() << " " << trialSpace.globalDofCount() << std::endl;
	for (int i=0; i < testSpace.globalDofCount(); ++i) {
		tmpErr += std::pow(std::abs(result[i,testRow]-oneRow[i,0]),2.0);
		globalNor += std::pow(std::abs(result[i,testRow]),2.0);
std::cout << result[i,testRow] << " =matr, row= " << oneRow[i,0] << " asdf " <<i << " oi " << std::abs(result[i,testRow]) << " asdf " << globalNor << " =gn, term= " << std::pow(std::abs(result[i,testRow]),2.0) << " rowtrans= " << oneRow[0,i] << " rowidx " << oneRow[i,testRow] << " rowtraidx " << oneRow[testRow,i] << "\n";
//		tmpErr = std::abs(result[0,i]);
//		if (std::abs(tmpErr) > std::abs(globalErr)) {
//			globalErr = tmpErr;
//			locmax = i;
//		}
	}
std::cout << std::sqrt(tmpErr) << " =err first row, nor= "  << std::sqrt(globalNor) <<" , testRow= " << testRow << " aspojfd "<<testSpace.globalDofCount()<<"\n";
*/
/*
	int nnz = 0;

//	arma::Mat<bool> qwer = result != 0;
	for (int j = 0; j < trialSpace.globalDofCount(); ++j) {
		globalNor = 0.0;
		for (int i =0; i < testSpace.globalDofCount(); ++i) {
			if (abs(result[i,j]-result(i,j))>0 ) {
				std::cerr << i << "=i,Error square and round indices, j=" << j << result[i,j] << "=sq, round=" << result(i,j) <<std::endl;
			}
			if (abs(result(i,j)) >0) {
//			if (abs(result[i,j]) >0) {
				nnz = nnz + 1;
//				if(qwer[i,j]) {
//				if(result[i,j] == 0) {
				if(result(i,j) == zeroR) {
//					std::cout << i << "=i,j=" << j << result[i,j] << "=result, qwer[i,j] = " << qwer[i,j] << std::endl;
					std::cout << i << "=i,j=" << j << result[i,j] << "=result" << abs(result(i,j)) << std::endl;
				}
			}
//			else if( (result[i,j] != 0)[0] == 1) {
			else if(result(i,j) != zeroR) {
					std::cout << i << "=i,j=" << j << result[i,j] << "=resultAbsZero" << abs(result(i,j)) << std::endl;
			}		
			tmpErr = std::abs(rois(i,j));

			if (j == 186) {
			}
			if (std::abs(tmpErr) > std::abs(globalNor)) {
				globalNor = tmpErr;
				locmaxs(j) = i;
			}
			else if ((j == 20) | (j == 0) ) {
			}
		}
		m << j << " " << locmaxs(j) << " / ";
	}
	std::cout << nnz << " = nnz, nbElem = " << testSpace.globalDofCount()*trialSpace.globalDofCount() << ", percentage = " << nnz/(0.0+testSpace.globalDofCount()*trialSpace.globalDofCount()) << std::endl;
*/
//	std::cout << " are corrs of testidx=0, idxmax = " << locmax << " " << globalErr << std::endl;
//    return std::unique_ptr<DiscreteBoundaryOperator<ResultType> >(
//                new DiscreteDenseBoundaryOperator<ResultType>(result));
}

};

} // namespace Bempp

#endif
