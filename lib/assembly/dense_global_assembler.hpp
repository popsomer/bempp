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
//#include "discrete_dense_boundary_operator.hpp"
#include "../fiber/default_local_assembler_for_integral_operators_on_surfaces.hpp"
//#include "../fiber/default_local_assembler_for_integral_operators_on_surfaces_imp.hpp"
#include "../fiber/quadrature_strategy.hpp"



namespace Fiber
{

/** \cond FORWARD_DECL */
template <typename ResultType> class LocalAssemblerForIntegralOperators;
template <typename ResultType> class LocalAssemblerForPotentialOperators;
//class GeometryFactory;//Peter
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
    DWFALBPeter(std::string strIn, const std::vector<int>& testIndices, const std::vector<std::vector<GlobalDofIndex> >& testGlobalDofs, const std::vector<std::vector<GlobalDofIndex> >& trialGlobalDofs, const std::vector<std::vector<BasisFunctionType> >& testLocalDofWeights, const std::vector<std::vector<BasisFunctionType> >& trialLocalDofWeights, Fiber::LocalAssemblerForIntegralOperators<ResultType>& assembler, arma::Mat<ResultType>& result, MutexType& mutex, arma::Col<ResultType> * solV, std::vector<ResultType> * rhsV, arma::Mat<ResultType> * wm) : //,std::unique_ptr<Fiber::LocalAssemblerForIntegralOperators<ResultType> > asmblr) :
        str(strIn), m_testIndices(testIndices), m_testGlobalDofs(testGlobalDofs), m_trialGlobalDofs(trialGlobalDofs), m_testLocalDofWeights(testLocalDofWeights), m_trialLocalDofWeights(trialLocalDofWeights), m_assembler(assembler), m_result(result), m_mutex(mutex), m_solV(solV), m_rhsV(rhsV), m_wm(wm) { }  //, p_asmblr(asmblr) {    }

//        DWFALBPeter(const DWFALBPeter&) = delete;


    void operator() (const tbb::blocked_range<size_t>& r) const {
        const int elementCount = m_testIndices.size();
        std::vector<arma::Mat<ResultType> > localResult;
	arma::Col<ResultType> solV = *m_solV;
	arma::Mat<ResultType> wm = *m_wm;
/*
std::stringstream s;
s << "testI = "; 
for(int i = 0; i < m_testIndices.size(); ++i){
	s << m_testIndices[i] << " ";
}
s << " = testIndx, size testGlobDofs = "<< m_testGlobalDofs.size()<< " and " << m_testGlobalDofs[0].size() << std::endl;
for(int i = 0; i < m_testGlobalDofs.size(); ++i){
	for(int j = 0; j < m_testGlobalDofs[i].size(); ++j){
		s << m_testGlobalDofs[i][j] << " ";
	}
	s << " / ";
}
s << std::endl;

s << "trialGlobalDofs = " << m_trialGlobalDofs.size()<< " and " << m_trialGlobalDofs[0].size() << std::endl;
for(int i = 0; i < m_trialGlobalDofs.size(); ++i){
	for(int j = 0; j < m_trialGlobalDofs[i].size(); ++j){
		s << m_trialGlobalDofs[i][j] << " ";
	}
	s << " / ";
}
s << std::endl;
std::cout << s.str();
*/
//	std::unique_ptr<Fiber::LocalAssemblerForIntegralOperators<ResultType> > p_asmblr = std::unique_ptr<Fiber::LocalAssemblerForIntegralOperators<ResultType> >(m_assembler);
//	std::unique_ptr<Fiber::LocalAssemblerForIntegralOperators<ResultType> > p_asmblr = &m_assembler;

//	std::cout << r.begin() << "iurogijaosijfd" << r.end() << std::endl;
        for (size_t trialIndex = r.begin(); trialIndex != r.end(); ++trialIndex) {
            // Evaluate integrals over pairs of the current trial element and all the test elements
//            m_assembler.evaluateLocalWeakForms(TEST_TRIAL, m_testIndices, trialIndex, ALL_DOFS, localResult);
//            m_assembler.evaluateLocalWeakFormsPeter(str,TEST_TRIAL, m_testIndices, trialIndex, ALL_DOFS, localResult); //Does not get overridden, use pointer
//            p_asmblr->evaluateLocalWeakFormsPeter(str,TEST_TRIAL, m_testIndices, trialIndex, ALL_DOFS, localResult);
/////		(&m_assembler)->evaluateLocalWeakFormsPeter(str,TEST_TRIAL, m_testIndices, trialIndex, ALL_DOFS, localResult);


//	const GeneralElementarySingularIntegralOperator<BFT,RT,RT> bla = dynamic_cast<const GeneralElementarySingularIntegralOperator<BFT,RT,RT>& > (*asdf);
//		DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<    BasisFunctionType, KernelType, ResultType, GeometryFactory> asCast = dynamic_cast<DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<    BasisFunctionType, KernelType, ResultType, GeometryFactory> > (m_assembler);
//		Fiber::DefaultLocalAssemblerForIntegralOperatorsOnSurfaces asCast = dynamic_cast<Fiber::DefaultLocalAssemblerForIntegralOperatorsOnSurfaces> (m_assembler);
//  typedef typename ScalarTraits<ResultType>::RealType KernelType;
//	typedef CoordinateType KernelType;
//	Fiber::DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<BasisFunctionType, KernelType, ResultType, GeometryFactory> asCast = dynamic_cast<Fiber::DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<BasisFunctionType, KernelType, ResultType, GeometryFactory> > (m_assembler);
//	const Fiber::DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<BasisFunctionType, KernelType, ResultType, GeometryFactory> asCast = dynamic_cast<const Fiber::DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<BasisFunctionType, KernelType, ResultType, GeometryFactory>& > (m_assembler);
	
            const int trialDofCount = m_trialGlobalDofs[trialIndex].size();
            // Global assembly
            {
                MutexType::scoped_lock lock(m_mutex);
                // Loop over test indices
                for (int testIndex = 0; testIndex < elementCount; ++testIndex) {
//			m_result(testIndex, 0) += m_wm[testIndex,trialIndex]*m_solV[trialIndex];
//			m_result(testIndex, 0) += m_wm(testIndex,trialIndex)*m_solV(trialIndex);
//			m_result(testIndex, 0) += wm(testIndex,trialIndex)*solV(trialIndex);
			m_result(trialIndex, 0) += wm(testIndex,trialIndex)*solV(testIndex);
//if (trialIndex == 1) {
//	std::cout << wm(testIndex,trialIndex) << solV(trialIndex) << testIndex << solV(testIndex) << std::endl;
//}
/*
                    const int testDofCount = m_testGlobalDofs[testIndex].size();
//		if ((r.begin() == 1) && (r.end() == 2) && (testIndex == 0) ) {
		if ((trialIndex == 0) && (testIndex == 0) ) {
			std::cout << "locres = " << localResult[testIndex] << std::endl; // Should be (+1.106e-03,+1.960e-04)
			std::cout << testDofCount << "iurogi" << trialIndex << "jaosijfd" << trialDofCount << " aosijfd" << elementCount << std::endl;	
			std::cout << "typeAssembler=" << typeid(m_assembler).name() << '\n';
		}
                    // Add the integrals to appropriate entries in the operator's matrix
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
                            m_result(testGlobalDof, trialGlobalDof) += conj(m_testLocalDofWeights[testIndex][testDof]) * m_trialLocalDofWeights[trialIndex][trialDof] *localResult[testIndex](testDof, trialDof);
                        }
                    }*/
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
//    std::unique_ptr<Fiber::LocalAssemblerForIntegralOperators<ResultType> > p_asmblr;
    // mutable OK because write access to this matrix is protected by a mutex
    arma::Mat<ResultType>& m_result;
    // mutex must be mutable because we need to lock and unlock it
    MutexType& m_mutex;
	arma::Col<ResultType> * m_solV;
	std::vector<ResultType> * m_rhsV;
	arma::Mat<ResultType> * m_wm;
};

// Build a list of lists of global DOF indices corresponding to the local DOFs on each element of space.grid().
//template <typename BasisFunctionType> void gatherGlobalDofs( const Space<BasisFunctionType>& space, std::vector<std::vector<GlobalDofIndex> >& globalDofs, std::vector<std::vector<BasisFunctionType> >& localDofWeights)

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
//    std::auto_ptr<EntityIterator<0> > it = view.entityIterator<0>();
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

/** \ingroup weak_form_assembly_internal
 *  \brief Dense-mode assembler.
 */
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

/*    typedef Fiber::LocalAssemblerForIntegralOperators<ResultType>
    LocalAssemblerForIntegralOperators;

    static std::auto_ptr<DiscreteBoundaryOperator<ResultType> >
    assembleDetachedWeakForm(
            const Space<BasisFunctionType>& testSpace,
            const Space<BasisFunctionType>& trialSpace,
            LocalAssemblerForIntegralOperators& assembler,
            const Context<BasisFunctionType, ResultType>& context); From v2.0*/
//    static std::auto_ptr<DiscreteBoundaryOperator<ResultType> > assembleDetachedWeakFormPeter(std::string str, const Space<BasisFunctionType>& testSpace, const Space<BasisFunctionType>& trialSpace, LocalAssemblerForIntegralOperators& assembler, const Context<BasisFunctionType, ResultType>& context); // Peter
//    static std::auto_ptr<DiscreteBoundaryOperator<ResultType> > assembleDetachedWeakFormPeter(float str, const Space<BasisFunctionType>& testSpace, const Space<BasisFunctionType>& trialSpace, LocalAssemblerForIntegralOperators& assembler, const Context<BasisFunctionType, ResultType>& context); // Peter
//   static std::auto_ptr<DiscreteBoundaryOperator<ResultType> > assembleDetachedWeakFormPeter(const Space<BasisFunctionType>& testSpace, const Space<BasisFunctionType>& trialSpace, LocalAssemblerForIntegralOperators& assembler, const Context<BasisFunctionType, ResultType>& context,std::string str); // Peter


//	template <typename BasisFunctionType, typename ResultType>
//static std::auto_ptr<DiscreteBoundaryOperator<ResultType> >
//DenseGlobalAssembler<BasisFunctionType, ResultType>::
//static std::unique_ptr<DiscreteBoundaryOperator<ResultType> > assembleDetachedWeakFormPeter(std::string str, const Space<BasisFunctionType>& testSpace, const Space<BasisFunctionType>& trialSpace, LocalAssemblerForIntegralOperators& assembler, const Context<BasisFunctionType, ResultType>& context)
//static std::unique_ptr<DiscreteBoundaryOperator<ResultType> > assembleDetachedWeakFormPeter(std::string str, const Space<BasisFunctionType>& testSpace, const Space<BasisFunctionType>& trialSpace, LocalAssemblerForIntegralOperators& assembler, const Context<BasisFunctionType, ResultType>& context, std::unique_ptr<arma::Col<ResultType> > solV, std::unique_ptr<std::vector<ResultType> > rhsV, std::unique_ptr<arma::Mat<ResultType> > wm)
static std::unique_ptr<DiscreteBoundaryOperator<ResultType> >
assembleDetachedWeakFormPeter(std::string str, const Space<BasisFunctionType>& testSpace, const Space<BasisFunctionType>& trialSpace, LocalAssemblerForIntegralOperators& assembler, const Context<BasisFunctionType, ResultType>& context, arma::Col<ResultType> * solV, std::vector<ResultType> * rhsV, arma::Mat<ResultType> * wm)
//static std::unique_ptr<DiscreteBoundaryOperator<ResultType> >
//assembleDetachedWeakFormPeter(std::string str, const Space<BasisFunctionType>& testSpace, const Space<BasisFunctionType>& trialSpace, std::unique_ptr<LocalAssemblerForIntegralOperators> asmblr, const Context<BasisFunctionType, ResultType>& context)
{
//	std::cout << "In denseglobalassembler: " << str << std::endl;
//	LocalAssemblerForIntegralOperators& assembler = *asmblr;
/*   
    // Make a vector of all element indices
    std::vector<int> testIndices(testElementCount);
    for (int i = 0; i < testElementCount; ++i)
        testIndices[i] = i;
*/
	const AssemblyOptions& options = context.assemblyOptions();
    // Global DOF indices corresponding to local DOFs on elements
    std::vector<std::vector<GlobalDofIndex> > testGlobalDofs, trialGlobalDofs;
    std::vector<std::vector<BasisFunctionType> > testLocalDofWeights,
        trialLocalDofWeights;
    ggdsPeter(testSpace, testGlobalDofs, testLocalDofWeights);
//    gatherGlobalDofs(testSpace, testGlobalDofs, testLocalDofWeights);
    if (&testSpace == &trialSpace) {
        trialGlobalDofs = testGlobalDofs;
        trialLocalDofWeights = testLocalDofWeights;
    } else
        ggdsPeter(trialSpace, trialGlobalDofs, trialLocalDofWeights);
//        gatherGlobalDofs(trialSpace, trialGlobalDofs, trialLocalDofWeights);
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

//	std::cout << "Iasdflakjsnldkfjnn denseglobalassembler: "  << std::endl;
//    arma::Mat<ResultType> result(testSpace.globalDofCount(), trialSpace.globalDofCount());
    arma::Mat<ResultType> result(testSpace.globalDofCount(), 1);
    result.fill(0.); // Create and fill the operator's matrix

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
//	std::cout << "asiodjfsapoijdfosaijf" << maxThreadCount << std::endl;
    tbb::task_scheduler_init scheduler(maxThreadCount);
/*
std::cout << "testI = "; // << std::copy(testIndices) 
for(int i = 0; i < testIndices.size(); ++i){
	std::cout << testIndices[i] << " ";
}
std::cout << " = testIndx, size testGlobDofs = "<<testGlobalDofs.size()<< " and " << testGlobalDofs[0].size() << std::endl;// << testGlobalDofs << " = testGlobDofs, trialGlobDofs = " << trialGlobalDofs << std::endl;

for(int i = 0; i < testGlobalDofs.size(); ++i){
//	std::cout << testGlobalDofs[i] << std::endl;
	for(int j = 0; j < testGlobalDofs[i].size(); ++j){
//		std::cout << i << "=i,j=" << j << " ";
		std::cout << testGlobalDofs[i][j] << " ";
	}
//	std::cout << " for i=" << i << std::endl;
//	std::cout << testGlobalDofs[i].size() << " ";
//	std::cout << testGlobalDofs[i][0] << " ";
	std::cout << " / ";
//	std::cout << std::endl;
}
//std::cout << testGlobalDofs[30][1];
std::cout << std::endl;

std::cout << "trialGlobalDofs = " << trialGlobalDofs.size()<< " and " << trialGlobalDofs[0].size() << std::endl;
for(int i = 0; i < trialGlobalDofs.size(); ++i){
	for(int j = 0; j < trialGlobalDofs[i].size(); ++j){
		std::cout << trialGlobalDofs[i][j] << " ";
	}
//	std::cout << testGlobalDofs[i].size() << " ";
	std::cout << " / ";
}
std::cout << std::endl;
*/


//std::cout << "apsojfd" << std::endl;
    {
        Fiber::SerialBlasRegion region;
        tbb::parallel_for(tbb::blocked_range<size_t>(0, trialElementCount), Body(str, testIndices, testGlobalDofs, trialGlobalDofs,testLocalDofWeights, trialLocalDofWeights, assembler, result, mutex, solV, rhsV, wm));
//        tbb::parallel_for(tbb::blocked_range<size_t>(0, trialElementCount), Body(str, testIndices, testGlobalDofs, trialGlobalDofs,testLocalDofWeights, trialLocalDofWeights, assembler, result, mutex,asmblr));
    }

	std::vector<ResultType> rhsVe = *rhsV;
	ResultType tmpErr, re, rh;
	ResultType globalNor = 0.0;
	ResultType globalNorRers = 0.0;
	ResultType globalErr = 0.0;
//std::cout << result << " = result, testp.globdofcoutn = " << testSpace.globalDofCount() << std::endl; //rhsVe << std::endl;

	for (int i=0; i < testSpace.globalDofCount(); ++i) {;
//		std::cout << i << " " << rhsV[i] << " "; // << result[i,0] << " " << std::endl;
		rh = rhsVe[i];
//		re = result[i];
//		re = result[i,0];
		re = result[0,i];
//		std::cout << i << " " << rh << " " << re << " " << std::endl;
		tmpErr = std::abs(re-rh);
		globalNor += std::abs(rh);
/*
		tmpErr = std::pow(std::abs(re-rh),2);
		globalNor += std::pow(std::abs(rh),2);*/
		globalNorRers += std::abs(re);
		globalErr += tmpErr;
	}
/*
//	std::cout << "tmpErrs = ";
//	for(int i=0; i < testSpace.globalDofCount(); ++i) {
	for (int idx = 0; idx < testSpace.globalDofCount(); ++idx ) {
//		tmpErr = std::abs(result[i,0] -rhsV[i]);
std::cout << "apsoijfdpsaoijfd ";
std::cout << idx << " ";
		re = result[idx,0];
//		re = result[i,0];
//		rh = rhsV[i];
//		rh = rhsVe[i];
		rh = rhsVe[idx];
//		tmpErr = std::abs(re-rh)^2;
		tmpErr = std::pow(std::abs(re-rh),2);
//		std::cout << tmpErr << " ";
		globalNor += std::pow(std::abs(rh),2);
//		globalNor += std::abs(rh)^2;
		globalNorRers += std::abs(re);
		globalErr += tmpErr;
	}*/

//	std::cout << std::endl << std::sqrt(globalErr)/std::sqrt(globalNor) <<
	std::cout << std::endl << globalErr/globalNor << " = relErr, GlobalErr = " << globalErr << " , " << globalNor << " = globalNor, globalNorRers = " << globalNorRers << std::endl;
//    return std::auto_ptr<DiscreteBoundaryOperator<ResultType> >(new DiscreteDenseBoundaryOperator<ResultType>(result));
    return std::unique_ptr<DiscreteBoundaryOperator<ResultType> >(
                new DiscreteDenseBoundaryOperator<ResultType>(result));
}

};

} // namespace Bempp

#endif
