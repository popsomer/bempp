// Copyright (C) 2011-2012 by the Bem++ Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#ifndef fiber_default_local_assembler_for_integral_operators_on_surfaces_hpp
#define fiber_default_local_assembler_for_integral_operators_on_surfaces_hpp

#include "../common/common.hpp"

#include "local_assembler_for_integral_operators.hpp"

#include "_2d_array.hpp"
#include "accuracy_options.hpp"
#include "default_local_assembler_for_operators_on_surfaces_utilities.hpp"
#include "element_pair_topology.hpp"
#include "numerical_quadrature.hpp"
#include "parallelization_options.hpp"
#include "shared_ptr.hpp"
#include "test_kernel_trial_integrator.hpp"
#include "verbosity_level.hpp"

#include <boost/static_assert.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/mutex.h>
#include <cstring>
#include <climits>
#include <set>
#include <utility>
#include <vector>

//Peter:
#include "nonseparable_numerical_test_kernel_trial_integrator.hpp" // For casting

namespace Fiber {

/** \cond FORWARD_DECL */
class OpenClHandler;
template <typename CoordinateType> class CollectionOfShapesetTransformations;
template <typename ValueType> class CollectionOfKernels;
template <typename BasisFunctionType, typename KernelType, typename ResultType>
class TestKernelTrialIntegral;

template <typename CoordinateType> class RawGridGeometry;

template <typename CoordinateType>
class QuadratureDescriptorSelectorForIntegralOperators;
template <typename CoordinateType> class DoubleQuadratureRuleFamily;
/** \endcond */

template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
class DefaultLocalAssemblerForIntegralOperatorsOnSurfaces
    : public LocalAssemblerForIntegralOperators<ResultType> {
public:
  typedef typename ScalarTraits<ResultType>::RealType CoordinateType;

  DefaultLocalAssemblerForIntegralOperatorsOnSurfaces(
      const shared_ptr<const GeometryFactory> &testGeometryFactory,
      const shared_ptr<const GeometryFactory> &trialGeometryFactory,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &testRawGeometry,
      const shared_ptr<const RawGridGeometry<CoordinateType>> &trialRawGeometry,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          testShapesets,
      const shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>> &
          trialShapesets,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &testTransformations,
      const shared_ptr<const CollectionOfKernels<KernelType>> &kernel,
      const shared_ptr<const CollectionOfShapesetTransformations<
          CoordinateType>> &trialTransformations,
      const shared_ptr<const TestKernelTrialIntegral<
          BasisFunctionType, KernelType, ResultType>> &integral,
      const shared_ptr<const OpenClHandler> &openClHandler,
      const ParallelizationOptions &parallelizationOptions,
      VerbosityLevel::Level verbosityLevel, bool cacheSingularIntegrals,
      const shared_ptr<const QuadratureDescriptorSelectorForIntegralOperators<
          CoordinateType>> &quadDescSelector,
      const shared_ptr<const DoubleQuadratureRuleFamily<CoordinateType>> &
          quadRuleFamily);
//template <typename BasisFunctionType, typename KernelType, typename ResultType,         typename GeometryFactory>
//DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<   BasisFunctionType, KernelType, ResultType, GeometryFactory>::
//    DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<BasisFunctionType, KernelType, ResultType, GeometryFactory>( DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<BasisFunctionType, KernelType, ResultType, GeometryFactory> ds) {
/*    DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<BasisFunctionType, KernelType, ResultType, GeometryFactory> (const DefaultLocalAssemblerForIntegralOperatorsOnSurfaces<BasisFunctionType, KernelType, ResultType, GeometryFactory>& ds) {
        m_testGeometryFactory = ds.m_testGeometryFactory;
      m_trialGeometryFactory = ds.m_trialGeometryFactory;
      m_testRawGeometry = ds.m_testRawGeometry;
	m_trialRawGeometry = ds.m_trialRawGeometry;
      m_testShapesets = ds.m_testShapesets;
	m_trialShapesets = ds.m_trialShapesets;
      m_testTransformations = ds.m_testTransformations;	
	m_kernels = ds.m_kernels;
      m_trialTransformations = ds.m_trialTransformations;
	m_integral = ds.m_integral;
      m_openClHandler = ds.m_openClHandler;
      m_parallelizationOptions = ds.m_parallelizationOptions;
      m_verbosityLevel = ds.m_verbosityLevel;
	m_quadDescSelector = ds.m_quadDescSelector;
      m_quadRuleFamily = ds.m_quadRuleFamily;
Utilities::checkConsistencyOfGeometryAndShapesets(*ds.m_testRawGeometry, *ds.m_testShapesets);
  Utilities::checkConsistencyOfGeometryAndShapesets(*ds.m_trialRawGeometry, *ds.m_trialShapesets);
//  if (cacheSingularIntegrals)
//    cacheSingularLocalWeakForms();
}*/
  virtual ~DefaultLocalAssemblerForIntegralOperatorsOnSurfaces();

public:
  virtual void
  evaluateLocalWeakForms(CallVariant callVariant,
                         const std::vector<int> &elementIndicesA,
                         int elementIndexB, LocalDofIndex localDofIndexB,
                         std::vector<arma::Mat<ResultType>> &result,
                         CoordinateType nominalDistance = -1.);



//Peter, in local_as...:
/*  void evaluateLocalWeakFormsPeter(std::string str, CallVariant callVariant,
                         const std::vector<int> &elementIndicesA,
                         int elementIndexB, LocalDofIndex localDofIndexB,
                         std::vector<arma::Mat<ResultType>> &result,
                         CoordinateType nominalDistance = -1.)  {
*/
	void evaluateLocalWeakFormsPeter(std::string str, CallVariant callVariant,
                           const std::vector<int> &elementIndicesA,
                           int elementIndexB, LocalDofIndex localDofIndexB,
                           std::vector<arma::Mat<ResultType>> &result,
                           CoordinateType nominalDistance) {
	if (elementIndexB == 0) {
//		std::cout << "djfoij " << str << std::endl;
	}
	typedef Shapeset<BasisFunctionType> Shapeset;

	const int elementACount = elementIndicesA.size();
	result.resize(elementACount);
	// TODO: remove this unnecessary copy
	const std::vector<const Shapeset *> &m_basesA =
		callVariant == TEST_TRIAL ? *m_testShapesets : *m_trialShapesets;
	const std::vector<const Shapeset *> &m_basesB =
		callVariant == TEST_TRIAL ? *m_trialShapesets : *m_testShapesets;
	std::vector<const Shapeset *> basesA(elementACount);
	for (int i = 0; i < elementACount; ++i)
		basesA[i] = m_basesA[elementIndicesA[i]];
	const Shapeset &basisB = *m_basesB[elementIndexB];


// Find cached matrices; select integrators to calculate non-cached ones
  typedef std::pair<const Integrator *, const Shapeset *> QuadVariant;
  const QuadVariant CACHED(0, 0);
  std::vector<QuadVariant> quadVariants(elementACount);
  for (int i = 0; i < elementACount; ++i) {
    // Try to find matrix in cache
    const arma::Mat<ResultType> *cachedLocalWeakForm = 0;
    if (callVariant == TEST_TRIAL) {
	if ((elementIndexB == 0) && (i==0)){
//		std::cout << "TestTrial" << &m_cache << std::endl;
	}
      const int testElementIndex = elementIndicesA[i];
      const int trialElementIndex = elementIndexB;
      for (size_t n = 0; n < m_cache.extent(0); ++n)
        if (m_cache(n, trialElementIndex).first == testElementIndex) {
          cachedLocalWeakForm = &m_cache(n, trialElementIndex).second;
          break;
        }
    } else {
	if ((elementIndexB == 0) && (i==0)){
//		std::cout << "No TestTrial" << std::endl;
	}
      const int testElementIndex = elementIndexB;
      const int trialElementIndex = elementIndicesA[i];
      for (size_t n = 0; n < m_cache.extent(0); ++n)
        if (m_cache(n, trialElementIndex).first == testElementIndex) {
          cachedLocalWeakForm = &m_cache(n, trialElementIndex).second;
          break;
        }
    }
/*
    if (cachedLocalWeakForm) { // Matrix found in cache
      quadVariants[i] = CACHED;
      if (localDofIndexB == ALL_DOFS) {
        result[i] = *cachedLocalWeakForm;
	if ((elementIndexB == 0) && (i==0)){
		std::cout << *cachedLocalWeakForm << "Matrix found in cache ALL_DOFS" << result[i] << std::endl;
	}
	}
      else {
        if (callVariant == TEST_TRIAL) {
	if ((elementIndexB == 0) && (i==0)){
		std::cout << "Matrix found in cache TEST_TRIAL" << std::endl;
	}
          result[i] = cachedLocalWeakForm->col(localDofIndexB);
	}
        else{
	if ((elementIndexB == 0) && (i==0)){
		std::cout << "Matrix found in cache NO test_trial" << std::endl;
	}
          result[i] = cachedLocalWeakForm->row(localDofIndexB);
	}
      }
    } else {
*/
	if ((elementIndexB == 0) && (i==0)){
//		std::cout << "No matrix found in cache" << std::endl;
//		std::cout << "Forced recomputation of cached matrix by commenting if(cache..)" << std::endl;
	}
      const Integrator *integrator =
          callVariant == TEST_TRIAL
              ? &selectIntegrator(elementIndicesA[i], elementIndexB,
                                  nominalDistance)
              : &selectIntegrator(elementIndexB, elementIndicesA[i],
                                  nominalDistance);
      quadVariants[i] = QuadVariant(integrator, basesA[i]);
//    }
  }


// Integration will proceed in batches of test elements having the same
  // "quadrature variant", i.e. integrator and shapeset

  // Find all the unique quadrature variants present
  typedef std::set<QuadVariant> QuadVariantSet;
  // Set of unique quadrature variants
  QuadVariantSet uniqueQuadVariants(quadVariants.begin(), quadVariants.end());

  std::vector<int> activeElementIndicesA;
  activeElementIndicesA.reserve(elementACount);
  std::vector<arma::Mat<ResultType> *> activeLocalResults;
  activeLocalResults.reserve(elementACount);

  // Now loop over unique quadrature variants
  for (typename QuadVariantSet::const_iterator it = uniqueQuadVariants.begin();
       it != uniqueQuadVariants.end(); ++it) {
    const QuadVariant activeQuadVariant = *it;
    if (activeQuadVariant == CACHED)
      continue;
    const Integrator &activeIntegrator = *it->first;

	if ((elementIndexB == 0) && (it == uniqueQuadVariants.begin()) ) {
//		std::cout << "type actInt = " << typeid(activeIntegrator).name() << std::endl;
//		const NonseparableNumericalTestKernelTrialIntegrator<BasisFunctionType, KernelType, ResultType, GeometryFactory> bla = dynamic_cast<const NonseparableNumericalTestKernelTrialIntegrator<BasisFunctionType, KernelType, ResultType, GeometryFactory> & > (activeIntegrator);
	}
    const Shapeset &activeBasisA = *it->second;

    // Find all the test elements for which quadrature should proceed
    // according to the current quadrature variant
    activeElementIndicesA.clear();
    activeLocalResults.clear();
    for (int indexA = 0; indexA < elementACount; ++indexA)
      if (quadVariants[indexA] == activeQuadVariant) {
        activeElementIndicesA.push_back(elementIndicesA[indexA]);
        activeLocalResults.push_back(&result[indexA]);
      }

	if ((elementIndexB == 0) && (it==uniqueQuadVariants.begin() )){
//		std::cout << "not activeQuadvariant=cached, actLocRes = "<< *activeLocalResults[0] << std::endl;
	}
    // Integrate!
//    activeIntegrator.integrate(callVariant, activeElementIndicesA,
//                               elementIndexB, activeBasisA, basisB,
//                               localDofIndexB, activeLocalResults);
    activeIntegrator.integratePeter(str,callVariant, activeElementIndicesA, elementIndexB, activeBasisA, basisB, localDofIndexB, activeLocalResults);

	if ((elementIndexB == 0) && (it==uniqueQuadVariants.begin() )){
//		std::cout << "result actLocRes = "<< *activeLocalResults[0] << std::endl;
	}
    // // Distribute the just calculated integrals into the result array
    // // that will be returned to caller
    // int i = 0;
    // for (int indexA = 0; indexA < elementACount; ++indexA)
    //     if (quadVariants[indexA] == activeQuadVariant)
    //         result[indexA] = localResult.slice(i++);
  }

}

  virtual void
  evaluateLocalWeakForms(const std::vector<int> &testElementIndices,
                         const std::vector<int> &trialElementIndices,
                         Fiber::_2dArray<arma::Mat<ResultType>> &result,
                         CoordinateType nominalDistance = -1.);

  virtual CoordinateType estimateRelativeScale(CoordinateType minDist) const;

/*
GeometricalData<CoordinateType> getTestGeomData() {
	return m_integral->getTestGeomData();
}
GeometricalData<CoordinateType> getTrialGeomData() {
	return m_integral->getTrialGeomData();
}*/

private:
  /** \cond PRIVATE */
  typedef TestKernelTrialIntegrator<BasisFunctionType, KernelType, ResultType>
  Integrator;
  typedef typename Integrator::ElementIndexPair ElementIndexPair;
  typedef DefaultLocalAssemblerForOperatorsOnSurfacesUtilities<
      BasisFunctionType> Utilities;

  /** \brief Alternative comparison functor for pairs.
   *
   *  This functor can be used to sort pairs according to the second member
   *  and then, in case of equality, according to the first member. */
  template <typename T1, typename T2> struct alternative_less {
    bool operator()(const std::pair<T1, T2> &a,
                    const std::pair<T1, T2> &b) const {
      return a.second < b.second ||
             (!(b.second < a.second) && a.first < b.first);
    }
  };

  /** \brief Comparison functor for element-index pairs.
   *
   *  This functor sorts element index pairs first after the trial element
   *  index (second member) and then, in case of equality, after the test
   *  element index (first member) */
  typedef alternative_less<typename ElementIndexPair::first_type,
                           typename ElementIndexPair::second_type>
  ElementIndexPairCompare;

  /** \brief Set of element index pairs.
   *
   *  The alternative sorting (first after the trial element index) is used
   *  because profiling has shown that evaluateLocalWeakForms is called more
   *  often in the TEST_TRIAL mode (with a single trial element index) than
   *  in the TRIAL_TEST mode. Therefore the singular integral cache is
   *  indexed with trial element index, and this sorting mode makes it easier
   *  to construct such cache. */
  typedef std::set<ElementIndexPair, ElementIndexPairCompare>
  ElementIndexPairSet;

  bool testAndTrialGridsAreIdentical() const;

  void cacheSingularLocalWeakForms();
  void findPairsOfAdjacentElements(ElementIndexPairSet &pairs) const;
  void cacheLocalWeakForms(const ElementIndexPairSet &elementIndexPairs);

  const Integrator &selectIntegrator(int testElementIndex,
                                     int trialElementIndex,
                                     CoordinateType nominalDistance = -1.);

  const Integrator &getIntegrator(const DoubleQuadratureDescriptor &index);

private:
  shared_ptr<const GeometryFactory> m_testGeometryFactory;
  shared_ptr<const GeometryFactory> m_trialGeometryFactory;
  shared_ptr<const RawGridGeometry<CoordinateType>> m_testRawGeometry;
  shared_ptr<const RawGridGeometry<CoordinateType>> m_trialRawGeometry;
  shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
  m_testShapesets;
  shared_ptr<const std::vector<const Shapeset<BasisFunctionType> *>>
  m_trialShapesets;
  shared_ptr<const CollectionOfShapesetTransformations<CoordinateType>>
  m_testTransformations;
  shared_ptr<const CollectionOfKernels<KernelType>> m_kernels;
  shared_ptr<const CollectionOfShapesetTransformations<CoordinateType>>
  m_trialTransformations;
  shared_ptr<const TestKernelTrialIntegral<BasisFunctionType, KernelType,
                                           ResultType>> m_integral;
  shared_ptr<const OpenClHandler> m_openClHandler;
  ParallelizationOptions m_parallelizationOptions;
  VerbosityLevel::Level m_verbosityLevel;
  shared_ptr<const QuadratureDescriptorSelectorForIntegralOperators<
      CoordinateType>> m_quadDescSelector;
  shared_ptr<const DoubleQuadratureRuleFamily<CoordinateType>> m_quadRuleFamily;

  typedef tbb::concurrent_unordered_map<DoubleQuadratureDescriptor,
                                        Integrator *> IntegratorMap;
  IntegratorMap m_testKernelTrialIntegrators;
  mutable tbb::mutex m_integratorCreationMutex;

  enum {
    INVALID_INDEX = INT_MAX
  };
  typedef _2dArray<std::pair<int, arma::Mat<ResultType>>> Cache;
  /** \brief Singular integral cache.
   *
   *  This cache stores the preevaluated local weak forms expressed by
   *  singular integrals. A particular item it stored in r'th row and c'th
   *  column stores, in its second member, the local weak form calculated for
   *  the test element with index it.first and the trial element with index
   *  c. In each column, the items are sorted after increasing test element
   *  index. At the end of each column there can be unused items with test
   *  element index set to INVALID_INDEX (= INT_MAX, so that the sorting is
   *  preserved). */
  Cache m_cache;
  /** \endcond */
};

} // namespace Fiber

#include "default_local_assembler_for_integral_operators_on_surfaces_imp.hpp"

#endif
