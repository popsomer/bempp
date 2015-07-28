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

#ifndef fiber_nonseparable_numerical_test_kernel_trial_integrator_hpp
#define fiber_nonseparable_numerical_test_kernel_trial_integrator_hpp

#include "../common/common.hpp"

#include "test_kernel_trial_integrator.hpp"

#include <tbb/concurrent_unordered_map.h>
#include <tbb/enumerable_thread_specific.h>

//Peter:
#include "geometrical_data.hpp" // For DOMAIN_INDEX


namespace Fiber {

/** \cond FORWARD_DECL */
class OpenClHandler;
template <typename CoordinateType> class CollectionOfShapesetTransformations;
template <typename ValueType> class CollectionOfKernels;
template <typename CoordinateType> class RawGridGeometry;
template <typename BasisFunctionType, typename KernelType, typename ResultType>
class TestKernelTrialIntegral;
/** \endcond */

/** \brief Integration over pairs of elements on non-tensor-product point grids.
 */
template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
class NonseparableNumericalTestKernelTrialIntegrator
    : public TestKernelTrialIntegrator<BasisFunctionType, KernelType,
                                       ResultType> {
public:
  typedef TestKernelTrialIntegrator<BasisFunctionType, KernelType, ResultType>
  Base;
  typedef typename Base::CoordinateType CoordinateType;
  typedef typename Base::ElementIndexPair ElementIndexPair;

  NonseparableNumericalTestKernelTrialIntegrator(
      const arma::Mat<CoordinateType> &localTestQuadPoints,
      const arma::Mat<CoordinateType> &localTrialQuadPoints,
      const std::vector<CoordinateType> quadWeights,
      const GeometryFactory &testGeometryFactory,
      const GeometryFactory &trialGgeometryFactory,
      const RawGridGeometry<CoordinateType> &testRawGeometry,
      const RawGridGeometry<CoordinateType> &trialRawGeometry,
      const CollectionOfShapesetTransformations<CoordinateType> &
          testTransformations,
      const CollectionOfKernels<KernelType> &kernel,
      const CollectionOfShapesetTransformations<CoordinateType> &
          trialTransformations,
      const TestKernelTrialIntegral<BasisFunctionType, KernelType, ResultType> &
          integral,
      const OpenClHandler &openClHandler);
  virtual ~NonseparableNumericalTestKernelTrialIntegrator();

  virtual void
  integrate(CallVariant callVariant, const std::vector<int> &elementIndicesA,
            int elementIndexB, const Shapeset<BasisFunctionType> &basisA,
            const Shapeset<BasisFunctionType> &basisB,
            LocalDofIndex localDofIndexB,
            const std::vector<arma::Mat<ResultType> *> &result) const;


//Peter
//template <typename BasisFunctionType, typename KernelType, typename ResultType,
//          typename GeometryFactory>
//void NonseparableNumericalTestKernelTrialIntegrator<
//    BasisFunctionType, KernelType, ResultType, GeometryFactory>::
    void integratePeter(std::string str, CallVariant callVariant, const std::vector<int> &elementIndicesA,
              int elementIndexB, const Shapeset<BasisFunctionType> &basisA,
              const Shapeset<BasisFunctionType> &basisB,
              LocalDofIndex localDofIndexB,
              const std::vector<arma::Mat<ResultType> *> &result) const {
  const int pointCount = m_quadWeights.size();
  const int elementACount = elementIndicesA.size();

  if (result.size() != elementIndicesA.size())
    throw std::invalid_argument(
        "NonseparableNumericalTestKernelTrialIntegrator::integrate(): "
        "arrays 'result' and 'elementIndicesA' must have the same number "
        "of elements");
  if (pointCount == 0 || elementACount == 0)
    return;
  // TODO: in the (pathological) case that pointCount == 0 but
  // geometryCount != 0, set elements of result to 0.

  // Evaluate constants

	//Peter
	if (elementIndexB == 0) {
//		std::cout << "nonseparable_n_.._imp.hpp integrate(CallVariant)" << std::endl;
	}
	else{
//		std::cout << "nonsepar" << elementIndexB << ", ";
	}
  const int dofCountA = basisA.size();
  const int dofCountB = localDofIndexB == ALL_DOFS ? basisB.size() : 1;
  const int testDofCount = callVariant == TEST_TRIAL ? dofCountA : dofCountB;
  const int trialDofCount = callVariant == TEST_TRIAL ? dofCountB : dofCountA;

  BasisData<BasisFunctionType> testBasisData, trialBasisData;
  GeometricalData<CoordinateType> &testGeomData = m_testGeomData.local();
  GeometricalData<CoordinateType> &trialGeomData = m_trialGeomData.local();

  size_t testBasisDeps = 0, trialBasisDeps = 0;
  size_t testGeomDeps = 0, trialGeomDeps = 0;

  m_testTransformations.addDependencies(testBasisDeps, testGeomDeps);
  m_trialTransformations.addDependencies(trialBasisDeps, trialGeomDeps);
  m_kernels.addGeometricalDependencies(testGeomDeps, trialGeomDeps);
  m_integral.addGeometricalDependencies(testGeomDeps, trialGeomDeps);

  typedef typename GeometryFactory::Geometry Geometry;

  std::unique_ptr<Geometry> geometryA, geometryB;
  const RawGridGeometry<CoordinateType> *rawGeometryA = 0, *rawGeometryB = 0;
  if (callVariant == TEST_TRIAL) {
    geometryA = m_testGeometryFactory.make();
    geometryB = m_trialGeometryFactory.make();
    rawGeometryA = &m_testRawGeometry;
    rawGeometryB = &m_trialRawGeometry;
  } else {
    geometryA = m_trialGeometryFactory.make();
    geometryB = m_testGeometryFactory.make();
    rawGeometryA = &m_trialRawGeometry;
    rawGeometryB = &m_testRawGeometry;
  }

  CollectionOf3dArrays<BasisFunctionType> testValues, trialValues;
  CollectionOf3dArrays<KernelType> kernelValues;

  for (size_t i = 0; i < result.size(); ++i) {
    assert(result[i]);
    result[i]->set_size(testDofCount, trialDofCount);
  }

  rawGeometryB->setupGeometry(elementIndexB, *geometryB);
  if (callVariant == TEST_TRIAL) {
    basisA.evaluate(testBasisDeps, m_localTestQuadPoints, ALL_DOFS,
                    testBasisData);
    basisB.evaluate(trialBasisDeps, m_localTrialQuadPoints, localDofIndexB,
                    trialBasisData);
    geometryB->getData(trialGeomDeps, m_localTrialQuadPoints, trialGeomData);
    if (trialGeomDeps & DOMAIN_INDEX)
      trialGeomData.domainIndex = rawGeometryB->domainIndex(elementIndexB);
    m_trialTransformations.evaluate(trialBasisData, trialGeomData, trialValues);
  } else {
    basisA.evaluate(trialBasisDeps, m_localTrialQuadPoints, ALL_DOFS,
                    trialBasisData);
    basisB.evaluate(testBasisDeps, m_localTestQuadPoints, localDofIndexB,
                    testBasisData);
    geometryB->getData(testGeomDeps, m_localTestQuadPoints, testGeomData);
    if (testGeomDeps & DOMAIN_INDEX)
      testGeomData.domainIndex = rawGeometryB->domainIndex(elementIndexB);
    m_testTransformations.evaluate(testBasisData, testGeomData, testValues);
  }

  // Iterate over the elements
  for (int indexA = 0; indexA < elementACount; ++indexA) {
    const int elementIndexA = elementIndicesA[indexA];
    rawGeometryA->setupGeometry(elementIndexA, *geometryA);
    if (callVariant == TEST_TRIAL) {
      geometryA->getData(testGeomDeps, m_localTestQuadPoints, testGeomData);
      if (testGeomDeps & DOMAIN_INDEX)
        testGeomData.domainIndex = rawGeometryA->domainIndex(elementIndexA);
      m_testTransformations.evaluate(testBasisData, testGeomData, testValues);
    } else {
      geometryA->getData(trialGeomDeps, m_localTrialQuadPoints, trialGeomData);
      if (trialGeomDeps & DOMAIN_INDEX)
        trialGeomData.domainIndex = rawGeometryA->domainIndex(elementIndexA);
      m_trialTransformations.evaluate(trialBasisData, trialGeomData,
                                      trialValues);
    }

//    m_kernels.evaluateAtPointPairs(testGeomData, trialGeomData, kernelValues);
    m_kernels.evaluateAtPointPairsPeter(str,testGeomData, trialGeomData, kernelValues);
    m_integral.evaluateWithNontensorQuadratureRule(
        testGeomData, trialGeomData, testValues, trialValues, kernelValues,
        m_quadWeights, *result[indexA]);
  }
}

  virtual void
  integrate(const std::vector<ElementIndexPair> &elementIndexPairs,
            const Shapeset<BasisFunctionType> &testShapeset,
            const Shapeset<BasisFunctionType> &trialShapeset,
            const std::vector<arma::Mat<ResultType> *> &result) const;

/*
GeometricalData<CoordinateType> getTestGeomData() {
	return &m_testGeomData.local();
}
GeometricalData<CoordinateType> getTrialGeomData() {
	return &m_trialGeomData.local();
}*/

private:
  enum ElementType {
    TEST,
    TRIAL
  };

  const BasisData<BasisFunctionType> &
  basisData(ElementType type,
            const Shapeset<BasisFunctionType> &shapeset) const;

  arma::Mat<CoordinateType> m_localTestQuadPoints;
  arma::Mat<CoordinateType> m_localTrialQuadPoints;
  std::vector<CoordinateType> m_quadWeights;

  const GeometryFactory &m_testGeometryFactory;
  const GeometryFactory &m_trialGeometryFactory;
  const RawGridGeometry<CoordinateType> &m_testRawGeometry;
  const RawGridGeometry<CoordinateType> &m_trialRawGeometry;

  const CollectionOfShapesetTransformations<CoordinateType> &
  m_testTransformations;
  const CollectionOfKernels<KernelType> &m_kernels;
  const CollectionOfShapesetTransformations<CoordinateType> &
  m_trialTransformations;
  const TestKernelTrialIntegral<BasisFunctionType, KernelType, ResultType> &
  m_integral;

  typedef tbb::concurrent_unordered_map<const Shapeset<BasisFunctionType> *,
                                        BasisData<BasisFunctionType> *>
  BasisDataCache;
  mutable BasisDataCache m_cachedTestBasisData;
  mutable BasisDataCache m_cachedTrialBasisData;

  const OpenClHandler &m_openClHandler;
  // thread-local static data for integrate() -- allocation and deallocation of
  // GeometricalData
  // is very time-consuming due to the presence of arma::Cube objects.
  mutable tbb::enumerable_thread_specific<GeometricalData<CoordinateType>>
  m_testGeomData, m_trialGeomData;
};

} // namespace Fiber

#include "nonseparable_numerical_test_kernel_trial_integrator_imp.hpp"

#endif
