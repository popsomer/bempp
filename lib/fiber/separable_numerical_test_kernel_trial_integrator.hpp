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

#ifndef fiber_separable_numerical_test_kernel_trial_integrator_hpp
#define fiber_separable_numerical_test_kernel_trial_integrator_hpp

#include "../common/common.hpp"

#include "bempp/common/config_opencl.hpp"

#include "test_kernel_trial_integrator.hpp"

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

/** \brief Integration over pairs of elements on tensor-product point grids. */
template <typename BasisFunctionType, typename KernelType, typename ResultType,
          typename GeometryFactory>
class SeparableNumericalTestKernelTrialIntegrator
    : public TestKernelTrialIntegrator<BasisFunctionType, KernelType,
                                       ResultType> {
public:
  typedef TestKernelTrialIntegrator<BasisFunctionType, KernelType, ResultType>
  Base;
  typedef typename Base::CoordinateType CoordinateType;
  typedef typename Base::ElementIndexPair ElementIndexPair;

  SeparableNumericalTestKernelTrialIntegrator(
      const arma::Mat<CoordinateType> &localTestQuadPoints,
      const arma::Mat<CoordinateType> &localTrialQuadPoints,
      const std::vector<CoordinateType> &testQuadWeights,
      const std::vector<CoordinateType> &trialQuadWeights,
      const GeometryFactory &testGeometryFactory,
      const GeometryFactory &trialGeometryFactory,
      const RawGridGeometry<CoordinateType> &testRawGeometry,
      const RawGridGeometry<CoordinateType> &trialRawGeometry,
      const CollectionOfShapesetTransformations<CoordinateType> &
          testTransformations,
      const CollectionOfKernels<KernelType> &kernels,
      const CollectionOfShapesetTransformations<CoordinateType> &
          trialTransformations,
      const TestKernelTrialIntegral<BasisFunctionType, KernelType, ResultType> &
          integral,
      const OpenClHandler &openClHandler, bool cacheGeometricalData = true);

  virtual ~SeparableNumericalTestKernelTrialIntegrator();

  virtual void
  integrate(CallVariant callVariant, const std::vector<int> &elementIndicesA,
            int elementIndexB, const Shapeset<BasisFunctionType> &basisA,
            const Shapeset<BasisFunctionType> &basisB,
            LocalDofIndex localDofIndexB,
            const std::vector<arma::Mat<ResultType> *> &result) const;



    void integratePeter(std::string str, CallVariant callVariant,
                 const std::vector<int> &elementIndicesA, int elementIndexB,
                 const Shapeset<BasisFunctionType> &basisA,
                 const Shapeset<BasisFunctionType> &basisB,
                 LocalDofIndex localDofIndexB,
                 const std::vector<arma::Mat<ResultType> *> &result) const {
//	if (! m_openClHandler.UseOpenCl() ) {
	if (m_openClHandler.UseOpenCl() ) {
		std::cerr << "ERRORpeter: separable_num...: integrateCl all commented out (in the CallVariant case)..." << std::endl;
		throw "sfd";
//		throw "separable_num...: integrateCpu all commented out...";
	}

	if (elementIndexB == 0) {
		std::cout << "separable_num..grator.hpp integrate(CallVariant) with integrateCpu iso integrateCl" << std::endl;
	}
	else{
//		std::cout << "separ" << elementIndexB << ", ";
	}

  const int testPointCount = m_localTestQuadPoints.n_cols;
  const int trialPointCount = m_localTrialQuadPoints.n_cols;
  const int elementACount = elementIndicesA.size();

  if (result.size() != elementIndicesA.size())
    throw std::invalid_argument(
        "SeparableNumericalTestKernelTrialIntegrator::integrate(): "
        "arrays 'result' and 'elementIndicesA' must have the same number "
        "of elements");
  if (testPointCount == 0 || trialPointCount == 0 || elementACount == 0)
    return;
  // TODO: in the (pathological) case that pointCount == 0 but
  // geometryCount != 0, set elements of result to 0.

  // Evaluate constants

  const int dofCountA = basisA.size();
  const int dofCountB = localDofIndexB == ALL_DOFS ? basisB.size() : 1;
  const int testDofCount = callVariant == TEST_TRIAL ? dofCountA : dofCountB;
  const int trialDofCount = callVariant == TEST_TRIAL ? dofCountB : dofCountA;

  BasisData<BasisFunctionType> testBasisData, trialBasisData;
  GeometricalData<CoordinateType> *testGeomData = &m_testGeomData.local();
  GeometricalData<CoordinateType> *trialGeomData = &m_trialGeomData.local();
  const GeometricalData<CoordinateType> *constTestGeomData = testGeomData;
  const GeometricalData<CoordinateType> *constTrialGeomData = trialGeomData;

  size_t testBasisDeps = 0, trialBasisDeps = 0;
  size_t testGeomDeps = 0, trialGeomDeps = 0;

  m_testTransformations.addDependencies(testBasisDeps, testGeomDeps);
  m_trialTransformations.addDependencies(trialBasisDeps, trialGeomDeps);
  m_kernels.addGeometricalDependencies(testGeomDeps, trialGeomDeps);
  m_integral.addGeometricalDependencies(testGeomDeps, trialGeomDeps);

  typedef typename GeometryFactory::Geometry Geometry;
  std::unique_ptr<Geometry> geometryA, geometryB;
  const RawGridGeometry<CoordinateType> *rawGeometryA = 0, *rawGeometryB = 0;
  if (!m_cacheGeometricalData) {
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
  }

  CollectionOf3dArrays<BasisFunctionType> testValues, trialValues;
  CollectionOf4dArrays<KernelType> kernelValues;

  for (size_t i = 0; i < result.size(); ++i) {
    assert(result[i]);
	if ((elementIndexB == 0) && (i == 0)) {
		std::cout << "before setsize : resi = " << *result[i] << std::endl;
	}
    result[i]->set_size(testDofCount, trialDofCount);
	if ((elementIndexB == 0) && (i == 0)) {
		std::cout << "after setsize: resi = " << *result[i] << std::endl;
	}
  }

  if (!m_cacheGeometricalData)
    rawGeometryB->setupGeometry(elementIndexB, *geometryB);
  if (callVariant == TEST_TRIAL) {
    basisA.evaluate(testBasisDeps, m_localTestQuadPoints, ALL_DOFS,
                    testBasisData);
    basisB.evaluate(trialBasisDeps, m_localTrialQuadPoints, localDofIndexB,
                    trialBasisData);
    if (m_cacheGeometricalData)
      constTrialGeomData = &m_cachedTrialGeomData[elementIndexB];
    else {
      geometryB->getData(trialGeomDeps, m_localTrialQuadPoints, *trialGeomData);
      if (trialGeomDeps & DOMAIN_INDEX)
        trialGeomData->domainIndex = rawGeometryB->domainIndex(elementIndexB);
    }
    m_trialTransformations.evaluate(trialBasisData, *constTrialGeomData,
                                    trialValues);
  } else {
    basisA.evaluate(trialBasisDeps, m_localTrialQuadPoints, ALL_DOFS,
                    trialBasisData);
    basisB.evaluate(testBasisDeps, m_localTestQuadPoints, localDofIndexB,
                    testBasisData);
    if (m_cacheGeometricalData)
      constTestGeomData = &m_cachedTestGeomData[elementIndexB];
    else {
      geometryB->getData(testGeomDeps, m_localTestQuadPoints, *testGeomData);
      if (testGeomDeps & DOMAIN_INDEX)
        testGeomData->domainIndex = rawGeometryB->domainIndex(elementIndexB);
    }
    m_testTransformations.evaluate(testBasisData, *constTestGeomData,
                                   testValues);
  }

  // Iterate over the elements
  for (int indexA = 0; indexA < elementACount; ++indexA) {
    const int elementIndexA = elementIndicesA[indexA];
    if (!m_cacheGeometricalData)
      rawGeometryA->setupGeometry(elementIndexA, *geometryA);
    if (callVariant == TEST_TRIAL) {
      if (m_cacheGeometricalData)
        constTestGeomData = &m_cachedTestGeomData[elementIndicesA[indexA]];
      else {
        geometryA->getData(testGeomDeps, m_localTestQuadPoints, *testGeomData);
        if (testGeomDeps & DOMAIN_INDEX)
          testGeomData->domainIndex = rawGeometryA->domainIndex(elementIndexA);
      }
      m_testTransformations.evaluate(testBasisData, *constTestGeomData,
                                     testValues);
    } else {
      if (m_cacheGeometricalData)
        constTrialGeomData = &m_cachedTrialGeomData[elementIndicesA[indexA]];
      else {
        geometryA->getData(trialGeomDeps, m_localTrialQuadPoints,
                           *trialGeomData);
        if (trialGeomDeps & DOMAIN_INDEX)
          trialGeomData->domainIndex = rawGeometryA->domainIndex(elementIndexA);
      }
      m_trialTransformations.evaluate(trialBasisData, *constTrialGeomData,
                                      trialValues);
    }

//    m_kernels.evaluateOnGrid(*constTestGeomData, *constTrialGeomData,kernelValues);
    m_kernels.evaluateOnGridPeter(str,*constTestGeomData, *constTrialGeomData,kernelValues);
	if ((elementIndexB == 0) && (indexA == 0)) {
		std::cout << "before evTenQuadRule = " << *result[indexA] << std::endl;
//		std::cout << testValues << " = test, trial= " << trialValues << std::endl;
		std::cout << testValues[0].extent(1) << " = testValExtent, trialValextent = " << trialValues[0].extent(1) << std::endl;
		std::cout << m_testQuadWeights.size() << " = testQuadSize, trialQuadSize = " << m_trialQuadWeights.size() << std::endl;
	}
//    m_integral.evaluateWithTensorQuadratureRule(
//        *constTestGeomData, *constTrialGeomData, testValues, trialValues,
//        kernelValues, m_testQuadWeights, m_trialQuadWeights, *result[indexA]);
    m_integral.evaluateWithTensorQuadratureRulePeter(str,
        *constTestGeomData, *constTrialGeomData, testValues, trialValues,
        kernelValues, m_testQuadWeights, m_trialQuadWeights, *result[indexA]);
	if ((elementIndexB == 0) && (indexA == 0)) {
		std::cout << "after evTenQuadRule = " << *result[indexA] << std::endl;
	}
  }
}
/*
GeometricalData<CoordinateType> getTestGeomData() {
	return &m_testGeomData.local();
}
GeometricalData<CoordinateType> getTrialGeomData() {
	return &m_trialGeomData.local();
}*/

  virtual void
  integrate(const std::vector<ElementIndexPair> &elementIndexPairs,
            const Shapeset<BasisFunctionType> &testShapeset,
            const Shapeset<BasisFunctionType> &trialShapeset,
            const std::vector<arma::Mat<ResultType> *> &result) const;

private:
  void integrateCpu(CallVariant callVariant,
                    const std::vector<int> &elementIndicesA, int elementIndexB,
                    const Shapeset<BasisFunctionType> &basisA,
                    const Shapeset<BasisFunctionType> &basisB,
                    LocalDofIndex localDofIndexB,
                    const std::vector<arma::Mat<ResultType> *> &result) const;

  void integrateCl(CallVariant callVariant,
                   const std::vector<int> &elementIndicesA, int elementIndexB,
                   const Shapeset<BasisFunctionType> &basisA,
                   const Shapeset<BasisFunctionType> &basisB,
                   LocalDofIndex localDofIndexB,
                   const std::vector<arma::Mat<ResultType> *> &result) const;

  void integrateCpu(const std::vector<ElementIndexPair> &elementIndexPairs,
                    const Shapeset<BasisFunctionType> &testShapeset,
                    const Shapeset<BasisFunctionType> &trialShapeset,
                    const std::vector<arma::Mat<ResultType> *> &result) const;

  void integrateCl(const std::vector<ElementIndexPair> &elementIndexPairs,
                   const Shapeset<BasisFunctionType> &testShapeset,
                   const Shapeset<BasisFunctionType> &trialShapeset,
                   const std::vector<arma::Mat<ResultType> *> &result) const;

  void precalculateGeometricalData();
  void precalculateGeometricalDataOnSingleGrid(
      const arma::Mat<CoordinateType> &localQuadPoints,
      const GeometryFactory &geometryFactory,
      const RawGridGeometry<CoordinateType> &rawGeometry, size_t geomDeps,
      std::vector<GeometricalData<CoordinateType>> &geomData);

  /**
   * \brief Returns an OpenCL code snippet containing the clIntegrate
   *   kernel function for integrating a single row or column
   */
  const std::pair<const char *, int> clStrIntegrateRowOrCol() const;

  arma::Mat<CoordinateType> m_localTestQuadPoints;
  arma::Mat<CoordinateType> m_localTrialQuadPoints;
  std::vector<CoordinateType> m_testQuadWeights;
  std::vector<CoordinateType> m_trialQuadWeights;

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

  const OpenClHandler &m_openClHandler;
  bool m_cacheGeometricalData;

  std::vector<GeometricalData<CoordinateType>> m_cachedTestGeomData;
  std::vector<GeometricalData<CoordinateType>> m_cachedTrialGeomData;
  mutable tbb::enumerable_thread_specific<GeometricalData<CoordinateType>>
  m_testGeomData, m_trialGeomData;

#ifdef WITH_OPENCL
  cl::Buffer *clTestQuadPoints;
  cl::Buffer *clTrialQuadPoints;
  cl::Buffer *clTestQuadWeights;
  cl::Buffer *clTrialQuadWeights;
#endif
};

} // namespace Fiber

#include "separable_numerical_test_kernel_trial_integrator_imp.hpp"

#endif
