// Copyright (C) 2011-2012 by the BEM++ Authors
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

#ifndef fiber_default_test_kernel_trial_integral_hpp
#define fiber_default_test_kernel_trial_integral_hpp

#include "test_kernel_trial_integral.hpp"

namespace Fiber {

/** \ingroup weak_form_elements
 *  \brief Default implementation of the TestKernelTrialIntegral interface.

  This class implements the interface defined by TestKernelTrialIntegral
  using a functor object to evaluate the integrand \f$I(x, y)\f$ of an integral
  of the form
  \f[ \int_\Gamma \int_\Sigma I(x, y)\, d\Gamma(x)\, d\Sigma(y), \f]
  where \f$\Gamma\f$ is a test element and \f$\Sigma\f$ a trial element,
  at individual pairs \f$(x, y\f$\f$) of test and trial points.

  \tparam Functor
    Type of the functor that will be passed to the constructor and used to
    evaluate th integrand at individual point pairs.

  The functor should provide the following interface:

  \code{.cpp}
class IntegrandFunctor
{
public:
    typedef ... BasisFunctionType;
    typedef ... KernelType;
    typedef ... ResultType;
    typedef ... CoordinateType;

    void addGeometricalDependencies(size_t& testGeomDeps, size_t& trialGeomDeps)
const;

    template <template <typename T> class CollectionOf2dSlicesOfConstNdArrays>
    ResultType evaluate(
            const ConstGeometricalDataSlice<CoordinateType>& testGeomData,
            const ConstGeometricalDataSlice<CoordinateType>& trialGeomData,
            const CollectionOf1dSlicesOfConst3dArrays<BasisFunctionType>&
testValues,
            const CollectionOf1dSlicesOfConst3dArrays<BasisFunctionType>&
trialValues,
            const CollectionOf2dSlicesOfConstNdArrays<KernelType>& kernelValues)
const;
};
  \endcode

  The addGeometricalDependencies() method should specify any geometrical data
  on which the integrand depends explicitly (not through kernels or shape
  function transformations). For example, if the integrand depends on the
  vectors normal to the surface at test and trial points, the function should
  have the form

  \code{.cpp}
void addGeometricalDependencies(size_t& testGeomDeps, size_t& trialGeomDeps)
const
{
    testGeomDeps |= NORMALS;
    trialGeomDeps |= NORMALS;
}
  \endcode

  The evaluate() method should compute the integrand -- excluding the quadrature
  weights! -- at a single (test point, trial point) pair. It is supplied with
the
  following parameters:

  \param[in] testGeomData
    Geometric data of a point located on the test element.
  \param[in] trialGeomData
    Geometric data of a point located on the trial element.
  \param[in] testValues
    Values of a collection of transformations of a test shape function at the
    test point. The number <tt>testValues[i](j)</tt> is the <em>j</em> component
value of the <em>i</em>th
    transformation of the TO BE CONTINUED
  \param[in] trialValues
    Values of a collection of transformations of a single trial function at the
trial point.
  \param[in] kernels
    Values of a collection of kernels at the (test point, trial point) pair.
 */
template <typename IntegrandFunctor>
class DefaultTestKernelTrialIntegral
    : public TestKernelTrialIntegral<
          typename IntegrandFunctor::BasisFunctionType,
          typename IntegrandFunctor::KernelType,
          typename IntegrandFunctor::ResultType> {
  typedef TestKernelTrialIntegral<typename IntegrandFunctor::BasisFunctionType,
                                  typename IntegrandFunctor::KernelType,
                                  typename IntegrandFunctor::ResultType> Base;

public:
  typedef typename Base::CoordinateType CoordinateType;
  typedef typename Base::BasisFunctionType BasisFunctionType;
  typedef typename Base::KernelType KernelType;
  typedef typename Base::ResultType ResultType;

  explicit DefaultTestKernelTrialIntegral(const IntegrandFunctor &functor)
      : m_functor(functor) {}

  virtual void addGeometricalDependencies(size_t &testGeomDeps,
                                          size_t &trialGeomDeps) const;

  virtual void evaluateWithTensorQuadratureRule(
      const GeometricalData<CoordinateType> &testGeomData,
      const GeometricalData<CoordinateType> &trialGeomData,
      const CollectionOf3dArrays<BasisFunctionType> &testValues,
      const CollectionOf3dArrays<BasisFunctionType> &trialValues,
      const CollectionOf4dArrays<KernelType> &kernelValues,
      const std::vector<CoordinateType> &testQuadWeights,
      const std::vector<CoordinateType> &trialQuadWeights,
      arma::Mat<ResultType> &result) const;

//Peter:
void evaluateWithTensorQuadratureRulePeter(std::string str,
        const GeometricalData<CoordinateType> &testGeomData,
        const GeometricalData<CoordinateType> &trialGeomData,
        const CollectionOf3dArrays<BasisFunctionType> &testValues,
        const CollectionOf3dArrays<BasisFunctionType> &trialValues,
        const CollectionOf4dArrays<KernelType> &kernelValues,
        const std::vector<CoordinateType> &testQuadWeights,
        const std::vector<CoordinateType> &trialQuadWeights,
        arma::Mat<ResultType> &result) const {
  // Evaluate constants
	std::cout << "evaluateWithTensorQuadratureRulePeter in default_test_kernel_trial_integral.hpp" << std::endl;
  const size_t testDofCount = testValues[0].extent(1);
  const size_t trialDofCount = trialValues[0].extent(1);

  const size_t testPointCount = testQuadWeights.size();
  const size_t trialPointCount = trialQuadWeights.size();

  // Assert that array dimensions are correct

  for (size_t i = 0; i < kernelValues.size(); ++i) {
    assert(kernelValues[i].extent(2) == testPointCount);
    assert(kernelValues[i].extent(3) == trialPointCount);
  }
  for (size_t i = 0; i < testValues.size(); ++i)
    assert(testValues[i].extent(2) == testPointCount);
  for (size_t i = 0; i < trialValues.size(); ++i)
    assert(trialValues[i].extent(2) == trialPointCount);

  assert(result.n_rows == testDofCount);
  assert(result.n_cols == trialDofCount);

  // Integrate

  for (size_t trialDof = 0; trialDof < trialDofCount; ++trialDof)
    for (size_t testDof = 0; testDof < testDofCount; ++testDof) {
      ResultType sum = 0.;
      for (size_t trialPoint = 0; trialPoint < trialPointCount; ++trialPoint) {
        const CoordinateType trialWeight =
            trialGeomData.integrationElements(trialPoint) *
            trialQuadWeights[trialPoint];
        ResultType partialSum = 0.;
        for (size_t testPoint = 0; testPoint < testPointCount; ++testPoint) {
          const CoordinateType testWeight =
              testGeomData.integrationElements(testPoint) *
              testQuadWeights[testPoint];
          partialSum += m_functor.evaluate(
                            testGeomData.const_slice(testPoint),
                            trialGeomData.const_slice(trialPoint),
                            testValues.const_slice(testDof, testPoint),
                            trialValues.const_slice(trialDof, trialPoint),
                            kernelValues.const_slice(testPoint, trialPoint)) *
                        testWeight;
        }
        sum += partialSum * trialWeight;
      }
      result(testDof, trialDof) = sum;
    }
}


  virtual void evaluateWithNontensorQuadratureRule(
      const GeometricalData<CoordinateType> &testGeomData,
      const GeometricalData<CoordinateType> &trialGeomData,
      const CollectionOf3dArrays<BasisFunctionType> &testValues,
      const CollectionOf3dArrays<BasisFunctionType> &trialValues,
      const CollectionOf3dArrays<KernelType> &kernelValues,
      const std::vector<CoordinateType> &quadWeights,
      arma::Mat<ResultType> &result) const;

private:
  IntegrandFunctor m_functor;
};

} // namespace Fiber

#include "default_test_kernel_trial_integral_imp.hpp"

#endif
