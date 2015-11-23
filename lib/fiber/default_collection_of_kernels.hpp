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

#ifndef fiber_default_collection_of_kernels_hpp
#define fiber_default_collection_of_kernels_hpp

#include "collection_of_kernels.hpp"

//Peter:
#include "../common/common.hpp"
#include "geometrical_data.hpp"
#include "scalar_traits.hpp"
#include "../common/complex_aux.hpp"
#include "collection_of_3d_arrays.hpp"
#include "collection_of_4d_arrays.hpp"
#include "modified_helmholtz_3d_single_layer_potential_kernel_functor.hpp"

namespace Fiber {

/** \ingroup weak_form_elements
 *  \brief Default implementation of a collection of kernels.

    This class implements the interface defined by CollectionOfKernels using
    a functor object to evaluate the kernels at specific point pairs.

    \tparam Functor
      Type of the functor that will be passed to the constructor and used to
      evaluate a number of kernels at individual point pairs.

    The Functor class should provide the following interface:

    \code
    class KernelCollectionFunctor
    {
    public:
        typedef ... ValueType; // can be float, double, std::complex<float>
                               // or std::complex<double>
        typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

        // Return the number of kernels in the collection
        int kernelCount() const;
        // Return the number of rows in the tensor being the value of i'th
 kernel
        int kernelRowCount(int i) const;
        // Return the number of columns in the tensor being the value of i'th
 kernel
        int kernelColCount(int i) const;

        // Specify types of geometrical data required by the kernels (see
        // documentation of CollectionOfKernels::addGeometricalDependencies
        // for details)
        void addGeometricalDependencies(size_t& testGeomDeps,
                                        size_t& trialGeomDeps) const;

        // Evaluate the kernels at the pair of points whose geometrical data
        // are provided in the testGeomData and trialGeomData arguments. The
        // (j, k)th element of the tensor being the value of i'th kernel should
        // be written to result[i](j, k).
        template <template <typename T> class CollectionOf2dSlicesOfNdArrays>
        void evaluate(
                const ConstGeometricalDataSlice<CoordinateType>& testGeomData,
                const ConstGeometricalDataSlice<CoordinateType>& trialGeomData,
                CollectionOf2dSlicesOfNdArrays<ValueType>& result) const;

        // (Optional)
        // Return an estimate of the magnitude of the kernel at test and trial
        // points lying in a given distance from each other. This estimate does
        // not need to be accurate; in practice, it is only useful to give it at
        // all for exponentially decaying kernels.  If this function is not
        // defined, the kernel behaves as if its estimated magnitude was 1
        // everywhere.
        CoordinateType estimateRelativeScale(CoordinateType distance) const;
    };
    \endcode

    See the Laplace3dSingleLayerPotentialKernelFunctor class for an example
    implementation of a (simple) kernel collection functor.
 */
template <typename Functor>
class DefaultCollectionOfKernels
    : public CollectionOfKernels<typename Functor::ValueType> {
  typedef CollectionOfKernels<typename Functor::ValueType> Base;

public:
  typedef typename Base::ValueType ValueType;
  typedef typename Base::CoordinateType CoordinateType;

  explicit DefaultCollectionOfKernels(const Functor &functor)
      : m_functor(functor) {}

  const Functor &functor() const { return m_functor; }

  Functor &functor() { return m_functor; }

  virtual void addGeometricalDependencies(size_t &testGeomDeps,
                                          size_t &trialGeomDeps) const;

  virtual void
  evaluateAtPointPairs(const GeometricalData<CoordinateType> &testGeomData,
                       const GeometricalData<CoordinateType> &trialGeomData,
                       CollectionOf3dArrays<ValueType> &result) const;


//Peter:
//template <typename Functor>
//template <template <typename T> class CollectionOf2dSlicesOfNdArrays>
void evaluateAtPointPairsPeter(std::string str, const GeometricalData<CoordinateType> &testGeomData, const GeometricalData<CoordinateType> &trialGeomData, CollectionOf3dArrays<ValueType> &result) const {
  assert(testGeomData.pointCount() == trialGeomData.pointCount());
  const size_t pointCount = testGeomData.pointCount();
  const size_t kernelCount = m_functor.kernelCount();
  result.set_size(kernelCount);
//std::cout << "evalAtPointPairsPeter\n";
//  CoordinateType m_waveNumber = m_functor.waveNumber();
//  ValueType m_waveNumber = m_functor.waveNumber();
//	std::string::size_type sz;     // alias of size_t
//  double earth = std::stod (orbits,&sz);
//  double moon = std::stod (orbits.substr(sz));
//	ValueType waveK = static_cast<CoordinateType>(std::stod(str,&sz) );
//	std::cerr << "In defColKern, str = " << str << std::endl;
//	ValueType waveK = static_cast<CoordinateType>(str);

  for (size_t k = 0; k < kernelCount; ++k)
    result[k].set_size(m_functor.kernelRowCount(k), m_functor.kernelColCount(k),
                       pointCount);

  for (size_t p = 0; p < pointCount; ++p) {
/*	if(std::is_same<Functor,ModifiedHelmholtz3dSingleLayerPotentialKernelFunctor<ValueType>>::value) {
		std::cout << "Right kernel" << std::endl;
		m_functor.evaluatePeter(str, testGeomData.const_slice(p), trialGeomData.const_slice(p), result.slice(p).self());		
	}
	else {
		std::cerr << "Wrong kernel" << std::endl;
		std::terminate();
	}*/
//		(m_functor*)->evaluatePeter(str, testGeomData.const_slice(p), trialGeomData.const_slice(p), result.slice(p).self());		
//std::unique_ptr<ModifiedHelmholtz3dSingleLayerPotentialKernelFunctor<ValueType> >(m_functor)->evaluatePeter(str, testGeomData.const_slice(p), trialGeomData.const_slice(p), result.slice(p).self());
//	const ModifiedHelmholtz3dSingleLayerPotentialKernelFunctor<ValueType>& tmp = dynamic_cast<const ModifiedHelmholtz3dSingleLayerPotentialKernelFunctor<ValueType>& >(m_functor);
//	ModifiedHelmholtz3dSingleLayerPotentialKernelFunctor<ValueType> tmp = static_cast<ModifiedHelmholtz3dSingleLayerPotentialKernelFunctor<ValueType> >(m_functor);
//		static_cast<ModifiedHelmholtz3dSingleLayerPotentialKernelFunctor<ValueType> >(m_functor).evaluatePeter(str, testGeomData.const_slice(p), trialGeomData.const_slice(p), result.slice(p).self());
		m_functor.evaluatePeter(str, testGeomData.const_slice(p), trialGeomData.const_slice(p), result.slice(p).self());

//	const ModifiedHelmholtz3dSingleLayerPotentialKernelFunctor<ValueType> bla = dynamic_cast<const ModifiedHelmholtz3dSingleLayerPotentialKernelFunctor<ValueType>& > (m_functor);
//	bla.evaluatePeter(str, testGeomData.const_slice(p), trialGeomData.const_slice(p), result.slice(p).self());

////    m_functor.evaluate(testGeomData.const_slice(p),
////                      trialGeomData.const_slice(p), result.slice(p).self());
/*
	const ConstGeometricalDataSlice<CoordinateType> &testGeomDataSl = testGeomData.const_slice(p);
	const ConstGeometricalDataSlice<CoordinateType> &trialGeomDataSl = trialGeomData.const_slice(p);
//	CollectionOf2dSlicesOfNdArrays<ValueType> &resultSl = result.slice(p).self();
	CollectionOf2dSlicesOf3dArrays<ValueType> &resultSl = result.slice(p).self();
	const int coordCount = 3;
	CoordinateType sum = 0;
	for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
		CoordinateType diff = testGeomDataSl.global(coordIndex) - trialGeomDataSl.global(coordIndex);
		sum += diff * diff;
	}
	CoordinateType distance = sqrt(sum);
	if (distance >= 2) {
		std::cerr << "Error, 2 <= distance =" << distance <<  std::endl;
	}
	CoordinateType wind = static_cast<CoordinateType>(1.0);
	CoordinateType cuto = static_cast<CoordinateType>(0.1);
	CoordinateType cutoSP = static_cast<CoordinateType>(0.3);
	if (false) {
		if (distance >= cuto*2) {
			wind = static_cast<CoordinateType>(0.0);
		}
		else if (distance >= cuto) {
			wind = exp(2*exp(-cuto/(distance-2*cuto) )/( (distance-2*cuto)/cuto-1) );
		} 
	}
//	result[0](0, 0) = static_cast<CoordinateType>(1.0 / (4.0 * M_PI)) / distance * exp(-m_waveNumber * distance)*wind;
//    result(0, 0) = static_cast<CoordinateType>(1.0 / (4.0 * M_PI)) / distance * exp(-m_waveNumber * distance)*wind;
	resultSl[0](0, 0) = static_cast<CoordinateType>(1.0 / (4.0 * M_PI)) / distance * exp(-waveK * distance)*wind;
*/
	}
}

void evaluateOnGridPeter(std::string str, const GeometricalData<CoordinateType> &testGeomData, const GeometricalData<CoordinateType> &trialGeomData, CollectionOf4dArrays<ValueType> &result) const {
  const size_t testPointCount = testGeomData.pointCount();
  const size_t trialPointCount = trialGeomData.pointCount();
  const size_t kernelCount = m_functor.kernelCount();
  result.set_size(kernelCount);
  for (size_t k = 0; k < kernelCount; ++k)
    result[k].set_size(m_functor.kernelRowCount(k), m_functor.kernelColCount(k),
                       testPointCount, trialPointCount);

#pragma ivdep
  for (size_t trialIndex = 0; trialIndex < trialPointCount; ++trialIndex)
    for (size_t testIndex = 0; testIndex < testPointCount; ++testIndex)
//      m_functor.evaluate(testGeomData.const_slice(testIndex), trialGeomData.const_slice(trialIndex), result.slice(testIndex, trialIndex).self());
      m_functor.evaluatePeter(str,testGeomData.const_slice(testIndex), trialGeomData.const_slice(trialIndex), result.slice(testIndex, trialIndex).self());
}


  virtual void
  evaluateOnGrid(const GeometricalData<CoordinateType> &testGeomData,
                 const GeometricalData<CoordinateType> &trialGeomData,
                 CollectionOf4dArrays<ValueType> &result) const;

  virtual std::pair<const char *, int> evaluateClCode() const;

  virtual CoordinateType estimateRelativeScale(CoordinateType distance) const;

private:
  Functor m_functor;
};

} // namespace Fiber

#include "default_collection_of_kernels_imp.hpp"

#endif
