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

#ifndef fiber_modified_helmholtz_3d_single_layer_potential_kernel_functor_hpp
#define fiber_modified_helmholtz_3d_single_layer_potential_kernel_functor_hpp

#include "../common/common.hpp"

#include "geometrical_data.hpp"
#include "scalar_traits.hpp"

#include "../common/complex_aux.hpp"

namespace Fiber {

/** \ingroup modified_helmholtz_3d
 *  \ingroup functors
 *  \brief Single-layer-potential kernel functor for the modified Helmholtz
 *  equation in 3D.
 *
 *  \tparam ValueType Type used to represent the values of the kernel. It can
 *  be one of: \c float, \c double, <tt>std::complex<float></tt> and
 *  <tt>std::complex<double></tt>. Note that setting \p ValueType to a real
 *  type implies that the wave number will also be purely real.
 *
 *  \see modified_helmholtz_3d
 */

template <typename ValueType_>
class ModifiedHelmholtz3dSingleLayerPotentialKernelFunctor {
public:
  typedef ValueType_ ValueType;
  typedef typename ScalarTraits<ValueType>::RealType CoordinateType;

  explicit ModifiedHelmholtz3dSingleLayerPotentialKernelFunctor(
      ValueType waveNumber)
      : m_waveNumber(waveNumber) {}

  int kernelCount() const { return 1; }
  int kernelRowCount(int /* kernelIndex */) const { return 1; }
  int kernelColCount(int /* kernelIndex */) const { return 1; }

  void addGeometricalDependencies(size_t &testGeomDeps,
                                  size_t &trialGeomDeps) const {
    testGeomDeps |= GLOBALS;
    trialGeomDeps |= GLOBALS;
  }

  ValueType waveNumber() const { return m_waveNumber; }

  template <template <typename T> class CollectionOf2dSlicesOfNdArrays>
  void evaluate(const ConstGeometricalDataSlice<CoordinateType> &testGeomData,
                const ConstGeometricalDataSlice<CoordinateType> &trialGeomData,
                CollectionOf2dSlicesOfNdArrays<ValueType> &result) const {
    const int coordCount = 3;
//std::cout << testGeomData.global(0) << " = test, trial = " << trialGeomData.global(0) << std::endl;

    CoordinateType sum = 0;
    for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
      CoordinateType diff =
          testGeomData.global(coordIndex) - trialGeomData.global(coordIndex);
      sum += diff * diff;
    }
    CoordinateType distance = sqrt(sum);
//Peter:
	if (distance >= 2) {
//		std::cerr << "Error, 2 <= distance =" << distance <<  std::endl;
	}
	CoordinateType wind = static_cast<CoordinateType>(1.0);
	CoordinateType cuto = static_cast<CoordinateType>(0.1);
	CoordinateType cutoSP = static_cast<CoordinateType>(0.3);
//if (testGeomData.global(0) < -7) { // for test without window
//	if (testGeomData.global(0) > 0.8) { // WRONG
//	if (testGeomData.global(0) < -0.8) { // Hardcoded
//	if (true) {
	if(false) {
		if (distance >= cuto*2) {
			wind = static_cast<CoordinateType>(0.0);
		}
		else if (distance >= cuto) {
			wind = exp(2*exp(-cuto/(distance-2*cuto) )/( (distance-2*cuto)/cuto-1) );
//			std::cout << wind;
		} 
	}
//throw std::logic_error("just for printing backtrace");
//std::terminate();
//	if(testGeomData.global(0) > 0
    result[0](0, 0) = static_cast<CoordinateType>(1.0 / (4.0 * M_PI)) / distance * exp(-m_waveNumber * distance)*wind;

  }


  template <template <typename T> class CollectionOf2dSlicesOfNdArrays>
  void evaluatePeter(std::string str, const ConstGeometricalDataSlice<CoordinateType> &testGeomData, const ConstGeometricalDataSlice<CoordinateType> &trialGeomData, CollectionOf2dSlicesOfNdArrays<ValueType> &result) const {
    const int coordCount = 3;
//std::cout << "Entered evaluatePeter, str=" << str << std::endl;
    CoordinateType sum = 0;
    for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
      CoordinateType diff =
          testGeomData.global(coordIndex) - trialGeomData.global(coordIndex);
      sum += diff * diff;
    }
    CoordinateType distance = sqrt(sum);
	if (distance >= 2) {
		std::cerr << "Error, 2 <= distance =" << distance <<  std::endl;
	}
	CoordinateType wind = static_cast<CoordinateType>(1.0);
	CoordinateType cuto = static_cast<CoordinateType>(0.1);
	CoordinateType cutoSP = static_cast<CoordinateType>(0.3);

	if (true) {
		CoordinateType percDecay = static_cast<CoordinateType>(0.8);
		std::string::size_type sz;   
		CoordinateType b = std::stof(str.substr(2),&sz);		
		CoordinateType a = (1-percDecay)*b;
		CoordinateType dist = std::sqrt( std::pow(testGeomData.global(1) - trialGeomData.global(1), 2.0) + std::pow(testGeomData.global(1) - trialGeomData.global(1), 2.0) );

//std::cout << dist << " = dist, b= " << b << "\n";
		// Distance without x to include stationary points (but also shadow when coll in illuminated...)
//	if ((m_testPos[testIndex].x < 0) || (m_trialPos[trialIndex].x > 0)) {
		if (dist <= a) {
			wind = static_cast<CoordinateType>(1.0);
//			std::cout << "Window 1 for fixedWindows with x=" << m_testPos[testIndex].x << std::endl;
		}
		else if (dist <= b) {
			wind = exp(2*exp(-(b-a)/(dist-a) )/((dist-a)/(b-a) -1));
		}
		else { 
			wind = static_cast<CoordinateType>(0.0); 
		}
	}
	if (false) {
		if (distance >= cuto*2) {
			wind = static_cast<CoordinateType>(0.0);
		}
		else if (distance >= cuto) {
			wind = exp(2*exp(-cuto/(distance-2*cuto) )/( (distance-2*cuto)/cuto-1) );
		} 
	}
	result[0](0, 0) = static_cast<CoordinateType>(1.0 / (4.0 * M_PI)) / distance * exp(-m_waveNumber * distance)*wind;
}


  CoordinateType estimateRelativeScale(CoordinateType distance) const {
    return exp(-realPart(m_waveNumber) * distance);
  }

private:
  ValueType m_waveNumber;
};

} // namespace Fiber

#endif
