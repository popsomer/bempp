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

    CoordinateType sum = 0;
    for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
      CoordinateType diff =
          testGeomData.global(coordIndex) - trialGeomData.global(coordIndex);
      sum += diff * diff;
    }
    CoordinateType distance = sqrt(sum);
    CoordinateType wind = static_cast<CoordinateType>(1.0);
    result[0](0, 0) = static_cast<CoordinateType>(1.0 / (4.0 * M_PI)) / distance * exp(-m_waveNumber * distance)*wind;

  }


  template <template <typename T> class CollectionOf2dSlicesOfNdArrays>
  void evaluatePeter(std::string str, const ConstGeometricalDataSlice<CoordinateType> &testGeomData, const ConstGeometricalDataSlice<CoordinateType> &trialGeomData, CollectionOf2dSlicesOfNdArrays<ValueType> &result) const {
    const int coordCount = 3;
    CoordinateType sum = 0;
    for (int coordIndex = 0; coordIndex < coordCount; ++coordIndex) {
      CoordinateType diff =
          testGeomData.global(coordIndex) - trialGeomData.global(coordIndex);
      sum += diff * diff;
    }
    CoordinateType percDecay = 0.65;
    CoordinateType distance = sqrt(sum);
	if (distance >= 2) {
		std::cerr << "Error, 2 <= distance =" << distance <<  std::endl;
	}
	CoordinateType wind = static_cast<CoordinateType>(1.0);

	if ((str.at(0) == 'k') | (str.at(0) == 'b') ) { // 'f' means fixed windows multiplying matrix entries, 'k' fixed windows multiplying kernel
		CoordinateType percDecay = static_cast<CoordinateType>(0.8);
		std::string::size_type sz;   
		CoordinateType b = std::stof(str.substr(2),&sz);		
		CoordinateType a = (1-percDecay)*b;
		CoordinateType dist = std::sqrt( std::pow(testGeomData.global(1) - trialGeomData.global(1), 2.0) + std::pow(testGeomData.global(2) - trialGeomData.global(2), 2.0) );

		if (dist <= a) {
			wind = static_cast<CoordinateType>(1.0);
		} else if (dist <= b) {
			wind = exp(2*exp(-(b-a)/(dist-a) )/((dist-a)/(b-a) -1));
		} else { 
			wind = static_cast<CoordinateType>(0.0); 
		}
	} else if (str.at(0) == 'j') { // only compression on the illuminated side
		std::string::size_type sz; // alias of size_t
		CoordinateType b = std::stof(str.substr(4),&sz);		
		CoordinateType a = (1-percDecay)*b;
		CoordinateType dist = std::sqrt( std::pow(testGeomData.global(0) - trialGeomData.global(0), 2.0) + std::pow(testGeomData.global(1) - trialGeomData.global(1), 2.0) + std::pow(testGeomData.global(2) - trialGeomData.global(2), 2.0) );
		if ( (dist <= a) || (testGeomData.global(0) > -0.15) ) {
		   wind = 1; // collocation point x negative is illuminated side, testIndex is row-index in dga.hpp
		} else if (dist <= b) {
		    wind = exp(2*exp(-(b-a)/(dist-a) )/((dist-a)/(b-a) -1));
		} else { 
			wind = static_cast<CoordinateType>(0.0); 
		}
		if(wind < std::abs(testGeomData.global(0) + 0.15) ) {
			wind = std::abs(testGeomData.global(0) + 0.15);
		}
		if ((testGeomData.global(0) < -0.95) && (std::abs(trialGeomData.global(0) - testGeomData.global(0)) < 1e-13) ) {
			std::cerr << "\n str is j: " << str << "\n\n";
		}
		if( (abs(wind) > 1e-13) && (trialGeomData.global(0) > b-0.15) && (testGeomData.global(0) <= -0.15) ) {
			std::cerr << "\n err, window should be zero: " << a << "=a, dist=" << dist << "=dist, tgd0=" << testGeomData.global(0) << "\n";
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
