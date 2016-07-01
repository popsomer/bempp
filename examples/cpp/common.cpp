
#include <complex>

#include "bempp/common/config_trilinos.hpp"
#include "meshes.hpp"
#include "meshes.cpp"

#include "assembly/assembly_options.hpp"
#include "assembly/abstract_boundary_operator_sum.hpp"
#include "assembly/boundary_operator.hpp"
#include "assembly/context.hpp"
#include "assembly/discrete_boundary_operator.hpp"
#include "assembly/evaluation_options.hpp"
#include "assembly/grid_function.hpp"
#include "assembly/interpolated_function.hpp"
#include "assembly/numerical_quadrature_strategy.hpp"
#include "assembly/surface_normal_independent_function.hpp"

#include "assembly/identity_operator.hpp"
#include "assembly/helmholtz_3d_single_layer_boundary_operator.hpp"
#include "assembly/modified_helmholtz_3d_single_layer_boundary_operator.hpp"
#include "assembly/helmholtz_3d_double_layer_boundary_operator.hpp"
#include "assembly/helmholtz_3d_single_layer_potential_operator.hpp"
#include "assembly/helmholtz_3d_double_layer_potential_operator.hpp"

#include "common/boost_make_shared_fwd.hpp"
#include "common/scalar_traits.hpp"
#include "grid/grid_factory.hpp"
#include "grid/grid.hpp"

#include "linalg/preconditioner.hpp"
#include "linalg/default_iterative_solver.hpp"
#include "linalg/default_direct_solver.hpp"
#include "space/piecewise_linear_continuous_scalar_space.hpp"
#include "space/piecewise_constant_scalar_space.hpp"

#include "common/armadillo_fwd.hpp"
#include <cmath>
#include <iostream>
#include <memory>
#include "assembly/general_elementary_singular_integral_operator.hpp"


#include <math.h>
#include "../../build/external/include/Trilinos/Thyra_BelosLinearOpWithSolve_def.hpp"
#include "space/piecewise_polynomial_continuous_scalar_space.hpp"
#include "../../build/external/include/boost/math/special_functions/bessel.hpp" // for sph_bessel

#define ARMA_USE_CXX11 // to be able to define arma::Mat asdf = { { ...}, {...} }

using namespace Bempp;

typedef double BFT; // basis function type
typedef std::complex<double> RT; // result type (type used to represent discrete operators)
typedef double CT; // coordinate type

typedef std::complex<long double> Hpc;
typedef long double Hp; // High precisions for calculations

RT waveNumber;
class MyFunctor
{
public:
    typedef RT ValueType;
    typedef ScalarTraits<RT>::RealType CoordinateType;
    // Number of components of the function's argument
    int argumentDimension() const { return 3; }
    int resultDimension() const { return 1; }
    // Evaluate the function at the point "point" and store result in the array "result"
    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
	RT imu = std::complex<double>(0, 1);
        result(0) = -exp(waveNumber*imu*point(0));
    }
};


