// gedit /opt/fb/bempp/lib/assembly/dense_global_assembler.hpp
// gedit /opt/fb/bempp/lib/assembly/general_elementary_singular_integral_operator.hpp

//// pushd /opt/bemppNew/bempp/build; make -j6 2>~/Desktop/Doctoraat/GreenBempp/brol; popd
// pushd ../..; make tutorial_dirichlet -j6 2>/opt/fb/bempp/build/examples/cpp/brol; popd

// cd  /opt/fb/bempp/build/examples/cpp/
// pushd ../..; make tutorial_dirichlet -j6; popd

//// compilecommand werkt niet:
//// g++ -I/opt/bemppNew/bempp/build/include -I/opt/bemppNew/bempp/build/include/bempp ~/Desktop/Doctoraat/GreenBempp/helmholtzRep.cpp -L/opt/bemppNew/bempp/build/lib -L/opt/bemppNew/bempp/build/external/lib -I/opt/bemppNew/bempp/build/external/include/ -I/opt/anaconda/inst/include/python2.7/ -I/opt/bemppNew/bempp/build/external/include/Trilinos -w -g -Wall -Werror -std=gnu++11 -lbempp -lteuchoscore -lpthread -Wl,-rpath,/opt/bemppNew/bempp/lib

//// opt/fb/bempp/build/examples/cpp/make

#include <complex>

//#include "bempp/common/config_alugrid.hpp"
#include "bempp/common/config_trilinos.hpp"//

#include "meshes.hpp"
#include "meshes.cpp"//

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
//Extra na 2.9.0
#include "assembly/general_elementary_singular_integral_operator.hpp"

using namespace Bempp;


typedef double BFT; // basis function type
typedef std::complex<double> RT; // result type (type used to represent discrete operators)
typedef double CT; // coordinate type

RT waveNumber = 1;

class MyFunctor
{
public:
    // Type of the function's values (e.g. float or std::complex<double>)
    typedef RT ValueType;
    // Type of coordinates (must be the "real part" of ValueType)
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

int main(int argc, char* argv[])
{
	std::cout << "iosauhdfiosaujhdfs" << std::endl;
//	argv[1] = "/home/peter/Desktop/Doctoraat/Bol/sphere0.msh";
//	shared_ptr<Grid> grid = loadTriangularMeshFromFile(argv[1]);
	shared_ptr<Grid> grid = loadTriangularMeshFromFile("/home/peter/Desktop/Doctoraat/Bol/sphere0.msh");
	std::cout << "palosidjfaoslidjnfs" << std::endl;
	PiecewiseLinearContinuousScalarSpace<BFT> HplusHalfSpace(grid);
	PiecewiseConstantScalarSpace<BFT> HminusHalfSpace(grid);
	AssemblyOptions assemblyOptions;

	std::cout << "iuyshbakjsndf" << std::endl;
	assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW); // Less junk
	// No ACA at first
	AccuracyOptions accuracyOptions;
	accuracyOptions.doubleRegular.setRelativeQuadratureOrder(1);
	NumericalQuadratureStrategy<BFT, RT> quadStrategy(accuracyOptions);
	Context<BFT, RT> context(make_shared_from_ref(quadStrategy), assemblyOptions);

	std::cout << "poasijdfoiajs" << std::endl;
	BoundaryOperator<BFT, RT> slpOp = helmholtz3dSingleLayerBoundaryOperator<BFT>(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HplusHalfSpace), make_shared_from_ref(HminusHalfSpace),waveNumber);

	boost::shared_ptr<const Bempp::AbstractBoundaryOperator<double, std::complex<double> > > asdf = slpOp.abstractOperator();
	std::cout << "iohussdifgoahs'sadf" << std::endl;
	const GeneralElementarySingularIntegralOperator<BFT,RT,RT> bla = dynamic_cast<const GeneralElementarySingularIntegralOperator<BFT,RT,RT>& > (*asdf);

	std::cout << "aowisuehfoasdhf" << std::endl;
	boost::shared_ptr<const Bempp::DiscreteBoundaryOperator<RT> > weak = bla.weakFormPeter("soifdj",context);
//	std::cerr << weak->asMatrix() << std::endl;

	std::cout << "Assemble rhs" << std::endl;
	GridFunction<BFT, RT> rhs(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HminusHalfSpace), // is this the right choice?
            surfaceNormalIndependentFunction(MyFunctor()));

	// Initialize the solver
#ifdef WITH_TRILINOS
	std::cout << "Initialize solver TRILINOS" << std::endl;
	DefaultIterativeSolver<BFT, RT> solver(slpOp,ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);
	solver.initializeSolver(defaultGmresParameterList(1e-5));
	Solution<BFT, RT> solution = solver.solve(rhs);
#else
	std::cout << "Initialize solver without trilinos" << std::endl;
	DefaultDirectSolver<BFT, RT> solver(slp, rhs);
	solver.solve();
#endif

	const GridFunction<BFT, RT>& solFun = solution.gridFunction();

	// Uncomment the block below if you are solving the problem on a sphere and
	// you want to compare the numerical and analytical solution.
	arma::Col<RT> solutionCoefficients = solFun.coefficients();
	std::cout << "solCoef(1) = " << solutionCoefficients(1) << std::endl;
	arma::Col<RT> deviation = solutionCoefficients - static_cast<RT>(-1.);
	// % in Armadillo -> elementwise multiplication
	RT stdDev = sqrt(arma::accu(deviation % deviation)/static_cast<RT>(solutionCoefficients.n_rows));
	std::cout << "Standard deviation = " << stdDev << std::endl;

	Helmholtz3dSingleLayerPotentialOperator<BFT> slPot (waveNumber);
	EvaluationOptions evalOptions = EvaluationOptions();
	int pointCount = 6;
	arma::Mat<CT> points(3,pointCount);
	points(0,0) = 1.0;	points(1,0) = 0.0;	points(2,0) = 0.0;
	points(0,1) = -1.0;	points(1,1) = 0.0;	points(2,1) = 0.0;

	points(0,2) = 0.0;	points(1,2) = 1.0;	points(2,2) = 0.0;
	points(0,3) = 0.0;	points(1,3) = -1.0;	points(2,3) = 0.0;

	points(0,4) = 0.0;	points(1,4) = 0.0;	points(2,4) = 1.0;
	points(0,5) = 0.0;	points(1,5) = 0.0;	points(2,5) = -1.0;

	std::cout << "pts: " << points << std::endl;
	arma::Mat<RT> potRes = slPot.evaluateAtPoints(solFun, points, quadStrategy, evalOptions);
	std::cout << "pot: " << potRes << std::endl;

	arma::Mat<RT> diri = potRes;
	MyFunctor tmp = MyFunctor();
	for (int i = 0; i < pointCount; ++i) {
		arma::Col<CT> pt(3);
		pt(0) = points(0,i);
		pt(1) = points(1,i);
		pt(2) = points(2,i);
		arma::Col<RT> t(1);
		tmp.evaluate(pt,t);
		diri(i) = t(0);
	}
	std::cout << diri << "=diri, potRes: " << potRes << std::endl;
	arma::Mat<RT> errBC = potRes-diri;
	std::cout << "errBC= " << errBC << std::endl;
}


