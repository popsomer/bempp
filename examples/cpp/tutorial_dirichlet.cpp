//gedit /opt/fb/bempp/lib/fiber/modified_helmholtz_3d_single_layer_potential_kernel_functor.hpp
// gedit /opt/fb/bempp/lib/linalg/default_iterative_solver.cpp

// gedit /opt/fb/bempp/lib/assembly/dense_global_assembler.hpp

// cd  /opt/fb/bempp/build/examples/cpp/
// pushd ../..; make tutorial_dirichlet -j6; popd
// ulimit -v 6000000
// ./tutorial_dirichlet >res 2>testgeom

//Simpson:
// cd build
// cmake -DCMAKE_BUILD_TYPE=Release -DWITH_FENICS=ON .. -DCMAKE_CXX_FLAGS:STRING=-lpthread
// cd /export/home1/NoCsBack/nines/fb/bempp/build/examples/cpp/
// vim /export/home1/NoCsBack/nines/fb/bempp/examples/cpp/tutorial_dirichlet.cpp 
// pushd ../..; make tutorial_dirichlet -j14; popd
// ulimit -v 62000000

#include <complex>

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
#include "assembly/general_elementary_singular_integral_operator.hpp"

using namespace Bempp;


typedef double BFT; // basis function type
typedef std::complex<double> RT; // result type (type used to represent discrete operators)
typedef double CT; // coordinate type

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
void fixedWindows() {
// Add code of fixed Windows here
std::cout << "Entered fixedWindows()\n";
}
int main(int argc, char* argv[])
{
fixedWindows();
arma::Mat<RT> ks = arma::exp2(arma::linspace<arma::Mat<RT>>(3,4,2));
const int kl = ks.size();

const int avm = 100;
arma::Mat<BFT> thetas = arma::zeros(avm,1);
arma::Mat<BFT> phis = arma::zeros(avm,1);
for(int av = 0; av < avm; av ++) {
	thetas(av) = M_PI*std::rand();
	phis(av) = M_PI*2*std::rand();
}
arma::Mat<CT> points(3,avm*avm);
// theta in [0,pi) en phi in [0, 2pi)
for (int thi =0; thi < avm; ++thi) {
	for(int phih =0; phih < avm; ++phih) {
		int idx = thi*avm+phih;
		points(0,idx) = cos(phis(phih))*sin(thetas(thi));
		points(1,idx) = sin(phis(phih))*sin(thetas(thi));
		points(2,idx) = cos(thetas(thi));
	}
}

arma::Mat<BFT> conds = arma::zeros(3,2);
arma::Mat<BFT> percs = arma::zeros(3,1);

arma::Mat<BFT> errBCavm = arma::zeros(3,2);
arma::Mat<BFT> errSol = arma::zeros(3,1);

arma::Mat<BFT> times = arma::zeros(3,2);

//for(int sim = 0; sim < 3; sim++) {
for(int sim = 0; sim < 1; sim++) {

	tbb::tick_count start = tbb::tick_count::now();

	shared_ptr<Grid> grid;
	if(sim == 0) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere1.msh");
		waveNumber = ks(0);
	} else if(sim == 1) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere2.msh");	
		waveNumber = ks(1);
	} else{
		grid = loadTriangularMeshFromFile("../../../meshes/man.msh");
		waveNumber = ks(0);
	}
	PiecewiseLinearContinuousScalarSpace<BFT> HplusHalfSpace(grid);
	PiecewiseConstantScalarSpace<BFT> HminusHalfSpace(grid);
	AssemblyOptions assemblyOptions;
	assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
	// No ACA (AcaOptions) 

	AccuracyOptions accuracyOptions;
	accuracyOptions.doubleRegular.setRelativeQuadratureOrder(1);
	NumericalQuadratureStrategy<BFT, RT> quadStrategy(accuracyOptions);
	Context<BFT, RT> context(make_shared_from_ref(quadStrategy), assemblyOptions);

	BoundaryOperator<BFT, RT> slpOp = helmholtz3dSingleLayerBoundaryOperator<BFT>(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HplusHalfSpace), make_shared_from_ref(HminusHalfSpace),waveNumber);

	arma::Mat<RT> wm = slpOp.weakForm()->asMatrix();

	std::cout << "Assemble rhs" << std::endl;
	GridFunction<BFT, RT> rhs(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HminusHalfSpace), // is this the right choice?
            surfaceNormalIndependentFunction(MyFunctor()));

	// Initialize the solver
#ifdef WITH_TRILINOS
	std::cout << "Initialize solver TRILINOS" << std::endl;
	DefaultIterativeSolver<BFT, RT> solver(slpOp,ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);
	solver.initializeSolver(defaultGmresParameterList(1e-5));

	std::cout << "TRILINOS" << std::endl;
	Solution<BFT, RT> solution = solver.solve(rhs);
#else
	std::cout << "Initialize solver without trilinos" << std::endl;
	DefaultDirectSolver<BFT, RT> solver(slp, rhs);
	solver.solve();
#endif
	tbb::tick_count end = tbb::tick_count::now();
        times(sim,0) = (end - start).seconds();

	std::cout << "ended std calc\n";

// ------------------------ Validation of original -----------------------
	conds(sim,0) = arma::cond(wm);

	const GridFunction<BFT, RT>& solFunOr = solution.gridFunction();
	arma::Col<RT> solutionCoefficientsOr = solFunOr.coefficients(); // Same as armaSolution from def_iter_solver.cpp

	std::ifstream input("/home/peter/Desktop/Doctoraat/GreenBempp/simpsonRes/rhsV");
	std::vector<RT> rhsVe{std::istream_iterator<RT>(input), std::istream_iterator<RT>() };
        input.close();

	Helmholtz3dSingleLayerPotentialOperator<BFT> slPot (waveNumber);
	EvaluationOptions evalOptions = EvaluationOptions();
	arma::Mat<RT> potResOr = slPot.evaluateAtPoints(solFunOr, points, quadStrategy, evalOptions);
	arma::Mat<RT> diri = potResOr;
	MyFunctor tmp = MyFunctor();
	for (int i = 0; i < avm*avm; ++i) {
		arma::Col<CT> pt(3);
		pt(0) = points(0,i);
		pt(1) = points(1,i);
		pt(2) = points(2,i);
		arma::Col<RT> t(1);
		tmp.evaluate(pt,t);
		diri(i) = t(0);
	}
	errBCavm(sim,0) = mean(mean(abs(potResOr - diri) ));


// -------------- Compression -----------------------------------------
	std::string str = "f 0.8"; //"f "+ std::to_string(Ts(ti)); //"c ";
	start = tbb::tick_count::now();

	BoundaryOperator<BFT, RT> slpOpCompr = helmholtz3dSingleLayerBoundaryOperator<BFT>(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HplusHalfSpace), make_shared_from_ref(HminusHalfSpace),waveNumber);
	boost::shared_ptr<const Bempp::AbstractBoundaryOperator<double, std::complex<double> > > asdf = slpOpCompr.abstractOperator();
	const GeneralElementarySingularIntegralOperator<BFT,RT,RT> bla = dynamic_cast<const GeneralElementarySingularIntegralOperator<BFT,RT,RT>& > (*asdf);

//	arma::Mat<RT> wmDummy(3,3);
//	wmDummy.fill(0.);
//	std::vector<RT> rhsVeDummy;
	boost::shared_ptr<const Bempp::DiscreteBoundaryOperator<RT> > weakCompr = bla.weakFormPeter(str,context,&solutionCoefficientsOr, &rhsVe, &wm);
	
	DefaultIterativeSolver<BFT, RT> solverCompr(weakCompr, str, slpOpCompr, ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);
	solverCompr.initializeSolver(defaultGmresParameterList(1e-5));
	Solution<BFT, RT> solutionCompr = solverCompr.solve(rhs);

	const GridFunction<BFT, RT>& solFunCompr = solutionCompr.gridFunction();
	arma::Col<RT> solutionCoefficientsCompr = solFunCompr.coefficients();

	end = tbb::tick_count::now();
        times(sim,1) = (end - start).seconds();

// ----------------   Validation of compressed  ------------------------------------

	arma::Mat<RT> wmC = weakCompr->asMatrix();
	conds(sim,1) = arma::cond(wmC);
	percs(sim,0) = arma::accu(wmC != 0)/(0.0+wmC.n_elem);
	errSol(sim,0) = arma::norm(solutionCoefficientsCompr -solutionCoefficientsOr)/arma::norm(solutionCoefficientsOr);
	arma::Mat<RT> potResCompr = slPot.evaluateAtPoints(solFunCompr, points, quadStrategy, evalOptions);
	errBCavm(sim,1) = mean(mean(abs(potResCompr - diri) ));
}

std::ofstream myfile;
myfile.open ("res");
myfile << "Output of tutorial_dirichlet.\n";
myfile << real(ks) << " = ks " << std::endl;
myfile << percs << " = perc, conds = " << std::endl << conds << std::endl;
myfile << times << " = times" << std::endl << std::endl;
myfile << errBCavm << " = errBCavm, errSol = " << std::endl << errSol << std::endl;
myfile.close();

std::cout << real(ks) << " = ks " << std::endl;
std::cout << percs << " = perc, conds = " << conds << std::endl;
std::cout << times << " = times, " << errBCavm << " = errBCavm, errSol = " << errSol << std::endl;
//correlations();
}
//void correlations() {
// add correlation code here later
//}

