//gedit /opt/fb/bempp/build/external/include/Trilinos/Thyra_BelosLinearOpWithSolve_def.hpp
//gedit /opt/fb/bempp/build/external/src/Trilinos/packages/stratimikos/adapters/belos/src/Thyra_BelosLinearOpWithSolve_def.hpp

//gedit /opt/fb/bempp/lib/fiber/modified_helmholtz_3d_single_layer_potential_kernel_functor.hpp
//gedit /opt/fb/bempp/lib/fiber/default_collection_of_kernels.hpp
//gedit /opt/fb/bempp/lib/fiber/default_test_kernel_trial_integral.hpp

//gedit /opt/fb/bempp/lib/fiber/separable_numerical_test_kernel_trial_integrator.hpp
// gedit /opt/fb/bempp/lib/linalg/default_iterative_solver.cpp

// gedit /opt/fb/bempp/lib/assembly/dense_global_assembler.hpp
// gedit /opt/fb/bempp/lib/assembly/general_elementary_singular_integral_operator.hpp
// gedit /opt/fb/bempp/lib/fiber/default_local_assembler_for_integral_operators_on_surfaces.hpp
// gedit /opt/fb/bempp/lib/fiber/local_assembler_for_integral_operators.hpp

// cd  /opt/fb/bempp/build/examples/cpp/
// pushd ../..; make tutorial_dirichlet -j6; popd
// ulimit -v 6000000
// ./tutorial_dirichlet >res 2>testgeom
//// pushd /opt/bemppNew/bempp/build; make -j6 2>~/Desktop/Doctoraat/GreenBempp/brol; popd
// pushd ../..; make tutorial_dirichlet -j6 2>/opt/fb/bempp/build/examples/cpp/brol; popd
// ./tutorial_dirichlet 2> ~/Desktop/Doctoraat/GreenBempp/compr

//// compilecommand werkt niet:
//// g++ -I/opt/bemppNew/bempp/build/include -I/opt/bemppNew/bempp/build/include/bempp ~/Desktop/Doctoraat/GreenBempp/helmholtzRep.cpp -L/opt/bemppNew/bempp/build/lib -L/opt/bemppNew/bempp/build/external/lib -I/opt/bemppNew/bempp/build/external/include/ -I/opt/anaconda/inst/include/python2.7/ -I/opt/bemppNew/bempp/build/external/include/Trilinos -w -g -Wall -Werror -std=gnu++11 -lbempp -lteuchoscore -lpthread -Wl,-rpath,/opt/bemppNew/bempp/lib

//// opt/fb/bempp/build/examples/cpp/make
// cd build
// cmake -DCMAKE_BUILD_TYPE=Release -DWITH_FENICS=ON .. -DCMAKE_CXX_FLAGS:STRING=-lpthread


//Simpson:
// cd /export/home1/NoCsBack/nines/fb/bempp/build/examples/cpp/
// vim ../../../examples/cpp/tutorial_dirichlet.cpp 
// pushd ../..; make tutorial_dirichlet -j14; popd
//// cd /export/home1/NoCsBack/nines/fb/bempp/build/examples/cpp
//// vim /export/home1/NoCsBack/nines/fb/bempp/examples/cpp/tutorial_dirichlet.cpp

//  PID USER      PR  NI    VIRT    RES    SHR S  %CPU %MEM     TIME+ COMMAND                                   
// 27994 peter     20   0 8702728 5.035g  27548 R  2339  8.0  68:58.30 tutorial_dirich  


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

RT waveNumber;
//RT waveNumber = 10;
//RT waveNumber = 5;
//RT waveNumber = 1;
//RT waveNumber = 100;
//RT waveNumber = 25;
//RT waveNumber = 40;
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
/*
'conds', zeros(Tl*kl,2), 'nbGm', nbGm, 'avm', avm, 'taus', taus, 'rtests', rtests, 'errBCavm', zeros(Tl*kl,2+nbGm),...
	'errTrueF', zeros(Tl*kl,2+nbGm), 'nnz', zeros(Tl*kl,2), 'perc', zeros(Tl*kl,2), 'errSol', zeros(Tl*kl,2+nbGm), 'errBCcol', ...
	zeros(Tl*kl,2+nbGm), 'compresErr', zeros(Tl*kl,2), 'timeSol', zeros(Tl*kl,2+nbGm), 'iterGm', zeros(Tl*kl,nbGm), 'timeA', ...
	zeros(Tl*kl,2), 'ks', ks, 'field'
*/
//arma::Mat<RT> ks = exp2(arma::linspace(3,3,1));
arma::Mat<RT> ks = arma::exp2(arma::linspace<arma::Mat<RT>>(3,4,2));
const int kl = ks.size();

//arma::Mat<BFT> Ts = arma::linspace<arma::Mat<BFT>>(0.6,0.8,2);
arma::Mat<BFT> Ts = arma::linspace<arma::Mat<BFT>>(0.03,2.0,10);

const int Tl = Ts.size();
//std::cout << ks(1) << "=ks1, length=" << kl << Ts << std::endl;

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

arma::Mat<BFT> conds = arma::zeros(kl,1+Tl);
arma::Mat<BFT> percs = arma::zeros(kl,Tl);

arma::Mat<BFT> errBCavm = arma::zeros(kl,1+Tl);
arma::Mat<BFT> errSol = arma::zeros(kl,Tl);

arma::Mat<BFT> times = arma::zeros(kl,1+Tl);

for(int ki = 0; ki < ks.size(); ki++) {

	tbb::tick_count start = tbb::tick_count::now();

	waveNumber = ks(ki);
//	std::string mfs = "/home/peter/Desktop/Doctoraat/Bol/sphere" + std::to_string(ki+1) + ".msh";
//	std::string mfs = "../../../meshes/sphere" + std::to_string(ki) + ".msh";
	std::string mfs = "../../../meshes/sphere" + std::to_string(ki+1) + ".msh";
//	std::string mfs = "../../../meshes/sphere1.msh"; // Others take too much memory on standard weakForm() without ACA etc
//std::cout << "meshString = " << mfs << std::endl;
	shared_ptr<Grid> grid = loadTriangularMeshFromFile(mfs.c_str());	
std::cout << "Loaded mesh\n";
//	std::cout << "iosauhdfiosaujhdfs" << std::endl;
//	argv[1] = "/home/peter/Desktop/Doctoraat/Bol/sphere0.msh";
//	shared_ptr<Grid> grid = loadTriangularMeshFromFile(argv[1]);
//	shared_ptr<Grid> grid = loadTriangularMeshFromFile("/home/peter/Desktop/Doctoraat/Bol/sphere0.msh");
//	shared_ptr<Grid> grid = loadTriangularMeshFromFile("/home/peter/Desktop/Doctoraat/Bol/sphere1.msh");
//	shared_ptr<Grid> grid = loadTriangularMeshFromFile("/home/peter/Desktop/Doctoraat/Bol/sphere2.msh");
//      shared_ptr<Grid> grid = loadTriangularMeshFromFile("/home/peter/Desktop/Doctoraat/Bol/sphere3.msh");
//	shared_ptr<Grid> grid = loadTriangularMeshFromFile("/home/peter/Desktop/Doctoraat/Man/manBetter.msh");
//	shared_ptr<Grid> grid = loadTriangularMeshFromFile("/home/peter/Desktop/Doctoraat/Tweebollen/tweebollen1.msh");
//	std::cout << "palosidjfaoslidjnfs" << std::endl;
	PiecewiseLinearContinuousScalarSpace<BFT> HplusHalfSpace(grid);
	PiecewiseConstantScalarSpace<BFT> HminusHalfSpace(grid);
	AssemblyOptions assemblyOptions;
//	assemblyOptions.enableSingularIntegralCaching(false);
	std::cout << "iuyshbakjsndf" << std::endl;
	assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW); // Less junk
//	assemblyOptions.setVerbosityLevel(VerbosityLevel::HIGH); // More info (progress % matrix)
	// No ACA at first
//	AcaOptions acaOptions; // Default parameters for ACA
//	acaOptions.eps = 1e-5;
//	assemblyOptions.switchToAcaMode(acaOptions); // But now do ACA to compute corr using H-m


	AccuracyOptions accuracyOptions;
	accuracyOptions.doubleRegular.setRelativeQuadratureOrder(1);
	NumericalQuadratureStrategy<BFT, RT> quadStrategy(accuracyOptions);
	Context<BFT, RT> context(make_shared_from_ref(quadStrategy), assemblyOptions);

//	std::cout << "asOptions.cache = " << assemblyOptions.isSingularIntegralCachingEnabled() << std::endl;
	BoundaryOperator<BFT, RT> slpOp = helmholtz3dSingleLayerBoundaryOperator<BFT>(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HplusHalfSpace), make_shared_from_ref(HminusHalfSpace),waveNumber);

//	DiscreteBoundaryOperator<RT> weak = slpOp.weakForm();
//	arma::Mat<RT> wm = weak->asMatrix();
/////	arma::Mat<RT> wm = slpOp.weakForm()->asMatrix();
//	std::cerr << wm << std::endl;
/*
	std::fstream myStream;
	std::cout << " Writing out matrix in tut_dir" << std::endl;
	myStream.open("/home/peter/Desktop/Doctoraat/GreenBempp/simpsonRes/A",std::ios::out);
	myStream << wm;*/

std::cout << "Making weakForm() \n";
boost::shared_ptr<const Bempp::DiscreteBoundaryOperator<RT> > weakTestMem = slpOp.weakForm();

	std::cout << "Assemble rhs" << std::endl;
	GridFunction<BFT, RT> rhs(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HminusHalfSpace), // is this the right choice?
            surfaceNormalIndependentFunction(MyFunctor()));



	// Initialize the solver
#ifdef WITH_TRILINOS
	std::cout << "Initialize solver TRILINOS" << std::endl;
	DefaultIterativeSolver<BFT, RT> solver(slpOp,ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);
	std::cout << " solver TRILINOS" << std::endl;
	solver.initializeSolver(defaultGmresParameterList(1e-5));

	std::cout << "TRILINOS" << std::endl;
	Solution<BFT, RT> solution = solver.solve(rhs);
#else
	std::cout << "Initialize solver without trilinos" << std::endl;
	DefaultDirectSolver<BFT, RT> solver(slp, rhs);
	solver.solve();
#endif

	std::cout << "TRILINOS done" << std::endl;
	tbb::tick_count end = tbb::tick_count::now();
        times(ki,0) = (end - start).seconds();

std::cout << "ended std calc\n";
//	conds(ki,0) = arma::cond(wm);

/////	const GridFunction<BFT, RT>& solFun = solution.gridFunction();
	const GridFunction<BFT, RT>& solFunOr = solution.gridFunction();

	// Uncomment the block below if you are solving the problem on a sphere and
	// you want to compare the numerical and analytical solution.
/////	arma::Col<RT> solutionCoefficients = solFun.coefficients();
	arma::Col<RT> solutionCoefficientsOr = solFunOr.coefficients();
//	std::cout << solutionCoefficients << std::endl;//Look whether the same as armaSolution from def_iter_solver.cpp: indeed

//	std::ifstream iStream("rhsV", std::ios::binary);
        /*for( size_t i = 0; i < appointments.size(); i++ ) {            
            appointments[i].read(iStream);
        }*/


//	std::ifstream input("/home/peter/Desktop/Doctoraat/GreenBempp/simpsonRes/rhsV");
//	arma::Col<RT> rhsVe{std::istream_iterator<RT>(input), std::istream_iterator<RT>() };
////	std::vector<RT> rhsVe{std::istream_iterator<RT>(input), std::istream_iterator<RT>() };
//        input.close();

//std::cerr << rhsVe[0]<< "= rhsV[0], rhsV[200] = " << rhsVe[200] << std::endl;


std::cout << "Got solution vector \n";

	Helmholtz3dSingleLayerPotentialOperator<BFT> slPot (waveNumber);
	EvaluationOptions evalOptions = EvaluationOptions();

	arma::Mat<RT> potResOr = slPot.evaluateAtPoints(solFunOr, points, quadStrategy, evalOptions);
//	arma::Mat<RT> potResCompr = slPot.evaluateAtPoints(solFunCompr, points, quadStrategy, evalOptions);


std::cout << "Rearranging potResOr\n";
	arma::Mat<RT> diri = potResOr;
	MyFunctor tmp = MyFunctor();
//	for (int i = 0; i < pointCount; ++i) {
	for (int i = 0; i < avm*avm; ++i) {
		arma::Col<CT> pt(3);
		pt(0) = points(0,i);
		pt(1) = points(1,i);
		pt(2) = points(2,i);
		arma::Col<RT> t(1);
		tmp.evaluate(pt,t);
		diri(i) = t(0);
	}
//	std::cout << diri << "=diri, potRes: " << potRes << std::endl;
	arma::Mat<RT> errBCOr = potResOr - diri;
	errBCavm(ki,0) = mean(mean(abs(errBCOr) ));



// -------------- Compression -----------------------------------------
   for(int ti = 0; ti < Ts.size(); ti ++) {
//   for(int ti = 0; ti < 0; ti ++) {

	std::cerr << ki << " = ki, ti = " << ti << ", perc = " << (0.0+ki*Ts.size() + ti)/(0.0+Ts.size()*ks.size()) <<std::endl;

//	std::cout << "aowisuehfoasdhf" << std::endl;
//	boost::shared_ptr<const Bempp::DiscreteBoundaryOperator<RT> > weak = bla.weakFormPeter(" Passed From tutorial_dirichlet.cpp ",context);
//	std::stringstream sstream;
//	sstream << "f, " << waveNumber << " =k " << std::endl;
//	sstream << "f 0.8" << std::endl;
//	std::string str = sstream.str();
	std::string str = "f " + std::to_string(Ts(ti));
	start = tbb::tick_count::now();

	BoundaryOperator<BFT, RT> slpOpCompr = helmholtz3dSingleLayerBoundaryOperator<BFT>(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HplusHalfSpace), make_shared_from_ref(HminusHalfSpace),waveNumber);
	boost::shared_ptr<const Bempp::AbstractBoundaryOperator<double, std::complex<double> > > asdf = slpOpCompr.abstractOperator();
//	std::cout << "iohussdifgoahs'sadf" << std::endl;
	const GeneralElementarySingularIntegralOperator<BFT,RT,RT> bla = dynamic_cast<const GeneralElementarySingularIntegralOperator<BFT,RT,RT>& > (*asdf);

//arma::Col<RT> * rsa = &solutionCoefficients;
//std::unique_ptr<arma::Col<RT>, void (*)(arma::Col<RT>*)> rsa(solutionCoefficients);
//std::unique_ptr<arma::Col<RT>, void (*)(arma::Col<RT>*)> rsa(solutionCoefficients, ::arma::Col<RT>::free);
//	std::unique_ptr<arma::Col<RT> > point {solutionCoefficients};
//	boost::shared_ptr<const Bempp::DiscreteBoundaryOperator<RT> > weakCompr = bla.weakFormPeter(str,context,std::auto_ptr<arma::Col<RT> > (solutionCoefficients), std::unique_ptr<std::vector<RT> > (rhsVe), std::unique_ptr<arma::Mat<RT> > (wm) );
//	boost::shared_ptr<const Bempp::DiscreteBoundaryOperator<RT> > weakCompr = bla.weakFormPeter(str,context,std::unique_ptr<arma::Col<RT> > (solutionCoefficients), std::unique_ptr<std::vector<RT> > (rhsVe), std::unique_ptr<arma::Mat<RT> > (wm) );
//*solutionCoefficients, *rhsVe, std::unique_ptr<> (wm) );

//	boost::shared_ptr<const Bempp::DiscreteBoundaryOperator<RT> > weakCompr = bla.weakFormPeter(str,context,&solutionCoefficientsOr, &rhsVe, &wm);
	arma::Mat<RT> wmDummy(3,3);
	wmDummy.fill(0.);
	std::vector<RT> rhsVeDummy;
	boost::shared_ptr<const Bempp::DiscreteBoundaryOperator<RT> > weakCompr = bla.weakFormPeter(str,context,&solutionCoefficientsOr, &rhsVeDummy, &wmDummy);
//	arma::Mat<RT> wmC = slpOpCompr.weakFormPeter()->asMatrix();
// Computing weak form should be done in (defit)solver, only here for printing matrix here-> wrong, need info rhsVe etc
//	arma::Mat<RT> wmC = weakCompr->asMatrix();
	
//	std::cerr << weak->asMatrix() << std::endl; // For comparisons and validation
//	std::cerr << wmC << std::endl; // For comparisons and validation

//	std::cerr << (weak->asMatrix())[0,0] << std::endl; // Should be(0.0022927,0.000174001)
//	std::cout << wm[0,0] << wm[0,1] << wm[2,0] << wm[50,66] << std::endl; // Should be (0.0022927,0.000174001)(0.00110578,0.000196013)(0.0022927,0.000174001)(0.000275424,0.000255719)
///////	std::cout << wmC[0,0] << wmC[0,1] << wmC[2,0] << wmC[50,66] << std::endl;


//	return 1;
//	weak = slpOp.weakForm();
//	wm = weak->asMatrix();
//	std::cerr << "Real result = " << wm[0,0] << wm[0,1] << wm[2,0] << wm[50,66] << std::endl;
//	return 1;


	DefaultIterativeSolver<BFT, RT> solverCompr(weakCompr, str, slpOpCompr, ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);
	solverCompr.initializeSolver(defaultGmresParameterList(1e-5));
	Solution<BFT, RT> solutionCompr = solverCompr.solve(rhs);

	const GridFunction<BFT, RT>& solFunCompr = solutionCompr.gridFunction();
	arma::Col<RT> solutionCoefficientsCompr = solFunCompr.coefficients();


	end = tbb::tick_count::now();
        times(ki,1+ti) = (end - start).seconds();



//////	arma::Mat<RT> wmC = weakCompr->asMatrix();
////////////////	conds(ki,1+ti) = arma::cond(wmC);
//	percs(ki,ti) = (arma::nonzeros(wmC)).size()/(0.0+wmC.size());
//	arma::Col<RT> nz = arma::nonzeros(wmC);
//	percs(ki,ti) = (nz.size())/(0.0+wmC.size());
//	percs(ki,ti) = arma::accu(wmC != 0)/(0.0+wmC.size());
/////	percs(ki,ti) = arma::accu(wmC != 0)/(0.0+wmC.n_elem);
//	std::cout << wmC.n_elem << "=nelem, acum = " << arma::accu(wmC != 0) << std::endl;

	
//	std::cout << "rel err sol coeffs = " << arma::norm(solutionCoefficientsCompr -solutionCoefficientsOr)/arma::norm(solutionCoefficientsOr) << std::endl;
	errSol(ki,ti) = arma::norm(solutionCoefficientsCompr -solutionCoefficientsOr)/arma::norm(solutionCoefficientsOr);


// ----------------   Validation  ------------------------------------
//	std::cout << "solCoef(1) = " << solutionCoefficients(1) << std::endl;
//	arma::Col<RT> deviation = solutionCoefficients - static_cast<RT>(-1.);
	// % in Armadillo -> elementwise multiplication
//	RT stdDev = sqrt(arma::accu(deviation % deviation)/static_cast<RT>(solutionCoefficients.n_rows));
//	std::cout << "Standard deviation = " << stdDev << std::endl;

//	Helmholtz3dSingleLayerPotentialOperator<BFT> slPot (waveNumber);
//	EvaluationOptions evalOptions = EvaluationOptions();
//	int nth = 5;
//	int nph = 4;
//	int pointCount = nth*nph;
//	int pointCount = 6;
//	arma::Mat<CT> points(3,pointCount);
/*
	points(0,0) = 1.0;	points(1,0) = 0.0;	points(2,0) = 0.0;
	points(0,1) = -1.0;	points(1,1) = 0.0;	points(2,1) = 0.0;

	points(0,2) = 0.0;	points(1,2) = 1.0;	points(2,2) = 0.0;
	points(0,3) = 0.0;	points(1,3) = -1.0;	points(2,3) = 0.0;

	points(0,4) = 0.0;	points(1,4) = 0.0;	points(2,4) = 1.0;
	points(0,5) = 0.0;	points(1,5) = 0.0;	points(2,5) = -1.0;
*/

/*
	for (int thi =0; thi < nth; ++thi) {
		CT theta = M_PI*2*thi/nth;
		for(int pih =0; pih < nph; ++pih) {
			CT phi = M_PI*2*pih/nph;
			int idx = thi*nph+pih;
			points(0,idx) = cos(phi)*sin(theta);
			points(1,idx) = sin(phi)*sin(theta);
			points(2,idx) = cos(theta);
		}
	}
*/	
//	std::cout << "pts: " << points << std::endl;
//	std::cout << waveNumber << " = waveNumber, pts: " << points.t() << std::endl;//Transpose
//	arma::Mat<RT> potRes = slPot.evaluateAtPoints(solFun, points, quadStrategy, evalOptions);
//	arma::Mat<RT> potResOr = slPot.evaluateAtPoints(solFunOr, points, quadStrategy, evalOptions);
	arma::Mat<RT> potResCompr = slPot.evaluateAtPoints(solFunCompr, points, quadStrategy, evalOptions);

//	arma::Mat<RT> diri = potResOr;
//	MyFunctor tmp = MyFunctor();
//	for (int i = 0; i < pointCount; ++i) {
//	for (int i = 0; i < avm*avm; ++i) {
//		arma::Col<CT> pt(3);
//		pt(0) = points(0,i);
//		pt(1) = points(1,i);
//		pt(2) = points(2,i);
//		arma::Col<RT> t(1);
//		tmp.evaluate(pt,t);
//		diri(i) = t(0);
//	}
//	std::cout << diri << "=diri, potRes: " << potRes << std::endl;
//	arma::Mat<RT> errBCOr = potResOr - diri;
	arma::Mat<RT> errBCCompr = potResCompr - diri;
//	std::cout << "errBC= " << errBC << std::endl;
//errBC=     (-1.522e-02,+4.633e-02)    (-5.599e-02,-1.191e-01)    (+6.368e-02,-4.378e-02)    (+6.293e-02,-4.343e-02)    (+6.436e-02,-4.436e-02)    (+6.602e-02,-4.498e-02) if no compr
//errBC=     (+5.403e-01,+8.415e-01)    (-5.702e-02,-1.183e-01)    (+6.487e-02,-6.107e-02)    (+6.355e-02,-6.238e-02)    (+6.567e-02,-6.163e-02)    (+6.654e-02,-6.501e-02) if -0.8 and cuto 0.1
//	std::cout << mean(mean(abs(errBC) )) << " = mean abs err BC, " << std::endl; //<< mean(abs(errBC) ) << std::endl;
//	std::cout << mean(mean(abs(errBCOr) )) << " = Original mean abs err BC, Compr err BC = " << mean(mean(abs(errBCCompr) )) << std::endl;
//	std::cout << mean(mean(abs(errBCOr) )) << " = Original mean abs err BC"  << std::endl;
	errBCavm(ki,1+ti) = mean(mean(abs(errBCCompr) ));
   }
}
std::ofstream myfile;
myfile.open ("res");
myfile << "Output of tutorial_dirichlet.\n";
myfile << real(ks) << " = ks, Ts = " << std::endl << Ts << std::endl;
myfile << percs << " = perc, conds = " << std::endl << conds << std::endl;
myfile << times << " = times" << std::endl << std::endl;
myfile << errBCavm << " = errBCavm, errSol = " << std::endl << errSol << std::endl;
myfile.close();

std::cout << real(ks) << " = ks, Ts = " << Ts << std::endl;
std::cout << percs << " = perc, conds = " << conds << std::endl;
std::cout << times << " = times, " << errBCavm << " = errBCavm, errSol = " << errSol << std::endl;
}

