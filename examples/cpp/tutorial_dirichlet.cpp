// gedit /opt/fb/bempp/lib/fiber/modified_helmholtz_3d_single_layer_potential_kernel_functor.hpp
// gedit /opt/fb/bempp/lib/linalg/default_iterative_solver.cpp

// gedit /opt/fb/bempp/lib/assembly/dense_global_assembler.hpp
// gedit /opt/fb/bempp/build/external/include/Trilinos/Thyra_BelosLinearOpWithSolve_def.hpp

// cd /opt/fb/bempp/build/examples/cpp/
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



#include <math.h> //Peter
#include "../../build/external/include/Trilinos/Thyra_BelosLinearOpWithSolve_def.hpp"//Peter

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
class solSphere
{
public:
    typedef RT ValueType;
    typedef ScalarTraits<RT>::RealType CoordinateType;
    int argumentDimension() const { return 3; }
    int resultDimension() const { return 1; }
    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
//	using namespace boost::math;
	RT imu = std::complex<double>(0, 1);
//        result(0) = 0.0;
	RT legnm1;
	RT legn;
	RT tmp;
	for(int n = 0; n < 1.4*std::abs(waveNumber)+40; n ++) {
		if(n == 0) {
			legn = 1;
		}
		else if(n == 1) {
			legn = point(0);
			legnm1 = 1;
		}
		else {
			tmp = legn;
//			legn = ((2*n+1)*point(0)*legn-n*legnm1)/(n+1);
			legn = (std::complex<double>(2*n+1,0)*point(0)*legn-std::complex<double>(n, 0)*legnm1)/std::complex<double>(n+1,0);
			legnm1 = tmp;
		}
//		result += (2.0*n+1.0)*imu^n;//*boost::cyl_bessel_j(n,wavek*1)
//		result += (2.0*n+1.0)*std::pow(imu,n)*boost::math::cyl_bessel_j(n,waveNumber*1.0);
//		result -= (2.0*n+1.0)*std::pow(imu,n)*jn(n,std::abs(waveNumber)*1.0)*legendre_pL(n,point(0));
		result -= (2.0*n+1.0)*std::pow(imu,n)*jn(n,std::abs(waveNumber)*1.0)*legn;
	}
    }
};

void fixedWindows() {
// Add code of fixed Windows here
std::cout << "Entered fixedWindows()\n";

arma::Mat<RT> ks = arma::exp2(arma::linspace<arma::Mat<RT>>(3,5,3));
const int kl = ks.size();
//const int kl = 1;

//arma::Mat<BFT> Ts = arma::linspace<arma::Mat<BFT>>(0.6,0.8,2);
arma::Mat<BFT> Ts = arma::linspace<arma::Mat<BFT>>(0.03,2.0,10);

const int Tl = Ts.size();

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

arma::Mat<RT> wm;
BoundaryOperator<BFT, RT> slpOp;
//DefaultIterativeSolver<BFT, RT> solver;
Solution<BFT, RT> solution;
//GridFunction<BFT, RT>& solFun;

for(int ki = 0; ki < kl; ki++) {
//for(int ki = 0; ki < 1; ki++) {

	tbb::tick_count start = tbb::tick_count::now();

	waveNumber = ks(ki);
	std::string mfs = "../../../meshes/sphere" + std::to_string(ki+1) + ".msh";
//	std::string mfs = "../../../meshes/sphere0.msh";
	shared_ptr<Grid> grid = loadTriangularMeshFromFile(mfs.c_str());

	PiecewiseLinearContinuousScalarSpace<BFT> HplusHalfSpace(grid);
	PiecewiseConstantScalarSpace<BFT> HminusHalfSpace(grid);
	AssemblyOptions assemblyOptions;
	assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW); // Less junk
//	assemblyOptions.setVerbosityLevel(VerbosityLevel::HIGH); // More info (progress % matrix)
	// No ACA (AcaOptions) 

	AccuracyOptions accuracyOptions;
	accuracyOptions.doubleRegular.setRelativeQuadratureOrder(1);
	NumericalQuadratureStrategy<BFT, RT> quadStrategy(accuracyOptions);
	Context<BFT, RT> context(make_shared_from_ref(quadStrategy), assemblyOptions);

	slpOp = helmholtz3dSingleLayerBoundaryOperator<BFT>(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HplusHalfSpace), make_shared_from_ref(HminusHalfSpace),waveNumber);

	wm = slpOp.weakForm()->asMatrix();

//	std::cout << "Assemble rhs" << std::endl;
	GridFunction<BFT, RT> rhs(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HminusHalfSpace), // is this the right choice?
            surfaceNormalIndependentFunction(MyFunctor()));

	// Initialize the solver
#ifdef WITH_TRILINOS
//	std::cout << "Initialize solver TRILINOS" << std::endl;
	DefaultIterativeSolver<BFT, RT>* solver = new DefaultIterativeSolver<BFT, RT>(slpOp,ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);
	solver->initializeSolver(defaultGmresParameterList(1e-5));

	std::cout << "TRILINOS" << std::endl;
	solution = solver->solve(rhs);
#else
//	std::cout << "Initialize solver without trilinos" << std::endl;
//	DefaultDirectSolver<BFT, RT> solver(slp, rhs);
//	solver.solve();
#endif
	tbb::tick_count end = tbb::tick_count::now();
        times(ki,0) = (end - start).seconds();

std::cout << "ended std calc\n";
	conds(ki,0) = arma::cond(wm);

//	const GridFunction<BFT, RT>& solFunOr = solution.gridFunction();
	GridFunction<BFT, RT> solFun = solution.gridFunction();
	arma::Col<RT> solutionCoefficientsOr = solFun.coefficients(); // Same as armaSolution from def_iter_solver.cpp

//	std::ifstream input("/home/peter/Desktop/Doctoraat/GreenBempp/simpsonRes/rhsV");
//	arma::Col<RT> rhsVe{std::istream_iterator<RT>(input), std::istream_iterator<RT>() };
////	std::vector<RT> rhsVe{std::istream_iterator<RT>(input), std::istream_iterator<RT>() };
//        input.close();

	Helmholtz3dSingleLayerPotentialOperator<BFT> slPot (waveNumber);
	EvaluationOptions evalOptions = EvaluationOptions();
	arma::Mat<RT> potResOr = slPot.evaluateAtPoints(solFun, points, quadStrategy, evalOptions);
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
	arma::Mat<RT> errBCOr = potResOr - diri;
	errBCavm(ki,0) = mean(mean(abs(errBCOr) ));


// -------------- Compression -----------------------------------------
   for(int ti = 0; ti < Tl; ti ++) {
	time_t now = time(0);
	std::cerr << "-----------------------------------\n" << ki << " = ki, ti = " << ti << ", perc = " << (0.0+ki*Tl + ti)/(0.0+Tl*kl) << " " << asctime(localtime(&now) ); // << "\n ----------------------------- \n";

//	std::cerr << "-----------------------------------\n" << ki << " = ki, ti = " << ti << ", perc = " << (0.0+ki*Ts.size() + ti)/(0.0+Ts.size()*ks.size()) << " " << asctime(localtime(&now) ) << "\n ----------------------------- \n";

	std::string str = "f " + std::to_string(Ts(ti));
	start = tbb::tick_count::now();

	slpOp = helmholtz3dSingleLayerBoundaryOperator<BFT>(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HplusHalfSpace), make_shared_from_ref(HminusHalfSpace),waveNumber);
	boost::shared_ptr<const Bempp::AbstractBoundaryOperator<double, std::complex<double> > > asdf = slpOp.abstractOperator();
	const GeneralElementarySingularIntegralOperator<BFT,RT,RT> bla = dynamic_cast<const GeneralElementarySingularIntegralOperator<BFT,RT,RT>& > (*asdf);

	arma::Mat<RT> wmDummy(3,3);
	wmDummy.fill(0.);
	std::vector<RT> rhsVeDummy;
	boost::shared_ptr<const Bempp::DiscreteBoundaryOperator<RT> > weakCompr = bla.weakFormPeter(str,context,&solutionCoefficientsOr, &rhsVeDummy, &wmDummy);
	
//	solver.~DefaultIterativeSolver<BFT, RT>(weakCompr, str, slpOp, ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE)();
//	new (solver) DefaultIterativeSolver<BFT, RT>(weakCompr, str, slpOp, ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);
	solver = new DefaultIterativeSolver<BFT, RT>(weakCompr, str, slpOp, ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);
	solver->initializeSolver(defaultGmresParameterList(1e-5));
	solution = solver->solve(rhs);

	solFun = solution.gridFunction();
	arma::Col<RT> solutionCoefficientsCompr = solFun.coefficients();

	end = tbb::tick_count::now();
        times(ki,1+ti) = (end - start).seconds();

// ----------------   Validation of compressed  ------------------------------------

	wm = weakCompr->asMatrix();
	conds(ki,1+ti) = arma::cond(wm);
	percs(ki,ti) = arma::accu(wm != 0)/(0.0+wm.n_elem);
	errSol(ki,ti) = arma::norm(solutionCoefficientsCompr -solutionCoefficientsOr)/arma::norm(solutionCoefficientsOr);
	arma::Mat<RT> potResCompr = slPot.evaluateAtPoints(solFun, points, quadStrategy, evalOptions);
	arma::Mat<RT> errBCCompr = potResCompr - diri;
	errBCavm(ki,1+ti) = mean(mean(abs(errBCCompr) ));
   }
std::ofstream myfile;
myfile.open ("res");
myfile << "Output of tutorial_dirichlet.\n";
myfile << real(ks) << " = ks, Ts = " << std::endl << Ts << std::endl;
myfile << percs << " = perc, conds = " << std::endl << conds << std::endl;
myfile << times << " = times" << std::endl << std::endl;
myfile << errBCavm << " = errBCavm, errSol = " << std::endl << errSol << std::endl;
myfile.close();
}

std::cout << real(ks) << " = ks, Ts = " << Ts << std::endl;
std::cout << percs << " = perc, conds = " << conds << std::endl;
std::cout << times << " = times, " << errBCavm << " = errBCavm, errSol = " << errSol << std::endl;

}


int main(int argc, char* argv[])
{
fixedWindows();
/*
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
//	assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
	assemblyOptions.setVerbosityLevel(VerbosityLevel::HIGH);
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
//	conds(sim,0) = arma::cond(wm);

	const GridFunction<BFT, RT>& solFunOr = solution.gridFunction();
	arma::Col<RT> solutionCoefficientsOr = solFunOr.coefficients(); // Same as armaSolution from def_iter_solver.cpp


	std::ifstream inputt("rhs");
	std::vector<RT> rhsVee{std::istream_iterator<RT>(inputt), std::istream_iterator<RT>() };
        inputt.close();
CT bnorr = 0.0;
CT berrr = 0.0;
for (int i=0; i < rhsVee.size(); ++i) {
	RT err = -rhsVee[i];
	for (int j=0; j < rhsVee.size(); ++j) {
		err += wm[i,j]*solutionCoefficientsOr[j];
	}
	bnorr += std::pow(std::abs(rhsVee[i]),2.0);
	berrr += std::pow(std::abs(err),2.0);
}
std::cout << "\n\n" << std::sqrt(berrr) << " =err,L2 nor= " << std::sqrt(bnorr) << "\n";



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


GridFunction<BFT, RT> projSol(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HminusHalfSpace), surfaceNormalIndependentFunction(solSphere()));

//arma::Col<RT> projVect(projSol.projections(boundaryOp->dualToRange()));
#ifdef WITH_TRILINOS
//	DefaultIterativeSolver<BFT, RT> solver(slpOp,ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);
	solver.saveProjections(projSol,"projVect");
	solver.saveProjections(rhs,"rhs");
#endif

	std::ifstream inputp("projVect");
	std::vector<RT> projVe{std::istream_iterator<RT>(inputp), std::istream_iterator<RT>() };
        inputp.close();

	std::ifstream input("rhs");
	std::vector<RT> rhsVe{std::istream_iterator<RT>(input), std::istream_iterator<RT>() };
        input.close();

//	std::vector<RT> diff = projVe-rhsVe;
	RT diff = 0.0;
	RT nor = 0.0;
	RT norV = 0.0;
	for(int i=0; i < projVe.size(); i++) {
//		diff = diff + (projVe[i]-rhsVe[i])^2.0;
//		diff = diff + std::abs(projVe[i]-rhsVe[i])*std::abs(projVe[i]-rhsVe[i]);
		diff = diff + std::abs(projVe[i]-solutionCoefficientsOr[i])*std::abs(projVe[i]-solutionCoefficientsOr[i]);
		norV = norV + std::abs(projVe[i])*std::abs(projVe[i]);
		nor = nor + std::abs(solutionCoefficientsOr[i])*std::abs(solutionCoefficientsOr[i]);
	}
	diff = std::sqrt(diff);
	nor = std::sqrt(nor);
	norV = std::sqrt(norV);
std::cout << projVe[0] << "asdf" << solutionCoefficientsOr[0] << "oiuh" << rhsVe[0] << "asdf" << projVe.size() << "\n";
	std::cout << errBCavm(sim,0) << "stopping" << diff << "apsoijfd" << nor << " , norV = " << norV << "\n";
return 0;
//break;
//continue;


// -------------- Compression -----------------------------------------
//	std::string str = "f 0.8"; //"f "+ std::to_string(Ts(ti)); //"c ";
	std::string str = "t 0"; //"t 2010";
	start = tbb::tick_count::now();

	BoundaryOperator<BFT, RT> slpOpCompr = helmholtz3dSingleLayerBoundaryOperator<BFT>(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HplusHalfSpace), make_shared_from_ref(HminusHalfSpace),waveNumber);
	boost::shared_ptr<const Bempp::AbstractBoundaryOperator<double, std::complex<double> > > asdf = slpOpCompr.abstractOperator();
	const GeneralElementarySingularIntegralOperator<BFT,RT,RT> bla = dynamic_cast<const GeneralElementarySingularIntegralOperator<BFT,RT,RT>& > (*asdf);

//	arma::Mat<RT> wmDummy(3,3);
//	wmDummy.fill(0.);
//	std::vector<RT> rhsVeDummy;
	boost::shared_ptr<const Bempp::DiscreteBoundaryOperator<RT> > weakCompr = bla.weakFormPeter(str,context,&solutionCoefficientsOr, &rhsVe, &wm);
	
	arma::Mat<RT> oneRow = weakCompr->asMatrix();


size_t testRow = 0;
std::string::size_type sz;
testRow = std::stof(str.substr(2),&sz);
RT tmpErr = 0.0;
RT globalNor = 0.0;
RT aprxb = 0.0;
	for (int i=0; i < oneRow.size(); ++i) {
		aprxb += oneRow[i,0]*projVe[i];
		tmpErr += std::pow(std::abs(wm[i,testRow]-oneRow[i,0]),2.0);
		globalNor += std::pow(std::abs(wm[i,testRow]),2.0);
if( (i == 0) || (i == oneRow.size() ) ) {std::cout << wm[i,testRow] << " =matr, row= " << oneRow[i,0] << " asdf " <<i << " oi "  << " rowtrans= " << oneRow[0,i] << " rowidx " << oneRow[i,testRow] << " rowtraidx " << oneRow[testRow,i] << "\n"; }
	}
std::cout << std::sqrt(tmpErr) << " =err first row, nor= "  << std::sqrt(globalNor) <<" , testRow= " << testRow << " aspojfd " << std::sqrt(tmpErr)/std::sqrt(globalNor) << "\n";

std::cout << oneRow.size() << " aos " << wm.size() << " s " << tmpErr << " e " << globalNor << " t " << tmpErr/globalNor << "\n\n";

std::cout << aprxb << " = aprxb, b = " << rhsVe[testRow] << " , testRow = " << testRow << "\n";

CT bnor = 0.0;
CT berr = 0.0;
for (int i=0; i < oneRow.size(); ++i) {
//	RT err = rhsVe[i];
	RT err = -rhsVe[i];
	for (int j=0; j < oneRow.size(); ++j) {
//		err += wm[i,j]*projVe[j];
//		err += wm[j,i]*solutionCoefficientsOr[j];
		err += wm[i,j]*solutionCoefficientsOr[j];
//		err += wm[j,i]*projVe[j];
	}
//	std::cout << " " << err;
	bnor += std::pow(std::abs(rhsVe[i]),2.0);
	berr += std::pow(std::abs(err),2.0);
}
std::cout << "\n\n" << std::sqrt(berr) << " =err,L2 nor= " << std::sqrt(bnor) << "\n";

break;

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
*/
//correlations();
}
//void correlations() {
// add correlation code here later
//}

