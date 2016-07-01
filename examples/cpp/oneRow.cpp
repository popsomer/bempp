// Validate the asymptotic compression scheme by comparing one row of the compressed BEM matrix.

// gedit /opt/fb/bempp/lib/fiber/modified_helmholtz_3d_single_layer_potential_kernel_functor.hpp
// gedit /opt/fb/bempp/lib/assembly/dense_global_assembler.hpp

// cd /opt/fb/bempp/build/examples/cpp/
// pushd ../..; make tutorial_dirichlet -j6 && popd && ./tutorial_dirichlet || popd

// str = "c 0.8" corr through dist to corr, "d  0.6" corr through phys dist
// "f  0.6" fixed windows elements, "k  0.6" fixed windows kernel
// " i 0.6" illuminated only element, "j  1.9" illuminated only kernel
// "n   " normal BEM, "t          0" only row 0 full
// "a    1.2   7570 " only row 7570 fixed windows elements, "b    0.6  7570 " only row 7570 fixed windows kernel
// "e    1.5   0  " only row 0 with distance to correlation threshold: still O(N^2) because of computation of correlations
// number appearing from position 4 = b = max correlation distance with nonzero weight
// number from 9 = which row

#include "common.cpp"

class CstFct
{
public:
    typedef RT ValueType;
    typedef ScalarTraits<RT>::RealType CoordinateType;
    RT val;
    CstFct(RT given) {
	val = given;
    }
    int argumentDimension() const { return 3; }
    int resultDimension() const { return 1; }
    inline void evaluate(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
        result(0) = val;
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
	Hpc imu = std::complex<double>(0, 1);
	Hpc legnm1;
	Hpc legn;
	result = 0.0;
	Hp r = sqrt(pow(point(0)+0.0L,2.0L) +std::pow(point(1)+0.0L,2.0L) +std::pow(point(2)+0.0L,2.0L));
	if(abs(r-1) > 1e-6) {
		std::cerr << "r not close to 1\n";
		exit(1);
	}
	Hpc hnm1 = exp(imu*(abs(waveNumber)+0.0L))/(abs(waveNumber)+0.0L); 
	for(int n = 0; n < 1.4*std::abs(waveNumber)+40; n ++) {
		if(n == 0) {
			legn = 1;
		}
		else if(n == 1) {
			legn = point(0)/r;
			legnm1 = 1;
		}
		else {
			Hpc tmp = legn;
			legn = (std::complex<Hp>((2*n-1)*point(0)/r,0)*legn-std::complex<Hp>(n-1, 0)*legnm1)/std::complex<Hp>(n,0); // l+1 = n in recurrence relation
			legnm1 = tmp;
		}
		Hpc besk = boost::math::sph_bessel(n, abs(waveNumber) );
		Hpc besyk = boost::math::sph_neumann(n, abs(waveNumber) );
		Hpc term = -(2.0L*n+1.0)*std::pow(imu,n+0.0L)*besk*(abs(waveNumber)+0.0L)*legn*(hnm1-(n+1)/(abs(waveNumber)+0.0L)*(besk+imu*besyk) )/(besk+imu*besyk);
		if (abs(result(0)) > 1e6 ) {
			std::cerr << n << "=n,legn=" << legn << "=legn, bes-n-1/2=" << besyk/std::pow(-1.0L,n+1.0L)/std::pow(acos(-1.0L)/2.0L/(abs(waveNumber)+0.0L), 0.5L) << "\n";
			std::cerr << "besk=" << besk << "=besk, besyk=" << besyk << "=besyk, hnm1= " << hnm1 << "=hnm1\n";
			std::cerr << term << " = term, result=" << result << "\n";
			std::cerr << waveNumber << "=wavenr, \n ERROR: normal derivative cannot be this high, printed info above\n";
			exit(1);
		}
		result(0) += term;
		hnm1 = besk+imu*besyk;
	} // Add the normal derivative of the incident field
	result(0) = -std::complex<Hp>(real(result(0)), imag(result(0))) - imu*(abs(waveNumber)+0.0L)*(point(0)+0.0L)*exp((abs(waveNumber)+0.0L)*imu*(point(0)+0.0L));
    }
};

void oneRow(int kl)
{
struct stat buffer;   
if (stat (("../../../meshes/sphere" + std::to_string(kl) + ".msh").c_str(), &buffer) != 0) {
	std::cerr << "Error: meshes/sphere" << std::to_string(kl) << ".msh does not exist.\n Please create all sphere*.msh by loading sphere.geo into gmsh, executing 2D meshing, refining (*-1) times and saving.\n Or use other arguments for tutorial_dirichlet.\n";
	exit(kl+1);
}
const int avm = 40;
arma::Mat<BFT> thetas = arma::zeros(avm,1);
arma::Mat<BFT> phis = arma::zeros(avm,1);
std::uniform_real_distribution<BFT> unif(0, 1);
std::random_device rd;
std::mt19937 gen(11); // Constant seed
//std::default_random_engine gen; // Random seed
for(int av = 0; av < avm; av ++) {
	thetas(av) = M_PI*unif(gen);
	phis(av) = M_PI*2*unif(gen);
}
arma::Mat<CT> points(3,avm*avm);
// theta in [0,pi) and phi in [0, 2pi)
for (int thi =0; thi < avm; ++thi) {
	for(int phih =0; phih < avm; ++phih) {
		int idx = thi*avm+phih;
		points(0,idx) = cos(phis(phih))*sin(thetas(thi));
		points(1,idx) = sin(phis(phih))*sin(thetas(thi));
		points(2,idx) = cos(thetas(thi));
	}
}
arma::Mat<RT> ks = arma::exp2(arma::linspace<arma::Mat<RT>>(3,kl+2,kl));
char typeSim[][10] = {"t        ", "b   0.6  ", "a   0.6  ", "a   0.4  ", "a   0.2  "};
int nrTpts = 4;
arma::Mat<CT> testPts(nrTpts,3);
testPts(0,0) = -1.0; testPts(0,1) = 0.0; testPts(0,2) = 0.0;
testPts(1,0) = -0.3; testPts(1,1) = std::sqrt(1.0-0.09-0.25); testPts(1,2) = -0.5;
testPts(2,0) = 0.0; testPts(2,1) = 1.0/std::sqrt(2.0); testPts(2,2) = -1.0/std::sqrt(2.0);
testPts(3,0) = 0.6; testPts(3,1) = -0.2; testPts(3,2) = std::sqrt(1.0-0.36-0.04); // Best case in the illuminated region, Illuminated region, Transition region and Shadow region
int nrTsim = sizeof(typeSim)/sizeof(typeSim[0]);

std::cout << kl << "=kl, nrTpts=" << nrTpts << "=nrTpts, nrTsim=" << nrTsim << "\n";

arma::Cube<BFT> percs(nrTsim, nrTpts, kl, arma::fill::none); percs.fill(-1.0);
arma::Mat<BFT> errBCpts(nrTpts, kl, arma::fill::none); errBCpts.fill(-1.0);
arma::Col<BFT> errBCproj(kl, arma::fill::none); errBCproj.fill(-1.0);
arma::Cube<BFT> errowAprojb(nrTsim, nrTpts, kl, arma::fill::none); errowAprojb.fill(-1.0);
arma::Cube<BFT> times(nrTsim, nrTpts, kl, arma::fill::none); times.fill(-1.0);
arma::Col<BFT> timesProj(kl, arma::fill::none); timesProj.fill(-1.0);

// Initialise unused Solution because GridFunction.setCoefficients(...) cannot be called on an uninitialised GridFunction
waveNumber = ks(0);
AssemblyOptions assemblyOptions;
assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
AccuracyOptions accuracyOptions;
accuracyOptions.doubleRegular.setRelativeQuadratureOrder(1);
NumericalQuadratureStrategy<BFT, RT> quadStrategy(accuracyOptions);
Context<BFT, RT> context(make_shared_from_ref(quadStrategy), assemblyOptions);

for(int ki = 0; ki < kl; ki++) {
	tbb::tick_count starttk = tbb::tick_count::now();
	waveNumber = ks(ki);
//	std::string mfs = "../../../meshes/sphere" + std::to_string(ki+3) + ".msh"; // Shift the meshes for higher accuracy, but lower maximal wavenumbers.
	std::string mfs = "../../../meshes/sphere" + std::to_string(ki+1) + ".msh"; // When k doubles, double the number of elements in each spatial direction (# elements *= 4)
	shared_ptr<Grid> grid = loadTriangularMeshFromFile(mfs.c_str());
	PiecewiseLinearContinuousScalarSpace<BFT> HplusHalfSpace(grid);
	PiecewiseConstantScalarSpace<BFT> HminusHalfSpace(grid); // Supertype ScalarSpace is abstract so cannot define Hplus/min HS as such to choose cst/lin/quadr depending on for example str.at(2) like:
//	PiecewisePolynomialContinuousScalarSpace<BFT> HplusHalfSpace(grid,2);
	GridFunction<BFT, RT> rhs(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HminusHalfSpace), surfaceNormalIndependentFunction(MyFunctor()));
	time_t now = time(0);
	std::cerr << asctime(localtime(&now) ) << " = time before making GridFunction projSol at k = " << waveNumber << ", time after making GridF projSol = ";
	GridFunction<BFT, RT> projSol(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HminusHalfSpace), surfaceNormalIndependentFunction(solSphere()));
	now = time(0);
	std::cerr << asctime(localtime(&now) );
	// Divide by the norms of basis functions
	GridFunction<BFT, RT> normbas(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HminusHalfSpace), surfaceNormalIndependentFunction(CstFct(1.0)));
	arma::Col<RT> rhsCol = rhs.projections(make_shared_from_ref(HminusHalfSpace));
	arma::Col<RT> normBasCol = normbas.projections(make_shared_from_ref(HminusHalfSpace));
	arma::Col<RT> projCol = projSol.projections(make_shared_from_ref(HminusHalfSpace));
	for(int i=0; i < projCol.n_rows; i++) {
		projCol(i) /= normBasCol(i);
	}
	projSol.setCoefficients(projCol);

	arma::Mat<RT> diri = arma::Mat<RT>(1,avm*avm);
	MyFunctor tmp = MyFunctor();
	arma::Col<CT> pt(3);
	arma::Col<RT> t(1);
	for (int i = 0; i < avm*avm; ++i) {
		pt(0) = points(0,i);
		pt(1) = points(1,i);
		pt(2) = points(2,i);
		tmp.evaluate(pt,t);
		diri(i) = t(0);
	}
	EvaluationOptions evalOptions = EvaluationOptions();
	Helmholtz3dSingleLayerPotentialOperator<BFT> slPot (waveNumber);
	errBCproj(ki) = arma::mean(arma::mean(abs(slPot.evaluateAtPoints(projSol, points, quadStrategy, evalOptions) - diri) ))/arma::mean(arma::mean(abs(diri)));

	BoundaryOperator<BFT, RT> slpOp = helmholtz3dSingleLayerBoundaryOperator<BFT>(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HplusHalfSpace), make_shared_from_ref(HminusHalfSpace),waveNumber);
	boost::shared_ptr<const Bempp::AbstractBoundaryOperator<double, std::complex<double> > > asdf = slpOp.abstractOperator();
	const GeneralElementarySingularIntegralOperator<BFT,RT,RT> bla = dynamic_cast<const GeneralElementarySingularIntegralOperator<BFT,RT,RT>& > (*asdf);
	arma::Mat<RT> wmDummy(3,3, arma::fill::zeros);
	std::vector<RT> rhsVeDummy;
	arma::Col<RT> scoDummy;

	timesProj(ki) = (tbb::tick_count::now() - starttk).seconds();

	std::unique_ptr<GridView> gv = grid->leafView();
	arma::Mat<CT> vertices;
	arma::Mat<int> corners;
	arma::Mat<char> aux;
	gv->getRawElementData(vertices, corners, aux);
	arma::Mat<CT> baryCenters(3,corners.n_cols);
	for(int elmt = 0; elmt < corners.n_cols; elmt ++) {
		baryCenters(0,elmt) = (vertices(0,corners(0,elmt)) + vertices(0,corners(1,elmt)) + vertices(0,corners(2,elmt)) )/3;
		baryCenters(1,elmt) = (vertices(1,corners(0,elmt)) + vertices(1,corners(1,elmt)) + vertices(1,corners(2,elmt)) )/3;
		baryCenters(2,elmt) = (vertices(2,corners(0,elmt)) + vertices(2,corners(1,elmt)) + vertices(2,corners(2,elmt)) )/3;
	}

	for(int pti = 0; pti < nrTpts; pti ++) {
		pt(0) = testPts(pti,0);	pt(1) = testPts(pti,1);	pt(2) = testPts(pti,2);
		tmp.evaluate(pt,t);
		arma::Mat<RT> evB = slPot.evaluateAtPoints(projSol, pt, quadStrategy, evalOptions);
		errBCpts(pti,ki) = abs(evB(0,0) - t(0))/abs(t(0));
		int rowPti = -1;
		CT minDist = 2.1;
		for(int elmt = 0; elmt < corners.n_cols; elmt ++) {
			CT dist = std::sqrt(std::pow(pt(0) - baryCenters(0,elmt), 2.0) + std::pow(pt(1) - baryCenters(1,elmt), 2.0) + std::pow(pt(2) -baryCenters(2,elmt), 2.0) );
			if(dist < minDist) {
				minDist = dist;
				rowPti = elmt;
			}
		}
		for(int tsi = 0; tsi < nrTsim; tsi ++) {
			tbb::tick_count start = tbb::tick_count::now();
			std::string str = typeSim[tsi] + std::to_string(rowPti);
			time_t now = time(0);
			std::cout << tsi << "= tsi, pti=" << pti << " =pti, waveNr= " << waveNumber << "=k, rowPti =" << rowPti << ", str=" << str << ", " << asctime(localtime(&now) ); //"-------\n\n";
			std::string::size_type sz;
			int rowCheck = std::stof(str.substr(9),&sz);
			if(rowCheck != rowPti) {
				std::cerr << rowCheck << " = rowCheck is not rowPti = " << rowPti << ", exiting\n";
				exit(rowCheck-rowPti);
			}
			boost::shared_ptr<const Bempp::DiscreteBoundaryOperator<RT> > weakCompr = bla.weakFormPeter(str,context,&scoDummy, &rhsVeDummy, &wmDummy);
			arma::Mat<RT> wm = weakCompr->asMatrix(); // wm now contains one row of the compressed matrix if str starts with t, a or b
			percs(tsi,pti,ki) = arma::accu(wm != 0)/(0.0+wm.n_elem);

			RT err = -rhsCol(rowPti);
			for (int j=0; j < projCol.n_rows; ++j) {
				err += wm(0,j)*projCol(j); // projCol is arma::col so () iso [] which gives a wrong result iso an error
			}
			errowAprojb(tsi, pti, ki) = std::abs(err)/std::abs(rhsCol(rowPti) );
        		times(tsi,pti,ki) = (tbb::tick_count::now() - start).seconds(); // Could add timeProj
		} // End of loop over types of simulation
		std::ofstream myfile;
		myfile.open("resOneRow");
		myfile << "Output of tutorial_dirichlet with one row.\n";
		myfile << real(ks) << " = ks\ntestPts=" << testPts << "\ntypeSim = ";
		for(int qwer = 0; qwer < nrTsim; qwer ++) {
			myfile << typeSim[qwer] + std::string(" \n");
		}
		myfile << times << " = times, timesProj = \n" << timesProj << "\n";
		myfile << errBCproj << "= errBCproj, errBCpts = " << errBCpts << "\n";
		myfile << percs << " = perc, errowAprojb=" << errowAprojb << "\n";
		myfile.close();
		errowAprojb.save("errowAprojb.dat", arma::raw_ascii);
	}
} // End of loop over wavenumbers

nrTpts = std::system("cat res");
}


