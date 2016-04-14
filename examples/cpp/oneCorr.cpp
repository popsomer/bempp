// gedit /opt/fb/bempp/lib/fiber/modified_helmholtz_3d_single_layer_potential_kernel_functor.hpp
// gedit /opt/fb/bempp/lib/assembly/dense_global_assembler.hpp

// cd /opt/fb/bempp/build/examples/cpp/
// pushd ../..; make tutorial_dirichlet -j6 && popd && ./tutorial_dirichlet || popd

//Simpson:
// cd build
// cmake -DCMAKE_BUILD_TYPE=Release -DWITH_FENICS=ON .. -DCMAKE_CXX_FLAGS:STRING=-lpthread
// cd /export/home1/NoCsBack/nines/fb/bempp/build/examples/cpp/
// vim /export/home1/NoCsBack/nines/fb/bempp/examples/cpp/tutorial_dirichlet.cpp 
// pushd ../..; make tutorial_dirichlet -j14; popd
// ulimit -v 62000000

// str = "c 0.8" corr through dist to corr, "d  0.6" corr through phys dist
// "f  0.6" fixed windows elements, "k  0.6" fixed windows kernel
// " i 0.6" illuminated only element, "j  1.9" illuminated only kernel
// "n   " normal BEM, "t          0" only row 0 full
// "a    1.2   7570 " only row 7570 fixed windows elements, "b    0.6  7570 " only row 7570 fixed windows kernel
// "e    1.5   0  " only row 0 with distance to correlation threshold: still O(N^2) because of computation of correlations
// number appearing from position 4 = b = max correlation distance with nonzero weight
// number from 9 = which row

#include "common.cpp"

void oneCorr()
{
waveNumber = 16;

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

int nrTpts = 4;
arma::Mat<CT> testPts(nrTpts,3);
testPts(0,0) = -1.0; testPts(0,1) = 0.0; testPts(0,2) = 0.0;
testPts(1,0) = -0.3; testPts(1,1) = std::sqrt(1.0-0.09-0.25); testPts(1,2) = -0.5;
testPts(2,0) = 0.0; testPts(2,1) = 1.0/std::sqrt(2.0); testPts(2,2) = -1.0/std::sqrt(2.0);
testPts(3,0) = 0.6; testPts(3,1) = -0.2; testPts(3,2) = std::sqrt(1.0-0.36-0.04); // Best case in the illuminated region, Illuminated region, Transition region and Shadow region
int nrTsim = sizeof(typeSim)/sizeof(typeSim[0]);

std::cout << kl << "=kl, nrTpts=" << nrTpts << "=nrTpts\n";

arma::Vec<BFT> errBCpts(nrTpts, arma::fill::none); errBCpts.fill(-1.0);
arma::Vec<BFT> times(nrTpts, arma::fill::none); times.fill(-1.0);

// Initialise unused Solution because GridFunction.setCoefficients(...) cannot be called on an uninitialised GridFunction
waveNumber = ks(0);
AssemblyOptions assemblyOptions;
assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
AccuracyOptions accuracyOptions;
accuracyOptions.doubleRegular.setRelativeQuadratureOrder(1);
NumericalQuadratureStrategy<BFT, RT> quadStrategy(accuracyOptions);
Context<BFT, RT> context(make_shared_from_ref(quadStrategy), assemblyOptions);

tbb::tick_count starttk = tbb::tick_count::now();
shared_ptr<Grid> grid = loadTriangularMeshFromFile("../../../meshes/sphere1.msh");
PiecewiseLinearContinuousScalarSpace<BFT> HplusHalfSpace(grid);
PiecewiseConstantScalarSpace<BFT> HminusHalfSpace(grid); 
GridFunction<BFT, RT> rhs(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HminusHalfSpace), surfaceNormalIndependentFunction(MyFunctor()));
time_t now = time(0);

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
/*
std::cout << vertices.n_cols << "=cols,rows=" << vertices.n_rows << ", first vert (x,y,z)= (" << vertices(0,0) << ", " << vertices(1,0) << ", " << vertices(2,0) << ")\n";
std::cout << corners.n_cols << "=corCols, corRows=" << corners.n_rows << " asdf " << aux.n_rows << "," << aux.n_cols << "\n";
std::cout << "first cor= " << corners(0,0) << ", " << corners(1,0) << "," << corners(2,0) << ", " << corners(3,0) << "\n";
std::cout << "second cor= " << corners(0,1) << ", " << corners(1,1) << "," << corners(2,1) << ", " << corners(3,1) << "\n";
for(int zx =0; zx<2; zx ++) {
	std::cout << zx << "=element has first corner (" << vertices(0,corners(0,zx)) << ", " << vertices(1,corners(0,zx)) << ", " << vertices(2,corners(0,zx)) << ")\n";
	std::cout << zx << "=element has second corner (" << vertices(0,corners(1,zx)) << ", " << vertices(1,corners(1,zx)) << ", " << vertices(2,corners(1,zx)) << ")\n";
	std::cout << zx << "=element has third corner (" << vertices(0,corners(2,zx)) << ", " << vertices(1,corners(2,zx)) << ", " << vertices(2,corners(2,zx)) << ")\n";
}
exit(1);
*/
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
//		arma::Mat<RT> evB = slPot.evaluateAtPoints(solFunNew, pt, quadStrategy, evalOptions);
		errBCpts(pti,ki) = abs(evB(0,0) - t(0))/abs(t(0));
//std::cout << evB << "=evB, t=" << t << "\n";
//		int rowPti = DenseGlobalAssembler<BFT,RT>::closestElement(HplusHalfSpace,pt);
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

//std::cerr << rhsCol.n_cols << "aposdijf\n";
			RT err = -rhsCol(rowPti);
//std::cerr << projCol.n_rows << "jaspojff\n";
			for (int j=0; j < projCol.n_rows; ++j) {
//				std::cerr << j << "aposdijf" << wm.n_cols << "\n";
				err += wm(0,j)*projCol(j); // projCol is arma::col so () iso [] which gives a wrong result iso an error
//				err += wm(rowPti,j)*projCol(j); // projCol is arma::col so () iso [] which gives a wrong result iso an error
			}
			errowAprojb(tsi, pti, ki) = std::abs(err)/std::abs(rhsCol(rowPti) );
        		times(tsi,pti,ki) = (tbb::tick_count::now() - start).seconds(); // Could add timeProj
		} // End of loop over types of simulation
		std::ofstream myfile;
		myfile.open("res");
		myfile << "Output of tutorial_dirichlet with one row.\n";
		myfile << real(ks) << " = ks\ntestPts=" << testPts << "\ntypeSim = ";// << std::str(typeSim) << "\n, testPts=" << testPts;
		for(int qwer = 0; qwer < nrTsim; qwer ++) {
			myfile << typeSim[qwer] + std::string(" \n");
		}
		myfile << times << " = times, timesProj = \n" << timesProj << "\n";
		myfile << errBCproj << "= errBCproj, errBCpts = " << errBCpts << "\n";
		myfile << percs << " = perc, errowAprojb=" << errowAprojb << "\n";
		myfile.close();
		errowAprojb.save("errowAprojb.dat", arma::raw_ascii); 
//		if(pti == 3) { exit(0); }
	}

//nrTpts = std::system("cat res");
}


