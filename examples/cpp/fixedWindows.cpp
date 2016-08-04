void fixedWindows() {

arma::Mat<RT> ks = arma::exp2(arma::linspace<arma::Mat<RT>>(3,4,2));
const int kl = ks.size();

arma::Mat<BFT> Ts = arma::linspace<arma::Mat<BFT>>(0.03,2.0,10);

const int Tl = Ts.size();

const int avm = 100;
arma::Mat<BFT> thetas = arma::zeros(avm,1);
arma::Mat<BFT> phis = arma::zeros(avm,1);
std::uniform_real_distribution<BFT> unif(0, 1);
std::random_device rd;
std::mt19937 gen(rd());
for(int av = 0; av < avm; av ++) {
	thetas(av) = M_PI*unif(gen);
	phis(av) = M_PI*2*unif(gen);
}
arma::Mat<CT> points(3,avm*avm);
arma::Mat<CT> pointsInt(3,avm*avm);
for (int thi =0; thi < avm; ++thi) {
	for(int phih =0; phih < avm; ++phih) {
		int idx = thi*avm+phih;
		points(0,idx) = cos(phis(phih))*sin(thetas(thi));
		points(1,idx) = sin(phis(phih))*sin(thetas(thi));
		points(2,idx) = cos(thetas(thi));
		BFT rtmp = unif(gen);
		pointsInt(0,idx) = rtmp*cos(phis(phih))*sin(thetas(thi));
		pointsInt(1,idx) = rtmp*sin(phis(phih))*sin(thetas(thi));
		pointsInt(2,idx) = rtmp*cos(thetas(thi));
	}
}

arma::Mat<BFT> conds = arma::zeros(kl,1+Tl);
arma::Mat<BFT> percs = arma::zeros(kl,Tl);

arma::Mat<BFT> errBCavm = arma::zeros(kl,1+Tl);
arma::Mat<BFT> errAxb = arma::zeros(kl,1+Tl);
arma::Mat<BFT> errInt = arma::zeros(kl,1+Tl);
arma::Mat<BFT> errSol = arma::zeros(kl,Tl);
arma::Mat<BFT> times = arma::zeros(kl,1+Tl);


for(int ki = 0; ki < kl; ki++) {
	tbb::tick_count start = tbb::tick_count::now();
	waveNumber = ks(ki);
	std::string mfs = "../../../meshes/sphere" + std::to_string(ki+1) + ".msh";
	shared_ptr<Grid> grid = loadTriangularMeshFromFile(mfs.c_str());

	PiecewiseLinearContinuousScalarSpace<BFT> HplusHalfSpace(grid);
	PiecewiseConstantScalarSpace<BFT> HminusHalfSpace(grid);
	AssemblyOptions assemblyOptions;
	assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
//	assemblyOptions.setVerbosityLevel(VerbosityLevel::HIGH); // More info (progress % matrix)
	// No ACA (AcaOptions) 
	AccuracyOptions accuracyOptions;
	accuracyOptions.doubleRegular.setRelativeQuadratureOrder(1);
	NumericalQuadratureStrategy<BFT, RT> quadStrategy(accuracyOptions);
	Context<BFT, RT> context(make_shared_from_ref(quadStrategy), assemblyOptions);

	BoundaryOperator<BFT, RT> slpOp = helmholtz3dSingleLayerBoundaryOperator<BFT>(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HplusHalfSpace), make_shared_from_ref(HminusHalfSpace),waveNumber);

	arma::Mat<RT> wm = slpOp.weakForm()->asMatrix();
	GridFunction<BFT, RT> rhs(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HminusHalfSpace), surfaceNormalIndependentFunction(MyFunctor()));

	// Initialize the solver
#ifdef WITH_TRILINOS
	DefaultIterativeSolver<BFT, RT>* solver = new DefaultIterativeSolver<BFT, RT>(slpOp,ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);
	solver->initializeSolver(defaultGmresParameterList(1e-5));
	solver->saveProjections(rhs,"rhsV");
	Solution<BFT, RT> solution = solver->solve(rhs);
#endif
	tbb::tick_count end = tbb::tick_count::now();
        times(ki,0) = (end - start).seconds();
	conds(ki,0) = arma::cond(wm);

	GridFunction<BFT, RT> solFun = solution.gridFunction();
	arma::Col<RT> solutionCoefficientsOr = solFun.coefficients(); // Same as armaSolution from def_iter_solver.cpp

	std::ifstream input("rhsV");
	std::vector<RT> rhsVe{std::istream_iterator<RT>(input), std::istream_iterator<RT>() };
        input.close();
	Helmholtz3dSingleLayerPotentialOperator<BFT> slPot (waveNumber);
	EvaluationOptions evalOptions = EvaluationOptions();
	arma::Mat<RT> potResOr = slPot.evaluateAtPoints(solFun, points, quadStrategy, evalOptions);
	arma::Mat<RT> diri = potResOr;
	arma::Mat<RT> dirInt = potResOr;
	MyFunctor tmp = MyFunctor();
	arma::Col<RT> t(1);
	for (int i = 0; i < avm*avm; ++i) {
		arma::Col<CT> pt(3);
		pt(0) = points(0,i);
		pt(1) = points(1,i);
		pt(2) = points(2,i);
		t.fill(0.);
		tmp.evaluate(pt,t);
		diri(i) = t(0);

		pt(0) = pointsInt(0,i);
		pt(1) = pointsInt(1,i);
		pt(2) = pointsInt(2,i);
		t.fill(0.);
		tmp.evaluate(pt,t);
		dirInt(i) = t(0);
	}
	arma::Mat<RT> errBCOr = potResOr - diri;
	errBCavm(ki,0) = mean(mean(abs(errBCOr) ))/mean(mean(abs(diri)));
	errInt(ki,0) = mean(mean(abs(slPot.evaluateAtPoints(solFun, pointsInt, quadStrategy, evalOptions)-dirInt) )); // solFun now has solution of original problem
	errInt(ki,0) = mean(mean(abs(slPot.evaluateAtPoints(solFun, pointsInt, quadStrategy, evalOptions)-dirInt) ))/mean(mean(abs(dirInt))); // solFun now has solution of original problem
	arma::Mat<RT> zxcv = slPot.evaluateAtPoints(solFun, pointsInt, quadStrategy, evalOptions);
	std::cout << dirInt(0) << " = dirint, slpint = " << zxcv(0) << "\n";

	CT bnorr = 0.0;
	BFT l1err = 0.0;
	BFT l1nor = 0.0;
	for (int i=0; i < wm.n_rows; ++i) {
	    RT err = -rhsVe[i];
	    for (int j=0; j < wm.n_cols; ++j) {
		err += wm(i,j)*solutionCoefficientsOr(j); // sco is arma::col so () iso []
	    }
	    bnorr += std::pow(std::abs(rhsVe[i]),2.0);
	    errAxb(ki,0) += std::pow(std::abs(err),2.0);
	    l1err += std::abs(err);
	    l1nor += std::abs(rhsVe[i]);
	}
	errAxb(ki,0) = std::sqrt(errAxb(ki,0)/bnorr);
	std::cout << l1nor << " =l1nor orig, l1err orig= " << l1err << "\n";

// -------------- Compression -----------------------------------------
   for(int ti = 0; ti < Tl; ti ++) {
	time_t now = time(0);
	std::cerr << "-----------------------------------\n" << ki << " = ki, ti = " << ti << ", perc = " << (0.0+ki*Tl + ti)/(0.0+Tl*kl) << " " << asctime(localtime(&now) ); 

	std::string str = "f   " + std::to_string(Ts(ti)); // Multiply element
//	std::string str = "k   " + std::to_string(Ts(ti)); // Multiply kernel
	start = tbb::tick_count::now();

	slpOp = helmholtz3dSingleLayerBoundaryOperator<BFT>(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HplusHalfSpace), make_shared_from_ref(HminusHalfSpace),waveNumber);
	boost::shared_ptr<const Bempp::AbstractBoundaryOperator<double, std::complex<double> > > asdf = slpOp.abstractOperator();
	const GeneralElementarySingularIntegralOperator<BFT,RT,RT> bla = dynamic_cast<const GeneralElementarySingularIntegralOperator<BFT,RT,RT>& > (*asdf);

	arma::Mat<RT> wmDummy(3,3);
	wmDummy.fill(0.);
	std::vector<RT> rhsVeDummy;
	boost::shared_ptr<const Bempp::DiscreteBoundaryOperator<RT> > weakCompr = bla.weakFormPeter(str,context,&solutionCoefficientsOr, &rhsVeDummy, &wmDummy);
	
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
	errBCavm(ki,1+ti) = mean(mean(abs(errBCCompr) ))/mean(mean(abs(diri)));
	errInt(ki,1+ti) = mean(mean(abs(slPot.evaluateAtPoints(solFun, pointsInt, quadStrategy, evalOptions)-dirInt) ))/mean(mean(abs(dirInt))); // solFun now has solution of compressed problem

	l1err = 0.0;
	l1nor = 0.0;
	for (int i=0; i < rhsVe.size(); ++i) {
	    RT err = -rhsVe[i];
	    for (int j=0; j < rhsVe.size(); ++j) {
		err += wm(i,j)*solutionCoefficientsOr(j); // sco is arma::col so () iso [] and wm is now compressed matrix
	    }
	    errAxb(ki,1+ti) += std::pow(std::abs(err),2.0);
	    if(std::abs(err)/std::sqrt(bnorr) > 0.025) {
		std::cout << err+rhsVe[i] << " =incorrect val row " << i << ", b_i= " << rhsVe[i] << ", relnorerr= " << std::abs(err)/std::sqrt(bnorr) << "\n";
	    }
	    l1nor += std::abs(rhsVe[i]);
	    l1err += std::abs(err);
	}
	std::cout << l1nor << " =l1nor, l1err= " << l1err << "\n";
	errAxb(ki,1+ti) = std::sqrt(errAxb(ki,1+ti)/bnorr);
   }
std::ofstream myfile;
myfile.open ("resFixWind");
myfile << "Output of tutorial_dirichlet.\n";
myfile << real(ks) << " = ks, Ts = " << std::endl << Ts << std::endl;
myfile << percs << " = perc, conds = " << std::endl << conds << std::endl;
myfile << times << " = times" << std::endl << std::endl;
myfile << errBCavm << " = errBCavm, errSol = " << std::endl << errSol << std::endl;
myfile << errInt << " = errInt, errAxb = \n" << errAxb << "\n";
myfile.close();
}

std::cout << real(ks) << " = ks, Ts = " << Ts << std::endl;
std::cout << percs << " = perc, conds = " << conds << std::endl;
std::cout << times << " = times, " << errBCavm << " = errBCavm, errSol = " << errSol << std::endl;
std::cout << errInt << " = errInt, errAxb = " << errAxb << "\n";
}
