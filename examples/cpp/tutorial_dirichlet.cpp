// gedit /opt/fb/bempp/lib/fiber/modified_helmholtz_3d_single_layer_potential_kernel_functor.hpp
// gedit /opt/fb/bempp/lib/linalg/default_iterative_solver.cpp

// gedit /opt/fb/bempp/lib/assembly/dense_global_assembler.hpp
// gedit /opt/fb/bempp/build/external/include/Trilinos/Thyra_BelosLinearOpWithSolve_def.hpp

// cd /opt/fb/bempp/build/examples/cpp/
// pushd ../..; make tutorial_dirichlet -j6 && popd && ./tutorial_dirichlet || popd

// pushd ../..; make tutorial_dirichlet -j6; popd
// ulimit -v 6000000
// ./tutorial_dirichlet >res 2>testgeom
// pushd ../..; make tutorial_dirichlet -j6; popd; ./tutorial_dirichlet 


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
// number appearing from position 4 = b = max correlation distance with nonzero weight
// number from 9 = which row

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
#include "space/piecewise_polynomial_continuous_scalar_space.hpp" // Peter
//#include <boost> //Peter: voor sph_bessel
#include "../../build/external/include/boost/math/special_functions/bessel.hpp" //Peter: voor sph_bessel
//#include <boost::math> //Peter: voor sph_bessel




using namespace Bempp;


typedef double BFT; // basis function type
typedef std::complex<double> RT; // result type (type used to represent discrete operators)
typedef double CT; // coordinate type


typedef std::complex<long double> Hpc;
typedef long double Hp; // High precisions for calculations

RT waveNumber;
//Hpc waveNumber;
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
class CstFct
{
public:
    typedef RT ValueType;
    typedef ScalarTraits<RT>::RealType CoordinateType;
    RT val;
    CstFct(RT given) {
	val = given;
    }
    // Number of components of the function's argument
    int argumentDimension() const { return 3; }
    int resultDimension() const { return 1; }
    // Evaluate the function at the point "point" and store result in the array "result"
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
	Hpc tmp;
	result = 0.0;
	Hp r = sqrt(pow(point(0)+0.0L,2.0L) +std::pow(point(1)+0.0L,2.0L) +std::pow(point(2)+0.0L,2.0L));
	if(abs(r-1) > 1e-6) {
		std::cerr << "r not close to 1\n";
		exit(1);
	}
	Hp ppvn = abs(waveNumber)+40;
	Hp prevNorm = abs(waveNumber)+4;
	Hpc hnm1 = exp(imu*(abs(waveNumber)+0.0L))/(abs(waveNumber)+0.0L); //*(1-2/imu/waveNumber);
	for(int n = 0; n < 1.4*std::abs(waveNumber)+40; n ++) {
		if(n == 0) {
			legn = 1;
		}
		else if(n == 1) {
			legn = point(0)/r;
			legnm1 = 1;
		}
		else {
			tmp = legn;
			legn = (std::complex<Hp>((2*n-1)*point(0)/r,0)*legn-std::complex<Hp>(n-1, 0)*legnm1)/std::complex<Hp>(n,0); // l+1 = n in recurrence relation
			legnm1 = tmp;
		}
/*
		Hp sinsumk = 0.0;
		Hp cossumk = 0.0;
		for(int k =0; k <= floor(n*0.5); k ++) {
//			sinsumk += pow(-1.0L,k+0.0L)*tgamma(1.0L+n+2*k)/pow(abs(waveNumber)+0.0L,(2*k+1.0L))/pow(2.0L,2.0L*k)/tgamma(1.0L+2*k)/tgamma(1.0L+n-2*k);
			sinsumk += pow(-1.0L,k+0.0L)*exp(lgamma(1.0L+n+2*k) - log(abs(waveNumber)+0.0L)*(2*k+1.0L) -log(2.0L)*2.0L*k -lgamma(1.0L+2*k) -lgamma(1.0L+n-2*k) );
		}
		for(int k =0; k <= floor((n-1)*0.5); k ++) {
//			cossumk += pow(-1.0L,k+0.0L)*tgamma(2.0L+n+2*k)/pow(abs(waveNumber)+0.0L,(2*k+2.0L))/pow(2.0L,2.0L*k+1)/tgamma(2.0L+2*k)/tgamma(n-2*k+0.0L);
			cossumk += pow(-1.0L,k+0.0L)*exp(lgamma(2.0L+n+2*k)-log(abs(waveNumber)+0.0L)*(2*k+2.0L) -log(2.0L)*(2.0L*k+1) -lgamma(2.0L+2*k) -lgamma(n-2*k+0.0L) );
		}
		Hpc besk = std::complex<Hp>(sin(abs(waveNumber)-n*acos(-1.0L)/2),0)*sinsumk + std::complex<Hp>(cos(abs(waveNumber)-n*acos(-1.0L)/2),0)*cossumk;
		Hpc besyk = -std::complex<Hp>(cos(abs(waveNumber)-n*acos(-1.0L)/2),0)*sinsumk + std::complex<Hp>(sin(abs(waveNumber)-n*acos(-1.0L)/2),0)*cossumk;
		
		besk = boost::cyl_bessel_j(n+0.5,abs(waveNumber));		
*/

/*
		Hpc besk = 0.0L;
		Hpc besyk = 0.0L;
//		for(int k=0; k < 1.5*n+20; k ++) {// Upper bound for k was guessed so check whether break-cond is met
		for(int k=0; k < n+20+abs(waveNumber); k ++) { // Upper bound for k was computed for when temporarily no osc
//		int k=0;
//		for(; abs(2.0L*k*log(abs(waveNumber)/2.0L) -lgamma(k+1.0L) -lgamma(k-n+0.5L) ) < -20; k ++) {// besyk term will be highest
//			besk += exp(4.0L*k*log(abs(waveNumber)/2.0L) -lgamma(2*k+1.0L) -lgamma(n+2*k+1.5L) )*(1.0L-abs(waveNumber)*abs(waveNumber)/4.0L/(2*k+1)/(n+2*k+1.5L));
			besk += exp(4.0L*k*log(abs(waveNumber)/2.0L) -lgamma(2*k+1.0L) -lgamma(n+2*k+1.5L) )*(1.0L-abs(waveNumber)*abs(waveNumber)/4.0L/(1.0L+2.0L*k)/(n+2*k+1.5L));
//			besk += pow(-1.0L,k)*exp(2.0L*k*log(abs(waveNumber)/2.0L) -lgamma(k+1.0L) -lgamma(n+k+1.5L) );
//			besyk += pow(-1.0L,k)*exp(2.0L*k*log(abs(waveNumber)/2.0L) -lgamma(k+1.0L) -lgamma(k-n+0.5L+0.L*imu) );
//			if( (k-n < 0) && (((n-k) % 2) == 1) ) { // because then gamma(k-n+0/5) < 0 and lgamma imag but double in cpp with only real component
			if( (k-n < 0) && (((n-2*k) % 2) == 1) ) { // because then gamma(k-n+0/5) < 0 and lgamma imag but double in cpp with only real component, changed meaning of k to 2k
//				besyk -= pow(-1.0L,k)*exp(2.0L*k*log(abs(waveNumber)/2.0L) -lgamma(k+1.0L) -lgamma(k-n+0.5L) ); // lgamma(k-n+0.5L) takes real part and imag part is pi
//				besyk -= exp(4.0L*k*log(abs(waveNumber)/2.0L) -lgamma(2*k+1.0L) -lgamma(2*k-n+0.5L) )*(1.0L-abs(waveNumber)*abs(waveNumber)/4.0L/std::max(1.0L,2.0L*k)/(2*k-n-0.5L));
				besyk -= exp(4.0L*k*log(abs(waveNumber)/2.0L) -lgamma(2*k+1.0L) -lgamma(2*k-n+0.5L) )*(1.0L-abs(waveNumber)*abs(waveNumber)/4.0L/(1.0L+2.0L*k)/(2*k-n+0.5L));
			} else {
//				besyk += pow(-1.0L,k)*exp(2.0L*k*log(abs(waveNumber)/2.0L) -lgamma(k+1.0L) -lgamma(k-n+0.5L) );
				besyk += exp(4.0L*k*log(abs(waveNumber)/2.0L) -lgamma(2*k+1.0L) -lgamma(2*k-n+0.5L) )*(1.0L-abs(waveNumber)*abs(waveNumber)/4.0L/(1.0L+2.0L*k)/(2*k-n+0.5L));
			}
//			if(abs(besk) > 1e6) { // besk = -inf, possibly because of 1/0
//			if(!std::isfinite(abs(besk)) ) { // besk = -inf, possibly because of 1/0
			if(abs(waveNumber-64.0) < 0.2) {
				std::cout << n << "=n,k=" << k << "=k, besk=" << besk << "=besk, exp=" << exp(4.0L*k*log(abs(waveNumber)/2.0L) -lgamma(2*k+1.0L) -lgamma(n+2*k+1.5L) ) << "=exp\n";
				std::cout << abs(waveNumber)*abs(waveNumber)/4.0L << "=fct, max=" << std::max(1.0L,2.0L*k) << "=max, n+2k+3/2=" << n+2*k+1.5L << "\n";
//				exit(1);
			}
//			std::cout << n << "=n, k=" << k << "=k, besk=" << besk << " = j, y= " << besyk << " = y, termy = " << pow(-1.0L,k)*exp(2.0L*k*log(abs(waveNumber)/2.0L) -lgamma(k+1.0L) -lgamma(k-n+0.5L) ) << "=term, cond1=" << (k-n < 0) << "=con1,cond2=" << (((n-k) % 2) == 1) <<"\n"; // wolfram alpha: sum( (-1)^k*exp(2*k*ln(8/2)-lgamma(k+1)-lgamma(k-10+0.5)), k=0..34)
//			std::cout << pow(-1.0L,k) << " = -1^k, zspk= " << exp(2.0L*k*log(abs(waveNumber)/2.0L) ) << "=zspk, fak1 = " << exp(-lgamma(k+1.0L)) << "=fak1, gam= " << exp(-lgamma(k-n+0.5L) ) << "=gamma, lgamma=" << -lgamma(k-n+0.5L) <<"\n";
		}
		besk *= sqrt(acos(-1.0L))/2.0L*std::pow(abs(waveNumber)/2.0L, n+0.0L);
//		besk *= sqrt(acos(-1.0L)/2.0L/abs(waveNumber))*std::pow(abs(waveNumber)/2.0L, n+0.5L);
		Hpc besMn = std::pow(abs(waveNumber)/2.0L, -0.5L-n)*besyk;
		besyk *= std::pow(-1.0L,n+1.0L)*sqrt(acos(-1.0L)/abs(waveNumber)/2.0L)*std::pow(abs(waveNumber)/2.0L, -n-0.5L);

*/
		Hpc besk = boost::math::sph_bessel(n, abs(waveNumber) );
		Hpc besyk = boost::math::sph_neumann(n, abs(waveNumber) );
//		besyk *= std::pow(-1.0L,n+1.0L)*sqrt(acos(-1.0L))/abs(waveNumber)*std::pow(2.0L/abs(waveNumber), n+0.0L);
		// Upper bound for k was guessed so check whether break-cond is met
		Hpc term = -(2.0L*n+1.0)*std::pow(imu,n+0.0L)*besk*(abs(waveNumber)+0.0L)*legn*(hnm1-(n+1)/(abs(waveNumber)+0.0L)*(besk+imu*besyk) )/(besk+imu*besyk);
		CoordinateType curNor = abs(term);

//		if(n==1) { 
//		if ((curNor > 3.1*prevNorm) && (n >9) ) {
//		if ((curNor > 3.1*prevNorm) && (prevNorm > 1.2*ppvn) && (n > 1.2*abs(waveNumber)+10) ) {
//		if( (!std::isnormal(1/3.0L+abs(2.0L*k*log(abs(waveNumber)/2.0L) -lgamma(k+1.0L) -lgamma(k-n+0.5L) ) ) ) || (abs(result(0)) > 1e6 ) ) {
//		if ((curNor > 3.1*prevNorm) && (n > 1.2*abs(waveNumber)+10) ) {
//		if (n == 10) {
		if (abs(result(0)) > 1e6 ) {
//			Hpc besMn = besyk/std::pow(-1.0L,n+1.0L)/sqrt(acos(-1.0L)/2.0L/(abs(waveNumber)+0.0L));
//			Hpc besMn = besyk/std::pow(-1.0L,n+1.0L)/std::pow(acos(-1.0L)/2.0L/(abs(waveNumber)+0.0L), 0.5L);
//std::cerr << k << "=k";
			std::cerr << n << "=n,legn=" << legn << "=legn, bes-n-1/2=" << besyk/std::pow(-1.0L,n+1.0L)/std::pow(acos(-1.0L)/2.0L/(abs(waveNumber)+0.0L), 0.5L) << "\n"; //, sinsumk=" << sinsumk << "=sinsumk\n";
			std::cerr << "besk=" << besk << "=besk, besyk=" << besyk << "=besyk, hnm1= " << hnm1 << "=hnm1\n";
			std::cerr << term << " = term, prevNorm = " << prevNorm << "=prevNorm, result=" << result << "\n";
			std::cerr << pow(-1.0L,1.5*n+20)*exp(2.0L*(1.5*n+20)*log(abs(waveNumber)/2.0L) -lgamma(1.5*n+20+1.0L) -lgamma(1.5*n+20-n+0.5L) ) << "=maxterm, mod-0.5=" << 0 % 2 << "=mod-0.5, mod-1.5=" << 1 % 2 << "=mod-1/5, mod-2.5=" << 2 % 2 << ", factbesk=" << sqrt(acos(-1.0L))/2.0L*std::pow(abs(waveNumber)/2.0L, n+0.0L) << "\n";
			std::cerr << waveNumber << "=wavenr,abs k-32=" << abs(waveNumber-32.0) << ", sphbes=" << boost::math::sph_bessel(n, abs(waveNumber) ) << " \n ERROR: normal derivative cannot be this high, printed info above\n";
			exit(1);
		}
		result(0) += term;
		ppvn = prevNorm;
		prevNorm = curNor;
		hnm1 = besk+imu*besyk;
	}
//	result(0) = result(0); // Do nothing
//	result(0) = -result(0);
//	result(0) += imu*(abs(waveNumber)+0.0L)*(point(0)+0.0L)*exp((abs(waveNumber)+0.0L)*imu*(point(0)+0.0L));
//	result(0) -= imu*(abs(waveNumber)+0.0L)*(point(0)+0.0L)*exp((abs(waveNumber)+0.0L)*imu*(point(0)+0.0L));
//	result(0) = -std::complex<Hp>(real(result(0)), imag(result(0))) + imu*(abs(waveNumber)+0.0L)*(point(0)+0.0L)*exp((abs(waveNumber)+0.0L)*imu*(point(0)+0.0L));
	result(0) = -std::complex<Hp>(real(result(0)), imag(result(0))) - imu*(abs(waveNumber)+0.0L)*(point(0)+0.0L)*exp((abs(waveNumber)+0.0L)*imu*(point(0)+0.0L));
//	result(0) = imu*(abs(waveNumber)+0.0L)*(point(0)+0.0L)*exp((abs(waveNumber)+0.0L)*imu*(point(0)+0.0L));
//	result(0) = -imu*(abs(waveNumber)+0.0L)*(point(0)+0.0L)*exp((abs(waveNumber)+0.0L)*imu*(point(0)+0.0L));
    }

    inline void evaluateOld(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
//	using namespace boost::math;
	RT imu = std::complex<double>(0, 1);
//        result(0) = 0.0;
	RT legnm1;
	RT legn;
	RT tmp;
	result = 0.0;
	CoordinateType r = std::sqrt(std::pow(point(0),2) +std::pow(point(1),2) +std::pow(point(2),2));
//	for(int n = 0; n < 1.4*std::abs(waveNumber)+40; n ++) {
	for(unsigned n = 0; n < 1.4*std::abs(waveNumber)+40; n ++) {
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
//			legn = (std::complex<double>(2*n+1,0)*point(0)*legn-std::complex<double>(n, 0)*legnm1)/std::complex<double>(n+1,0);
			legn = (std::complex<double>(2*n-1,0)*point(0)/r*legn-std::complex<double>(n-1, 0)*legnm1)/std::complex<double>(n,0); // l+1 = n in recurrence relation and maybe /r is useful though r approx 1
			legnm1 = tmp;
		}
//		result += (2.0*n+1.0)*imu^n;//*boost::cyl_bessel_j(n,wavek*1)
//		result += (2.0*n+1.0)*std::pow(imu,n)*boost::math::cyl_bessel_j(n,waveNumber*1.0);
//		result -= (2.0*n+1.0)*std::pow(imu,n)*jn(n,std::abs(waveNumber)*1.0)*legendre_pL(n,point(0));
//		result -= (2.0*n+1.0)*std::pow(imu,n)*jn(n,std::abs(waveNumber)*1.0)*legn;
//		result -= (2.0*n+1.0)*std::pow(imu,n)*math::sph_bessel(n,std::abs(waveNumber)*1.0)*legn;
//		result -= (2.0*n+1.0)*std::pow(imu,n)*sph_jn(n,std::abs(waveNumber)*1.0)*legn;
//		result -= (2.0*n+1.0)*std::pow(imu,n)*sqrt(acos(-1.0)/2/std::abs(waveNumber) )*jn(n+0.5,std::abs(waveNumber)*1.0)*legn;
		result -= waveNumber*(2.0*n+1.0)*std::pow(imu,n)*sqrt(acos(-1.0)/2/std::abs(waveNumber) )*jn(n+0.5,std::abs(waveNumber)*1.0)*legn;
	}
//	CoordinateType r = std::sqrt(std::pow(point(0),2) +std::pow(point(1),2) +std::pow(point(2),2));
	if(result.n_elem != 1) {
		std::cerr << " Error nelem not 1. , " << result(1);
//	} else if((r > 1+1e-7) || (r < 1-1e-7)) {
	} else if((r > 1+1e-2) || (r < 1-1e-2)) {
		std::cerr << " r out of bounds=" << r;
	}
	result -= imu*waveNumber*point(0)*exp(waveNumber*imu*point(0)); // r=1
//	result += imu*waveNumber*point(0)*exp(waveNumber*imu*point(0));
//	result = -result + imu*waveNumber*point(0)*exp(waveNumber*imu*point(0));
//	result = -result - imu*waveNumber*point(0)*exp(waveNumber*imu*point(0));
//	result = imu*waveNumber*point(0)*exp(waveNumber*imu*point(0));
//	result = -imu*waveNumber*point(0)*exp(waveNumber*imu*point(0));
//	result = -result;
//	result = +result;
    }

    inline double factorial(int n)
    {
	return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
    }
    inline void evaluateF(const arma::Col<CoordinateType>& point,
                         arma::Col<ValueType>& result) const {
//	RT imu = std::complex<double>(0, 1);
//	RT legnm1;
//	RT legn;
//	RT tmp;
	Hpc imu = std::complex<double>(0, 1);
	Hpc legnm1;
	Hpc legn;
	Hpc tmp;
	result = 0.0;
//	CoordinateType r = std::sqrt(std::pow(point(0),2) +std::pow(point(1),2) +std::pow(point(2),2));
//	CoordinateType prevNorm = abs(waveNumber)+4;
	Hp r = sqrt(pow(point(0)+0.0L,2.0L) +std::pow(point(1)+0.0L,2.0L) +std::pow(point(2)+0.0L,2.0L));
	Hp prevNorm = abs(waveNumber)+4;
	for(int n = 0; n < 1.4*std::abs(waveNumber)+40; n ++) {
		if(n == 0) {
			legn = 1;
		}
		else if(n == 1) {
			legn = point(0)/r;
			legnm1 = 1;
		}
		else {
			tmp = legn;
//			legn = (std::complex<double>(2*n+1,0)*point(0)/r*legn-std::complex<double>(n, 0)*legnm1)/std::complex<double>(n+1,0);
//			legn = (std::complex<double>(2*n-1,0)*point(0)/r*legn-std::complex<double>(n-1, 0)*legnm1)/std::complex<double>(n,0); // l+1 = n in recurrence relation
			legn = (std::complex<Hp>((2*n-1)*point(0)/r,0)*legn-std::complex<Hp>(n-1, 0)*legnm1)/std::complex<Hp>(n,0); // l+1 = n in recurrence relation
			legnm1 = tmp;
		}
//		std::cout << n << " =n, x= " << point(0)/r << ", legn = " << legn;
//		result += (2.0*n+1.0)*std::pow(imu,n)*sqrt(acos(-1.0)/2/std::abs(waveNumber*r) )*jn(n+0.5,std::abs(waveNumber)*r)*legn*(jn(n+0.5,abs(waveNumber)*r) + imu*yn(n+0.5,abs(waveNumber)*r) )/(jn(n+0.5,abs(waveNumber)) + imu*jn(n+0.5,abs(waveNumber)))/sqrt(r);
//		result -= (2.0*n+1.0)*std::pow(imu,n)*sqrt(acos(-1.0)/2/std::abs(waveNumber) )*jn(n+0.5,std::abs(waveNumber))*legn*(jn(n+0.5,abs(waveNumber)*r) + imu*yn(n+0.5,abs(waveNumber)*r) )/(jn(n+0.5,abs(waveNumber)) + imu*yn(n+0.5,abs(waveNumber)))/sqrt(r);

//		ValueType sinsumk = 0.0;
		Hp sinsumk = 0.0;
		Hp cossumk = 0.0;
		Hp sinsumkr = 0.0;
		Hp cossumkr = 0.0;
		for(int k =0; k <= floor(n*0.5); k ++) {
//			sinsumk += pow(-1.0,k)*this->factorial(n+2*k)/pow(waveNumber,(2*k+1))/pow(2.0,2*k)/factorial(2*k)*1.0/factorial(n-2*k);
			sinsumk += pow(-1.0L,k+0.0L)*tgamma(1.0L+n+2*k)/pow(abs(waveNumber)+0.0L,(2*k+1.0L))/pow(2.0L,2.0L*k)/tgamma(1.0L+2*k)/tgamma(1.0L+n-2*k);
			sinsumkr += pow(-1.0L,k+0.0L)*tgamma(1.0L+n+2*k)/pow(abs(waveNumber)*r+0.0L,(2*k+1.0L))/pow(2.0L,2.0L*k)/tgamma(1.0L+2*k)/tgamma(1.0L+n-2*k);
		}
		for(int k =0; k <= floor((n-1)*0.5); k ++) {
			cossumk += pow(-1.0L,k+0.0L)*tgamma(2.0L+n+2*k)/pow(abs(waveNumber)+0.0L,(2*k+2.0L))/pow(2.0L,2.0L*k+1)/tgamma(2.0L+2*k)/tgamma(n-2*k+0.0L);
			cossumkr += pow(-1.0L,k+0.0L)*tgamma(2.0L+n+2*k)/pow(abs(waveNumber)*r+0.0L,(2*k+2.0L))/pow(2.0L,2*k+1.0L)/tgamma(2.0L+2*k)/tgamma(n-2*k+0.0L);
//			cossumkr += pow(-1.0,k)*tgamma(1+n+2*k+1)/pow(waveNumber*r,(2*k+2))/pow(2.0,2*k+1)/tgamma(1+2*k+1)*1.0/tgamma(n-2*k);
		}
//		ValueType besk = sin(waveNumber-n*acos(-1.0)/2)*sinsumk + cos(waveNumber-n*acos(-1.0)/2)*cossumk;
		Hpc besk = std::complex<Hp>(sin(abs(waveNumber)-n*acos(-1.0L)/2),0)*sinsumk + std::complex<Hp>(cos(abs(waveNumber)-n*acos(-1.0L)/2),0)*cossumk;
		Hpc beskr = std::complex<Hp>(sin(abs(waveNumber)*r-n*acos(-1.0L)/2),0)*sinsumkr + std::complex<Hp>(cos(abs(waveNumber)*r-n*acos(-1.0L)/2),0)*cossumkr;
		Hpc besyk = -std::complex<Hp>(cos(abs(waveNumber)-n*acos(-1.0L)/2),0)*sinsumk + std::complex<Hp>(sin(abs(waveNumber)-n*acos(-1.0L)/2),0)*cossumk;
		Hpc besykr = -std::complex<Hp>(cos(abs(waveNumber)*r-n*acos(-1.0L)/2),0)*sinsumkr + std::complex<Hp>(sin(abs(waveNumber)*r-n*acos(-1.0L)/2),0)*cossumkr;
		
		Hpc term = (2.0L*n+1.0)*std::pow(imu,n+0.0L)*besk*legn*(beskr+imu*besykr)/(besk+imu*besyk);
		CoordinateType curNor = abs(term);//norm(term);
//		if(curNor > prevNorm*3.1) {
//		if(curNor > 1) {
		if((curNor > 3.1*prevNorm) && (n >9) ) {
//			std::cerr << curNor << "=curnor, prevnor=" << prevNorm << "\n";
			break;
		}
		result(0) -= term;
		prevNorm = curNor;
		if(false) {
		    std::cerr << n << "=n,term=" << (2.0L*n+1.0)*std::pow(imu,n)*besk*legn*(beskr+imu*besykr)/(besk+imu*besyk) << "=term,res=" << result << r << "=r, costheta=" << point(0)/r << "=coshteta, jn=" << besk << "\n";
//		    std::cerr << n << "=n,term=" << (2.0*n+1.0)*std::pow(imu,n)*sqrt(acos(-1.0)/2/std::abs(waveNumber) )*jn(n+0.5,std::abs(waveNumber))*legn*(jn(n+0.5,abs(waveNumber)*r) + imu*yn(n+0.5,abs(waveNumber)*r) )/(jn(n+0.5,abs(waveNumber)) + imu*yn(n+0.5,abs(waveNumber)))/sqrt(r) << "=term,res=" << result << r << "=r, costheta=" << point(0)/r << "=coshteta, jn=" << jn(n,abs(waveNumber)) << " " << jn(n+0.5,std::abs(waveNumber)) << " " << jn(0.5+n,std::abs(waveNumber)) << "\n";
//		    std::cerr//< sph_bessel(n,waveNumber) << "\n";
		    std::cerr << legn << "=legn, nr=" << (2.0L*n+1.0)*std::pow(imu,n) << "=nr, bes=" << besk << "\n"; //sqrt(acos(-1.0)/2/std::abs(waveNumber) )*jn(n+0.5,std::abs(waveNumber)) << "\n";
		    std::cerr << beskr+imu*besykr << "=hankel kr, hankel k=" << besk + imu*besyk << "\n";
//		    std::cerr << sqrt(acos(-1.0)/2/std::abs(waveNumber)*r)*(jn(n+0.5,abs(waveNumber)*r) + imu*yn(n+0.5,abs(waveNumber)*r) ) << "=hankel krm hankel k=" << sqrt(acos(-1.0)/2/std::abs(waveNumber) )*(jn(n+0.5,abs(waveNumber)) + imu*yn(n+0.5,abs(waveNumber))) << "\n";
		}
//		std::cout << n << " =n, costheta= " << point(0)/r << ", r= " << r << ", k=" << waveNumber << ", result = " << result << "\n";
//wolfram alpha:-sum( sqrt(-1)^n*(2*n+1)*j_n(8)*(h_n^(1)(0.1))/(h_n^(1)(8))*p(n,0.0646678), n=0..0)
// sphbess=
	}

//exit(0);
/*
	if(result.n_elem != 1) {
		std::cerr << " Error nelem not 1. , " << result(1);
	}
	if(abs(result(0)) > 3) {
//	if(true) {
		for(int n = 0; n < 1.4*std::abs(waveNumber)+40; n ++) {
		    std::cerr << n << "=n, " << (2.0*n+1.0)*std::pow(imu,n)*sqrt(acos(-1.0)/2/std::abs(waveNumber*r) )*jn(n+0.5,std::abs(waveNumber)*r)*legn*(jn(n+0.5,abs(waveNumber)*r) + imu*yn(n+0.5,abs(waveNumber)*r) )/(jn(n+0.5,abs(waveNumber)) + imu*yn(n+0.5,abs(waveNumber)))/sqrt(r) << ", r=" << r << ", waveNumber = " << waveNumber << ", costheta" << point(0)/r << "\n";
//		    std::cerr << (2.0*n+1.0)*std::pow(imu,n)*sqrt(acos(-1.0)/2/std::abs(waveNumber*r) ) << ", jn=" << jn(n+0.5,std::abs(waveNumber)*r) << ", legn = " << legn << ", ratio sph hankel = " << (jn(n+0.5,abs(waveNumber)*r) + imu*yn(n+0.5,abs(waveNumber)*r) )/(jn(n+0.5,abs(waveNumber)) + imu*yn(n+0.5,abs(waveNumber)))/sqrt(r) << "\n";
		    std::cerr << (2.0*n+1.0)*std::pow(imu,n)*sqrt(acos(-1.0)/2/std::abs(waveNumber) ) << ", jn=" << jn(n+0.5,std::abs(waveNumber) ) << ", legn = " << legn << ", ratio sph hankel = " << (jn(n+0.5,abs(waveNumber)*r ) + imu*yn(n+0.5,abs(waveNumber)*r) )/(jn(n+0.5,abs(waveNumber)) + imu*yn(n+0.5,abs(waveNumber)))/sqrt(r) << "\n";
		    std::cerr << jn(n+0.5,abs(waveNumber)*r) << " = jnkr, iynkr= " << imu*yn(n+0.5,abs(waveNumber)*r) << ", jnk = " << jn(n+0.5,abs(waveNumber)) << ", ynk = " << imu*yn(n+0.5,abs(waveNumber)) << "\n";
		}
		exit(0);
	}
*/
    }
};

void oneRow()
{
//arma::Mat<RT> ks = arma::exp2(arma::linspace<arma::Mat<RT>>(2,3,2));
//arma::Mat<RT> ks = arma::exp2(arma::linspace<arma::Mat<RT>>(3,4,2));
//arma::Mat<RT> ks = arma::exp2(arma::linspace<arma::Mat<RT>>(3,5,3));
arma::Mat<RT> ks = arma::exp2(arma::linspace<arma::Mat<RT>>(3,7,5));
const int kl = ks.size();

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
//int maxss = 6;
int nrSimWoWind = 8;
int maxss = 3*nrSimWoWind;
arma::Mat<BFT> percs = arma::zeros(maxss,1);
arma::Mat<BFT> errBCavm = arma::zeros(maxss,1);
errBCavm.fill(-1.0); // No meaningful values for when only making one row because no exact solution vector to compare with
arma::Mat<BFT> errBCproj = arma::zeros(maxss,1);
arma::Mat<BFT> errAxb = arma::zeros(maxss,1);
errAxb.fill(-1.0); // No meaningful values for when only making one row because no exact solution vector to compare with
arma::Mat<BFT> errAprojb = arma::zeros(maxss,1);
arma::Mat<BFT> errProj = arma::zeros(maxss,1);
errProj.fill(-1.0); // No meaningful values for when only making one row because no exact solution vector to compare with
arma::Mat<BFT> times = arma::zeros(maxss,1);

Solution<BFT, RT> solution; // Reuse the Solution<BFT,RT> from first simulation and overwrite the coefficients because GridFunction.setCoefficients(...) cannot be called on an uninitialised GridFunction

for(int sim = 0; sim < maxss; sim++) {
//for(int sim = nrSimWoWind; sim < maxss; sim++) {
	tbb::tick_count start = tbb::tick_count::now();
	shared_ptr<Grid> grid;
	std::string str;
//	ScalarSpace<BFT> HplusHalfSpace = PiecewiseConstantScalarSpace<BFT>(loadTriangularMeshFromFile("../../../meshes/sphere0.msh"));
//	ScalarSpace<BFT> HminusHalfSpace = PiecewiseConstantScalarSpace<BFT>(loadTriangularMeshFromFile("../../../meshes/sphere0.msh"));
//	ScalarSpace<BFT> HplusHalfSpace = NULL;
//	ScalarSpace<BFT> HminusHalfSpace = &NULL;
//	int minOrd = 0;
	if(sim == 0) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere1.msh");
		waveNumber = ks(0);
//		PiecewiseLinearContinuousScalarSpace<BFT> HplusHalfSpace(grid);
//		PiecewiseConstantScalarSpace<BFT> HminusHalfSpace(grid);
		str = "n               "; // need at least one space after 'n' for str(2) == 'l'-test in dga.hpp and 9 to be sure for when thinking t
	} else if(sim == 1) {
		if(false){
			grid = loadTriangularMeshFromFile("../../../meshes/sphere2.msh");
			waveNumber = ks(1);
			str = "n    ";
		} else if (false) {
			grid = loadTriangularMeshFromFile("../../../meshes/sphere4.msh");	
			waveNumber = ks(3);
			str = "t         0";
		} else {
			grid = loadTriangularMeshFromFile("../../../meshes/sphere1.msh");	
			waveNumber = ks(0);
			str = "t         0";
		}
	} else if(sim == 2) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere1.msh");
//		grid = loadTriangularMeshFromFile("../../../meshes/sphere2.msh");	
		waveNumber = ks(1);
//		minOrd = 1;
		str = "t         0 ";
	} else if(sim == 3) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere1.msh");	
		waveNumber = ks(1);
		str = "n           ";
	} else if(sim == 4) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere2.msh");	
		waveNumber = ks(1);
//		PiecewisePolynomialContinuousScalarSpace<BFT> HplusHalfSpace(grid,2);
//		PiecewiseLinearContinuousScalarSpace<BFT> HminusHalfSpace(grid);
//		minOrd = 1;
		str = "t         0 ";
	} else if(sim == 5) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere3.msh");	
		waveNumber = ks(2);
		str = "t         0 ";
	} else if(sim == 6) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere4.msh");	
		waveNumber = ks(3);
		str = "t         0 ";
	} else if(sim == 7) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere5.msh");	
		waveNumber = ks(4);
		str = "t         0 ";
	// Start using windows, now a = multiply elements
	} else if(sim == nrSimWoWind + 0) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere1.msh");
		waveNumber = ks(0);
		str = "f   0.8     ";
	} else if(sim == nrSimWoWind + 1) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere1.msh");	
		waveNumber = ks(0);
		str = "a   0.8   0 ";
	} else if(sim == nrSimWoWind + 2) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere1.msh");	
		waveNumber = ks(1);
		str = "a   0.8   0 ";
	} else if(sim == nrSimWoWind + 3) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere1.msh");	
		waveNumber = ks(1);
		str = "f   0.8     ";
	} else if(sim == nrSimWoWind + 4) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere2.msh");	
		waveNumber = ks(1);
		str = "a   0.8   0 ";
	} else if(sim == nrSimWoWind + 5) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere3.msh");	
		waveNumber = ks(2);
		str = "a   0.8   0 ";
	} else if(sim == nrSimWoWind + 6) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere4.msh");	
		waveNumber = ks(3);
		str = "a   0.8   0 ";
	} else if(sim == nrSimWoWind + 7) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere5.msh");	
		waveNumber = ks(4);
		str = "a   0.8   0 ";
	// Now b = multiply kernel
	} else if(sim == 2*nrSimWoWind + 0) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere1.msh");
		waveNumber = ks(0);
		str = "k   0.8      ";
	} else if(sim == 2*nrSimWoWind + 1) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere1.msh");	
		waveNumber = ks(0);
		str = "b   0.8   0 ";
	} else if(sim == 2*nrSimWoWind + 2) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere1.msh");	
		waveNumber = ks(1);
		str = "b   0.8    0 ";
	} else if(sim == 2*nrSimWoWind + 3) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere1.msh");	
		waveNumber = ks(1);
		str = "k   0.8       ";
	} else if(sim == 2*nrSimWoWind + 4) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere2.msh");	
		waveNumber = ks(1);
		str = "b   0.8   0 ";
	} else if(sim == 2*nrSimWoWind + 5) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere3.msh");	
		waveNumber = ks(2);
		str = "b   0.8   0 ";
	} else if(sim == 2*nrSimWoWind + 6) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere4.msh");	
		waveNumber = ks(3);
		str = "b   0.8   0 ";
	} else if(sim == 2*nrSimWoWind + 7) {
		grid = loadTriangularMeshFromFile("../../../meshes/sphere5.msh");	
		waveNumber = ks(4);
		str = "b   0.8   0 ";
	} else{
		grid = loadTriangularMeshFromFile("../../../meshes/sphere2.msh");
		waveNumber = ks(1);
		str = "n               ";
//		minOrd = 1;
	}
//	CT factorProj = 410.36/1.74179; // Probably dependent on k, mesh and choice of Piecewise ScalarSpaces
	PiecewiseLinearContinuousScalarSpace<BFT> HplusHalfSpace(grid);
	PiecewiseConstantScalarSpace<BFT> HminusHalfSpace(grid); // cannot do this in the ifblocks above because supertype ScalarSpace is abstract
//	PiecewisePolynomialContinuousScalarSpace<BFT> HplusHalfSpace(grid,minOrd+1);
//	PiecewisePolynomialContinuousScalarSpace<BFT> HminusHalfSpace(grid,minOrd);

	std::cout << "\n----  sim = " << sim << ", waveNr= " << waveNumber << "-------\n\n";
	AssemblyOptions assemblyOptions;
	assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);

	AccuracyOptions accuracyOptions;
	accuracyOptions.doubleRegular.setRelativeQuadratureOrder(1);
	NumericalQuadratureStrategy<BFT, RT> quadStrategy(accuracyOptions);
	Context<BFT, RT> context(make_shared_from_ref(quadStrategy), assemblyOptions);
	
	GridFunction<BFT, RT> rhs(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HminusHalfSpace),  surfaceNormalIndependentFunction(MyFunctor()));

	time_t now = time(0);
	std::cerr << asctime(localtime(&now) ) << "= time before making GridFunction projSol \n";
	GridFunction<BFT, RT> projSol(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HminusHalfSpace), surfaceNormalIndependentFunction(solSphere()));
	now = time(0);
	std::cerr << asctime(localtime(&now) ) << "= time after making GridF projSol \n";
	// Divide by the norms of basis functions
	GridFunction<BFT, RT> normbas(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HminusHalfSpace), surfaceNormalIndependentFunction(CstFct(1.0)));
//	now = time(0);
//	std::cerr << asctime(localtime(&now) ) << "= time after making GridFunction normBas \n";
	
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
//	pt(0) = 0.0;
//	pt(1) = 1.0/std::sqrt(2.0);
//	pt(2) = 1.0/std::sqrt(2.0);
//	SolSphere ssph = solSphere();
//	ssph.evaluate(pt,t);
	// but not clear how to evaluate projSol at pt...
	

	BoundaryOperator<BFT, RT> slpOp = helmholtz3dSingleLayerBoundaryOperator<BFT>(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HplusHalfSpace), make_shared_from_ref(HminusHalfSpace),waveNumber);
	boost::shared_ptr<const Bempp::AbstractBoundaryOperator<double, std::complex<double> > > asdf = slpOp.abstractOperator();
	const GeneralElementarySingularIntegralOperator<BFT,RT,RT> bla = dynamic_cast<const GeneralElementarySingularIntegralOperator<BFT,RT,RT>& > (*asdf);

	arma::Mat<RT> wmDummy(3,3);
	wmDummy.fill(0.);
	std::vector<RT> rhsVeDummy;
	arma::Col<RT> scoDummy;
	boost::shared_ptr<const Bempp::DiscreteBoundaryOperator<RT> > weakCompr = bla.weakFormPeter(str,context,&scoDummy, &rhsVeDummy, &wmDummy);
	
	now = time(0);
//	std::cerr << asctime(localtime(&now) ) << "= time before converting wm to matrix\n";
	arma::Mat<RT> wm = weakCompr->asMatrix(); // Warning: wm now contains one row of the compressed matrix if str starts with t or ..?.., else use 'n'
	now = time(0);
	std::cerr << asctime(localtime(&now) ) << " = time after conv wm\n";

//	DefaultIterativeSolver<BFT, RT> solver(slpOp,ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);
	DefaultIterativeSolver<BFT, RT> solver(weakCompr, str, slpOp,ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);

	tbb::tick_count end = tbb::tick_count::now();
        times(sim,0) = (end - start).seconds(); // Overwrite with solution time added, when str.at(0) == n
	now = time(0);
	std::cerr << asctime(localtime(&now) ) << "\n";

//	solver.saveProjections(normbas,"normbas");	
//	Vector<ResultType> norVect(normbas.projections(slpOp->dualToRange()));
	std::fstream myStream;
	myStream.open("normbas",std::ios::out);
	myStream << normbas.projections(slpOp.dualToRange());
	myStream.close();

	now = time(0);
//	std::cerr << asctime(localtime(&now) ) << ": starting computing and saving projections\n";
	solver.saveProjections(projSol,"projVectOld");
	now = time(0);
//	std::cerr << asctime(localtime(&now) ) << "=time after saving proj\n";
	solver.saveProjections(rhs,"rhs");
//std::cerr << "aposidjf\n";
	std::ifstream inputn("normbas");
	std::vector<RT> normBasVe{std::istream_iterator<RT>(inputn), std::istream_iterator<RT>() };
        inputn.close();
//std::cerr << "aasdff\n";

	std::ifstream inputp("projVectOld");
	std::vector<RT> projVeOld{std::istream_iterator<RT>(inputp), std::istream_iterator<RT>() };
        inputp.close();

//std::cerr << "aposufgha\n";
	std::ifstream inputt("rhs");
	std::vector<RT> rhsVe{std::istream_iterator<RT>(inputt), std::istream_iterator<RT>() };
        inputt.close();
//std::cout << "pasoijfd\n";
	std::cout << "\n norbas(0)=" << normBasVe[0] << "=nb[0],nb[1]=" << normBasVe[1] << "=nb1, pv0=" << projVeOld[0] << "=pv0, pv1=" << projVeOld[1] << "\n";

	// Change projSol to divide by the norms of the basis functions
	arma::Col<RT> correctedProj = arma::Col<RT>(normBasVe.size());
	for(int i=0; i < projVeOld.size(); i++) {
		correctedProj(i) = projVeOld[i]/normBasVe[i];
//		newNor += std::pow(std::abs(solutionCoefficientsNew(i)),2.0);
	}
//	projSol.setCoefficients(correctedProj);	
	projSol.setProjections(make_shared_from_ref(HminusHalfSpace), correctedProj);

	solver.saveProjections(projSol,"projVect");
	std::ifstream inputpn("projVect");
	std::vector<RT> projVe{std::istream_iterator<RT>(inputpn), std::istream_iterator<RT>() };
        inputpn.close();

	std::cout << projVe[0] << "=pv0, pv1=" << projVe[1] << "=pv1, cp0=" << correctedProj(0) << "\n\n";

	int startI = 0;
	int endI = rhsVe.size();
	EvaluationOptions evalOptions = EvaluationOptions();
	Helmholtz3dSingleLayerPotentialOperator<BFT> slPot (waveNumber);

//	if(sim == 0) {
	if ((str.at(0) == 'n') | (str.at(0) == 'f')  | (str.at(0) == 'k') ) {
		solver.initializeSolver(defaultGmresParameterList(1e-5));
		solution = solver.solve(rhs);
		end = tbb::tick_count::now();
	        times(sim,0) = (end - start).seconds();

		std::cout << "ended std calc\n";
		const GridFunction<BFT, RT>& solFunOr = solution.gridFunction();
		arma::Col<RT> solutionCoefficientsOr = solFunOr.coefficients(); // Same as armaSolution from def_iter_solver.cpp
		arma::Mat<RT> potResOr = slPot.evaluateAtPoints(solFunOr, points, quadStrategy, evalOptions);
		errBCavm(sim,0) = arma::mean(arma::mean(abs(potResOr - diri) ))/arma::mean(arma::mean(abs(diri)));
		RT diff = 0.0;
		RT nor = 0.0;
		RT norV = 0.0;
//		RT errFact = 0.0;
		for(int i=0; i < projVe.size(); i++) {
			diff = diff + std::abs(projVe[i]-solutionCoefficientsOr(i))*std::abs(projVe[i]-solutionCoefficientsOr(i)); // sco is arma:col so () iso []
			norV = norV + std::abs(projVe[i])*std::abs(projVe[i]);
			nor = nor + std::abs(solutionCoefficientsOr(i))*std::abs(solutionCoefficientsOr(i));
//			errFact = errFact + std::abs(projVe[i]*factorProj -solutionCoefficientsOr(i))*std::abs(projVe[i]*factorProj -solutionCoefficientsOr(i)); // sco is arma:col so () iso []
		}
		diff = std::sqrt(diff);
		nor = std::sqrt(nor);
		norV = std::sqrt(norV);
//		errProj(sim,0) = std::abs(std::sqrt(errFact)/nor);
		errProj(sim,0) = std::abs(diff/nor);
		std::cout << errBCavm(sim,0) << "=errBCavm, diff =" << diff << "=diff, nor=" << nor << " =norm solCoeffs, norV = " << norV << "=norV\n";

		CT bnorr = 0.0;
		CT berrr = 0.0;
		for (int i=0; i < rhsVe.size(); ++i) {
			RT err = -rhsVe[i];//rhsVee is vector so [] iso ()
			for (int j=0; j < rhsVe.size(); ++j) {
				err += wm(i,j)*solutionCoefficientsOr(j); // sco is arma::col so () iso []
			}
			bnorr += std::pow(std::abs(rhsVe[i]),2.0);
			berrr += std::pow(std::abs(err),2.0);
		}
		errAxb(sim, 0) = sqrt(berrr/bnorr);
		RT corrDiff = 0.0;
		for(int i=0; i < projVe.size(); i++) {
			corrDiff += std::pow(std::abs(projVe[i]*nor/norV -solutionCoefficientsOr(i)), 2); 
		}
		std::cout << std::sqrt(corrDiff) << "=corrDiff, errAxb=" << errAxb(sim,0) << "\n";
	} else {
		std::string::size_type sz;
//		startI = std::stof(str.substr(4),&sz);
		startI = std::stof(str.substr(9),&sz);
		endI = startI+1;
	}
	percs(sim,0) = arma::accu(wm != 0)/(0.0+wm.n_elem);
	CT newNor = 0.0;
//	std::vector<RT> scn = projVe;
//	std::transform(scn.begin(),scn.end(),scn.begin(), std::bind1st(std::multiplies<RT>(),factorProj));
//	arma::Col<RT> solutionCoefficientsNew = arma::vec<RT>(scn);
	arma::Col<RT> solutionCoefficientsNew = arma::Col<RT>(projVe.size());
	for(int i=0; i < projVe.size(); i++) {
		solutionCoefficientsNew(i) = projVe[i]; //projVeOld[i]*factorProj;
		newNor += std::pow(std::abs(solutionCoefficientsNew(i)),2.0);
	}
	
	CT bnorr = 0.0;
	CT berrr = 0.0;
	std::cout << std::sqrt(newNor) << "=newNor, startI=" << startI << "apoapsodjf" << endI << "=endI, mrow=" << wm.n_rows << "=mrows, mcols=" << wm.n_cols << "\n";
	for (int i=startI; i < endI; ++i) {
		RT err = -rhsVe[i];//rhsVee is vector so [] iso ()
		for (int j=0; j < rhsVe.size(); ++j) {
			err += wm(i,j)*solutionCoefficientsNew(j); // sco is arma::col so () iso []
		}
		bnorr += std::pow(std::abs(rhsVe[i]),2.0);
		berrr += std::pow(std::abs(err),2.0);
	}
	errAprojb(sim, 0) = sqrt(berrr/bnorr);
	GridFunction<BFT, RT> solFunNew = solution.gridFunction();
//	const GridFunction<BFT, RT> solFunNew();
	solFunNew.setCoefficients(solutionCoefficientsNew);	
	errBCproj(sim,0) = arma::mean(arma::mean(abs(slPot.evaluateAtPoints(solFunNew, points, quadStrategy, evalOptions) - diri) ))/arma::mean(arma::mean(abs(diri)));
	std::cout << errBCproj(sim,0) << " = err BC when using scaled projection coeffs\n";
	std::cout << percs(sim,0) << "perc, errAprojb=" << errAprojb(sim,0) << "\n";

	std::ofstream myfile;
	myfile.open("res");
	myfile << "Output of tutorial_dirichlet with one row.\n";
	myfile << real(ks) << " = ks " << std::endl;
	myfile << percs << " = perc, errAxb = \n" << errAxb << "=errAxb, errAprojb=" << errAprojb << "\n";
	myfile << times << " = times, errProj = \n" << errProj << "\n";
	myfile << errBCavm << " = errBCavm, errBCproj=" << errBCproj << "\n";
	myfile.close();
}

maxss = std::system("cat res");
//std::cout << real(ks) << " = ks, errAprojb= " << errAprojb << "\n";
//std::cout << percs << " = perc, errAxb = " << errAxb << "\n";
//std::cout << times << " = times, errProj" << errProj << "\n";
//std::cout << errBCavm << "=errBCavm, errBCproj=" << errBCproj << "\n";
}


void fixedWindows() {

// Add code of fixed Windows here
std::cout << "Entered fixedWindows() " <<  "\n";
/*
arma::Mat<RT> ks = arma::exp2(arma::linspace<arma::Mat<RT>>(3,5,3));
//const int kl = ks.size();
const int kl = 1;

arma::Mat<BFT> Ts = arma::linspace<arma::Mat<BFT>>(0.6,0.8,2);
//arma::Mat<BFT> Ts = arma::linspace<arma::Mat<BFT>>(0.03,2.0,10);

const int Tl = Ts.size();

const int avm = 100;
arma::Mat<BFT> thetas = arma::zeros(avm,1);
arma::Mat<BFT> phis = arma::zeros(avm,1);
std::uniform_real_distribution<BFT> unif(0, 1);
std::random_device rd;
std::mt19937 gen(rd());
for(int av = 0; av < avm; av ++) {
//	thetas(av) = M_PI*std::rand();
//	phis(av) = M_PI*2*std::rand();
	thetas(av) = M_PI*unif(gen);
	phis(av) = M_PI*2*unif(gen);
}
arma::Mat<CT> points(3,avm*avm);
arma::Mat<CT> pointsInt(3,avm*avm);
// theta in [0,pi) en phi in [0, 2pi)
for (int thi =0; thi < avm; ++thi) {
	for(int phih =0; phih < avm; ++phih) {
		int idx = thi*avm+phih;
		points(0,idx) = cos(phis(phih))*sin(thetas(thi));
		points(1,idx) = sin(phis(phih))*sin(thetas(thi));
		points(2,idx) = cos(thetas(thi));
//		BFT rtmp = std::rand();
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

//arma::Mat<RT> wm;
//BoundaryOperator<BFT, RT> slpOp;
//DefaultIterativeSolver<BFT, RT> solver;
//Solution<BFT, RT> solution;
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
//	PiecewisePolynomialContinuousScalarSpace<BFT> HplusHalfSpace(grid,2);
//	PiecewiseLinearContinuousScalarSpace<BFT> HminusHalfSpace(grid);
//	PiecewisePolynomialContinuousScalarSpace<BFT> HplusHalfSpace(grid,3);
//	PiecewisePolynomialContinuousScalarSpace<BFT> HminusHalfSpace(grid,2);


	AssemblyOptions assemblyOptions;
	assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW); // Less junk
//	assemblyOptions.setVerbosityLevel(VerbosityLevel::HIGH); // More info (progress % matrix)
	// No ACA (AcaOptions) 

	AccuracyOptions accuracyOptions;
	accuracyOptions.doubleRegular.setRelativeQuadratureOrder(1);
	NumericalQuadratureStrategy<BFT, RT> quadStrategy(accuracyOptions);
	Context<BFT, RT> context(make_shared_from_ref(quadStrategy), assemblyOptions);

//	slpOp = helmholtz3dSingleLayerBoundaryOperator<BFT>(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HplusHalfSpace), make_shared_from_ref(HminusHalfSpace),waveNumber);
	BoundaryOperator<BFT, RT> slpOp = helmholtz3dSingleLayerBoundaryOperator<BFT>(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HplusHalfSpace), make_shared_from_ref(HminusHalfSpace),waveNumber);

//	std::vector< Point3D<CoordinateType> > testPos;
//	testSpace.getGlobalDofPositions(testPos);
//	std::vector< Point3D<CoordinateType> > trialPos;
//	trialSpace.getGlobalDofPositions(trialPos);

//	wm = slpOp.weakForm()->asMatrix();
	arma::Mat<RT> wm = slpOp.weakForm()->asMatrix();
//	std::cout << "Assemble rhs" << std::endl;
	GridFunction<BFT, RT> rhs(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HminusHalfSpace), // is this the right choice?
            surfaceNormalIndependentFunction(MyFunctor()));

	// Initialize the solver
#ifdef WITH_TRILINOS
//	std::cout << "Initialize solver TRILINOS" << std::endl;
	DefaultIterativeSolver<BFT, RT>* solver = new DefaultIterativeSolver<BFT, RT>(slpOp,ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);
	solver->initializeSolver(defaultGmresParameterList(1e-5));
	solver->saveProjections(rhs,"rhsV");
	std::cout << "TRILINOS" << std::endl;
//	solution = solver->solve(rhs);
	Solution<BFT, RT> solution = solver->solve(rhs);
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

//	arma::Col<RT> rhsVe{std::istream_iterator<RT>(input), std::istream_iterator<RT>() };
//std::cout << "apsoijfd\n";
//	std::ifstream input("/home/peter/Desktop/Doctoraat/GreenBempp/simpsonRes/rhsV")
	std::ifstream input("rhsV");
	std::vector<RT> rhsVe{std::istream_iterator<RT>(input), std::istream_iterator<RT>() };
        input.close();
//std::cout << "soghaosiufhd\n";
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
		t.fill(0.);//t = 0;
		tmp.evaluate(pt,t);
		diri(i) = t(0);

		pt(0) = pointsInt(0,i);
		pt(1) = pointsInt(1,i);
		pt(2) = pointsInt(2,i);
		t.fill(0.);
		tmp.evaluate(pt,t);
		dirInt(i) = t(0);
	}
//std::cout << "oiuhaoisdf\n";
	arma::Mat<RT> errBCOr = potResOr - diri;
	errBCavm(ki,0) = mean(mean(abs(errBCOr) ))/mean(mean(abs(diri)));
	errInt(ki,0) = mean(mean(abs(slPot.evaluateAtPoints(solFun, pointsInt, quadStrategy, evalOptions)-dirInt) )); // solFun now has solution of original problem
	errInt(ki,0) = mean(mean(abs(slPot.evaluateAtPoints(solFun, pointsInt, quadStrategy, evalOptions)-dirInt) ))/mean(mean(abs(dirInt))); // solFun now has solution of original problem
	
	arma::Mat<RT> zxcv = slPot.evaluateAtPoints(solFun, pointsInt, quadStrategy, evalOptions);
	std::cout << dirInt(0) << " = dirint, slpint = " << zxcv(0) << "\n";
	std::cout << wm.n_rows << " = wmnr, rhsSiz = " << rhsVe.size() << "\n";


//std::cout << "apsoijfd\n";
	CT bnorr = 0.0;
	BFT l1err = 0.0;
	BFT l1nor = 0.0;
//	for (int i=0; i < rhsVe.size(); ++i) {
	for (int i=0; i < wm.n_rows; ++i) {
	    RT err = -rhsVe[i];
//	    for (int j=0; j < rhsVe.size(); ++j) {
	    for (int j=0; j < wm.n_cols; ++j) {
		err += wm(i,j)*solutionCoefficientsOr(j); // sco is arma::col so () iso []
	    }
	    if(i % (wm.n_rows/10) == 0) {
//		std::cout << err << " = err, bnor = " << bnorr << ", eaxb= " << errAxb(ki,0) << "\n";
	    }
	    bnorr += std::pow(std::abs(rhsVe[i]),2.0);
	    errAxb(ki,0) += std::pow(std::abs(err),2.0);
	    l1err += std::abs(err);
	    l1nor += std::abs(rhsVe[i]);
	}
	errAxb(ki,0) = std::sqrt(errAxb(ki,0)/bnorr);
	std::cout << l1nor << " =l1nor orig, l1err orig= " << l1err << "\n";

//CT qwer = 0.0;
//CT berrr = 0.0;
//(1.4838e-07,2.27629e-07) = err, bnor = 0, berrr = 0
//(-2.3406e-07,-2.94841e-07) = err, bnor = 0.00538152, berrr = 4.30607e-11
//(3.29072e-07,-2.21606e-08) = err, bnor = 0.0108281, berrr = 8.46947e-11
//(1.40028e-07,-3.20568e-07) = err, bnor = 0.0161018, berrr = 1.2864e-10
//(-3.75027e-07,-2.12503e-07) = err, bnor = 0.021414, berrr = 1.72332e-10
//(1.41165e-07,-2.91829e-07) = err, bnor = 0.0267195, berrr = 2.16245e-10
//(-3.66345e-07,1.34379e-07) = err, bnor = 0.0322728, berrr = 2.56041e-10
//(-3.99113e-07,-3.60391e-07) = err, bnor = 0.0378162, berrr = 2.97497e-10
//(-2.11767e-07,-6.79234e-09) = err, bnor = 0.0436223, berrr = 3.44035e-10
//(2.21862e-07,3.64388e-07) = err, bnor = 0.0490523, berrr = 3.83923e-10
//for (int i=0; i < rhsVe.size(); ++i) {
//	RT err = -rhsVe[i];
//	for (int j=0; j < rhsVe.size(); ++j) {
//		err += wm(i,j)*solutionCoefficientsOr(j); // sco is arma::col so () iso []
//	}
//	if(i % (rhsVe.size()/10) == 0) {
//	    std::cout << err << " = err, qwer = " << qwer << ", berrr = " << berrr << "\n";
//	}
//	qwer += std::pow(std::abs(rhsVe[i]),2.0);
//	berrr += std::pow(std::abs(err),2.0);
//}
//std::cout << std::sqrt(berrr) << " =err,L2 nor= " << std::sqrt(qwer) << "\n";



// -------------- Compression -----------------------------------------
   for(int ti = 0; ti < Tl; ti ++) {
	time_t now = time(0);
	std::cerr << "-----------------------------------\n" << ki << " = ki, ti = " << ti << ", perc = " << (0.0+ki*Tl + ti)/(0.0+Tl*kl) << " " << asctime(localtime(&now) ); // << "\n ----------------------------- \n";

//	std::cerr << "-----------------------------------\n" << ki << " = ki, ti = " << ti << ", perc = " << (0.0+ki*Ts.size() + ti)/(0.0+Ts.size()*ks.size()) << " " << asctime(localtime(&now) ) << "\n ----------------------------- \n";

//	std::string str = "f   " + std::to_string(Ts(ti));
	std::string str = "k   " + std::to_string(Ts(ti));
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
	    if(std::abs(err)/std::sqrt(bnorr) > 0.01) {
//		std::cout << err+rhsVe[i] << " =incorrect val row " << i << ", b_i= " << rhsVe[i] << ", relnorerr= " << std::abs(err)/std::sqrt(bnorr) << "\n";
	    }
	    if(i % (rhsVe.size()/10) == 0) {
//		std::cout << err << " = error for row " << i << ", b_i = " << rhsVe[i] << "\n";
		std::cout << err+rhsVe[i] << " =val row " << i << ", b_i= " << rhsVe[i] << ", relnorerr= " << std::abs(err)/std::sqrt(bnorr) << "\n";
	    }
	    l1nor += std::abs(rhsVe[i]);
	    l1err += std::abs(err);
	}
	std::cout << l1nor << " =l1nor, l1err= " << l1err << "\n";
	errAxb(ki,1+ti) = std::sqrt(errAxb(ki,1+ti)/bnorr);
   }
std::ofstream myfile;
myfile.open ("res");
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
*/
}


int main(int argc, char* argv[])
{
fixedWindows();
oneRow();
/*

arma::Mat<RT> ks = arma::exp2(arma::linspace<arma::Mat<RT>>(3,4,2));
const int kl = ks.size();

const int avm = 20; //100;
arma::Mat<BFT> thetas = arma::zeros(avm,1);
arma::Mat<BFT> phis = arma::zeros(avm,1);
std::uniform_real_distribution<BFT> unif(0, 1);
std::random_device rd;
//std::mt19937 gen(rd());
std::mt19937 gen(11); // Constant seed??
//std::default_random_engine gen;
for(int av = 0; av < avm; av ++) {
//	thetas(av) = M_PI*std::rand();
//	phis(av) = M_PI*2*std::rand();
	thetas(av) = M_PI*unif(gen);
	phis(av) = M_PI*2*unif(gen);
}
arma::Mat<CT> points(3,avm*avm);
arma::Mat<CT> pointsInt(3,avm*avm);
arma::Mat<CT> pointsExt(3,avm*avm);
// theta in [0,pi) en phi in [0, 2pi)
for (int thi =0; thi < avm; ++thi) {
	for(int phih =0; phih < avm; ++phih) {
		int idx = thi*avm+phih;
		points(0,idx) = cos(phis(phih))*sin(thetas(thi));
		points(1,idx) = sin(phis(phih))*sin(thetas(thi));
		points(2,idx) = cos(thetas(thi));
//		BFT rtmp = std::rand();
		BFT rtmp = unif(gen);
		pointsInt(0,idx) = rtmp*cos(phis(phih))*sin(thetas(thi));
		pointsInt(1,idx) = rtmp*sin(phis(phih))*sin(thetas(thi));
		pointsInt(2,idx) = rtmp*cos(thetas(thi));
//		pointsExt(0,idx) = (rtmp+11.1)*cos(phis(phih))*sin(thetas(thi));
//		pointsExt(1,idx) = (rtmp+11.1)*sin(phis(phih))*sin(thetas(thi));
//		pointsExt(2,idx) = (rtmp+11.1)*cos(thetas(thi));
		pointsExt(0,idx) = (rtmp+1.1)*cos(phis(phih))*sin(thetas(thi));
		pointsExt(1,idx) = (rtmp+1.1)*sin(phis(phih))*sin(thetas(thi));
		pointsExt(2,idx) = (rtmp+1.1)*cos(thetas(thi));
	}
}
int maxss = 1;
arma::Mat<BFT> conds = arma::zeros(maxss,2);
arma::Mat<BFT> percs = arma::zeros(maxss,1);

arma::Mat<BFT> errBCavm = arma::zeros(maxss,2);
arma::Mat<BFT> errInt = arma::zeros(maxss,2);
arma::Mat<BFT> errExt = arma::zeros(maxss,2);
arma::Mat<BFT> errAxb = arma::zeros(maxss,2);
arma::Mat<BFT> errSol = arma::zeros(maxss,1);

arma::Mat<BFT> times = arma::zeros(maxss,2);

//for(int sim = 0; sim < 3; sim++) {
//for(int sim = 0; sim < 1; sim++) {
for(int sim = 0; sim < maxss; sim++) {

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
//	PiecewisePolynomialContinuousScalarSpace<BFT> HplusHalfSpace(grid,2);
//	PiecewiseLinearContinuousScalarSpace<BFT> HminusHalfSpace(grid);
	AssemblyOptions assemblyOptions;
	assemblyOptions.setVerbosityLevel(VerbosityLevel::LOW);
//	assemblyOptions.setVerbosityLevel(VerbosityLevel::HIGH);
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
//	std::cout << "Initialize solver TRILINOS" << std::endl;
	DefaultIterativeSolver<BFT, RT> solver(slpOp,ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);
	solver.initializeSolver(defaultGmresParameterList(1e-5));
//	std::cout << "TRILINOS" << std::endl;
	Solution<BFT, RT> solution = solver.solve(rhs);
#endif
	tbb::tick_count end = tbb::tick_count::now();
        times(sim,0) = (end - start).seconds();

	std::cout << "ended std calc\n";

// ------------------------ Validation of original -----------------------
	conds(sim,0) = arma::cond(wm);
	const GridFunction<BFT, RT>& solFunOr = solution.gridFunction();
	arma::Col<RT> solutionCoefficientsOr = solFunOr.coefficients(); // Same as armaSolution from def_iter_solver.cpp

	Helmholtz3dSingleLayerPotentialOperator<BFT> slPot (waveNumber);
	EvaluationOptions evalOptions = EvaluationOptions();
	arma::Mat<RT> potResOr = slPot.evaluateAtPoints(solFunOr, points, quadStrategy, evalOptions);
	arma::Mat<RT> potInt = slPot.evaluateAtPoints(solFunOr, pointsInt, quadStrategy, evalOptions);
	arma::Mat<RT> potExt = slPot.evaluateAtPoints(solFunOr, pointsExt, quadStrategy, evalOptions);
	arma::Mat<RT> diri = potResOr;
	arma::Mat<RT> dirInt = potResOr;
	arma::Mat<RT> dirExt = potResOr;
	MyFunctor tmp = MyFunctor();
	solSphere temss = solSphere();
	arma::Col<CT> pt(3);
	for (int i = 0; i < avm*avm; ++i) {
		pt(0) = points(0,i);
		pt(1) = points(1,i);
		pt(2) = points(2,i);
		arma::Col<RT> t(1);
		tmp.evaluate(pt,t);
		diri(i) = t(0);

		pt(0) = pointsInt(0,i);
		pt(1) = pointsInt(1,i);
		pt(2) = pointsInt(2,i);
		t.fill(0.);
		tmp.evaluate(pt,t);
		dirInt(i) = t(0);

		pt(0) = pointsExt(0,i);
		pt(1) = pointsExt(1,i);
		pt(2) = pointsExt(2,i);
		t.fill(0.);
		temss.evaluateF(pt,t);
		dirExt(i) = t(0);
		if (i % (avm*avm/10) == 0) {
//		    std::cout << i << " =i, pt= " << pt << "\n";
//		    std::cout << dirInt(i) << " =dint, potint = " << potInt(i) << " " << dirInt(i,0) << " =dint, potint = " << potInt(i,0) << "\n";
//		    std::cout << dirInt(i) << " =dint, potint = " << potInt(i) << dirInt(0,i) << " =dint, potint = " << potInt(0,i) << "\n";
//		    std::cout << dirExt(i) << " =dExt, potExt = " << potExt(i) <<"\n"; //<< dirExt(0,i) << " =dExt, potExt = " << potExt(0,i) << "\n";
		}
	}
	errBCavm(sim,0) = mean(mean(abs(potResOr - diri) ))/mean(mean(abs(diri)));
	errInt(sim,0) = mean(mean(abs(slPot.evaluateAtPoints(solFunOr, pointsInt, quadStrategy, evalOptions) - dirInt) ))/mean(mean(abs(dirInt)));
	errExt(sim,0) = mean(mean(abs(slPot.evaluateAtPoints(solFunOr, pointsExt, quadStrategy, evalOptions) - dirExt) ))/mean(mean(abs(dirExt)));
//	errExt(sim,0) = mean(mean(abs(slPot.evaluateAtPoints(solFunOr, pointsExt, quadStrategy, evalOptions) + dirExt) ))/mean(mean(abs(dirExt)));
	pt(0) = pointsExt(0,0);
	pt(1) = pointsExt(1,0);
	pt(2) = pointsExt(2,0);
	std::cout << potExt(0,0) << " =potExt,dirExt = " << dirExt(0,0) << ", mean = " << mean(mean(abs(dirExt))) << ", pt = " << pt << "\n";
//	std::cout << potExt(1,0) << " =potExt,dirExt = " << dirExt(1,0) << "\n"; Gives index out of bounds because points*(3,avm*avm) 
	std::cout << potExt(0,1) << " =potExt,dirExt = " << dirExt(0,1) << "\n";



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

//	std::ifstream input("rhs");
//	std::vector<RT> rhsVe{std::istream_iterator<RT>(input), std::istream_iterator<RT>() };
//        input.close();

	std::ifstream inputt("rhs");
	std::vector<RT> rhsVe{std::istream_iterator<RT>(inputt), std::istream_iterator<RT>() };
        inputt.close();
CT bnorr = 0.0;
CT berrr = 0.0;
for (int i=0; i < rhsVe.size(); ++i) {
	RT err = -rhsVe[i];//rhsVee is vector so [] iso ()
//	RT err = -rhsVee(i);
	for (int j=0; j < rhsVe.size(); ++j) {
//		err += wm[i,j]*solutionCoefficientsOr[j];
		err += wm(i,j)*solutionCoefficientsOr(j); // sco is arma::col so () iso []
	}
	if(i % (rhsVe.size()/10) == 0) {
//	    std::cout << err << " = err, bnor = " << bnorr << ", berrr = " << berrr << ", b_i = " << rhsVe[i] << "\n";
	}
	bnorr += std::pow(std::abs(rhsVe[i]),2.0);
//	bnorr += std::pow(std::abs(rhsVee(i)),2.0);
	berrr += std::pow(std::abs(err),2.0);
}
//std::cout << "\n\n" << std::sqrt(berrr) << " =err,L2 nor= " << std::sqrt(bnorr) << "\n";
	errAxb(sim, 0) = sqrt(berrr/bnorr);

//	std::vector<RT> diff = projVe-rhsVe;
	RT diff = 0.0;
	RT nor = 0.0;
	RT norV = 0.0;
	for(int i=0; i < projVe.size(); i++) {
//		diff = diff + (projVe[i]-rhsVe[i])^2.0;
//		diff = diff + std::abs(projVe[i]-rhsVe[i])*std::abs(projVe[i]-rhsVe[i]);
//		diff = diff + std::abs(projVe[i]-solutionCoefficientsOr[i])*std::abs(projVe[i]-solutionCoefficientsOr[i]);
//		norV = norV + std::abs(projVe[i])*std::abs(projVe[i]);
//		nor = nor + std::abs(solutionCoefficientsOr[i])*std::abs(solutionCoefficientsOr[i]);
		diff = diff + std::abs(projVe[i]-solutionCoefficientsOr(i))*std::abs(projVe[i]-solutionCoefficientsOr(i)); // sco is arma:col so () iso []
		norV = norV + std::abs(projVe[i])*std::abs(projVe[i]);
		nor = nor + std::abs(solutionCoefficientsOr(i))*std::abs(solutionCoefficientsOr(i));
//		diff = diff + std::abs(projVe(i)-solutionCoefficientsOr(i))*std::abs(projVe(i)-solutionCoefficientsOr(i));
//		norV = norV + std::abs(projVe(i))*std::abs(projVe(i));
//		nor = nor + std::abs(solutionCoefficientsOr(i))*std::abs(solutionCoefficientsOr(i));
//		std::cout << projVe[i] << " ";
	}
//	std::cout << "\n";
	diff = std::sqrt(diff);
	nor = std::sqrt(nor);
	norV = std::sqrt(norV);
std::cout << projVe[0] << "asdf" << solutionCoefficientsOr(0) << "oiuh" << rhsVe[0] << "asdf" << projVe.size() << "\n";
//std::cout << projVe(0) << "asdf" << solutionCoefficientsOr(0) << "oiuh" << rhsVe(0) << "asdf" << projVe.size() << "\n";
	std::cout << errBCavm(sim,0) << "stopping" << diff << "apsoijfd" << nor << " , norV = " << norV << "\n";
	
	RT corrDiff = 0.0;
	for(int i=0; i < projVe.size(); i++) {
		corrDiff += std::pow(std::abs(projVe[i]*nor/norV -solutionCoefficientsOr(i)), 2); 
	}
	std::cout << jn(0,0.9) << " asodufh " << jn(3,9.7) << " paosjfhd " << sqrt(acos(-1.0)/2/std::abs(waveNumber) ) << " cordiff = " << sqrt(corrDiff) << "\n";

	arma::Col<RT> solutionCoefficientsNew = solutionCoefficientsOr;	
	for(int i=0; i < projVe.size(); i++) {
//		solutionCoefficientsNew(i) = projVe[i];
		solutionCoefficientsNew(i) = projVe[i]*nor/norV;
	}
	GridFunction<BFT, RT> solFunNew = solution.gridFunction();
//	solFunOr.setCoefficients(solutionCoefficientsNew);
	solFunNew.setCoefficients(solutionCoefficientsNew);
	
	std::cout << mean(mean(abs(slPot.evaluateAtPoints(solFunNew, points, quadStrategy, evalOptions) - diri) ))/mean(mean(abs(diri))) << " = err BC when using projection coeffs\n";

//return 0;
//break;
//continue;


// -------------- Compression -----------------------------------------
//	std::string str = "c   0.8"; //"f   "+ std::to_string(Ts(ti)); //"c   ";
//	std::string str = "t   0"; //"t   2010";
//	std::string str = "d   0.1";
//	std::string str = "k     0.6";
	std::string str = "k  1.6";
//	std::string str = "f  1.6";
//	std::string str = "i  1.6";
//	std::string str = "j  1.6";
	start = tbb::tick_count::now();

	slpOp = helmholtz3dSingleLayerBoundaryOperator<BFT>(make_shared_from_ref(context), make_shared_from_ref(HminusHalfSpace), make_shared_from_ref(HplusHalfSpace), make_shared_from_ref(HminusHalfSpace),waveNumber);
	boost::shared_ptr<const Bempp::AbstractBoundaryOperator<double, std::complex<double> > > asdf = slpOp.abstractOperator();
	const GeneralElementarySingularIntegralOperator<BFT,RT,RT> bla = dynamic_cast<const GeneralElementarySingularIntegralOperator<BFT,RT,RT>& > (*asdf);

//	arma::Mat<RT> wmDummy(3,3);
//	wmDummy.fill(0.);
//	std::vector<RT> rhsVeDummy;
	boost::shared_ptr<const Bempp::DiscreteBoundaryOperator<RT> > weakCompr = bla.weakFormPeter(str,context,&solutionCoefficientsOr, &rhsVe, &wm);
	
	wm = weakCompr->asMatrix(); // Warning: wm now contains compressed matrix
//std::cout << "apsoijfd\n";

if (str.at(0) == 't') {
arma::Mat<RT> oneRow = weakCompr->asMatrix();
size_t testRow = 0;
std::string::size_type sz;
testRow = std::stof(str.substr(2),&sz);
RT tmpErr = 0.0;
RT globalNor = 0.0;
RT aprxb = 0.0;
//	for (int i=0; i < rhsVee.size(); ++i) { // oneRow.n_cols
	for (int i=0; i < oneRow.size(); ++i) {
		aprxb += oneRow(i,0)*projVe[i]; // projVe is vector so square brackets
//		aprxb += oneRow[i,0]*projVe[i];
//		tmpErr += std::pow(std::abs(wm[i,testRow]-oneRow[i,0]),2.0);
		tmpErr += std::pow(std::abs(wm(i,testRow)-oneRow(i,0)),2.0);
//		globalNor += std::pow(std::abs(wm[i,testRow]),2.0);
		globalNor += std::pow(std::abs(wm(i,testRow)),2.0);
if( (i == 0) || (i == oneRow.size() ) ) {
	std::cout << wm[i,testRow] << " =matr, row= " << oneRow[i,0] << " asdf " <<i << " oi "  << " rowtrans= " << oneRow[0,i] << " rowidx " << oneRow[i,testRow] << " rowtraidx " << oneRow[testRow,i] << "\n"; 
	std::cout << wm(i,testRow) << " =matr, row= " << oneRow(i,0) << " asdf " <<i << " oi "  << " rowtrans= " << oneRow(0,i) << " rowidx " << oneRow(i,testRow) << " rowtraidx " << oneRow(testRow,i) << "\n"; 

}
	}
std::cout << std::sqrt(tmpErr) << " =err first row, nor= "  << std::sqrt(globalNor) <<" , testRow= " << testRow << " aspojfd " << std::sqrt(tmpErr)/std::sqrt(globalNor) << "\n\n";

std::cout << oneRow.size() << " aos " << wm.size() << " s " << tmpErr << " e " << globalNor << " t " << tmpErr/globalNor << "\n\n";

std::cout << aprxb << " = aprxb, b = " << rhsVe[testRow] << " , testRow = " << testRow << "\n";
}


CT bnor = 0.0;
CT berr = 0.0;
for (int i=0; i < wm.n_rows; ++i) {
//	RT err = rhsVe[i];
	RT err = -rhsVe[i]; //rhsVe is vector so [] iso ()
//	RT err = -rhsVe(i);
//	for (int j=0; j < oneRow.size(); ++j) {
	for (int j=0; j < wm.n_cols; ++j) {
//		err += wm[i,j]*projVe[j];
//		err += wm[j,i]*solutionCoefficientsOr[j];
//		err += wm[i,j]*solutionCoefficientsOr[j];
		err += wm(i,j)*solutionCoefficientsOr(i); //sco is arma:col so () iso []
//		err += wm[j,i]*projVe[j];
	}
	if(i % (wm.n_rows/10) == 0) {
//	    std::cout << err << " = err, bnor = " << bnor << ", b_i = " << rhsVe[i] << "\n";
	}
//	std::cout << " " << err;
	bnor += std::pow(std::abs(rhsVe[i]),2.0);
//	bnor += std::pow(std::abs(rhsVe(i)),2.0);
	berr += std::pow(std::abs(err),2.0);
}
//std::cout << std::sqrt(berr) << " =err,L2 nor= " << std::sqrt(bnor) << "\n";
	errAxb(sim,1) = sqrt(berr/bnor);
//break;

	DefaultIterativeSolver<BFT, RT> solverCompr(weakCompr, str, slpOp, ConvergenceTestMode::TEST_CONVERGENCE_IN_DUAL_TO_RANGE);
	solverCompr.initializeSolver(defaultGmresParameterList(1e-5));
	Solution<BFT, RT> solutionCompr = solverCompr.solve(rhs);

	const GridFunction<BFT, RT>& solFunCompr = solutionCompr.gridFunction();
	arma::Col<RT> solutionCoefficientsCompr = solFunCompr.coefficients();

	end = tbb::tick_count::now();
        times(sim,1) = (end - start).seconds();

// ----------------   Validation of compressed  ------------------------------------
//	std::cout << "valcompr\n";
	conds(sim,1) = arma::cond(wm);
	percs(sim,0) = arma::accu(wm != 0)/(0.0+wm.n_elem);
	errSol(sim,0) = arma::norm(solutionCoefficientsCompr -solutionCoefficientsOr)/arma::norm(solutionCoefficientsOr);
	arma::Mat<RT> potResCompr = slPot.evaluateAtPoints(solFunCompr, points, quadStrategy, evalOptions);
	errBCavm(sim,1) = mean(mean(abs(potResCompr - diri) ))/mean(mean(abs(diri)));
	errInt(sim,1) = mean(mean(abs(slPot.evaluateAtPoints(solFunCompr, pointsInt, quadStrategy, evalOptions) - dirInt) ))/mean(mean(abs(dirInt)));
	errExt(sim,1) = mean(mean(abs(slPot.evaluateAtPoints(solFunCompr, pointsExt, quadStrategy, evalOptions) - dirExt) ))/mean(mean(abs(dirExt)));
//	errExt(sim,1) = mean(mean(abs(slPot.evaluateAtPoints(solFunCompr, pointsExt, quadStrategy, evalOptions) + dirExt) ))/mean(mean(abs(dirExt)));

	std::ofstream myfile;
	myfile.open ("res");
	myfile << "Output of tutorial_dirichlet.\n";
	myfile << real(ks) << " = ks " << std::endl;
	myfile << percs << " = perc, conds = " << std::endl << conds << std::endl;
	myfile << times << " = times" << std::endl << std::endl;
	myfile << errBCavm << " = errBCavm, errSol = " << std::endl << errSol << std::endl;
	myfile << errInt << " = errInt, errSol = \n" << errAxb << "\nerrExt=\n" << errExt <<"\n";
	myfile.close();
}


std::cout << real(ks) << " = ks " << std::endl;
std::cout << percs << " = perc, conds = " << conds << std::endl;
std::cout << times << " = times, " << errBCavm << " = errBCavm, errSol = " << errSol << std::endl;
std::cout << errInt << " = errInt, errAxb = " << errAxb << "errExt=" << errExt <<"\n";
*/
//correlations();
}
//void correlations() {
// add correlation code here later
//}

