// gedit /opt/fb/bempp/lib/fiber/modified_helmholtz_3d_single_layer_potential_kernel_functor.hpp
// gedit /opt/fb/bempp/lib/assembly/dense_global_assembler.hpp
// gedit /opt/fb/bempp/examples/cpp/oneRow.cpp

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

#include "oneRow.cpp"
#include "fixedWindows.cpp"

int main(int argc, char* argv[])
{

if(argc == 1) {
	fixedWindows();
	oneRow(7);
}
if (argv[1] == "help") {
	return std::system("cat README");
} else if (argv[1] == "oneRow") {
	oneRow(std::atoi(argv[2]) );
}
return 0;

}


