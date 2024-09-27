#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <omp.h>
#include <fftw3.h>
#include "fielddata.hpp"

int main(int argc, char **argv){
    std::string FileName = argv[1];

    int ng(2000);
    double boxsize(4000.);

    fftwf_init_threads();
    fftwf_plan_with_nthreads(omp_get_max_threads());
    std::cout << omp_get_max_threads() << " threads will be used in FFT." << std::endl;
    //fftwf_import_wisdom_from_filename(wisdom_file.c_str());

    FieldData overdensity(ng,boxsize,false);
    overdensity.load(FileName);
    //fftwf_export_wisdom_to_filename(wisdom_file.c_str());
    std::cout << "File loaded successfully." << std::endl;

    exit(0);
}
