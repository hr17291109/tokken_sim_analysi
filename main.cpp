#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <omp.h>
#include <fftw3.h>
#include "binneddata.hpp"
#include "fielddata.hpp"
#include "fileloader.hpp"
#include "cosmo.hpp"
#include "util.hpp"
#include "binning.hpp"
#include "param.hpp"
#include "filenames.hpp"

int main(int argc, char **argv){
    std::string FileBase = argv[1];
    int snapnum(atoi(argv[2]));

    // Cosmological parameters
    double Omega_cb; // This is Omega_c + Omega_b, not Omega_m
    double Omega_m;  // Omega_m (=Omega_cb + Omega_{massive-neutrino})
    double h0;       // Dimensionless hubble parameter (should be around 0.7)
    double As;       // Amplitude of the primordial fluctuations
    double ns;       // Spectral index
    double k0;       // The "pivot" wavenumber at which As is given

    // Simulation parameters
    double Box;      // in Mpc/h
    double redshift; // The output redshift
    int ng(2000); // Number of grid points per dim for FFT

    switch(snapnum){
        case 0:
            redshift = 0.61;
            break;
        case 1:
            redshift = 0.51;
            break;
        case 2:
            redshift = 0.38;
            break;
        default:
            redshift = 0;
            break;
    }
    std::cerr << "Redshift: " << redshift << std::endl;
    // paramfile1 = FileBase + "/params/output_sample.txt"; // parameter file for the initial-condition generator
    // paramfile2 = FileBase + "/params/input_param.ini"; // parameter file for CLASS Boltzmann solver

    param::parameter p1(FileBase+"/"+params_dir+"/"+nugenic_param_file);
    Omega_cb = p1.get<double>("Omega");
    Box = p1.get<double>("Box");
    // Npart1d = (long long int)p1.get<int>("Npart");
    As = p1.get<double>("As1");
    ns = p1.get<double>("ns1");
    h0 = p1.get<double>("HubbleParam");
    k0 = p1.get<double>("kpivot")/h0;

    std::cout << "Info from " << FileBase+"/"+params_dir+"/"+nugenic_param_file << std::endl;
    std::cout << "Box: " << Box << std::endl;
    // std::cerr << "Npart1d: " << Npart1d << std::endl;
    std::cout << "Omega_cb: " << Omega_cb << std::endl;
    std::cout << "As: " << As << std::endl;
    std::cout << "ns: " << ns << std::endl;
    std::cout << "h: " << h0 << std::endl;
    std::cout << "kpivot: " << k0 << std::endl;

    param::parameter p2(FileBase+"/"+params_dir+"/"+class_param_file);

    Omega_m = p2.get<double>("Omega_m");

    std::cout << "Info from " << FileBase+"/"+params_dir+"/"+class_param_file << std::endl;
    std::cout << "Omega_m: " << Omega_m << std::endl << std::endl;

    fftwf_init_threads();
    fftwf_plan_with_nthreads(omp_get_max_threads());
    std::cout << omp_get_max_threads() << " threads will be used in FFT." << std::endl;
    fftwf_import_wisdom_from_filename(wisdom_file.c_str());

    FieldData Df_lin(ng,Box,true);
    fftwf_export_wisdom_to_filename(wisdom_file.c_str());
    Df_lin.load(FileBase+"/"+white_dir+"/"+white_file);
    exit(0);
}