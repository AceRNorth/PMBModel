#include <iostream>
#include <string>
#include <filesystem>
#include "Simulation.h"
#include "inputval.h"
#include "constants.h"
#include "Params.h"
#include "globals.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

	int start_t = 0;

int main() {
    //std::string params_filename; 
    std::cout << "hello" <<std::endl;

    char boundary_type;
    char disp_type;
    std::string coords_filename; 
    std::string rainfall_filename; 
    std::string rel_times_filename;
    std::string mu_a_filename;

    // input parameters
	// progression parameters
	int num_runs;
	int max_t; 

	// area parameters
	int num_pat;  
	
	// model parameters
	double mu_j;
	double mu_a;
	double beta;
	double theta;
	double comp_power;
	int min_dev;

	// gene drive inheritance parameters
	double bias;
	double omega_M;
	double omega_F;

	// gene drive release parameters
	int driver_start;
	int num_driver_M;
	int num_driver_sites; 

	// dispersal parameters 
	double disp_rate; 
	double max_disp; 

	// aestivation parameters
	double psi;
	double mu_aes;
	int t_hide1;
	int t_hide2;
	int t_wake1; 
	int t_wake2;

	// seasonality parameters
	double alpha0_mean;
	double alpha0_variance;
	double alpha1;
	double amp;
	double resp;

	// data-recording parameters
	int rec_start; 
	int rec_end;
	int rec_interval_global;
	int rec_interval_local;
	int rec_sites_freq; 
	int set_label;
    std::cout << "hello2" <<std::endl;

    std::cin >> num_runs;
    std::cin >> max_t;
    std::cin >> num_pat;
    std::cin >> mu_j;
    std::cin >> mu_a;
    std::cin >> beta;
    std::cin >> theta;
    std::cin >> comp_power;
    std::cin >> min_dev;
    std::cin >> bias;
    std::cin >> omega_M;
    std::cin >> omega_F;
    std::cin >> driver_start;
    std::cin >> num_driver_M;
    std::cin >> num_driver_sites;
    std::cin >> disp_rate;
    std::cin >> max_disp;
    std::cin >> psi;
    std::cin >> mu_aes;
    std::cin >> t_hide1;
    std::cin >> t_hide2;
    std::cin >> t_wake1;
    std::cin >> t_wake2;
    std::cin >> alpha0_mean;
    std::cin >> alpha0_variance;
    std::cin >> alpha1;
    std::cin >> amp;
    std::cin >> resp;
    std::cin >> rec_start;
    std::cin >> rec_end;
    std::cin >> rec_interval_global;
    std::cin >> rec_interval_local;
    std::cin >> rec_sites_freq;
    std::cin >> set_label;
    std::cin >> boundary_type;
    std::cin >> disp_type;
    std::cin >> start_t;
    std::cin >> coords_filename;
    std::cin >> rainfall_filename; 
    std::cin >> rel_times_filename;
    std::cin >> mu_a_filename;
    std::cin >> OutputType;

    std::cout << "hello3" <<"   "<<OutputType<<std::endl;

    InputParams params;
    params.num_runs = num_runs;
    params.max_t = max_t;
    params.num_pat = num_pat;
    params.mu_j = mu_j;
    params.mu_a = mu_a;
    params.beta = beta;
    params.theta = theta;
    params.comp_power = comp_power;
    params.min_dev = min_dev;
    params.bias = bias;
    params.omega_M = omega_M;
    params.omega_F = omega_F;
    params.driver_start = driver_start;
    params.num_driver_M = num_driver_M;
    params.num_driver_sites = num_driver_sites;
    params.disp_rate = disp_rate;
    params.max_disp = max_disp;
    params.psi = psi;
    params.mu_aes = mu_aes;
    params.t_hide1 = t_hide1;
    params.t_hide2 = t_hide2;
    params.t_wake1 = t_wake1;
    params.t_wake2 = t_wake2;
    params.alpha0_mean = alpha0_mean;
    params.alpha0_variance = alpha0_variance;
    params.alpha1 = alpha1;
    params.amp = amp;
    params.resp = resp;
    params.rec_start = rec_start;
    params.rec_end = rec_end;
    params.rec_interval_global = rec_interval_global;
    params.rec_interval_local = rec_interval_local;
    params.rec_sites_freq = rec_sites_freq;
    params.set_label = set_label;
   // num_rel_sites=num_driver_sites;
    std::cout << "hello4" <<"  "<<set_label<<std::endl;

    BoundaryType boundary;
    if (boundary_type == 't') {
        boundary = Toroid;
    }
    else if (boundary_type == 'e') {
        boundary = Edge;
    }

    DispersalType disp;
    if (disp_type == 'd') {
        disp = DistanceKernel;
    }
    else if (disp_type == 'r') {
        disp = Radial;
    }

    std::cout << "hello5" <<"  "<<set_label<<std::endl;
    Simulation simulation(params);
    std::cout << "hello6" <<"  "<<set_label<<std::endl;
    InheritanceParams inher;
    inher.bias = params.bias;
    inher.omega_M = params.omega_M;
    inher.omega_F = params.omega_F;
    simulation.set_inheritance(inher);
    simulation.set_boundary_type(boundary);
    simulation.set_dispersal_type(disp);
    std::cout << "hello7" <<"  "<<set_label<<std::endl;

    if (coords_filename != "none") {
        auto coords_filepath = std::filesystem::path(coords_filename);
        simulation.set_coords(coords_filepath);
    }
    if (rainfall_filename != "none") {
        auto rainfall_filepath = std::filesystem::path(rainfall_filename);
        simulation.set_rainfall(rainfall_filepath);
    }
    if (rel_times_filename != "none") {
        auto rel_times_filepath = std::filesystem::path(rel_times_filename);
        simulation.set_release_times(rel_times_filepath);
    }


	/*-----------------input file with adult female mortality data--------------------------*/
    /*-------------this is for running CRT against Imperial malaria model with variable female mortality input---------*/
    
    if (mu_a_filename != "none") {
		std::ostringstream ggg;
		ggg.str(mu_a_filename);
		std::ifstream mortdata(ggg.str().c_str()); 
		double rr;
		std::string line ="";
		std::string cell;
		while(std::getline(mortdata, line ))
		{       
			std::stringstream lineStream(line);
			std::getline(lineStream,cell,',');
			rr=std::strtod(cell.c_str(),NULL);
			mu_a_list.push_back(rr);
		};
	    }
	
    /*--------------------------------------------------------------------------------------*/

    simulation.run_reps();

    return 0;
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
