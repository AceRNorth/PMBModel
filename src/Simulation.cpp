#include <array>
#include <vector>
#include <iostream> // for error messages
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include "Simulation.h"
#include "constants.h"
#include "Model.h"
#include "Record.h"
#include "inputval.h"
#include "globals.h"

#include <random>
#include <algorithm>  // Required for std::shuffle

/**
 * Simulation constructor.
 * @param[in] input simulation parameters
 * @see InputParams
 */
Simulation::Simulation(InputParams input)
{ 
	num_runs = input.num_runs;
	max_t = input.max_t;

	model_params = new ModelParams;
	model_params->area = new AreaParams;
	model_params->area->num_pat = input.num_pat;
	model_params->life = new LifeParams;
	model_params->life->mu_j = input.mu_j;
	model_params->life->mu_a = input.mu_a;
	model_params->life->beta = input.beta;
	model_params->life->theta = input.theta;
	model_params->life->comp_power = input.comp_power;
	model_params->life->min_dev = input.min_dev;
	model_params->rel = new ReleaseParams;
	model_params->rel->release_times.push_back(input.driver_start);
	model_params->rel->num_driver_M = input.num_driver_M;
	model_params->rel->num_driver_sites = input.num_driver_sites;
	model_params->disp = new DispersalParams;
	model_params->disp->disp_rate = input.disp_rate;
	model_params->disp->max_disp = input.max_disp;
	model_params->aes = new AestivationParams;
	model_params->aes->psi = input.psi;
	model_params->aes->mu_aes = input.mu_aes;
	model_params->aes->t_hide1 = input.t_hide1;
	model_params->aes->t_hide2 = input.t_hide2;
	model_params->aes->t_wake1 = input.t_wake1;
	model_params->aes->t_wake2 = input.t_wake2;
	model_params->initial = new InitialPopsParams;
	model_params->initial->initial_WJ = input.initial_WJ;
	model_params->initial->initial_WM = input.initial_WM;
	model_params->initial->initial_WV = input.initial_WV;
	model_params->initial->initial_WF = input.initial_WF;

	sine_rainfall_params = new SineRainfallParams;
	sine_rainfall_params->alpha1 = input.alpha1;
	sine_rainfall_params->amp = input.amp;
	input_rainfall_params = new InputRainfallParams;
	input_rainfall_params->alpha1 = input.alpha1;
	input_rainfall_params->resp = input.resp;
	(input_rainfall_params->rainfall).clear();
	alpha0_mean = input.alpha0_mean;
	alpha0_variance = input.alpha0_variance;

	rec_params = new RecordParams;
	rec_params->rec_start = input.rec_start;
	rec_params->rec_end = input.rec_end;
	rec_params->rec_interval_global = input.rec_interval_global;
	rec_params->rec_interval_local = input.rec_interval_local;
	rec_params->rec_sites_freq = input.rec_sites_freq;
	rec_params->set_label = input.set_label;
	
	for (int i=0; i < constants::num_gen; ++i) {
		for (int j=0; j < constants::num_gen; ++j) {
			for (int k=0; k < 2*constants::num_gen; ++k) {
				inher_fraction[i][j][k] = 0;
			}
		}
	}
	InheritanceParams mendelian = {0, 0, 0};
	set_inheritance(mendelian);

	sites_coords.clear();
	humans.clear();
	release_sites.clear();
	boundary_type = BoundaryType::Toroid;
	disp_type = DispersalType::DistanceKernel;
}

/**
 * Simulation destructor.
 */
Simulation::~Simulation() 
{
	delete model_params->area;
	delete model_params->life;
	delete model_params->rel;
	delete model_params->disp;
	delete model_params->aes;
	delete model_params->initial;
	delete model_params;
	delete sine_rainfall_params;
	delete input_rainfall_params;
	delete rec_params;
}

/**
 * Sets the sites' coordinates and the release sites from a file, unless any errors are thrown.
 * The file should be structured in three columns, with x and y coordinates and whether the site is a release site (y/n) respectively. Values should be delimited by white space, and site rows should be delimited by new lines. 
 * @warning Setting specific release sites will change the mode of gene drive release from random to scheduled. The specified release sites will be used for each of the releases (if using multiple release times). Mosquitoes will be released at all release sites for each release time.
 * @param[in] filepath coordinates filepath, can be a relative or absolute filepath, or a filename (only if the file is in the build directory)
 */
void Simulation::set_coords(const std::filesystem::path& filepath) 
{
	sites_coords.clear();
	humans.clear();
	release_sites.clear();

	if (!std::filesystem::exists(filepath) || !std::filesystem::is_regular_file(filepath)) {
		std::cerr << "Invalid filename. To enter a filename, the file should be in the build directory. Otherwise, the filepath should be provided (either relative to 'build' or absolute)." << std::endl;
	}
	else {
		std::ifstream file(filepath);
		std::string line;
		std::vector<Point> temp_coords;
		std::vector<int> temp_rel_sites;
		if (file.is_open()) {
			for(int i=0; std::getline(file, line); ++i) {
				std::stringstream linestream(line);
				if (line.size() == 0) break;

				double x, y,bd,rand_effect;
				int type,ass,index,relnum,tert;
				char is_rel_site;
				int err = 0;
				if (!read_and_validate_type(linestream, x, "x" + std::to_string(i+1), "double")) err++;
				if (!read_and_validate_type(linestream, y, "y" + std::to_string(i+1), "double")) err++;
				if (!read_and_validate_type(linestream, bd, "build dens" + std::to_string(i+1), "double")) err++;
				if (!read_and_validate_type(linestream, type, "patch type" + std::to_string(i+1), "int")) err++;
				if (!read_and_validate_type(linestream, rand_effect, "random effect" + std::to_string(i+1), "double")) err++;
				if (!read_and_validate_type(linestream, tert, "prevalence tertile" + std::to_string(i+1), "int")) err++;
				
				if (err == 0) {
					temp_coords.push_back({x, y});
					humans.push_back(bd*rand_effect);
					patch_type.push_back(type);
					patch_tertile.push_back(tert);
					if (type == 5)  release_sites.push_back(i); 
				}
			}
		}
		file.close();
		sites_coords = temp_coords;
	}
		std::cerr<<" coords length  "<<sites_coords.size()<<std::endl;
		std::cerr<<" num release sites  "<<release_sites.size()<<std::endl;

	//	for(int i=0;i<release_sites.size();i++)std::cout<<release_sites[i]<<std::endl;
}

/** 
 * Sets the boundary type for the model.
 * @param[in] boundary boundary type
 * @see BoundaryType, BoundaryStrategy
 */
void Simulation::set_boundary_type(BoundaryType boundary) 
{
	boundary_type = boundary;
}

/**
 * Sets the dispersal type for the model.
 * @param[in] disp dispersal type
 * @see DispersalType, Dispersal
 */
void Simulation::set_dispersal_type(DispersalType disp)
{
	disp_type = disp;
}

/**
 * Sets the daily rainfall values from a file, unless any errors are thrown.
 * These can be daily values for a year cycle (365 days), or daily values for all the simulated days (max_t). 
 * The values should be delimited by new lines.
 * @note This option will use the input parameters alpha1 and resp previously provided and ignore amp.
 * @param[in] filepath rainfall filename, can be a relative or absolute filepath, or a filename (only if the file is in the build directory)
 * @see Seasonality
 */
void Simulation::set_rainfall(const std::filesystem::path& filepath)
{
	input_rainfall_params->rainfall.clear();

	if (!std::filesystem::exists(filepath) || !std::filesystem::is_regular_file(filepath)) {
		std::cerr << "Invalid filename. To enter a filename, the file should be in the build directory. Otherwise, the filepath should be provided (either relative to 'build' or absolute)." << std::endl;
	}
	else {
		std::ifstream file(filepath);
		std::string line;
		std::vector<double> temp;
		if (file.is_open()) {
			for(int i=0; std::getline(file, line); ++i) {
				std::stringstream linestream(line);
				if (line.size() == 0) break;

				double r_d;
				int err = 0;
				if (!read_and_validate_type(linestream, r_d, "rainfall_day" + std::to_string(i+1), "double")) err++;
				if (!check_bounds("rainfall_day" + std::to_string(i+1), r_d, 0.0, true)) err++;

				if (err == 0) {
					temp.push_back(r_d);
				}
			}
		}
		file.close();

		if (temp.size() == 365 || temp.size() == max_t) {
			input_rainfall_params->rainfall = temp;
		}
		else {
			std::cerr << "Error: the number of valid daily rainfall values in the file is not 365 or max_t." << std::endl;
		}	
	}
}

/**
 * Sets the release times for gene drive release from a file, unless any errors are thrown.
 * Release times should be simulation day numbers within the maximum time of simulation.
 * Values should be delimited by new lines.
 * @note When using this option, the input parameter driver_start previously entered will be ignored.
 * @param[in] filepath release times filename, can be a relative or absolute filepath, or a filename (only if the file is in the build directory)
 * @see GDRelease
 */
void Simulation::set_release_times(const std::filesystem::path& filepath) 
{
	if (!std::filesystem::exists(filepath) || !std::filesystem::is_regular_file(filepath)) {
		std::cerr << "Invalid filename. To enter a filename, the file should be in the build directory. Otherwise, the filepath should be provided (either relative to 'build' or absolute)." << std::endl;
	}

	else {
		std::ifstream file(filepath);
		std::string line;
		std::vector<int> temp;
		int tot_err = 0;
		if (file.is_open()) {
			for(int i=0; std::getline(file, line); ++i) {
				std::stringstream linestream(line);
				if (line.size() == 0) break;

				int r_d;
				int err = 0;
				if (!read_and_validate_type(linestream, r_d, "release_day" + std::to_string(i+1), "int")) err++;
				if (!check_bounds("release_day" + std::to_string(i+1), r_d, 0, true, max_t, true)) err++;

				if (err == 0) {
					temp.push_back(r_d);
				}
				else tot_err += 1;
			}
		}
		file.close();

		if (tot_err == 0) {
			model_params->rel->release_times = temp;
		}
		else {
			std::cerr << "There were errors in the file. The simulation will run with rel_times = driver_start" << std::endl;
		}	
	}
}

/**
 * Sets the values of the f_{ijk} inheritance fraction for the gene drive considering r2 non-functional resistance alleles. f_{ijk} denotes the fraction of genotype k offspring from mother with genotype i mated to father with genotype j.
 * @details The order of elements in each matrix axis is the following genotype order: WW, WD, DD, WR, RR and DR, composed of wild-type (W), drive-type (D) and resistant-type (R) alleles. 
 * @note Six genotypes are counted and not nine because WD and DW genotypes are counted together, and likewise for the other heterozygous genotypes.
 * @note  DD, RR and DR females are assumed to be sterile as they don't possess one functional copy of the dsx gene and thus produce no offspring.
 * @param[in] inher_params inheritance parameters
 * @see InheritanceParams, InputParams:gamma, InputParams:xi, InputParams::e
 */
void Simulation::set_inheritance(InheritanceParams inher_params)
{
	double omega_M = inher_params.omega_M;
	double omega_F = inher_params.omega_F;
	double b=(1+inher_params.bias)*0.5;
	std::cout<<"Sim1  "<<b<<std::endl;
//juvenile types: (wwf,wzf,zzf,wwm,wzm,zzm), where w is wildtype, z is transgenic, f is female and m is male
	std::array<double, 6> f_ww_ww = {0.5, 0, 0, 0.5, 0, 0};
	std::array<double, 6> f_ww_wz = {(1-b)*0.5,(1-b)*0.5 , 0,b*0.5 , b*0.5, 0};
	std::array<double, 6> f_ww_zz = {0,1-b , 0, 0, b, 0};
	std::array<double, 6> f_wz_ww = {0.25, 0.25, 0, 0.25,0.25, 0};
	std::array<double, 6> f_wz_wz = {(1-b)*0.25,(1-b)*0.5,(1-b)*0.25,b*0.25,b*0.5,b*0.25};
	std::array<double, 6> f_wz_zz = {0,(1-b)*0.5,(1-b)*0.5,0,b*0.5,b*0.5};
	std::array<double, 6> f_zz_ww = {0,0.5,0,0,0.5,0};
	std::array<double, 6> f_zz_wz = {0,(1-b)*0.5,(1-b)*0.5,0,b*0.5,b*0.5};
	std::array<double, 6> f_zz_zz = {0, 0,1-b , 0, 0, b};

	std::cout<<"Sim2  "<<b<<std::endl;
	std::cout << std::extent<decltype(inher_fraction), 0>::value << " x "
          << std::extent<decltype(inher_fraction), 1>::value << " x "
          << std::extent<decltype(inher_fraction), 2>::value << "\n";
	std::cout<<"Sim2  "<<b<<std::endl;

// need to update inheritance code below!!

	for (int k=0; k<6; ++k) {
		for (int i=0; i<3;++i) {
			for (int j=0; j<3; ++j) {
				if (i==0) {
					if (j==0) inher_fraction[i][j][k] = f_ww_ww[k];
					else if (j==1) inher_fraction[i][j][k] = omega_M*f_ww_wz[k];
					else if (j==2) inher_fraction[i][j][k] = omega_M*f_ww_zz[k];
				}
				else if (i==1) {
					if (j==0) inher_fraction[i][j][k] = omega_F*f_wz_ww[k];
					else if (j==1) inher_fraction[i][j][k] = omega_F*omega_M*f_wz_wz[k];
					else if (j==2) inher_fraction[i][j][k] = omega_F*omega_M*f_wz_zz[k];
				}
				else if (i==2) {
					if (j==0) inher_fraction[i][j][k] = omega_F*f_zz_ww[k];
					else if (j==1) inher_fraction[i][j][k] = omega_F*omega_M*f_zz_wz[k];
					else if (j==2) inher_fraction[i][j][k] = omega_F*omega_M*f_zz_zz[k];
				}
std::cout<<i<<"  "<<j<<"  "<<k<<"   "<<inher_fraction[i][j][k]<<std::endl;
			}
		}
	}	
	std::cout<<"Sim3  "<<b<<std::endl;
}

/**
 * Runs the simulation num_runs times, recording data in output files.
 * @see InputParams::num_runs, Record
 */
void Simulation::run_reps() 
{
	for (int rep=1; rep <= num_runs; ++rep) {
		Model* model;
		if (!((input_rainfall_params->rainfall).empty())) {
			model = new Model(model_params, inher_fraction, input_rainfall_params, alpha0_mean, alpha0_variance, release_sites, boundary_type, disp_type, sites_coords,humans);
		}
		else {
			model = new Model(model_params, inher_fraction, sine_rainfall_params, alpha0_mean, alpha0_variance, release_sites, boundary_type, disp_type, sites_coords,humans);
		}
		Record data(rec_params, rep);
		model->initiate();
		data.record_coords(model->get_sites());

		for (int tt=start_t; tt <= max_t; ++tt) { // current day of the simulation 
			today=tt;
			model->run(tt);

			if (data.is_rec_global_time(tt)) {
				data.output_totals(tt, model->calculate_tot_J(), model->calculate_tot_M(), model->calculate_tot_V(),
				 model->calculate_tot_F());
				data.record_global(tt, model->calculate_tot_M_gen());
			}
			if (data.is_rec_local_time(tt)) {
				data.record_local(tt, model->get_sites());
			}
		}

		delete model;
	}
}
