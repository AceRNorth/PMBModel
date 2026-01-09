#include <vector>
#include <cassert>
#include <algorithm>
#include "Model.h"
#include "random.h"
#include "constants.h"
#include "Patch.h"
#include "Dispersal.h"
#include "GDRelease.h"
#include "Aestivation.h"
#include "globals.h"

using namespace constants;

/**
 * Model constructor.
 * Builds the components of the model for sinusoid rainfall-type seasonality.
 * @param[in] params 		model parameters
 * @param[in] inher_frac 	inheritance fraction
 * @param[in] season 		seasonality parameters for sinusoid rainfall
 * @param[in] a0_mean 		alpha0 mean (mean of the baseline contribution to the carrying capacity)
 * @param[in] a0_var 		alpha0 variance (variance of the baseline contribution to the carrying capacity)
 * @param[in] rel_sites		release sites (indices relative to the coords vector)
 * @param[in] boundary		boundary type
 * @param[in] disp_type		dispersal type
 * @param[in] coords		site coordinates vector
 */
Model::Model(ModelParams* params, const std::array<std::array<std::array <double, 2*constants::num_gen>, constants::num_gen>, constants::num_gen> &inher_frac, SineRainfallParams* season,
 double a0_mean, double a0_var, std::vector<int> rel_sites, BoundaryType boundary, DispersalType disp_type, std::vector<Point> coords,std::vector<double>humandens)
{
	num_pat = params->area->num_pat;
	initial_pops = params->initial;
	min_dev = params->life->min_dev;
	alpha0_mean = a0_mean;
	alpha0_variance = a0_var;
	dev_duration_probs.fill(0);
	inher_fraction = inher_frac;

	day_sim = 0;
	sites.clear();
	side_x = 0;
	side_y = 0;
	if (!coords.empty()) {
		assert(coords.size() == num_pat);

		// calculate sides for toroidal boundary
		auto compare_x = [](const Point& a, const Point& b){return a.x < b.x;}; 
		auto x_min = (*std::min_element(coords.begin(), coords.end(), compare_x)).x;
		auto x_max = (*std::max_element(coords.begin(), coords.end(), compare_x)).x;
		auto compare_y = [](const Point& a, const Point& b){return a.y < b.y;};
		auto y_min = (*std::min_element(coords.begin(), coords.end(), compare_y)).y;
		auto y_max = (*std::max_element(coords.begin(), coords.end(), compare_y)).y;
		side_x = x_max - x_min;
		side_y = y_max - y_min;
	
		for (int i=0; i < num_pat; ++i) {
			Patch* pp = new Patch(this, params->life, alpha0(), coords[i],humandens[i]);
			sites.push_back(pp);
		}
	}
	else {
		side_x = 1;
		side_y = 1;
		for (int i=0; i < num_pat; ++i) {
			Patch* pp = new Patch(this, params->life, alpha0(), side_x, side_y);
			sites.push_back(pp);
		}
	}

	Aestivation* new_aestivation = new Aestivation(params->aes, sites.size());
	aestivation = new_aestivation;

	Dispersal* new_disp;
	if (disp_type == DistanceKernel) {
		new_disp = new DistanceKernelDispersal(params->disp, boundary, side_x, side_y);
	}
	else if (disp_type == Radial) {
		new_disp = new RadialDispersal(params->disp, boundary, side_x, side_y);
	}
	else {
		new_disp = new DistanceKernelDispersal(params->disp, boundary, side_x, side_y);
	}
	dispersal = new_disp;

	GDRelease* new_gd_release;
	if (coords.empty()) {
		new_gd_release = new RandomGDRelease(params->rel);
	}
	else {
		new_gd_release = new SchedGDRelease(params->rel, rel_sites, sites);
	}
	gd_release = new_gd_release;

	Seasonality* new_seasonality = new SineRainfall(season);
	seasonality = new_seasonality;
}

/**
 * Model constructor.
 * Builds the components of the model for input rainfall-type seasonality.
 * @param[in] params 		model parameters
 * @param[in] inher_frac 	inheritance fraction
 * @param[in] season 		seasonality parameters for input rainfall data
 * @param[in] a0_mean 		alpha0 mean (mean of the baseline contribution to the carrying capacity)
 * @param[in] a0_var 		alpha0 variance (variance of the baseline contribution to the carrying capacity)
 * @param[in] rel_sites		release sites (indices relative to the coords vector)
 * @param[in] boundary		boundary type
 * @param[in] disp_type		dispersal type
 * @param[in] coords		site coordinates vector
 */
Model::Model(ModelParams* params, const std::array<std::array<std::array <double, 2*constants::num_gen>, constants::num_gen>, constants::num_gen> &inher_frac, InputRainfallParams *season,
 double a0_mean, double a0_var, std::vector<int> rel_sites, BoundaryType boundary, DispersalType disp_type, std::vector<Point> coords,std::vector<double>humandens)
{
	num_pat = params->area->num_pat;
	initial_pops = params->initial;
	min_dev = params->life->min_dev;
	alpha0_mean = a0_mean;
	alpha0_variance = a0_var;
	dev_duration_probs.fill(0);
	inher_fraction = inher_frac;

	day_sim = 0;
	sites.clear();

	if (!coords.empty()) {
		assert(coords.size() == num_pat);

		// calculate sides for toroidal boundary
		auto compare_x = [](const Point& a, const Point& b){return a.x < b.x;};
		auto x_min = (*std::min_element(coords.begin(), coords.end(), compare_x)).x;
		auto x_max = (*std::max_element(coords.begin(), coords.end(), compare_x)).x;
		auto compare_y = [](const Point& a, const Point& b){return a.y < b.y;};
		auto y_min = (*std::min_element(coords.begin(), coords.end(), compare_y)).y;
		auto y_max = (*std::max_element(coords.begin(), coords.end(), compare_y)).y;
		side_x = x_max - x_min;
		side_y = y_max - y_min;

		for (int i=0; i < num_pat; ++i) {
			Patch* pp = new Patch(this, params->life, alpha0(), coords[i],humandens[i]);
			sites.push_back(pp);
		}
	}
	else {
		side_x = 1;
		side_y = 1;
		for (int i=0; i < num_pat; ++i) {
			Patch* pp = new Patch(this, params->life, alpha0(), side_x, side_y);
			sites.push_back(pp);
		}
	}

	Aestivation* new_aestivation = new Aestivation(params->aes, sites.size());
	aestivation = new_aestivation;

	Dispersal* new_disp;
	if (disp_type == DistanceKernel) {
		new_disp = new DistanceKernelDispersal(params->disp, boundary, side_x, side_y);
	}
	else if (disp_type == Radial) {
		new_disp = new RadialDispersal(params->disp, boundary, side_x, side_y);
	}
	else {
		new_disp = new DistanceKernelDispersal(params->disp, boundary, side_x, side_y);
	}
	dispersal = new_disp;

	GDRelease* new_gd_release;
	if (!rel_sites.empty()) {
		new_gd_release = new SchedGDRelease(params->rel, rel_sites, sites);
	}
	else {
		new_gd_release = new RandomGDRelease(params->rel);
	}
	gd_release = new_gd_release;

	Seasonality* new_seasonality = new InputRainfall(season);
	seasonality = new_seasonality;
}

/**
 * Model destructor.
 */
Model::~Model() 
{
	delete aestivation;
	delete dispersal;
	delete gd_release;
	delete seasonality;
	
	for (auto pat : sites) {
		delete pat;
	}
}

/**
 * Returns a random value for the baseline contribution to the carrying capacity ``alpha0``. 
 * @details The random value is obtained as a random draw from a lognormal distribution with desired mean ``alpha0_mean`` and desired variance ``alpha0_variance``.
 * @return A random draw of ``alpha0``.
 * @see InputParams::alpha0_mean, InputParams::alpha0_variance, random_lognormal()
 */
double Model::alpha0() 
{
   return random_lognormal(alpha0_mean, alpha0_variance);
}

/**
 * Sets up the model architecture. Populates sites and sets dispersal connections. 
 */
void Model::initiate()
{
	populate_sites();
	set_dev_duration_probs(min_dev, constants::max_dev);
	dispersal->set_connecs(sites); 
}

/**
 * Populates all sites with a wild mosquito population of different types (age and sex).
 */
void Model::populate_sites() 
{
	for (auto pat : sites) {
		pat->populate(initial_pops->initial_WJ, initial_pops->initial_WM, initial_pops->initial_WV, initial_pops->initial_WF);
	}
}
 
/**
 * Sets probabilities of juvenile eclosion for different age groups.
 * @param[in] min_time minimum development time (in days)
 * @param[in] max_time maximum development time (in days)
 */
void Model::set_dev_duration_probs(int min_time, int max_time) 
{
	for (int a=0; a < max_time + 1; ++a) {
        if (a >= min_time) {
            dev_duration_probs[a] = 1.0 / (max_time - min_time);
        }
        else {
            dev_duration_probs[a] = 0;
        }
    }
}

/**
 * Runs the daily step of the model.
 * @details The life processes do not run during the initialisation day (day 0). However, gene drive release can be carried out on day 0.
 * @param[in] day simulation day
 */ 
void Model::run(int day)
{
	day_sim = day; // used later for seasonality
	gd_release->release_gene_drive(day, sites);
	if (day > 0) {
		run_step(day);
	}
}

/**
 * Runs the daily mosquito life-processes for all sites, including dispersal and aestivation.
 * @param[in] day simulation day
 */
void Model::run_step(int day) 
{
	juv_get_older();
	adults_die();
	virgins_mate();
	dispersal->adults_disperse(sites);
	lay_eggs();
	juv_eclose();
	if (aestivation->is_hide_time(day)) aestivation->hide(sites);
	if (aestivation->is_wake_time(day)) aestivation->wake(day, sites);
}

/**
 * Returns the total number of juveniles across all ages and genotypes and across all patches.
 * @return The total number of juveniles in the model run. 
 * @see Patch
 */
long long int Model::calculate_tot_J()
{
	long long int tot_J = 0;
	for (auto pat : sites) {
		tot_J += pat->calculate_tot_J();
	}
	return tot_J;
}
 
/**
 * Returns the total number of males across all genotypes and across all patches.
 * @return The total number of males in the model run.
 * @see Patch::get_M()
 */
long long int Model::calculate_tot_M()
{
	long long int tot_M = 0;
	for (auto pat : sites) {
		tot_M += pat->calculate_tot_M();
	}
	return tot_M;
}

/**
 * Returns the total number of virgin (unmated) females across all genotypes and across all patches.
 * @return The total number of virgin (unmated) females in the model run.
 * @see Patch
 */
long long int Model::calculate_tot_V()
{
	long long int tot_V = 0;
	for (auto pat : sites) {
		tot_V += pat->calculate_tot_V();
	}
	return tot_V;
}

/**
 * Returns the total number of mated females across all female and male sperm genotypes and across all patches.
 * @return The total number of mated females in the model run.
 * @see Patch::get_F()
 */
long long int Model::calculate_tot_F()
{
	long long int tot_F = 0;
	for (auto pat : sites) {
		tot_F += pat->calculate_tot_F();
	}
	return tot_F;
}

/**
 * Returns the total number of males of each genotype across all patches.
 * @return The total number of males in the model run, divided by genotype.
 * @see Patch::get_M()
 */
std::array<long long int, constants::num_gen> Model::calculate_tot_M_gen() 
{
	std::array<long long int, constants::num_gen> tot_M_gen;
	tot_M_gen.fill(0);
	for (auto pat : sites) {
		std::array<long long int, constants::num_gen> m_pat = pat->get_M();
		for (int i = 0; i < constants::num_gen; ++i) {
			tot_M_gen[i] += m_pat[i];
		}
	}
	return tot_M_gen;
}



/**
 * Returns the total number of females of each genotype across all patches.
 * @return The total number of females in the model run, divided by genotype.
 * @see Patch::get_M()
 */
std::array<long long int, constants::num_gen> Model::calculate_tot_F_gen() 
{
	std::array<long long int, constants::num_gen> tot_F_gen;
	tot_F_gen.fill(0);
	for (auto pat : sites) {
		std::array<std::array<long long int, constants::num_gen>, constants::num_gen> f_pat = pat->get_F();
		for (int i = 0; i < constants::num_gen; ++i) {
		for (int j = 0; j < constants::num_gen; ++j) {
			tot_F_gen[i] += f_pat[i][j];
		}
		}
	}
	return tot_F_gen;
}











/**
 * Returns the sites vector.
 * @details The patches within the vector are ordered according to order of creation. If specific coordinates have been set for patches, the sites vector will follow the same order as the coordinates vector of the Model constructor.
 * @return The sites vector, containing all Patch objects. 
 */
std::vector<Patch*> Model::get_sites() const
{
	return sites;
}

/**
 * Returns the current simulation day. 
 * @details Simulation days run from 0 (initialisation day) up to ``max_t`` (inclusive). 
 * @return The simulation day.
 */
int Model::get_day() const
{
	return day_sim;
}

/**
 * Returns the carrying capacity for the given alpha0 (the baseline contribution to carrying capacity) and for the current day.
 * @param[in] alpha0 baseline contribution to the carrying capacity
 * @see Seasonality::alpha()
 */
double Model::get_alpha(double alpha0)
{
	double alpha = seasonality->alpha(day_sim, alpha0);
	return alpha;
}

/**
 * Carries out juvenile aging across all patches.
 * @see Patch::juv_get_older()
 */
void Model::juv_get_older() 
{
	for (auto pat : sites) {
		pat->juv_get_older();
	}
}

/**
 * Carries out adult mortality across all patches.
 * @see Patch::adults_die()
 */
void Model::adults_die()
{
	for (auto pat : sites) {
		pat->adults_die();
	}
}

/**
 * Carries out mating across all patches.
 * @see Patch::virgins_mate()
 */
void Model::virgins_mate() 
{
	for (auto pat : sites) {
		pat->virgins_mate();
	}
}

/**
 * Carries out egg-laying across all patches.
 * @see Patch::lay_eggs()
 */
void Model::lay_eggs()
{
	for (auto pat : sites) {
		pat->lay_eggs(inher_fraction, dev_duration_probs);
	}
}

/**
 * Carries out juvenile eclosion across all patches.
 * @see Patch::juv_eclose()
 */ 

void Model::juv_eclose()
{
	for (auto pat : sites) {
		pat->juv_eclose();
	}
}
