#ifndef GENERALMETAPOP_MODEL_H
#define GENERALMETAPOP_MODEL_H

#include <vector>
#include <array>
#include "constants.h"
#include "Params.h"
#include "Patch.h"
#include "Dispersal.h"
#include "Aestivation.h"
#include "GDRelease.h"
#include "Seasonality.h"
#include "Point.h"

class Patch;
class Aestivation;
class Dispersal;
class GDRelease;
class Seasonality;

/**
 * Sets up and runs the model.
 * @details In addition to the set-up of the model, Model coordinates the different classes to run daily steps of mosquito life processes and dispersal as well as other events like gene drive release, aestivation and seasonality. It also contains methods to interface with the Record class.
 * @see AreaParams, ModelParams
 */ 
class Model {
public:
	Model(ModelParams *params, const std::array<std::array<std::array <double, constants::num_gen>, constants::num_gen>, 2*constants::num_gen> &inher_frac, SineRainfallParams *season, double a0_mean, double a0_var, std::vector<int> rel_sites = {}, BoundaryType boundary = BoundaryType::Toroid, DispersalType disp_type = DispersalType::DistanceKernel, std::vector<Point> coords = {},std::vector<double> humans = {});
	Model(ModelParams *params, const std::array<std::array<std::array <double, constants::num_gen>, constants::num_gen>, 2*constants::num_gen> &inher_frac, InputRainfallParams *season, double a0_mean, double a0_var, std::vector<int> rel_sites = {}, BoundaryType boundary = BoundaryType::Toroid, DispersalType disp_type = DispersalType::DistanceKernel, std::vector<Point> coords = {},std::vector<double> humans = {});
	
	~Model();
	void initiate();
	void run(int day);

	long long int calculate_tot_J(); 
	long long int calculate_tot_M();
	long long int calculate_tot_V();
	long long int calculate_tot_F();
	std::array<long long int, constants::num_gen> calculate_tot_M_gen();
	std::vector<Patch*> get_sites() const;
	int get_day() const;
	double get_alpha(double alpha0);

private:
	std::vector<Patch*> sites; /**< Vector of all Patch objects. Contains the population sites over the whole simulation area. */
	Dispersal* dispersal;
	Aestivation* aestivation;
	GDRelease* gd_release;
	Seasonality* seasonality;

	int day_sim; /**< Current day of the simulation. */ 

	// simulation area parameters
	int num_pat; /**< Number of population sites (patches) chosen for the simulation. */ 
	double side_x; /**< Size of one side of the simulation area (x-axis). */
	double side_y; /**< Size of one side of the simulation area (y-axis). */

	InitialPopsParams *initial_pops; /**< Initial population values - common for all Patches. */

	// juvenile development parameters - common for all Patches
	int min_dev; /**< Minimum development time for a juvenile (in days). */ 
	std::array<double, constants::max_dev+1> dev_duration_probs; /**< Array of probabilities of juvenile development duration for a new juvenile (index indicates the number of days to develop or, equivalently, the age class the new juvenile starts at). */ 

	// inheritance
	std::array<std::array<std::array <double, constants::num_gen>, constants::num_gen>, 2*constants::num_gen> inher_fraction; /**< Inheritance fraction. f_ijk is the fraction of genotype k offspring from mother with genotype i mated to father with genotype j. */

	// Patch construction parameters
	double alpha0_mean; /**< Mean of the baseline contribution to the carrying capacity. */ 
    double alpha0_variance; /**< Variance of the baseline contribution to the carrying capacity. */ 
	double alpha0();

	// initiation methods
	void populate_sites();
	void set_dev_duration_probs(int min_time, int max_time);

	// life-processes - interface with Patch
	void run_step(int day);
	void juv_get_older();
	void adults_die();
	void virgins_mate();
	void lay_eggs();
	void juv_eclose();
};

#endif //GENERALMETAPOP_MODEL_H
