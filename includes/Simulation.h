#ifndef GENERALMETAPOP_SIMULATION_H
#define GENERALMETAPOP_SIMULATION_H

#include <array>
#include <vector>
#include <string>
#include <filesystem>
#include "constants.h"
#include "Params.h"
#include "InputParams.h"
#include "Point.h"

/**
 * Sets up and controls the flow of the simulation.
 * @details Running set_inheritance() and run_reps() in this order is essential to run a simulation. The class also contains setters for advanced options, which should be run before run_reps(). 
 * @note Advanced options that have been set apply to all repetitions (runs) of a single simulation.
 * @warning Advanced options will typically give an error message if the setting process has been unsuccessful. If they are not re-set successfully before running the simulation the option will maintain its default behaviour.
 * @see InputParams::num_runs, InputParams::max_t
 */
class Simulation {
public:
	Simulation(InputParams input);
	~Simulation();
	void set_coords(const std::filesystem::path& filepath);
	void set_boundary_type(BoundaryType boundary);
	void set_dispersal_type(DispersalType disp);
	void set_rainfall(const std::filesystem::path& filepath);
	void set_release_times(const std::filesystem::path& filepath);
	void set_inheritance(InheritanceParams inher_params); 
	void run_reps();

private:
	int num_runs; /**< Number of simulation replicates to run. */ 
	int max_t; /**< Maximum simulated time (in days). */ 

	ModelParams *model_params; /**< Model parameters. */ 
	RecordParams *rec_params; /**< Data-recording parameters. */ 

	SineRainfallParams *sine_rainfall_params; /** Sinusoid rainfall seasonality parameters.  */
	InputRainfallParams *input_rainfall_params; /** Input data rainfall seasonality parameters. */
	double alpha0_mean; /**< Seasonality parameter. */
	double alpha0_variance; /**< Seasonality parameter. */ 

	// additional parameter options
	std::vector<Point> sites_coords; /**< 2D coordinates for the sites on the simulated square. */ 
	std::vector<double> humans; /**< 2D coordinates for the sites on the simulated square. */ 
	std::vector<int> release_sites; /**< Indices relative to coords vector for the gene drive release sites. */ 
//	double building_dens;
	BoundaryType boundary_type;
	DispersalType disp_type;

	// inheritance
	std::array<std::array<std::array <double, constants::num_gen>, constants::num_gen>, 2*constants::num_gen> inher_fraction; /**< Inheritance fraction. f_ijk is the fraction of genotype k offspring from mother with genotype i mated to father with genotype j. */
};

#endif //GENERALMETAPOP_SIMULATION_H
