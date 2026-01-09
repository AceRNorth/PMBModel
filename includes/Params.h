#ifndef GENERALMETAPOP_PARAMS_H
#define GENERALMETAPOP_PARAMS_H

#include <vector>

// ** These parameters should not be modified after being passed to the Simulation or Model **

/**
 * Simulation progression parameters.
 */ 
struct ProgressionParams {
	int num_runs; /**< Number of simulation replicates to run. */ 
	int max_t; /**< Maximum simulated time (in days). */ 
};

/**
 *  Model area parameters.
 */ 
struct AreaParams {
	int num_pat; /**< Number of population sites chosen for the simulation. */ 
};

/**
 * Model life-process parameters.
 */ 
struct LifeParams {
	double mu_j; /**< Juvenile density independent mortality rate per day. */ 
	double mu_a; /**< Adult mortality rate per day. */ 
	double beta; /**< Parameter that controls mating rate. */ 
	double theta; /**< Average egg laying rate of wildtype females (eggs per day). */ 
	double comp_power; /**< Parameter that controls the juvenile survival probability. */
	int min_dev; /**< Minimum development time for a juvenile (in days). */ 
}; 

/**
 * Gene drive inheritance parameters.
 */ 
struct InheritanceParams {
//	double gamma; 
//	double xi; 
//	double e; 
	double bias;
	double omega_M;
	double omega_F;
};

/**
 * Gene drive release model parameters.
 */ 
struct ReleaseParams {
	int num_driver_M; /**< Number of drive heterozygous (WD) male mosquitoes per release. */ 
	int num_driver_sites; /**< Number of gene drive release sites per year. */ 
	std::vector<int> release_times; /**< Days on which the gene drive mosquitoes will be released. */ 
};

/**
 * Dispersal model parameters.
 */ 
struct DispersalParams {
	double disp_rate; /**< Adult dispersal rate. */ 
	double max_disp; /**< Maximum distance at which two sites are connected. */ 
};

/**
 * Aestivation model parameters.
 */ 
struct AestivationParams {
	double psi; /**< Aestivation rate. */ 
	double mu_aes; /**< Aestivation mortality. */ 
	int t_hide1; /**< Start day of aestivation-hiding period (exclusive). */ 
	int t_hide2; /**< End day of aestivation-entering period (inclusive). */ 
	int t_wake1; /**< Start day of aestivation-waking period (exclusive). */ 
	int t_wake2; /**< End day of aestivation-waking period (inclusive). */ 
};

/**
 * Seasonality model parameters for sinusoid rainfall contribution to carrying capacity.
 */ 
struct SineRainfallParams {
    double alpha1; /**< Rainfall contribution factor to carrying capacity. */ 
    double amp; /**< Amplitude of rainfall fluctuations. */ 
};

/**
 * Seasonality model parameters for rainfall contribution to carrying capacity from rainfall data. 
 */ 
struct InputRainfallParams {
	double alpha1; /**< Rainfall contribution factor to carrying capacity. */ 
	double resp; /**< Carrying capacity's responsiveness to rainfall contribution. */ 
    std::vector<double> rainfall; /**< Daily rainfall data. */
};

/**
 * Initial population values for the model.
 */
struct InitialPopsParams {
	int initial_WJ; /**< Number of initial juvenile mosquitoes with wild homozygous (WW) genotype for each age group. */ 
	int initial_WM; /**< Number of initial adult male mosquitoes with wild homozygous (WW) genotype. */ 
	int initial_WV; /**< Number of initial adult unmated female (virgin) mosquitoes with wild homozygous (WW) genotype. */ 
	int initial_WF; /**< Number of initial adult mated female mosquitoes with wild homozygous (WW) genotype. */ 
};

/**
 * Model parameters. 
 */
struct ModelParams {
	AreaParams *area;
	InitialPopsParams *initial;
	LifeParams *life;
	AestivationParams *aes;
	DispersalParams *disp; 
	ReleaseParams *rel;
};

/**
 * Data-recording parameters.
 */ 
struct RecordParams {
	// recording window and intervals
	int rec_start; /**< Start time for the data recording window (in days) (inclusive). */ 
	int rec_end; /**< End time for the data recording window (in days) (inclusive). */ 
	int rec_interval_global; /**< Time interval for global data recording/output. */ 
	int rec_interval_local; /**< Time interval at which to collect/record local data (in days). */ 
	int rec_sites_freq; /**< Fraction of sites to collect local data for (1 is all sites, 10 is 1 in 10 etc). */ 

	// output filename labels
	int set_label; /**< 'Set of repetitions' index label for output files. */ 
};

#endif //GENERALMETAPOP_PARAMS_H
