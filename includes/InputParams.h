#ifndef GENERALMETAPOP_INPUTPARAMS_H
#define GENERALMETAPOP_INPUTPARAMS_H

/**
 * Input parameters for the simulation.
 */
struct InputParams {
	int num_runs; /**< Number of simulation replicates to run. */ 
	int max_t; /**<  Maximum simulated time (in days). */
	int num_pat; /**< Number of population sites chosen for the simulation. */ 

	double mu_j; /**< Juvenile density independent mortality rate per day. */ 
	double mu_a; /**< Adult mortality rate per day. */ 
	double beta; /**< Parameter that controls mating rate. */ 
	double theta; /**< Average egg laying rate of wildtype females (eggs per day). */ 
	double comp_power; /**< Parameter that controls the juvenile survival probability. */
	int min_dev; /**< Minimum development time for a juvenile (in days). */ 

//	double gamma; 
//	double xi; 
//	double e; 

	double bias;
	double omega_M;
	double omega_F;

	int driver_start; /**< Time to start releasing drive alleles into the mosquito population. */ 
	int num_driver_M; /**< Number of drive heterozygous (WD) male mosquitoes per release. */ 
	int num_driver_sites; /**< Number of gene drive release sites per year. */ 

	double disp_rate; /**< Adult dispersal rate. */ 
	double max_disp; /**< Maximum dispersal distance at which two sites are connected. */ 

	double psi; /**< Aestivation rate. */ 
	double mu_aes; /**< Aestivation mortality rate. */ 
	int t_hide1; /**< Start day of aestivation-hiding period (exclusive). */ 
	int t_hide2; /**< End day of aestivation-hiding period (inclusive). */ 
	int t_wake1; /**< Start day of aestivation-waking period (exclusive). */ 
	int t_wake2; /**< End day of aestivation-waking period (inclusive). */ 
    double alpha0_mean; /**< Mean of the baseline contribution to the carrying capacity. */ 
    double alpha0_variance; /**< Variance of the baseline contribution to the carrying capacity. */ 
    double alpha1; /**< Rainfall contribution factor to carrying capacity. */ 
    double amp; /**< Amplitude of rainfall fluctuations. */ 
    double resp; /**< Carrying capacity's responsiveness to rainfall contribution. */ 

	int initial_WJ = 10000; /**< Number of initial juvenile mosquitoes with wild homozygous (WW) genotype for each age group. */ 
	int initial_WM = 50000; /**< Number of initial adult male mosquitoes with wild homozygous (WW) genotype. */ 
	int initial_WV = 10000; /**< Number of initial adult unmated female (virgin) mosquitoes with wild homozygous (WW) genotype. */ 
	int initial_WF = 40000; /**< Number of initial adult mated female mosquitoes with wild homozygous (WW) genotype. */ 

	int rec_start; /**< Start time for the data recording window (inclusive). */ 
	int rec_end; /**< End time for the data recording window (inclusive). */ 
	int rec_interval_global; /**< Time interval for global data recording/output. */ 
	int rec_interval_local; /**< Time interval at which to collect/record local data (in days). */ 
	int rec_sites_freq; /**< Fraction of sites to collect local data for (1 is all sites, 10 is 1 in 10 etc). */ 
	int set_label; /**< 'Set of repetitions' index label for output files. */ 
};

#endif //GENERALMETAPOP_INPUTPARAMS_H
