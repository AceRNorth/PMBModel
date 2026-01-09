#ifndef GENERALMETAPOP_PATCH_H
#define GENERALMETAPOP_PATCH_H

#include <array>
#include "Patch.h"
#include "constants.h"
#include "Model.h"
#include "Params.h"
#include "Point.h"

class Model;

/**
 * Contains the information of a local mosquito population. The population is divided into four types: juveniles (J), males (M), virgin (unmated) females (V) and mated females (F). These are then subdivided into six genotypes (in this order): WW, WD, DD, DR, WR and RR composed of a wild-type allele (W), a drive-type allele (D) and a resistant-type allele (R). Juveniles are also subdivided into age groups, ordered from oldest (0 days left to eclosion) to youngest (``max_dev`` - 1 days left).
 * The population can carry out life-processes and interface with other classes to introduce gene drive mosquitoes into the Patch, carry out dispersal in and out of the Patch and aestivate. 
 * @image xml patch_diagram.png "Patch behaviour diagram"
 */
class Patch {
public:
	Patch(Model* mod, LifeParams* par, double a0, double side_x, double side_y); // *
	Patch(Model* mod, LifeParams* par, double a0, Point point, double hum);
	void populate(int initial_WJ, int initial_WM, int initial_WV, int initial_WF);

	Point get_coords() const;
	std::array<long long int, constants::num_gen> get_M() const;
	std::array<std::array<long long int, constants::num_gen>, constants::num_gen> get_F() const;

	long long int calculate_tot_J();
	long long int calculate_tot_M();
	long long int calculate_tot_V();
	long long int calculate_tot_F();

	// life-processes for the local site

	void juv_get_older();
	void adults_die();
	void virgins_mate();
	void lay_eggs(const std::array<std::array<std::array <double,2*constants::num_gen>, constants::num_gen>, constants::num_gen> &f,
	 const std::array<double, constants::max_dev+1> &dev_duration_probs);
	void juv_eclose();
	void update_comp();
	void update_mate();
	
	// interface to Dispersal

	void M_disperse_out(const std::array<long long int, constants::num_gen> &m_out);
	void F_disperse_out(const std::array<std::array<long long int, constants::num_gen>, constants::num_gen> &f_out);
	void M_disperse_in(int gen, long long int m_in);
	void F_disperse_in(int f_gen, int m_gen, long long int f_disp);

	// interface to Aestivation

	void F_hide(const std::array<std::array<long long int, constants::num_gen>, constants::num_gen> &f_try);
	void F_wake(const std::array<std::array<long long int, constants::num_gen>, constants::num_gen> &f_wake);

	// interface to GDRelease

	void add_driver_M(int num_driver_M);

private:
	Model* model; /**< Pointer to Model object */
	LifeParams* params;  /**< Pointer to LifeParams object */
	double alpha0; /**< Baseline contribution to the carrying capacity */ 

	Point coords; /**< Patch coordinates */ 

	double humandens;

	std::array<std::array<long long int, constants::max_dev+1>, 2*constants::num_gen> J; /**< Number of juvenile mosquitoes in the patch, divided by genotype and age group. @note The age groups are ordered from oldest (0 days left to eclosion) to youngest (max_dev - 1 days left). */
	std::array<long long int, constants::num_gen> M; /**< Number of male mosquitoes in the patch, divided by genotype. */ 
	std::array<long long int, constants::num_gen> V; /**< Number of virgin (unmated) female mosquitoes, divided by genotype. */ 
	std::array<std::array<long long int, constants::num_gen>, constants::num_gen> F; /**< Number of mated female mosquitoes F_{ij}, divided by female genotype i and male sperm genotype j. */ 

	long double comp; /**< Survival probability per juvenile per day (both density-dependent and independent factors). */ 
	long double mate_rate; /**< Probability of a virgin (unmated) female mating on a given day. */
};

#endif //GENERALMETAPOP_PATCH_H
