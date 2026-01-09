#include <cmath>
#include <array>
#include <vector>
#include "Patch.h"
#include "random.h"
#include "constants.h"
#include "globals.h"
#include <iostream> // for error messages

/**
 * @brief Patch constructor for randomly generated coordinates.
 * @details Sets coordinates randomly on the square simulation space of size side x side. 
 * @note Coordinates range from 0 to side for bound checking purposes. 
 * @param[in] mod 		Model pointer
 * @param[in] par 		life parameters
 * @param[in] a0 		alpha0 carrying capacity baseline
 * @param[in] side_x 	size of one side of the simulation area (x-axis)
 * @param[in] side_y 	size of one side of the simulation area (y-axis)
 */
Patch::Patch(Model* mod, LifeParams* par, double a0, double side_x, double side_y) 
{
	model = mod;
	params = par;
	alpha0 = a0;
	
	double x = random_real() * side_x;
	double y = random_real() * side_y;
	coords = {x, y};
	humandens=1;

	for (int i=0; i < 2*constants::num_gen; ++i) {
		for (int a=0; a < constants::max_dev + 1; ++a) {
			J[i][a] = 0; 
		}
	}
	for (int i=0; i < constants::num_gen; ++i) {
		M[i] = 0;
		V[i] = 0;
		for (int j=0; j < constants::num_gen; ++j) {
			F[i][j] = 0;
		}
	}
	comp = 0;
	mate_rate = 0;
}

/**
 * @brief Patch constructor for custom coordinates.
 * @param[in] mod 	Model pointer
 * @param[in] par 	life parameters
 * @param[in] a0 	alpha0 carrying capacity baseline
 * @param[in] point	patch coordinates
 */
Patch::Patch(Model* mod, LifeParams* par, double a0, Point point, double humans) 
{
	model = mod;
	params = par;
	alpha0 = a0;
	coords = point;
	humandens=humans;

	// include to be able to compare data to test data when testing
	double x = random_real();
	double y = random_real();
	for (int i=0; i < 2*constants::num_gen; ++i) {
		for (int a=0; a < constants::max_dev + 1; ++a) {
			J[i][a] = 0; 
		}
	}

	for (int i=0; i < constants::num_gen; ++i) {
		M[i] = 0;
		V[i] = 0;
		for (int j=0; j < constants::num_gen; ++j) {
			F[i][j] = 0;
		}
	}
	comp = 0;
	mate_rate = 0;
}

/**
 * @brief Populates the local site with a wild mosquito population.
 * @param[in] initial_WJ The number of initial wild juveniles.
 * @param[in] initial_WM The number of initial wild males.
 * @param[in] initial_WV The number of initial wild virgin (unmated) females.
 * @param[in] initial_WF The number of initial wild mated females.
 */
void Patch::populate(int initial_WJ, int initial_WM, int initial_WV, int initial_WF) 
{
	for (int a=0; a < constants::max_dev + 1; ++a) {
		J[0][a] += initial_WJ;
		J[3][a] += initial_WJ;
	}	
	M[0] = initial_WM;
//	V[0] = initial_WV;
	F[0][0] = initial_WF;

	update_comp();
	update_mate();
}

/**
 * @brief Returns the coordinates of the patch. 
 * @see Point
 */
Point Patch::get_coords() const
{
	return coords;
}

/**
 * @brief Returns the number of males in the patch, divided by genotype. 
 */
std::array<long long int, constants::num_gen> Patch::get_M() const
{
	return M;
}

/**
 * @brief Returns the number of females in the patch, divided by female genotype and male sperm genotype.
 */
std::array<std::array<long long int, constants::num_gen>, constants::num_gen> Patch::get_F() const
{
	return F;
}

/**
 * @brief Calculates the total number of juveniles in the patch.
 * @return The total number of juveniles in the patch across all ages and genotypes. 
 */
long long int Patch::calculate_tot_J() 
{
	long long int tot_J = 0;
	for (int i = 0; i < 2*constants::num_gen; ++i) {
		for (int a = 0; a < constants::max_dev+1; ++a) {
			tot_J += J[i][a];
		}
	}
	return tot_J;
}

/**
 * @brief Calculates the total number of males in the patch.
 * @return The total number of males in the patch across all genotypes.
 */
long long int Patch::calculate_tot_M() 
{
	long long int tot_M = 0;
	for (const auto& m : M) {
		tot_M += m;
	}
	return tot_M;
}

/**
 * @brief Calculates the total number of virgin (unmated) females in the patch.
 * @return The total number of virgin (unmated) females in the patch across all genotypes.
 */
long long int Patch::calculate_tot_V()
{
	long long int tot_V = 0;
	for (const auto& v : V) {
		tot_V += v;
	}
	return tot_V;
}

/**
 * @brief Calculates the total number of mated females in the patch.
 * @return The total number of mated females in the patch across all female and male sperm genotypes.
 */
long long int Patch::calculate_tot_F()
{
	long long int tot_F = 0;
	for (int i = 0; i < constants::num_gen; ++i) {
		for (int j = 0; j < constants::num_gen; ++j) {
			tot_F += F[i][j];
		}
	}
	return tot_F;
}

/**
 * @brief Removes males from the patch as they disperse out. 
 * @param[in] m_out The number of males dispersing out. 
 */
void Patch::M_disperse_out(const std::array<long long int, constants::num_gen> &m_out) 
{
	for (std::size_t i = 0; i < m_out.size(); ++i) {
		M[i] -= m_out[i];
	}
	update_mate();
}

/**
 * @brief Removes mated females from the patch as they disperse out. 
 * @param[in] f_out The number of mated females dispersing out. 
 */
void Patch::F_disperse_out(const std::array<std::array<long long int, constants::num_gen>, constants::num_gen> &f_out) 
{
	for (std::size_t i = 0; i < f_out.size(); ++i) {
		for (std::size_t j = 0; j < f_out[0].size(); ++j) {
			F[i][j] -= f_out[i][j];
		}
	}
}

/**
 * @brief Introduces new males into the patch as they disperse in. 
 * @param[in] gen	The genotype index for the males to disperse into M.
 * @param[in] m_in 	The number of males dispersing in for the given genotype.
 */
void Patch::M_disperse_in(int gen, long long int m_in) 
{
	M[gen] += m_in;
	update_mate();
}

/**
 * @brief Introduces new mated females into the patch as they disperse in. 
 * @param[in] f_gen		The female genotype index for the females to disperse into F.
 * @param[in] m_gen 	The male sperm genotype index for the females to disperse into F.
 * @param[in] f_disp 	The number of females dispersing in for the given genotype combination. 
 */
void Patch::F_disperse_in(int f_gen, int m_gen, long long int f_disp) 
{
	F[f_gen][m_gen] += f_disp;
}

/**
 * @brief Removes those females from the patch that attempt to go into aestivation.
 * @param[in] f_try The number of mated females that attempt to go into aestivation.
 */
void Patch::F_hide(const std::array<std::array<long long int, constants::num_gen>, constants::num_gen> &f_try)
{
	for (std::size_t i = 0; i < F.size(); ++i) {
		for (std::size_t j = 0; j < F[0].size(); ++j) {
			F[i][j] -= f_try[i][j];
		}
	}
}

/**
 * @brief Introduces back into the patch those females that wake up from aestivation.
 * @param[in] f_wake The number of mated females from the given patch that wake up from aestivation.
 */
void Patch::F_wake(const std::array<std::array<long long int, constants::num_gen>, constants::num_gen> &f_wake)
{
	for (std::size_t i = 0; i < F.size(); ++i) {
		for (std::size_t j = 0; j < F[0].size(); ++j) {
			F[i][j] += f_wake[i][j];
		}
	}
}

/**
 * @brief Introduces driver heterozygous males into the patch.
 * @param[in] num_driver_M The number of driver heterozygous males to introduce.
 */
void Patch::add_driver_M(int num_driver_M) 
{
	M[1] += num_driver_M;
	update_mate();
}

/**
 * @brief Ages the juveniles of different age groups in the patch by a day.
 * @note Not all juveniles survive the aging process.
 * @details The number of surviving individuals in an age group (for a given genotype) is determined by a binomial distibution of juvenile survival probability. 
 */
void Patch::juv_get_older() 
{
	for (int i=0; i < 2*constants::num_gen; ++i) {
		for (int a=0; a < constants::max_dev; ++a) {
			// number of juveniles that survive aging by a day are placed into the new older age group	
			J[i][a] = random_binomial(J[i][a+1], comp);
		}
		J[i][constants::max_dev] = 0; // youngest ones have all aged by one day so none left in this age group
	}
}

/**
 * @brief Removes dying adults from the patch.
 * @details Determines the number of adults that die in the given day and removes them from the patch. The number of males that die (for a given genotype) is determined by a binomial distribution of adult mortality, and similarly for virgin and mated females.
 * @see InputParams::mu_a
 */
void Patch::adults_die()
{
	double mu_a = params->mu_a;
	for (int i=0; i < constants::num_gen; ++i) {
		long long int m = random_binomial(M[i], mu_a); // number of males that die
		M[i] -= m;

		long long int v = random_binomial(V[i], mu_a);
		V[i] -= v;	

		for (int j=0; j < constants::num_gen; ++j) {
			long long int f = random_binomial(F[i][j], mu_a);
			F[i][j] -= f;
		}
	}

	update_mate();
}

/**
 * @brief Mates a fraction of the virgin females in the patch.
 * @details Determines the number of virgin females that mate in the given day (for a given female genotype) with a male of genotype j by using a binomial distribution of the mating rate. Then, tranforms the virgin females into mated females carrying male sperm of genotype j. Mating is carried out within the patch, and females only mate once. 
 */
void Patch::virgins_mate() 
{
	std::array<long long int, constants::num_gen> v;
	std::vector<long long int> v_c;
	for (int i=0; i < constants::num_gen; ++i) {
		v[i] = random_binomial(V[i], mate_rate); // how many V will mate
		if (v[i] > 0) {
			v_c = random_multinomial(v[i], M); // how many V with given genotype will carry each of the male genotypes
			for (int j=0; j < constants::num_gen; j++) {
				F[i][j] += v_c[j];
			}
			V[i] -= v[i];
		}
	}
	
}

// 
/**
 * @brief Females lay eggs, creating new juveniles. 
 * @details Determines the number of eggs laid with genotype k on the given day by using a Poisson distribution. Other relevant parameters include the egg laying rate. Determines the development duration of these new juveniles using a multinomial distribution of development duration probabilities. 
 * @param[in] inher_fraction 		inheritance fraction for new offspring
 * @param[in] dev_duration_probs 	probabilities for juvenile development duration of new offspring
 * @see Simulation::set_inheritance(), InputParams::theta
 */
void Patch::lay_eggs(const std::array<std::array<std::array <double, 2*constants::num_gen>, constants::num_gen>, constants::num_gen> &inher_fraction,
 const std::array<double, constants::max_dev+1> &dev_duration_probs)
{
	std::vector<long long int> j_new;
	for (int i=0; i < constants::num_gen; ++i) {
		for (int j=0; j < constants::num_gen; ++j) {
			for (int k=0; k < 2*constants::num_gen; ++k) {
				double num = (params->theta) * F[i][j] * inher_fraction[i][j][k]; // expected number of eggs laid with k genotype
				long long int eggs = random_poisson(num); // actual number of eggs laid sampled from random distribution

				j_new = random_multinomial(eggs, dev_duration_probs); // number of eggs that start in each different age class (according to different juvenile development times)
//if(i==0&&j==0)std::cout<<" lay eggs   "<<i<<"  "<<j<<"  "<<k<<"   "<<inher_fraction[i][j][k]<<"     "<<eggs<<std::endl;
				for (int t=0; t < constants::max_dev + 1; ++t) { // juveniles created with assigned remaining time to develop
					J[k][t] += j_new[t];
				}
			}
		}
	}

	update_comp();
}

/**
 * @brief Turns the oldest juveniles into adults.
 * @note Not all juveniles survive eclosion. 
 * @details The number of survivors is determined by a binomial distribution of the juvenile survival probability. Sex determination upon eclosion is determined by a binomial distribution with 0.5 probability. 
 */
void Patch::juv_eclose() 
{
	for (int i=0; i < constants::num_gen; ++i) {
		long long int surv = random_binomial(J[i][0], comp); // number of juveniles that survive eclosion
		J[i][0] = 0; // all the oldest juveniles either successfully eclose or die
		if (surv > 0) {	
//if(i==0)std::cout<<" females eclose   "<<i<<"     "<<surv<<std::endl;

			// roughly half of the juveniles become male and half female following a distribution
			V[i] += surv;
		}
	}
	for (int i=constants::num_gen; i < 2*constants::num_gen; ++i) {
		long long int surv = random_binomial(J[i][0], comp); // number of juveniles that survive eclosion
		J[i][0] = 0; // all the oldest juveniles either successfully eclose or die
		if (surv > 0) {	
//if(i==3)std::cout<<" males eclose   "<<i<<"     "<<surv<<std::endl;
			// roughly half of the juveniles become male and half female following a distribution
			M[i-constants::num_gen] += surv; 
		}
	}
	update_comp();
	update_mate();
}

/**
 * @brief Updates the juvenile survival probability of the patch. 
 * @details Relevant parameters include the juvenile mortality rate, the juvenile survival probability power and the carrying capacity. The survival probability is computed as 
 * \f[
 *  \textrm{comp} = (1 - \mu_j)    \sqrt[\textrm{comp_power}]{\frac{\alpha}{\alpha + \textrm{J}_{\textrm{tot}}}},
 * \f]
 * with \f$ \mu_j \f$ juvenile mortality rate and \f$ \alpha \f$ carrying capacity.
 * @see Seasonality::alpha(), InputParams::mu_j, InputParams::comp_power
 */
void Patch::update_comp()
{
	int d = model->get_day();
	double alpha = model->get_alpha(alpha0);
	alpha*=humandens;
//	std::cout<<humandens<<std::endl;
	long long int tot_J = calculate_tot_J();
	comp = (1 - (params->mu_j)) * std::pow(alpha / (alpha + tot_J), params->comp_power);
}

/**
 * @brief Updates the mating rate of the patch.
 * @details Relevant parameters include the beta parameter. The mating rate is computed as
 * \f[ 
 *   \textrm{mate_rate} =  \frac{\textrm{M}_{\textrm{tot}}}{\beta + \textrm{M}_{\textrm{tot}}}.
 * \f]
 * @see InputParams::beta
 */
void Patch::update_mate()
{
	long long int tot_M = calculate_tot_M();
	mate_rate = tot_M / ((params->beta) + tot_M);
}
