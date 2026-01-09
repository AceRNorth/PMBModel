#include <filesystem>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "Record.h"
#include "globals.h"

/**
 * Record constructor.
 * Creates LocalData, Totals and CoordinateList output .txt files.
 * @details Creates a subdirectory for output files in the current directory.
 * @param[in] rec_params 	recording parameters
 * @param[in] rep 			initial repetition label for the given set of runs
 */
Record::Record(RecordParams* rec_params, int rep) 
{
	rec_start = rec_params->rec_start;
	rec_end = rec_params->rec_end;
	rec_interval_global = rec_params->rec_interval_global;
	rec_interval_local = rec_params->rec_interval_local;
	rec_sites_freq = rec_params->rec_sites_freq;
	set_label = rec_params->set_label;
	rep_label = rep;

	// create folder for destination of output files 
	if (!std::filesystem::exists("./output_files")) {
		std::filesystem::create_directory("output_files");
	}
	std::filesystem::current_path("./output_files");
	
	os1 << "LocalData" << set_label << "run" << rep_label << ".txt"; 
	local_data.open(os1.str());
	os2 << "Totals" << set_label << "run" << rep_label << ".txt";
	global_data.open(os2.str());
	os3 << "CoordinateList" << set_label << "run" << rep_label << ".txt";
	coord_list.open(os3.str());

		//local_data << day << "\t" << pat+1 << "\t" << patch_type[pat] << "\t" << count << "\t" << freq<<std::endl;
	local_data << "Population info at each site\n";
	local_data << "Day" << "\t" << "Site" << "\t" << "patch type" << "\t" << "Fww" << "\t" << "Fwz" << "\t" << "Fzz"<< "\t" << "Mww" << "\t" << "Mwz" << "\t" << "Mzz"<<std::endl;

	global_data << "Total number of adults in each class\n";
	global_data << "Day" << "\t" << "Fww" << "\t" << "Fwz" << "\t" << "Fzz"<< "\t" << "Mww" << "\t" << "Mwz" << "\t" << "Mzz"<<std::endl;

	coord_list << "Coordinate list of the sites\n";
	coord_list << "Site" << "\t" << "x" << "\t" << "y" << std::endl;
}

/**
 * Record destructor.
 * Resets the current filepath so the output_files directory can be found in the next set of runs. 
 */
Record::~Record()
{
	std::filesystem::current_path("..");
}

/**
 * @brief Records the coordinates of the population sites. 
 * @details Relevant parameters include the fraction of sites to collect data for.
 * @param[in] sites vector of all Patch objects
 * @see InputParams::rec_sites_freq
 */
void Record::record_coords(const std::vector<Patch*> &sites) 
{
	const auto default_precision{std::cout.precision()};
	constexpr auto max_precision{std::numeric_limits<double>::digits10 + 1};
	for (int pat=0; pat < sites.size(); pat += rec_sites_freq) {
		auto coords = sites[pat]->get_coords();
		coord_list << pat+1 << "\t" << std::setprecision(max_precision) << coords.x << "\t" << coords.y << std::endl;
	}
	coord_list << std::setprecision(default_precision);
}

/**
 * @brief Records the total numbers of male mosquitoes for the given day, divided by genotype. 
 * @details The totals are assumed to be across all sites. 
 * @param[in] day 		simulation day
 * @param[in] tot_M_gen total number of males divided by genotype
 * @see Model::calculate_tot_M_gen(), Patch::get_M()
 */
void Record::record_global(int day, const std::array<long long int, constants::num_gen> &tot_M_gen,const std::array<long long int, constants::num_gen> &tot_F_gen)
{
	global_data << day;
	for (const auto& f_gen : tot_F_gen) {
		global_data << "\t" << f_gen;
	}
	for (const auto& m_gen : tot_M_gen) {
		global_data << "\t" << m_gen;
	}
	global_data << std::endl;
}
/*void Record::record_global(int day, const std::array<long long int, constants::num_gen> &tot_M_gen)
{
	global_data << day;
	for (const auto& m_gen : tot_M_gen) {
		global_data << "\t" << m_gen;
	}
	global_data << std::endl;
}*/

/**
 * @brief Outputs the total numbers of juvenile (J), male (M), virgin female (V) and mated female (F) mosquitoes for the given day.
 * @details The totals are assumed to be across all sites, and over all genotypes and age groups. 
 * @param[in] day 	simulation day
 * @param[in] tot_J	total number of juveniles
 * @param[in] tot_M	total number of males
 * @param[in] tot_V	total number of virgin (unmated) females
 * @param[in] tot_F	total number of mated females
 * @see Patch
 */
void Record::output_totals(int day, long long int tot_J, long long int tot_M, long long int tot_V, long long int tot_F)
{
	if (day == 0) {
		std::cout << "\n" << "rep " << rep_label << "\n";
		std::cout << "day" << "   " << "total J" << "   " << "total M" << "   " << "total V" << "   " << "total F" << "\n";
	}
	std::cout << day << "     " << tot_J << "   " << tot_M << "   " << tot_V << "   " << tot_F << std::endl;
}

/**
 * @brief Records the number of males at each site for the given day.
 * @details The number of males at each site is divided by genotype. Relevant parameters include the fraction of sites to collect data for.
 * @param[in] day 	simulation day
 * @param[in] sites vector of all Patch objects
 * @see InputParams::rec_sites_freq
 */
void Record::record_local(int day, const std::vector<Patch*> &sites) 
{
	int count,tot,gd,ind;
	int f0,f1,f2;
	for (int pat=0; pat < sites.size(); pat ++) 
	{
		if(patch_type[pat]>OutputType)
		{
		const auto& F = sites[pat]->get_F();
		const auto& M = sites[pat]->get_M();

		f0= F[0][0]+ F[0][1]+ F[0][2];
		f1= F[1][0]+ F[1][1]+ F[1][2];
		f2= F[2][0]+ F[2][1]+ F[2][2];
		local_data << day << "\t" << pat+1 << "\t" << patch_type[pat] << "\t" << f0 << "\t" << f1 <<  "\t" <<f2<< "\t" <<M[0]<< "\t" <<M[1]<< "\t" <<M[2]<<std::endl;
		}
	}
}

/**
 * @brief Determines if it is time to record global data.
 * @details The number of males at each site is divided by genotype.  
 * @param[in] day 	simulation day
 * @return As you would expect.
 * @see Record::record_global()
 */
bool Record::is_rec_global_time(int day)
{
	return day % rec_interval_global == 0;
}

/**
 * @brief Determines if it is time to record local data.
 * @note The initialisation day (day 0) will always be recorded, and the recording window will be inclusive of the start and end times. 
 * @details Other relevant parameters include the local recording interval. 
 * @param[in] day 	simulation day
 * @return As you would expect.
 * @see Record::record_local(), InputParams::rec_start, InputParams::rec_end, InputParams::rec_interval_local
 */
bool Record::is_rec_local_time(int day) 
{
	return (day == 0) || (day >= rec_start && day <= rec_end && day % rec_interval_local == 0);
}
