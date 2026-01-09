#ifndef GENERALMETAPOP_RECORD_H
#define GENERALMETAPOP_RECORD_H

#include <sstream>
#include <fstream>
#include <vector>
#include <array>
#include "constants.h"
#include "Patch.h"
#include "Params.h"

/** 
 * Records model data.
 * Creates CoordinateList, Totals and LocalData output .txt files for each repetition of the model.
 * CoordinateList: contains the coordinates of all recorded sites.
 * Totals: contains the total numbers of males of each genotype across all sites, for each recorded day of the simulation.
 * LocalData: contains the number of males of each genotype for each site on each simulated day that is recorded. 
 * @see RecordParams
 */ 
class Record {
public:
	Record(RecordParams *rec_params, int rep);
	~Record();
	void record_coords(const std::vector<Patch*> &sites);
	void record_global(int day, const std::array<long long int, constants::num_gen> &tot_M_gen, const std::array<long long int, constants::num_gen> &tot_F_gen);
	void output_totals(int day, long long int tot_J, long long int tot_M, long long int tot_V, long long int tot_F);
	void record_local(int day, const std::vector<Patch*> &sites);

	bool is_rec_local_time(int day);
	bool is_rec_global_time(int day);

private:
	// recording window and intervals

	int rec_start; /**< Start time for the data recording window (in days) (inclusive) */ 
	int rec_end; /**< End time for the data recording window (in days) (inclusive) */
	int rec_interval_global; /**< Time interval for global data recording (in days) */
	int rec_interval_local; /**< Time interval for local data recording (in days) */
	int rec_sites_freq; /**< Fraction of sites to collect local data for (1 is all sites, 10 is 1 in 10 etc.). */

	// output filename labels

	int set_label; /**< 'Set of repetitions' index label for output files. */
	int rep_label; /**< 'Repetition' index label for output files (within given set). */

	std::ostringstream os1; /**< Local data output filename */
	std::ostringstream os2; /**< Global data output filename */
	std::ostringstream os3; /**< Coordinate data output filename */
	std::ofstream local_data; /**< Local data filestream object */
	std::ofstream global_data; /**< Global data filestream object */
	std::ofstream coord_list; /**< Coordinate data filestream object */
};

#endif //GENERALMETAPOP_RECORD_H
