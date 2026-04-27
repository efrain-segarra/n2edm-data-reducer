#pragma once
#ifndef __ReduceTemperature_H__

#include <filesystem>
#include "constants.h"

namespace fs = std::filesystem;
using namespace std;

// ---------------------------------------------------------
// Temperature reducer
// ---------------------------------------------------------
inline double ReduceTemperature(const std::string& filepath ){

	// missing file flag
	if (!fs::exists(filepath)) {
		cerr << "Could not find file " << filepath << "\n";
		return DUMMY_VAL; 
	}

	// Select which temperature sensor to read
	int TEMPERATURE_COLUMN = 34;

	// Initialize and check temperature data
	tN2data N2data = {0};
	N2_ReadFile(filepath.c_str(), &N2data);
	if( N2data.NbRow == 0 ){
		cerr << "Unexpected temperature file size of " << N2data.NbRow << "\n";
		return DUMMY_VAL;
	}

	// Calculate average temperature in the cycle and return
	double temp_sum = 0.0;
	int it = 0;
	for( int r=0; r<N2data.NbRow; r++ ){
		double timestamp = ((double**)N2data.Data)[r][0];

		if( timestamp > UCN_Free_Precession_Start && timestamp < UCN_Free_Precession_Stop ){
			temp_sum += ((double**)N2data.Data)[r][TEMPERATURE_COLUMN];
			it++;
		}
	}
	temp_sum /= (double)it;
	N2_ClearConfig(&N2data);

	return temp_sum;
}

#endif
