#pragma once
#ifndef __ReduceHg_H__

#include <filesystem>
#include "constants.h"

namespace fs = std::filesystem;
using namespace std;

// ---------------------------------------------------------
// Hg reducer
// ---------------------------------------------------------
struct HgResult {
	double B_Hg_Top 	= DUMMY_VAL;
	double B_Hg_Bot 	= DUMMY_VAL;
	double B_Hg_Top_Err 	= DUMMY_VAL;
	double B_Hg_Bot_Err	= DUMMY_VAL;
};
inline HgResult ReduceHg(const std::string& filepath ){
	HgResult thisHgEvent;

	// missing file flag
	if (!fs::exists(filepath)) {
		cerr << "Could not find file " << filepath << "\n";
		return thisHgEvent; 
	}


	// Initialize and check hg data
	tN2data N2data = {0};
	N2_ReadFile(filepath.c_str(), &N2data);
	if( N2data.NbRow != 1 ){
		cerr << "Unexpected onlineAna_hg file size of " << N2data.NbRow << "\n";
		return thisHgEvent;
	}

	// Grab Hg frequencies
	double fHg_Top 		= ((double**)N2data.Data)[0][9];
	double fHg_Bot 		= ((double**)N2data.Data)[0][10];
	double fHg_Top_Err 	= ((double**)N2data.Data)[0][11];
	double fHg_Bot_Err 	= ((double**)N2data.Data)[0][12];

	// convert to field (pT)
	fHg_Top *= 1000000. / 7.5901152;
	fHg_Bot *= 1000000. / 7.5901152;
	fHg_Top_Err *= 1000000. / 7.5901152;
	fHg_Bot_Err *= 1000000. / 7.5901152;

	// Check for NaN
	if( fHg_Top != fHg_Top ){
		fHg_Top 	= DUMMY_VAL;
		fHg_Top_Err 	= DUMMY_VAL;
	}
	if( fHg_Bot != fHg_Bot ){
		fHg_Bot 	= DUMMY_VAL;
		fHg_Bot_Err 	= DUMMY_VAL;
	}
	thisHgEvent.B_Hg_Top = fHg_Top;
	thisHgEvent.B_Hg_Bot = fHg_Bot;
	thisHgEvent.B_Hg_Top_Err = fHg_Top_Err;
	thisHgEvent.B_Hg_Bot_Err = fHg_Bot_Err;

	return thisHgEvent;
}



#endif
