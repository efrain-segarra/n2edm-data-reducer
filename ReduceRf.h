#pragma once
#ifndef __ReduceRf_H__

#include <filesystem>
#include "constants.h"

namespace fs = std::filesystem;
using namespace std;

// Constants to be set from the RF file and then used by subsequent analyses
inline double UCN_Free_Precession_Start 	= DUMMY_VAL;
inline double UCN_Free_Precession_Stop   	= DUMMY_VAL;

// ---------------------------------------------------------
// Rf reducer
// ---------------------------------------------------------
struct RfResult{
	double Rf_Hg_Start 		= DUMMY_VAL;
	double Rf_Hg_Duration 		= DUMMY_VAL;
	double Rf_Hg_Freq		= DUMMY_VAL;
	double Rf_Ucn1_Start		= DUMMY_VAL;
	double Rf_Ucn1_Duration		= DUMMY_VAL;
	double Rf_Ucn1_Freq		= DUMMY_VAL;
	double Rf_Ucn2_Start		= DUMMY_VAL;
	double Rf_Ucn2_Duration		= DUMMY_VAL;
	double Rf_Ucn2_Freq		= DUMMY_VAL;
};
inline RfResult ReduceRf(const std::string& filepath ){
	RfResult thisRfEvent;

	// missing file flag
	if (!fs::exists(filepath)) {
		cerr << "Could not find file " << filepath << "\n";
		return thisRfEvent; 
	}


	// Initialize and check temperature data
	tN2data N2data = {0};
	N2_ReadFile(filepath.c_str(), &N2data);
	if( N2data.NbRow == 0 ){
		cerr << "Unexpected rf file size of " << N2data.NbRow << "\n";
		return thisRfEvent;
	}
	
	thisRfEvent.Rf_Hg_Start 	= ((double**)N2data.Data)[0][1];
	thisRfEvent.Rf_Hg_Duration 	= ((double**)N2data.Data)[0][2];
	thisRfEvent.Rf_Hg_Freq 		= ((double**)N2data.Data)[0][3];

	thisRfEvent.Rf_Ucn1_Start 	= ((double**)N2data.Data)[1][1];
	thisRfEvent.Rf_Ucn1_Duration 	= ((double**)N2data.Data)[1][2];
	thisRfEvent.Rf_Ucn1_Freq 	= ((double**)N2data.Data)[1][3];

	thisRfEvent.Rf_Ucn2_Start 	= ((double**)N2data.Data)[2][1];
	thisRfEvent.Rf_Ucn2_Duration 	= ((double**)N2data.Data)[2][2];
	thisRfEvent.Rf_Ucn2_Freq 	= ((double**)N2data.Data)[2][3];
	
	UCN_Free_Precession_Start 	= thisRfEvent.Rf_Ucn1_Start + thisRfEvent.Rf_Ucn1_Duration; // after ucn1 rf pulse
	UCN_Free_Precession_Stop   	= thisRfEvent.Rf_Ucn2_Start;				    // at start of ucn2 rf pulse
												    //
	N2_ClearConfig(&N2data);
	return thisRfEvent;
}


#endif
