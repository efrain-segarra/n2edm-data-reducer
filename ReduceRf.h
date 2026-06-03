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

	for( int r=0; r<N2data.NbRow; r++ ){
		double start = ((double**)N2data.Data)[r][1];
		double dur = ((double**)N2data.Data)[r][2];
		double freq = ((double**)N2data.Data)[r][3];
		
		// Assume Hg frequencies are below 10 Hz and only 1 of these pulses
		if( freq < 10 ){
			thisRfEvent.Rf_Hg_Start = start;
			thisRfEvent.Rf_Hg_Duration = dur;
			thisRfEvent.Rf_Hg_Freq = freq;
		}
		else if( freq > 20 && freq < 35 ){
			// If we do not yet have a UCN pulse saved,
			// put it in the first pulse channel:
			if( thisRfEvent.Rf_Ucn1_Start == DUMMY_VAL ){
				thisRfEvent.Rf_Ucn1_Start = start;
				thisRfEvent.Rf_Ucn1_Duration = dur;
				thisRfEvent.Rf_Ucn1_Freq = freq;
			}
			else{ // otherwise put it in the second pulse channel:
				thisRfEvent.Rf_Ucn2_Start = start;
				thisRfEvent.Rf_Ucn2_Duration = dur;
				thisRfEvent.Rf_Ucn2_Freq = freq;
				UCN_Free_Precession_Start 	= thisRfEvent.Rf_Ucn1_Start + thisRfEvent.Rf_Ucn1_Duration; // after ucn1 rf pulse
				UCN_Free_Precession_Stop   	= thisRfEvent.Rf_Ucn2_Start;				    // at start of ucn2 rf pulse
			}
		}
	}
	
												    //
	N2_ClearConfig(&N2data);
	return thisRfEvent;
}


#endif
