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
	double Hg_Delta_Top	= DUMMY_VAL;
	double Hg_Delta_Bot	= DUMMY_VAL;
	double Hg_Delta_Top_Err	= DUMMY_VAL;
	double Hg_Delta_Bot_Err	= DUMMY_VAL;
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

	// Create the detuning phase using Hg for UCN analysis
	double Delta_Top	= fHg_Top * gammaN/gammaHg;
	double Delta_Bot	= fHg_Bot * gammaN/gammaHg;
	double Delta_Top_Err	= fHg_Top_Err * gammaN/gammaHg;
	double Delta_Bot_Err	= fHg_Bot_Err * gammaN/gammaHg;

	// convert frequency to field (pT)
	fHg_Top *= 1000000. / gammaHg;		// Hz/muT --> pT
	fHg_Bot *= 1000000. / gammaHg;		// Hz/muT --> pT
	fHg_Top_Err *= 1000000. / gammaHg;	// Hz/muT --> pT
	fHg_Bot_Err *= 1000000. / gammaHg;	// Hz/muT --> pT

	// Check for NaN
	if( fHg_Top != fHg_Top ){
		fHg_Top 	= DUMMY_VAL;
		fHg_Top_Err 	= DUMMY_VAL;
		Delta_Top	= DUMMY_VAL;
		Delta_Top_Err	= DUMMY_VAL;
	}
	if( fHg_Bot != fHg_Bot ){
		fHg_Bot 	= DUMMY_VAL;
		fHg_Bot_Err 	= DUMMY_VAL;
		Delta_Bot	= DUMMY_VAL;
		Delta_Bot_Err	= DUMMY_VAL;
	}
	thisHgEvent.B_Hg_Top = fHg_Top;
	thisHgEvent.B_Hg_Bot = fHg_Bot;
	thisHgEvent.B_Hg_Top_Err = fHg_Top_Err;
	thisHgEvent.B_Hg_Bot_Err = fHg_Bot_Err;
	thisHgEvent.Hg_Delta_Top 	= Delta_Top;
	thisHgEvent.Hg_Delta_Bot 	= Delta_Bot;
	thisHgEvent.Hg_Delta_Top_Err 	= Delta_Top_Err;
	thisHgEvent.Hg_Delta_Bot_Err 	= Delta_Bot_Err;

	N2_ClearConfig(&N2data);
	return thisHgEvent;
}



#endif
