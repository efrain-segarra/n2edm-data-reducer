#pragma once
#ifndef __ReduceUcn_H__

#include <filesystem>
#include "constants.h"

namespace fs = std::filesystem;
using namespace std;

// ---------------------------------------------------------
// Ucn reducer
// ---------------------------------------------------------
struct UcnResult {
	double Ucn_Top = DUMMY_VAL;
	double Ucn_Bot = DUMMY_VAL;
	double A_Top = DUMMY_VAL;
	double A_Bot = DUMMY_VAL;
};
UcnResult ReduceUcn(const std::string& filepath ){
	UcnResult thisUcnEvent;

	// missing file flag
	if (!fs::exists(filepath)) {
		cerr << "Could not find file " << filepath << "\n";
		return thisUcnEvent; 
	}


	// Initialize and check ucn data
	tN2data N2data = {0};
	N2_ReadFile(filepath.c_str(), &N2data);
	if( N2data.NbRow != 1 ){
		cerr << "Unexpected onlineAna_ucn file size of " << N2data.NbRow << "\n";
		return thisUcnEvent;
	}

	// Grab ucn counts and asymmetry
	double Ucn_Top 		= ((uint64_t**)N2data.Data)[0][6];
	double Ucn_Bot 		= ((uint64_t**)N2data.Data)[0][7];
	double A_Top 		= ((double**)N2data.Data)[0][9];
	double A_Bot 		= ((double**)N2data.Data)[0][10];
	thisUcnEvent.Ucn_Top	= Ucn_Top;
	thisUcnEvent.Ucn_Bot	= Ucn_Bot;
	thisUcnEvent.A_Top	= A_Top;
	thisUcnEvent.A_Bot	= A_Bot;

	return thisUcnEvent;
}

#endif
