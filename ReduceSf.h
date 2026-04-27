#pragma once
#ifndef __ReduceSf_H__

#include <filesystem>
#include "constants.h"

namespace fs = std::filesystem;
using namespace std;

// ---------------------------------------------------------
// Sf reducer
// ---------------------------------------------------------
struct SfResult {
	double Sf_1		= DUMMY_VAL; // det1 == bot red
	double Sf_2		= DUMMY_VAL; // det2 == bot blue
	double Sf_3		= DUMMY_VAL; // det3 == top yellow
	double Sf_4		= DUMMY_VAL; // det4 == top white
	double Sf_Top 		= DUMMY_VAL; // 1==det3(ON)+det4(OFF) | 0==det3(OFF)+det4(ON)
	double Sf_Bot 		= DUMMY_VAL; // 1==det1(ON)+det2(OFF) | 0==det1(OFF)+det2(ON)
};
inline SfResult ReduceSf(const std::string& filepath ){
	SfResult thisSfEvent;

	// missing file flag
	if (!fs::exists(filepath)) {
		cerr << "Could not find file " << filepath << "\n";
		return thisSfEvent; 
	}

	// Initialize and check sf data
	tN2data N2data = {0};
	N2_ReadFile(filepath.c_str(), &N2data);
	if( N2data.NbRow == 0 ){
		cerr << "Unexpected sf file size of " << N2data.NbRow << "\n";
		return thisSfEvent; 
	}

	// Calculate average sf state in the cycle and return
	double bot_red = 0;
	double bot_blue = 0;
	double top_yellow = 0;
	double top_white = 0;
	int it = 0;
	for( int r=0; r<N2data.NbRow; r++ ){
		double timestamp = ((double**)N2data.Data)[r][0];

		if( timestamp > UCN_Free_Precession_Start && timestamp < UCN_Free_Precession_Stop ){
			bot_red 	+= ((double**)N2data.Data)[r][1];	// count1
			bot_blue 	+= ((double**)N2data.Data)[r][2];	// count2
			top_yellow 	+= ((double**)N2data.Data)[r][3];	// count3
			top_white	+= ((double**)N2data.Data)[r][4];	// count4
			it+=1;
		}
	}
	N2_ClearConfig(&N2data);
	bot_red /= it;
	bot_blue /= it;
	top_yellow /= it;
	top_white /= it;

	// BOTTOM 1(red)=ON + 2(blue)=OFF
	if( 	bot_red > 0.01 && bot_blue < 0.01){
		thisSfEvent.Sf_1 = 1;
		thisSfEvent.Sf_2 = -1;
		thisSfEvent.Sf_Bot = 1;
	}
	// BOTTOM 1(red)=OFF + 2(blue)=ON
	else if(bot_red < 0.01 && bot_blue > 0.01){
		thisSfEvent.Sf_1 = -1;
		thisSfEvent.Sf_2 = 1;
		thisSfEvent.Sf_Bot = -1;
	}
	// TOP 3(yellow)=ON + 4(white)=OFF
	if(	top_yellow > 0.01 && top_white < 0.01){
		thisSfEvent.Sf_3 = 1;
		thisSfEvent.Sf_4 = -1;
		thisSfEvent.Sf_Top = 1;
	}
	// TOP 3(yellow)=OFF + 4(white)=ON
	else if(top_yellow < 0.01 && top_white > 0.01){
		thisSfEvent.Sf_3 = -1;
		thisSfEvent.Sf_4 = 1;
		thisSfEvent.Sf_Top = -1;
	}

	return thisSfEvent;
}


#endif
