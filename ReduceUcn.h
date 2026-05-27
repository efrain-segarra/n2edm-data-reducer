#pragma once
#ifndef __ReduceUcn_H__

#include <filesystem>
#include "constants.h"

namespace fs = std::filesystem;
using namespace std;

// ---------------------------------------------------------
// Ucn reducer
// ---------------------------------------------------------
struct Detector {
	double counts		= DUMMY_VAL; 	// counts with bg subtracted
	double error		= DUMMY_VAL; 	// error with bg subtraction
	double bg_integral	= 0;		// bg integral counts
	double bg_slope 	= 0;		// bg slope
	double bg_intercept 	= 0;		// bg intercept
	double rss 		= 0;		// residual sum of squares
	double fit_quality	= 0;		// root-mean-square-error / expected poisson scatter
};
struct UcnResult {
	Detector Det[4];
	double Ucn_Top = DUMMY_VAL;
	double Ucn_Bot = DUMMY_VAL;
	double A_Top = DUMMY_VAL;
	double A_Bot = DUMMY_VAL;
	double A_Top_Err = DUMMY_VAL;
	double A_Bot_Err = DUMMY_VAL;
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
	if( N2data.NbRow == 0 ){
		cerr << "Unexpected livecount_ucn file size of " << N2data.NbRow << "\n";
		return thisUcnEvent;
	}


	double Sx = 0;
	double Sxx = 0;
	double Sy[4] = {0}, Sxy[4] = {0}, Signal[4] = {0};
	int N_Bg = 0;
	double T_bg_start = UCN_Free_Precession_Stop - 20;
	double T_bg_stop = UCN_Free_Precession_Stop  - 5; 
	double T_cnt_start = UCN_Free_Precession_Stop;  // Tfill + Tstore
	double T_cnt_stop = UCN_Counting_Stop;		// Tfill + Tstore + Tcount
	if( UCN_Free_Precession_Stop == 0 || UCN_Free_Precession_Stop == DUMMY_VAL  ){ 
		// direct shots so analyze different
		T_bg_start = 0;
		T_bg_stop = 0;
		T_cnt_start = 12;
		T_cnt_stop = UCN_Counting_Stop;
	}
	// Gather the full cycle information
	for( int r=0; r<N2data.NbRow; r++ ){
		double timestamp = ((double**)N2data.Data)[r][0];
		// Integrate background before stored ucn are released:
		if( timestamp > T_bg_start && timestamp < T_bg_stop ){
			N_Bg++;
			Sx 	+= timestamp;
			Sxx 	+= timestamp * timestamp;
			for( int det = 0; det < 4; det++){
				double c = ((uint64_t**)N2data.Data)[r][det+2];
				Sy[det] 	+= c;
				Sxy[det]	+= c * timestamp;
			}
		}
		// Integrate stored UCN counts
		if( timestamp > T_cnt_start && timestamp < T_cnt_stop ){
			for( int det = 0; det < 4; det++){
				double c = ((uint64_t**)N2data.Data)[r][det+2];
				Signal[det]	+= c;
			}
		}
        }


	double denom = (N_Bg * Sxx - Sx*Sx);
	for( int det = 0; det < 4; det++ ){
		thisUcnEvent.Det[det].bg_slope = (N_Bg * Sxy[det] - Sx * Sy[det]) / denom;
           	thisUcnEvent.Det[det].bg_intercept = (Sy[det] - thisUcnEvent.Det[det].bg_slope * Sx) / N_Bg;

            	// Background integral over signal window (Integral of At + B)
            	thisUcnEvent.Det[det].bg_integral = thisUcnEvent.Det[det].bg_slope * 0.5 * (T_cnt_stop*T_cnt_stop - T_cnt_start*T_cnt_start) 
                                 + thisUcnEvent.Det[det].bg_intercept * (T_cnt_stop - T_cnt_start);
		if( std::isnan(thisUcnEvent.Det[det].bg_integral) ){
			thisUcnEvent.Det[det].counts = Signal[det];
		}
		else{
			thisUcnEvent.Det[det].counts = Signal[det] - thisUcnEvent.Det[det].bg_integral;
		}
		

            	// Error Propagation (Approximate)
            	// Error = sqrt( SignalSum + (TimeRatio^2 * BackgroundSum) )
            	double time_ratio = (T_cnt_stop-T_cnt_start) / (T_bg_stop - T_bg_start);
            	thisUcnEvent.Det[det].error = std::sqrt(Signal[det] + std::pow(time_ratio, 2) * Sy[det]);
	}


	// Calculate residual squared sum
	for( int r=0; r<N2data.NbRow; r++ ){
		double timestamp = ((double**)N2data.Data)[r][0];
		if( timestamp > T_bg_start && timestamp < T_bg_stop ){
			for(int det = 0; det < 4; det++) {
				double expected = thisUcnEvent.Det[det].bg_slope * timestamp + thisUcnEvent.Det[det].bg_intercept;
				double actual = ((uint64_t**)N2data.Data)[r][det+2];
				thisUcnEvent.Det[det].rss += std::pow(actual - expected, 2);
			}
		}
	}
	for( int det = 0; det < 4; det++ ){
		double avg_bg_counts = Sy[det] / N_Bg;
		double expected_scatter = (avg_bg_counts > 0) ? std::sqrt(avg_bg_counts) : 1.0;
		double rmse = std::sqrt( thisUcnEvent.Det[det].rss / (N_Bg - 2 ) );
		thisUcnEvent.Det[det].fit_quality = rmse / expected_scatter;
	}


	// Now form asymmetries
	double Ucn_Top = thisUcnEvent.Det[2].counts + thisUcnEvent.Det[3].counts;
	double Ucn_Bot = thisUcnEvent.Det[0].counts + thisUcnEvent.Det[1].counts;
	double A_Top = DUMMY_VAL;
	double A_Bot = DUMMY_VAL;

	if( UCN_SF_3 == -1 && UCN_SF_4 == 1 ){
		A_Top = (thisUcnEvent.Det[2].counts - thisUcnEvent.Det[3].counts)/Ucn_Top;
	}
	else if( UCN_SF_3 == 1 && UCN_SF_4 == -1 ){
		A_Top = (thisUcnEvent.Det[3].counts - thisUcnEvent.Det[2].counts)/Ucn_Top;
	}

	if( UCN_SF_1 == 1 && UCN_SF_2 == -1 ){
		A_Bot = (thisUcnEvent.Det[1].counts - thisUcnEvent.Det[0].counts)/Ucn_Bot;
	}
	else if( UCN_SF_1 == -1 && UCN_SF_2 == 1 ){
		A_Bot = (thisUcnEvent.Det[0].counts - thisUcnEvent.Det[1].counts)/Ucn_Bot;
	}
	

	thisUcnEvent.Ucn_Top = Ucn_Top;
	thisUcnEvent.Ucn_Bot = Ucn_Bot;
	thisUcnEvent.A_Top = A_Top;
	thisUcnEvent.A_Bot = A_Bot;

	thisUcnEvent.A_Top_Err = (2.0 / std::pow(Ucn_Top,2) ) * 
		std::sqrt( std::pow( thisUcnEvent.Det[2].counts * thisUcnEvent.Det[3].error , 2 ) +
			   std::pow( thisUcnEvent.Det[3].counts * thisUcnEvent.Det[2].error , 2 ) 	);
	thisUcnEvent.A_Bot_Err = (2.0 / std::pow(Ucn_Bot,2) ) * 
		std::sqrt( std::pow( thisUcnEvent.Det[0].counts * thisUcnEvent.Det[1].error , 2 ) +
			   std::pow( thisUcnEvent.Det[1].counts * thisUcnEvent.Det[0].error , 2 ) 	);
	


	// Debugging
	/*
	for( int det = 0; det < 4; det++ ){
		cout << "Det " << det << "\n";
		cout << "\tCounts after bkg: " << thisUcnEvent.Det[det].counts << "\n";
		cout << "\tBkg integral: " << thisUcnEvent.Det[det].bg_integral << "\n";
		cout << "\tCounts error: " << thisUcnEvent.Det[det].error << "\n";
		cout << "\tBkg RSS: " << thisUcnEvent.Det[det].rss << "\n";
		cout << "\tBkg quality: " << thisUcnEvent.Det[det].fit_quality << "\n";
	}
	cout << "Sum of bot: " << thisUcnEvent.Det[0].counts + thisUcnEvent.Det[1].counts << "\n";
	cout << "Sum of top: " << thisUcnEvent.Det[2].counts + thisUcnEvent.Det[3].counts << "\n";
	*/

	N2_ClearConfig(&N2data);
	return thisUcnEvent;
}

#endif
