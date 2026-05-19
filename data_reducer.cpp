#include "data_reducer.h"

// TODO:
// - add hv
// - add csm
// - add trim coil information


int main(int argc, char** argv ){
	// ---------------------------------------------------------
	// Check arguments
	// ---------------------------------------------------------
	if( argc != 2 ){
		cerr << "Unexpected number of arguments. Instead use:" << endl;
		cerr << "\t./data_reducer [RUN NUMBER]" << endl;
		return 1;
	}

	// ---------------------------------------------------------
	// Read config file to figure out which subsystems we analyze
	// ---------------------------------------------------------
	TEnv config;
	if( config.ReadFile("configuration.cfg", kEnvUser) < 0) {
		cerr << "Could not find or read the configuration file in current directory!" << endl;
		cerr << "Make sure you are running the executable from the directory containing the config file" << endl;
		return 1;
	}
	bool do_rf	= config.GetValue("Process.Rf", 0);
	bool do_sf	= config.GetValue("Process.Sf", 0);
	bool do_hg	= config.GetValue("Process.Hg", 0);
	bool do_ucn	= config.GetValue("Process.Ucn", 0);
	bool do_csm	= config.GetValue("Process.Csm", 0);
	bool do_temp	= config.GetValue("Process.Temperature", 0);


	// ---------------------------------------------------------
	// Set flags for the logger for error/output messages
	// ---------------------------------------------------------
	SimpleLog_Setup(NULL, NULL, 0, 0, 0, "\t");
	SimpleLog_FilterLevel(SL_ERROR/*|SL_NOTICE SL_ALL*/); // Default is SL_ALL
	
	// ---------------------------------------------------------
	// Build the B0 field using gradient expansion and map coefficients
	// ---------------------------------------------------------
	// Load the gradient field map
	/*
	map<string,GradientInfo> FieldMap = Gradients("../../dataset/fieldmap/001_updated_optimized_all_B0_2022_2025.csv");
	// Calculate the fields at the Cs sensor locations
	for( int ch=0; ch < out_CsFieldCalc.size(); ch++){
		Bvector pos = CsCellPosition( ch+1 ); // daq channel starts at 1, not 0
		Bvector B0 = CalculateField( FieldMap, pos.rho, pos.phi, pos.z, "up");
	
		out_CsX.at(ch) = pos.rho * std::cos( pos.phi );
		out_CsY.at(ch) = pos.rho * std::sin( pos.phi );
		out_CsZ.at(ch) = pos.z;
		out_CsFieldCalc.at(ch) = sqrt(B0.rho*B0.rho + B0.phi*B0.phi + B0.z*B0.z);
	}
	*/


	// ---------------------------------------------------------
	// Load specific run from user
	// ---------------------------------------------------------
	int this_run = atoi(argv[1]);


	// ---------------------------------------------------------
	// Temporary overload for Hg analysis due to issues with online
	// analysis files during Dec data-taking
	// ---------------------------------------------------------
	if( this_run < 8009 ){
		LoadHgFallbackCsv("/Users/efrainsegarra/work/n2EDM/projects/dec2025_analysis/dataset/multi_run_plot_data.csv", this_run);
	}

	// ---------------------------------------------------------
	// Format strings and get base directory of run:
	// ---------------------------------------------------------
	std::string run_str = formatNumber(this_run, 6);	// "008003"
	std::string dir_part1 = run_str.substr(0, 3);		// "008"
	std::string dir_part2 = run_str.substr(3, 3);		// "003"
	//std::string base_run_dir = "/xdata/n2edmdata/" + dir_part1 + "/" + dir_part2 + "/";
	std::string base_run_dir = "/Users/efrainsegarra/work/n2EDM/projects/dec2025_analysis/dataset/" + string(argv[1]) + "/";
	std::string output_filename = "run_" + run_str + "_reduced.root";


	// ---------------------------------------------------------
	// Check for base directory:
	// ---------------------------------------------------------
	if (!fs::exists(base_run_dir)) {
		cerr << "Error: Run directory does not exist" << endl;
		return 1;
	}


	// ---------------------------------------------------------
	// Scan directory to find unique cycles:
	// ---------------------------------------------------------
	std::set<int> available_cycles;
	std::string run_prefix = run_str + "_";			// "008003_"

	for (const auto& entry : fs::directory_iterator(base_run_dir)) {
		if (entry.is_regular_file()) {
			std::string fname = entry.path().filename().string();

			// Check if file starts with our run prefix and is long enough
			if (fname.find(run_prefix) == 0 && fname.length() >= 13) {
				// Extract the cycle part (e.g., from "008003_000020...", extract "000020")
				std::string cycle_str = fname.substr(7, 6);
				try {
					available_cycles.insert(std::stoi(cycle_str));
				} catch (...) {
					// Ignore files that don't match the strict integer naming convention
				}
			}
		}
	}
	cout << "Found " << available_cycles.size() << " unique cycles in run " << this_run << endl;


	// ---------------------------------------------------------
	// Loop over the cycles and perform reductions but do not fill tree yet because
	// we need to fit the full cycle data to do visibility fit and THEN go back to 
	// calculate UCN frequency
	// ---------------------------------------------------------
	//	Cycle data in memory
	std::vector<CycleData> runbuffer;
	// 	containers for visibility fit later
	std::vector<double> xT;
	std::vector<double> xTerr;
	std::vector<double> yT;
	std::vector<double> yTerr;
	std::vector<double> zT;
	std::vector<double> zTerr;
	std::vector<double> xB;
	std::vector<double> xBerr;
	std::vector<double> yB;
	std::vector<double> yBerr;
	std::vector<double> zB;
	std::vector<double> zBerr;
	for( int current_cycle : available_cycles ){
		// Initialize struct for cycle data
		CycleData cycledata;

		// Format string for cycle reading
		std::string cycle_str = base_run_dir + run_str + "_" + formatNumber(current_cycle, 6); // "000020"
		
		// ---------------------------------------------------------
		// RF :
		// ---------------------------------------------------------
		if( do_rf ){
			// Construct file paths based on naming convention
			std::string rf_file		= cycle_str + "_000_rfgen_000.hd";
			// Run reduction
			RfResult out_rf			= ReduceRf( 		rf_file 		);
			cycledata.Rf_Hg_Start		= out_rf.Rf_Hg_Start;
			cycledata.Rf_Hg_Duration	= out_rf.Rf_Hg_Duration;
			cycledata.Rf_Hg_Freq		= out_rf.Rf_Hg_Freq;
			cycledata.Rf_Ucn1_Start	= out_rf.Rf_Ucn1_Start;
			cycledata.Rf_Ucn1_Duration	= out_rf.Rf_Ucn1_Duration;
			cycledata.Rf_Ucn1_Freq		= out_rf.Rf_Ucn1_Freq;
			cycledata.Rf_Ucn2_Start	= out_rf.Rf_Ucn2_Start;
			cycledata.Rf_Ucn2_Duration	= out_rf.Rf_Ucn2_Duration;
			cycledata.Rf_Ucn2_Freq		= out_rf.Rf_Ucn2_Freq;
		}

		// Flag for requiring RF due to dependencies
		if( UCN_Free_Precession_Start == DUMMY_VAL || UCN_Free_Precession_Stop == DUMMY_VAL ){
			cout << "No RF file found. Cannot continue with reducer due to later dependencies. Skipping cycle...\n";
			continue;
		}

		// ---------------------------------------------------------
		// Temp :
		// ---------------------------------------------------------
		if( do_temp ){
			// Construct file paths based on naming convention
			std::string temperature_file 	= cycle_str + "_000_temperature_000.hd";
			// Run reduction
			cycledata.Temperature		= ReduceTemperature(	temperature_file	);
		}

		// ---------------------------------------------------------
		// SF :
		// ---------------------------------------------------------
		if( do_sf ){
			// Construct file paths based on naming convention
			std::string sf_file	= cycle_str + "_000_sflipper_000.hd";
			// Run reduction
			SfResult out_sf		= ReduceSf(		sf_file			);
			cycledata.Sf_Top	= out_sf.Sf_Top;
			cycledata.Sf_Bot	= out_sf.Sf_Bot;
			cycledata.Sf_1		= out_sf.Sf_1;
			cycledata.Sf_2		= out_sf.Sf_2;
			cycledata.Sf_3		= out_sf.Sf_3;
			cycledata.Sf_4		= out_sf.Sf_4;
		}

		// ---------------------------------------------------------
		// Hg :
		// ---------------------------------------------------------
		if( do_hg ){
			HgResult out_hg;
			if( this_run < 8009 ){
				out_hg	= ReduceHgCsv( current_cycle );
			}
			else{
				// Construct file paths based on naming convention
				std::string hg_file   		= cycle_str + "_000_onlineAna_hgm_000.hd";
				// Run reduction
				out_hg			= ReduceHg(		hg_file			);
			}
			cycledata.B_Hg_Top		= out_hg.B_Hg_Top;
			cycledata.B_Hg_Bot		= out_hg.B_Hg_Bot;
			cycledata.B_Hg_Top_Err		= out_hg.B_Hg_Top_Err;
			cycledata.B_Hg_Bot_Err		= out_hg.B_Hg_Bot_Err;
			cycledata.Hg_Delta_Top		= out_hg.Hg_Delta_Top;
			cycledata.Hg_Delta_Bot		= out_hg.Hg_Delta_Bot;
			cycledata.Hg_Delta_Top_Err	= out_hg.Hg_Delta_Top_Err;
			cycledata.Hg_Delta_Bot_Err	= out_hg.Hg_Delta_Bot_Err;

		}

		// ---------------------------------------------------------
		// Ucn :
		// ---------------------------------------------------------
		if( do_ucn ){
			// Construct file paths based on naming convention
			std::string ucn_file  	= cycle_str + "_000_livecount_000.hd";
			// Run reduction
			UcnResult out_ucn 	= ReduceUcn(		ucn_file	);
			cycledata.Ucn_Top	= out_ucn.Ucn_Top;
			cycledata.Ucn_Bot	= out_ucn.Ucn_Bot;
			cycledata.A_Top	= out_ucn.A_Top;
			cycledata.A_Bot	= out_ucn.A_Bot;
			cycledata.A_Top_Err	= out_ucn.A_Top_Err;
			cycledata.A_Bot_Err	= out_ucn.A_Bot_Err;
		}

		// ---------------------------------------------------------
		// Combined Ucn and Hg :
		// ---------------------------------------------------------
		if( do_ucn && do_hg ){
			if( cycledata.Rf_Ucn1_Freq 	== DUMMY_VAL ) continue;
			if( cycledata.Hg_Delta_Top 	== DUMMY_VAL ) continue;
			if( cycledata.Hg_Delta_Top_Err 	== DUMMY_VAL ) continue;
			if( cycledata.Sf_Top 	   	== DUMMY_VAL ) continue;
			if( cycledata.A_Top 	  	== DUMMY_VAL ) continue;
			if( cycledata.A_Top_Err	   	== DUMMY_VAL ) continue;
			if( cycledata.B_Hg_Top		< 900000     ) continue; // safety flag on Hg for now

			xT	 	.push_back( cycledata.Rf_Ucn1_Freq - cycledata.Hg_Delta_Top );
			xTerr		.push_back( cycledata.Hg_Delta_Top_Err );
			yT		.push_back( cycledata.Sf_Top );
			yTerr		.push_back( 0. );
			zT	 	.push_back( cycledata.A_Top );
			zTerr	 	.push_back( cycledata.A_Top_Err );

			xB	 	.push_back( cycledata.Rf_Ucn1_Freq - cycledata.Hg_Delta_Bot );
			xBerr		.push_back( cycledata.Hg_Delta_Bot_Err );
			yB		.push_back( cycledata.Sf_Bot );
			yBerr		.push_back( 0. );
			zB	 	.push_back( cycledata.A_Bot );
			zBerr	 	.push_back( cycledata.A_Bot_Err );
		}

		// Save cycle data in memory
		cycledata.Cycle = current_cycle;
		cycledata.Run	 = this_run;
		runbuffer.push_back(cycledata);
	}

	// ---------------------------------------------------------
	// Perform UCN visibility and frequency analysis
	// ---------------------------------------------------------
	
	// Define graphs for the combined fit
	TGraph2DErrors * gr2d_top = new TGraph2DErrors( xT.size(), xT.data(), yT.data(), zT.data(), 
							   xTerr.data(), yTerr.data(), zTerr.data() );
	TGraph2DErrors * gr2d_bot = new TGraph2DErrors( xB.size(), xB.data(), yB.data(), zB.data(), 
							   xBerr.data(), yBerr.data(), zBerr.data() );
	
	// Define function to use
	TF2* ramsey_fit = new TF2("ramsey_fit","-[0] * cos( (TMath::Pi() / [4]) * (x - [1]) ) + [2]*((1 - y)/2.0) + [3]*((1 + y)/2.0)",-0.005, 0.005,-1.5, 1.5);

	ramsey_fit->SetParNames("alpha", "Phi", "delta_01", "delta_10", "delta_nu");
	double guess_alpha = 0.85;
	double guess_fn = 0.;
	double guess_phi_m1 = 0.0;
	double guess_phi_p1 = 0.0;
	double delta_nu = 1./(2*180. + 8*2/TMath::Pi());

	// Fit top chamber and get values
	ramsey_fit->SetParameters(guess_alpha, guess_fn, guess_phi_m1, guess_phi_p1, delta_nu);
	ramsey_fit->FixParameter(4, delta_nu);
	gr2d_top->Fit(ramsey_fit);
	double Vis_Top 		= ramsey_fit->GetParameter(0);
	double Vis_Top_Err 	= ramsey_fit->GetParError(0);
	double Phi_M1_Top 	= ramsey_fit->GetParameter(2);
	double Phi_M1_Top_Err 	= ramsey_fit->GetParError(2);
	double Phi_P1_Top 	= ramsey_fit->GetParameter(3);
	double Phi_P1_Top_Err 	= ramsey_fit->GetParError(3);

	// Fit bottom chamber and get values
	ramsey_fit->SetParameters(guess_alpha, guess_fn, guess_phi_m1, guess_phi_p1, delta_nu);
	ramsey_fit->FixParameter(4, delta_nu);
	gr2d_bot->Fit(ramsey_fit);
	double Vis_Bot 		= ramsey_fit->GetParameter(0);
	double Vis_Bot_Err 	= ramsey_fit->GetParError(0);
	double Phi_M1_Bot 	= ramsey_fit->GetParameter(2);
	double Phi_M1_Bot_Err 	= ramsey_fit->GetParError(2);
	double Phi_P1_Bot 	= ramsey_fit->GetParameter(3);
	double Phi_P1_Bot_Err 	= ramsey_fit->GetParError(3);


	// ---------------------------------------------------------
	// Fill the TTree from the memory struct
	// ---------------------------------------------------------
	TFile * outFile = new TFile(output_filename.c_str(), "RECREATE");
	TTree * outTree = new TTree("reduced_data","Reduced data from EDM run");
	int out_Run 		= this_run;
	int out_Cycle 		= 0;
	double out_Rf_Hg_Start		= DUMMY_VAL;
	double out_Rf_Hg_Duration	= DUMMY_VAL;
	double out_Rf_Hg_Freq		= DUMMY_VAL;
	double out_Rf_Ucn1_Start	= DUMMY_VAL;
	double out_Rf_Ucn1_Duration	= DUMMY_VAL;
	double out_Rf_Ucn1_Freq		= DUMMY_VAL;
	double out_Rf_Ucn2_Start	= DUMMY_VAL;
	double out_Rf_Ucn2_Duration	= DUMMY_VAL;
	double out_Rf_Ucn2_Freq		= DUMMY_VAL;
	double out_Temperature 		= DUMMY_VAL;
	double out_Sf_Top		= DUMMY_VAL;
	double out_Sf_Bot		= DUMMY_VAL;
	double out_Sf_1			= DUMMY_VAL;
	double out_Sf_2			= DUMMY_VAL;
	double out_Sf_3			= DUMMY_VAL;
	double out_Sf_4			= DUMMY_VAL;
	double out_B_Hg_Top		= DUMMY_VAL;
	double out_B_Hg_Bot		= DUMMY_VAL;
	double out_B_Hg_Top_Err		= DUMMY_VAL;
	double out_B_Hg_Bot_Err		= DUMMY_VAL;
	double out_Hg_Delta_Top		= DUMMY_VAL;
	double out_Hg_Delta_Bot		= DUMMY_VAL;
	double out_Hg_Delta_Top_Err	= DUMMY_VAL;
	double out_Hg_Delta_Bot_Err	= DUMMY_VAL;
	double out_Ucn_Top		= DUMMY_VAL;
	double out_Ucn_Bot		= DUMMY_VAL;
	double out_A_Top		= DUMMY_VAL;
	double out_A_Bot		= DUMMY_VAL;
	double out_A_Top_Err		= DUMMY_VAL;
	double out_A_Bot_Err		= DUMMY_VAL;
	double out_Fn_Top		= DUMMY_VAL;
	double out_Fn_Bot		= DUMMY_VAL;
	double out_Fn_Top_Err		= DUMMY_VAL;
	double out_Fn_Bot_Err		= DUMMY_VAL;
	outTree->Branch("Cycle"			,&out_Cycle		,"Cycle/I"		);
	outTree->Branch("Run"			,&out_Run		,"Run/I"		);
	if( do_rf ){
		outTree->Branch("Rf_Hg_Start"		,&out_Rf_Hg_Start	,"Rf_Hg_Start/D"	);
		outTree->Branch("Rf_Hg_Duration"	,&out_Rf_Hg_Duration	,"Rf_Hg_Duration/D"	);
		outTree->Branch("Rf_Hg_Freq"		,&out_Rf_Hg_Freq	,"Rf_Hg_Freq/D"		);
		outTree->Branch("Rf_Ucn1_Start"		,&out_Rf_Ucn1_Start	,"Rf_Ucn1_Start/D"	);
		outTree->Branch("Rf_Ucn1_Duration"	,&out_Rf_Ucn1_Duration	,"Rf_Ucn1_Duration/D"	);
		outTree->Branch("Rf_Ucn1_Freq"		,&out_Rf_Ucn1_Freq	,"Rf_Ucn1_Freq/D"	);
		outTree->Branch("Rf_Ucn2_Start"		,&out_Rf_Ucn2_Start	,"Rf_Ucn2_Start/D"	);
		outTree->Branch("Rf_Ucn2_Duration"	,&out_Rf_Ucn2_Duration	,"Rf_Ucn2_Duration/D"	);
		outTree->Branch("Rf_Ucn2_Freq"		,&out_Rf_Ucn2_Freq	,"Rf_Ucn2_Freq/D"	);
	}
	if( do_temp ){
		outTree->Branch("Temperature"		,&out_Temperature	,"Temperature/D"	);
	}
	if( do_sf ){
		outTree->Branch("Sf_Top"		,&out_Sf_Top		,"Sf_Top/D"		);
		outTree->Branch("Sf_Bot"		,&out_Sf_Bot		,"Sf_Bot/D"		);
		outTree->Branch("Sf_1"			,&out_Sf_1		,"Sf_1/D"		);
		outTree->Branch("Sf_2"			,&out_Sf_2		,"Sf_2/D"		);
		outTree->Branch("Sf_3"			,&out_Sf_3		,"Sf_3/D"		);
		outTree->Branch("Sf_4"			,&out_Sf_4		,"Sf_4/D"		);
	}
	if( do_hg ){
		outTree->Branch("B_Hg_Top"		,&out_B_Hg_Top		,"B_Hg_Top/D"		);
		outTree->Branch("B_Hg_Bot"		,&out_B_Hg_Bot		,"B_Hg_Bot/D"		);
		outTree->Branch("B_Hg_Top_Err"		,&out_B_Hg_Top_Err	,"B_Hg_Top_Err/D"	);
		outTree->Branch("B_Hg_Bot_Err"		,&out_B_Hg_Bot_Err	,"B_Hg_Bot_Err/D"	);
		outTree->Branch("Hg_Delta_Top"		,&out_Hg_Delta_Top	,"Hg_Delta_Top/D"	);
		outTree->Branch("Hg_Delta_Bot"		,&out_Hg_Delta_Bot	,"Hg_Delta_Bot/D"	);
		outTree->Branch("Hg_Delta_Top_Err"	,&out_Hg_Delta_Top_Err	,"Hg_Delta_Top_Err/D"	);
		outTree->Branch("Hg_Delta_Bot_Err"	,&out_Hg_Delta_Bot_Err	,"Hg_Delta_Bot_Err/D"	);
	}
	if( do_ucn ){
		outTree->Branch("Ucn_Top"		,&out_Ucn_Top		,"Ucn_Top/D"		);
		outTree->Branch("Ucn_Bot"		,&out_Ucn_Bot		,"Ucn_Bot/D"		);
		outTree->Branch("A_Top"			,&out_A_Top		,"A_Top/D"		);
		outTree->Branch("A_Bot"			,&out_A_Bot		,"A_Bot/D"		);
		outTree->Branch("A_Top_Err"		,&out_A_Top_Err		,"A_Top_Err/D"		);
		outTree->Branch("A_Bot_Err"		,&out_A_Bot_Err		,"A_Bot_Err/D"		);
	}
	if( do_hg && do_ucn ){
		outTree->Branch("Fn_Top"		,&out_Fn_Top		,"Fn_Top/D"		);
		outTree->Branch("Fn_Bot"		,&out_Fn_Bot		,"Fn_Bot/D"		);
		outTree->Branch("Fn_Top_Err"		,&out_Fn_Top_Err	,"Fn_Top_Err/D"		);
		outTree->Branch("Fn_Bot_Err"		,&out_Fn_Bot_Err	,"Fn_Bot_Err/D"		);
	}

	// Loop over all the cycle data to put into ROOT file
	for (const auto& cycledata : runbuffer) {
		// ---------------------------------------------------------
		// RF :
		// ---------------------------------------------------------
		if( do_rf ){
			out_Rf_Hg_Start		= cycledata.Rf_Hg_Start;
			out_Rf_Hg_Duration	= cycledata.Rf_Hg_Duration;
			out_Rf_Hg_Freq		= cycledata.Rf_Hg_Freq;
			out_Rf_Ucn1_Start	= cycledata.Rf_Ucn1_Start;
			out_Rf_Ucn1_Duration	= cycledata.Rf_Ucn1_Duration;
			out_Rf_Ucn1_Freq	= cycledata.Rf_Ucn1_Freq;
			out_Rf_Ucn2_Start	= cycledata.Rf_Ucn2_Start;
			out_Rf_Ucn2_Duration	= cycledata.Rf_Ucn2_Duration;
			out_Rf_Ucn2_Freq	= cycledata.Rf_Ucn2_Freq;
		}

		// ---------------------------------------------------------
		// Temp :
		// ---------------------------------------------------------
		if( do_temp ){
			out_Temperature		= cycledata.Temperature; 
		}

		// ---------------------------------------------------------
		// SF :
		// ---------------------------------------------------------
		if( do_sf ){
			out_Sf_Top		= cycledata.Sf_Top;
			out_Sf_Bot		= cycledata.Sf_Bot;
			out_Sf_1		= cycledata.Sf_1;
			out_Sf_2		= cycledata.Sf_2;
			out_Sf_3		= cycledata.Sf_3;
			out_Sf_4		= cycledata.Sf_4;
		}

		// ---------------------------------------------------------
		// Hg :
		// ---------------------------------------------------------
		if( do_hg ){
			out_B_Hg_Top		= cycledata.B_Hg_Top;
			out_B_Hg_Bot		= cycledata.B_Hg_Bot;
			out_B_Hg_Top_Err	= cycledata.B_Hg_Top_Err;
			out_B_Hg_Bot_Err	= cycledata.B_Hg_Bot_Err;
			out_Hg_Delta_Top	= cycledata.Hg_Delta_Top;
			out_Hg_Delta_Bot	= cycledata.Hg_Delta_Bot;
			out_Hg_Delta_Top_Err	= cycledata.Hg_Delta_Top_Err;
			out_Hg_Delta_Bot_Err	= cycledata.Hg_Delta_Bot_Err;
		}

		// ---------------------------------------------------------
		// Ucn :
		// ---------------------------------------------------------
		if( do_ucn ){
			out_Ucn_Top		= cycledata.Ucn_Top;
			out_Ucn_Bot		= cycledata.Ucn_Bot;
			out_A_Top		= cycledata.A_Top;
			out_A_Bot		= cycledata.A_Bot;
			out_A_Top_Err		= cycledata.A_Top_Err;
			out_A_Bot_Err		= cycledata.A_Bot_Err;
		}

		// ---------------------------------------------------------
		// Additional UCN+Hg manipulation to do UCN frequency extraction
		// ---------------------------------------------------------
		if( do_hg && do_ucn ){
			// Figure out which side of the ramsey fringe we are to do the 
			// correct sign [RF +/- delta_nu * acos(...)]
			double sign_top = 1;
			if( cycledata.Rf_Ucn1_Freq - cycledata.Hg_Delta_Top < 0 ) sign_top = 1;
			if( cycledata.Rf_Ucn1_Freq - cycledata.Hg_Delta_Top > 0 ) sign_top = -1;
			double sign_bot = 1;
			if( cycledata.Rf_Ucn1_Freq - cycledata.Hg_Delta_Bot < 0 ) sign_bot = 1;
			if( cycledata.Rf_Ucn1_Freq - cycledata.Hg_Delta_Bot > 0 ) sign_bot = -1;
			out_Fn_Top = cycledata.Rf_Ucn1_Freq + sign_top * delta_nu/TMath::Pi() 
				* std::acos( 
					( Phi_M1_Top*(1-cycledata.Sf_Top)/2. + Phi_P1_Top*(1+cycledata.Sf_Top)/2. 
					  - cycledata.A_Top )
					/ Vis_Top 
				);
			out_Fn_Bot = cycledata.Rf_Ucn1_Freq + sign_top * delta_nu/TMath::Pi() 
				* std::acos( 
					( Phi_M1_Bot*(1-cycledata.Sf_Bot)/2. + Phi_P1_Bot*(1+cycledata.Sf_Bot)/2. 
					  - cycledata.A_Bot )
					/ Vis_Bot 
				);
		}

		// Save to tree
		out_Cycle 	= cycledata.Cycle;
		out_Run		= cycledata.Run;
		outTree->Fill();
	}
	



	// ---------------------------------------------------------
	// Write and cleanup
	// ---------------------------------------------------------
	outFile->cd();
	gr2d_top->Write("Top_Ramsey");
	gr2d_bot->Write("Bot_Ramsey");
	outTree->Write(); // also deletes outTree
	outFile->Close();
	delete outFile;

	cout << "Reduction done, save output to " << output_filename << endl;

	return 0;
}



