#include "data_reducer.h"

// TODO:
// - add hv
// - add csm
// - add livecount
// - perform joint fit to extract fn per cycle

int main(int argc, char** argv ){
	// ---------------------------------------------------------
	// Check arguments
	// ---------------------------------------------------------
	if( argc != 2 ){
		cerr << "Unexpected number of arguments. Instead use:" << endl;
		cerr << "\t./data_reducer [RUN NUMBER]" << endl;
		return 1;
	}

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
	// Build the B0 field using gradient expansion and map coefficients
	// ---------------------------------------------------------
	// Load the gradient field map
	//map<string,GradientInfo> FieldMap = Gradients("../../dataset/fieldmap/001_updated_optimized_all_B0_2022_2025.csv");
	/*
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
	// Format strings and get base directory of run:
	// ---------------------------------------------------------
	std::string run_str = formatNumber(this_run, 6);	// "008003"
	std::string dir_part1 = run_str.substr(0, 3);		// "008"
	std::string dir_part2 = run_str.substr(3, 3);		// "003"
	std::string base_run_dir = "/xdata/n2edmdata/" + dir_part1 + "/" + dir_part2 + "/";
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
	// Initialize ROOT file
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
	double out_Ucn_Top		= DUMMY_VAL;
	double out_Ucn_Bot		= DUMMY_VAL;
	double out_A_Top		= DUMMY_VAL;
	double out_A_Bot		= DUMMY_VAL;
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
	}
	if( do_ucn ){
		outTree->Branch("Ucn_Top"		,&out_Ucn_Top		,"Ucn_Top/D"		);
		outTree->Branch("Ucn_Bot"		,&out_Ucn_Bot		,"Ucn_Bot/D"		);
		outTree->Branch("A_Top"			,&out_A_Top		,"A_Top/D"		);
		outTree->Branch("A_Bot"			,&out_A_Bot		,"A_Bot/D"		);
	}


	// ---------------------------------------------------------
	// Loop over the cycles:
	// ---------------------------------------------------------
	for( int current_cycle : available_cycles ){
		std::string cycle_str = base_run_dir + run_str + "_" + formatNumber(current_cycle, 6); // "000020"

		// ---------------------------------------------------------
		// RF :
		// ---------------------------------------------------------
		if( do_rf ){
			out_Rf_Hg_Start		= DUMMY_VAL;
			out_Rf_Hg_Duration	= DUMMY_VAL;
			out_Rf_Hg_Freq		= DUMMY_VAL;
			out_Rf_Ucn1_Start	= DUMMY_VAL;
			out_Rf_Ucn1_Duration	= DUMMY_VAL;
			out_Rf_Ucn1_Freq	= DUMMY_VAL;
			out_Rf_Ucn2_Start	= DUMMY_VAL;
			out_Rf_Ucn2_Duration	= DUMMY_VAL;
			out_Rf_Ucn2_Freq	= DUMMY_VAL;
			// Construct file paths based on naming convention
			std::string rf_file		= cycle_str + "_000_rfgen_000.hd";
			// Run reduction
			RfResult out_rf		= ReduceRf( 		rf_file 		);
			out_Rf_Hg_Start		= out_rf.Rf_Hg_Start;
			out_Rf_Hg_Duration	= out_rf.Rf_Hg_Duration;
			out_Rf_Hg_Freq		= out_rf.Rf_Hg_Freq;
			out_Rf_Ucn1_Start	= out_rf.Rf_Ucn1_Start;
			out_Rf_Ucn1_Duration	= out_rf.Rf_Ucn1_Duration;
			out_Rf_Ucn1_Freq	= out_rf.Rf_Ucn1_Freq;
			out_Rf_Ucn2_Start	= out_rf.Rf_Ucn2_Start;
			out_Rf_Ucn2_Duration	= out_rf.Rf_Ucn2_Duration;
			out_Rf_Ucn2_Freq	= out_rf.Rf_Ucn2_Freq;
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
			out_Temperature 	= DUMMY_VAL;
			// Construct file paths based on naming convention
			std::string temperature_file 	= cycle_str + "_000_temperature_000.hd";
			// Run reduction
			out_Temperature		= ReduceTemperature(	temperature_file	);
		}

		// ---------------------------------------------------------
		// SF :
		// ---------------------------------------------------------
		if( do_sf ){
			out_Sf_Top		= DUMMY_VAL;
			out_Sf_Bot		= DUMMY_VAL;
			out_Sf_1		= DUMMY_VAL;
			out_Sf_2		= DUMMY_VAL;
			out_Sf_3		= DUMMY_VAL;
			out_Sf_4		= DUMMY_VAL;
			// Construct file paths based on naming convention
			std::string sf_file		= cycle_str + "_000_sflipper_000.hd";
			// Run reduction
			SfResult out_sf		= ReduceSf(		sf_file			);
			out_Sf_Top		= out_sf.Sf_Top;
			out_Sf_Bot		= out_sf.Sf_Bot;
			out_Sf_1		= out_sf.Sf_1;
			out_Sf_2		= out_sf.Sf_2;
			out_Sf_3		= out_sf.Sf_3;
			out_Sf_4		= out_sf.Sf_4;
		}

		// ---------------------------------------------------------
		// Hg :
		// ---------------------------------------------------------
		if( do_hg ){
			out_B_Hg_Top		= DUMMY_VAL;
			out_B_Hg_Bot		= DUMMY_VAL;
			out_B_Hg_Top_Err	= DUMMY_VAL;
			out_B_Hg_Bot_Err	= DUMMY_VAL;
			// Construct file paths based on naming convention
			std::string hg_file   		= cycle_str + "_000_onlineAna_hgm_000.hd";
			// Run reduction
			HgResult out_hg		= ReduceHg(		hg_file			);
			out_B_Hg_Top		= out_hg.B_Hg_Top;
			out_B_Hg_Bot		= out_hg.B_Hg_Bot;
			out_B_Hg_Top_Err	= out_hg.B_Hg_Top_Err;
			out_B_Hg_Bot_Err	= out_hg.B_Hg_Bot_Err;
		}

		// ---------------------------------------------------------
		// Ucn :
		// ---------------------------------------------------------
		if( do_ucn ){
			out_Ucn_Top		= DUMMY_VAL;
			out_Ucn_Bot		= DUMMY_VAL;
			out_A_Top		= DUMMY_VAL;
			out_A_Bot		= DUMMY_VAL;
			// Construct file paths based on naming convention
			std::string ucn_file  		= cycle_str + "_000_onlineAna_ucn_000.hd";
			// Run reduction
			UcnResult out_ucn 	= ReduceUcn(		ucn_file		);
			out_Ucn_Top		= out_ucn.Ucn_Top;
			out_Ucn_Bot		= out_ucn.Ucn_Bot;
			out_A_Top		= out_ucn.A_Top;
			out_A_Bot		= out_ucn.A_Bot;
		}
		

		// Save to tree
		out_Cycle 	= current_cycle;
		out_Run		= this_run;
		outTree->Fill();
	}


	// ---------------------------------------------------------
	// Write and cleanup
	// ---------------------------------------------------------
	outFile->cd();
	outTree->Write(); // also deletes outTree
	outFile->Close();
	delete outFile;

	cout << "Reduction done, save output to " << output_filename << endl;

	return 0;
}



