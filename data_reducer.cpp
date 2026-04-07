#include "data_reducer.h"

// TODO:
// - only look at temperature during precession time (need relevant timestamp)
// - add rf, sf, hv, files

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
	// Set flags for the logger for error/output messages
	// ---------------------------------------------------------
	SimpleLog_Setup(NULL, NULL, 0, 0, 0, "\t");
	SimpleLog_FilterLevel(SL_ERROR/*|SL_NOTICE SL_ALL*/); // Default is SL_ALL
	
	// ---------------------------------------------------------
	// Build the B0 field using gradient expansion and map coefficients
	// ---------------------------------------------------------
	// Load the gradient field map
	map<string,GradientInfo> FieldMap = Gradients("../../dataset/fieldmap/001_updated_optimized_all_B0_2022_2025.csv");
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
	// Initialize ROOT file
	// ---------------------------------------------------------
	TFile * outFile = new TFile(output_filename.c_str(), "RECREATE");
	TTree * outTree = new TTree("reduced_data","Reduced data from EDM run");
	int out_Run 		= this_run;
	int out_Cycle 		= 0;
	double out_Temperature 	= DUMMY_VAL;
	double out_B_Hg_Top	= DUMMY_VAL;
	double out_B_Hg_Bot	= DUMMY_VAL;
	double out_B_Hg_Top_Err	= DUMMY_VAL;
	double out_B_Hg_Bot_Err	= DUMMY_VAL;
	double out_Ucn_Top	= DUMMY_VAL;
	double out_Ucn_Bot	= DUMMY_VAL;
	double out_A_Top	= DUMMY_VAL;
	double out_A_Bot	= DUMMY_VAL;

	outTree->Branch("Cycle"		,&out_Cycle		,"Cycle/I"		);
	outTree->Branch("Run"		,&out_Run		,"Run/I"		);
	outTree->Branch("Temperature"	,&out_Temperature	,"Temperature/D"	);
	outTree->Branch("B_Hg_Top"	,&out_B_Hg_Top		,"B_Hg_Top/D"		);
	outTree->Branch("B_Hg_Bot"	,&out_B_Hg_Bot		,"B_Hg_Bot/D"		);
	outTree->Branch("B_Hg_Top_Err"	,&out_B_Hg_Top_Err	,"B_Hg_Top_Err/D"	);
	outTree->Branch("B_Hg_Bot_Err"	,&out_B_Hg_Bot_Err	,"B_Hg_Bot_Err/D"	);
	outTree->Branch("Ucn_Top"	,&out_Ucn_Top		,"Ucn_Top/D"		);
	outTree->Branch("Ucn_Bot"	,&out_Ucn_Bot		,"Ucn_Bot/D"		);
	outTree->Branch("A_Top"		,&out_A_Top		,"A_Top/D"		);
	outTree->Branch("A_Bot"		,&out_A_Bot		,"A_Bot/D"		);




	// ---------------------------------------------------------
	// Loop over the cycles:
	// ---------------------------------------------------------
	for( int current_cycle : available_cycles ){
		std::string cycle_str = base_run_dir + run_str + "_" + formatNumber(current_cycle, 6); // "000020"

		// Reset branch variables to DUMMY_VAL for new event
		out_Temperature 	= DUMMY_VAL;
		out_B_Hg_Top		= DUMMY_VAL;
		out_B_Hg_Bot		= DUMMY_VAL;
		out_B_Hg_Top_Err	= DUMMY_VAL;
		out_B_Hg_Bot_Err	= DUMMY_VAL;
		out_Ucn_Top		= DUMMY_VAL;
		out_Ucn_Bot		= DUMMY_VAL;
		out_A_Top		= DUMMY_VAL;
		out_A_Bot		= DUMMY_VAL;

		// Construct file paths based on naming convention
		std::string temperature_file 	= cycle_str + "_000_temperature_000.hd";
		std::string hg_file   		= cycle_str + "_000_onlineAna_hgm_000.hd";
		std::string ucn_file  		= cycle_str + "_000_onlineAna_ucn_000.hd";
		std::string csm_file  		= cycle_str + "_000_csm_000.hd";
		
		// Run reductions 
		out_Temperature		= ReduceTemperature(	temperature_file	);

		HgResult out_hg		= ReduceHg(		hg_file			);
		out_B_Hg_Top		= out_hg.B_Hg_Top;
		out_B_Hg_Bot		= out_hg.B_Hg_Bot;
		out_B_Hg_Top_Err	= out_hg.B_Hg_Top_Err;
		out_B_Hg_Bot_Err	= out_hg.B_Hg_Bot_Err;

		UcnResult out_ucn 	= ReduceUcn(		ucn_file		);
		out_Ucn_Top		= out_ucn.Ucn_Top;
		out_Ucn_Bot		= out_ucn.Ucn_Bot;
		out_A_Top		= out_ucn.A_Top;
		out_A_Bot		= out_ucn.A_Bot;
		ReduceCsm(	csm_file	);
		



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



