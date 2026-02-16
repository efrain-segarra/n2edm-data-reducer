#include "data_reducer.h"


int main(int argc, char** argv ){
	// Set flags for the logger for error/output messages
	SimpleLog_Setup(NULL, NULL, 0, 0, 0, "\t");
	SimpleLog_FilterLevel(SL_ERROR/*|SL_NOTICE SL_ALL*/); // Default is SL_ALL

	// Create data structure to be used for EDM datafiles
	tN2data N2data = {0};

	// Initialize ROOT file
	std::unique_ptr<TFile> outFile( TFile::Open("test.root","RECREATE") );
	TTree * outTree = new TTree("reduced_data","Reduced data from EDM run");
	int out_Cycle;
	int out_Run;
	vector<double> out_CsField(16); 
	vector<double> out_CsX(16); 
	vector<double> out_CsY(16); 
	vector<double> out_CsZ(16); 
	vector<double> out_CsFieldCalc(16); 

	outTree->Branch("Cycle",	&out_Cycle	);
	outTree->Branch("Run",		&out_Run	);
	outTree->Branch("CsField",	&out_CsField	);
	outTree->Branch("CsX",		&out_CsX	);
	outTree->Branch("CsY",		&out_CsY	);
	outTree->Branch("CsZ",		&out_CsZ	);
	outTree->Branch("CsFieldCalc",	&out_CsFieldCalc);

	// Load the gradient field map
	map<string,GradientInfo> FieldMap = Gradients("../../dataset/fieldmap/001_updated_optimized_all_B0_2022_2025.csv");


	// Load specific run
	int this_run = 8003;

	// Calculate the Cs fields just once assuming mapping
	for( int ch=0; ch < out_CsFieldCalc.size(); ch++){
		Bvector pos = CsCellPosition( ch+1 ); // daq channel starts at 1, not 0
		Bvector B0 = CalculateField( FieldMap, pos.rho, pos.phi, pos.z, "up");
	
		out_CsX.at(ch) = pos.rho * std::cos( pos.phi );
		out_CsY.at(ch) = pos.rho * std::sin( pos.phi );
		out_CsZ.at(ch) = pos.z;
		out_CsFieldCalc.at(ch) = sqrt(B0.rho*B0.rho + B0.phi*B0.phi + B0.z*B0.z);
	}
	


	// Loop over all cycles in a run
	int start_cycle = 0;
	int max_cycle = 208;
	for( int this_cycle = start_cycle ; this_cycle <= max_cycle ; this_cycle++ ){

		// Create filename string given the cycle:
		string filename = std::format("../../dataset/8003/008003_{:06}_000_csm_000.hd",this_cycle);

		// Load data into structure
		N2_ReadFile(filename.c_str(), &N2data);

		// Create average array for frequencies
		vector<double> avg_ch_field(16);
		int it=0;

		// Loop over events/rows
		for( int r=0; r<N2data.NbRow; r++ ){
			double csm_timestamp = ((double**)N2data.Data)[r][0];

			// Hardcoded free precession window
			if( csm_timestamp > 47.2102 && csm_timestamp < 227.21 ){
			//if( csm_timestamp > 55. && csm_timestamp < 220. ){
				// Grab the frequencies for the 16 cells
				avg_ch_field.at(0)  += ((double**)N2data.Data)[r][34];
				avg_ch_field.at(1)  += ((double**)N2data.Data)[r][35];
				avg_ch_field.at(2)  += ((double**)N2data.Data)[r][36];
				avg_ch_field.at(3)  += ((double**)N2data.Data)[r][37];
				avg_ch_field.at(4)  += ((double**)N2data.Data)[r][38];
				avg_ch_field.at(5)  += ((double**)N2data.Data)[r][39];
				avg_ch_field.at(6)  += ((double**)N2data.Data)[r][40];
				avg_ch_field.at(7)  += ((double**)N2data.Data)[r][41];
				avg_ch_field.at(8)  += ((double**)N2data.Data)[r][42];
				avg_ch_field.at(9)  += ((double**)N2data.Data)[r][43];
				avg_ch_field.at(10) += ((double**)N2data.Data)[r][44];
				avg_ch_field.at(11) += ((double**)N2data.Data)[r][45];
				avg_ch_field.at(12) += ((double**)N2data.Data)[r][46];
				avg_ch_field.at(13) += ((double**)N2data.Data)[r][47];
				avg_ch_field.at(14) += ((double**)N2data.Data)[r][48];
				avg_ch_field.at(15) += ((double**)N2data.Data)[r][49];
				it+=1;
			}
		}	
		N2_ClearConfig(&N2data);

		// Scale the frequencies to Bfield
		for( int ch=0; ch< avg_ch_field.size(); ch++ ){
			avg_ch_field.at(ch) /= (double(it) * 1000 * 2 * gammaF / 1e6 );

		}

		// Save to tree
		out_Cycle	= this_cycle;
		out_Run		= this_run;
		out_CsField	= avg_ch_field;
		outTree->Fill();
	}
	Bvector test = CalculateField( FieldMap, 17.5, 270, -40., "down");
	cout << test.rho << " " << test.phi << " " << test.z << "\n";
	outFile->cd();
	outTree->Write(); // also deletes outTree
	outFile->Close();




	return 1;
}



