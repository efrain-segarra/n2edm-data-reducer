#pragma once
#ifndef __ReduceSummary_H__

namespace fs = std::filesystem;
using namespace std;

inline double UCN_Counting_Stop = DUMMY_VAL;

// ---------------------------------------------------------
// Summary reducer
// ---------------------------------------------------------
struct SummaryResult{
	double Tfill 	= DUMMY_VAL;
	double Tstore 	= DUMMY_VAL;
	double Tcount 	= DUMMY_VAL;
	double Tpump 	= DUMMY_VAL;
};

std::map<int, SummaryResult> summaryfile_cache;
inline SummaryResult ReduceSummary(const int cycle){
	return summaryfile_cache[cycle];
}

inline void LoadSummaryFile(const std::string& filepath ){
	SummaryResult thisSummaryEvent;

	std::ifstream file(filepath);
	std::string line;

	// missing file flag
	if (!file.is_open()) {
		cerr << "Could not find file " << filepath << "\n";
		return; 
	}

	while (std::getline(file, line)) {
		std::stringstream ss(line);
		std::string val;
		std::vector<std::string> row;

		while (std::getline(ss, val, '\t')) {
			row.push_back(val);
		}

		// Check expected size
		if( row.size() == 17 ){
			int this_cycle = std::stoi(row.at(0));
			thisSummaryEvent.Tfill 	= std::stod(row.at(3));
			thisSummaryEvent.Tstore = std::stod(row.at(4));
			thisSummaryEvent.Tpump 	= std::stod(row.at(5));
			thisSummaryEvent.Tcount = std::stod(row.at(6));
			summaryfile_cache[this_cycle] = thisSummaryEvent;
		}	
		else{
			cerr << "Incorrect summary file size. Skipping...\n";
			return;
		}
	}

	return;
}



#endif
