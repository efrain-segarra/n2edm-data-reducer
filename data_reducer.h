#pragma once
#ifndef __data_reducer_H__

#include <iostream>
#include <string>
#include <filesystem>
#include <set>
#include <iomanip>
#include <sstream>

//#include <limits>
#include <fstream>
#include <format>
#include <vector>
#include <map>
#include <cstdio>

// ROOT headers (needs to be first due to conflict with n2dataread function defs)
#include "TFile.h"
#include "TTree.h"
#include "TVector.h"
#include "TEnv.h"

// n2dataread headers
#include "SimpleLog/SimpleLog.h"
#include "ListHash/ListHash.h"
#include "N2readData.h"

#include "constants.h"
// reduction headers
#include "ReduceRf.h"
#include "ReduceSf.h"
#include "ReduceHg.h"
#include "ReduceUcn.h"
#include "ReduceTemperature.h"
#include "ReduceCsm.h"

namespace fs = std::filesystem;

using namespace std;



// ---------------------------------------------------------
// Helper: Zero-pad integers (e.g., 8003 -> "008003")
// ---------------------------------------------------------
std::string formatNumber(int num, int width = 6) {
	std::ostringstream oss;
	oss << std::setw(width) << std::setfill('0') << num;
	return oss.str();
}



// ---------------------------------------------------------
// Functions for field map
// ---------------------------------------------------------

// Helper function for reading CSV of gradient table
std::vector<std::string> parseCSVRow(const std::string& line) {
	std::vector<std::string> result;
	std::string cell;
	bool inside_quotes = false;

	for (char c : line) {
		if (c == '"') {
			inside_quotes = !inside_quotes; // Toggle state
		} 
		else if (c == ',' && !inside_quotes) {
			result.push_back(cell);         // End of cell
			cell.clear();
		} 
		else {
			cell += c;                      // Add character to current cell
		}
	}
	result.push_back(cell); // Add the last cell
	return result;
}

// Struct for reading the gradient table
struct GradientInfo{
	// mu = average
	// sigma = reproducibility
	// tau = repeatibility
	double mu_up_2025;
	double mu_up_2022;
	double sigma_up_2025;
	double sigma_up_2022;
	double mu_down_2025;
	double mu_down_2022;
	double sigma_down_2025;
	double sigma_down_2022;
	double tau_2025;
	double tau_2022;
};

map<string,GradientInfo> Gradients(string path_to_csv){
	// What we will return
	map<string,GradientInfo> gradientmap;

	ifstream file(path_to_csv);
	string line;
	getline(file,line); // read first info line and throw away
	while( getline(file,line)){
		if( line.empty()) continue;
		vector<string> rowVec = parseCSVRow(line);
		if (rowVec.size() >= 11) {
			try {
				GradientInfo thisrow;
				thisrow.mu_up_2025		= std::stod(rowVec[1]);
				thisrow.mu_up_2022		= std::stod(rowVec[2]);
				thisrow.sigma_up_2025		= std::stod(rowVec[3]);
				thisrow.sigma_up_2022		= std::stod(rowVec[4]);
				thisrow.mu_down_2025		= std::stod(rowVec[5]);
				thisrow.mu_down_2022		= std::stod(rowVec[6]);
				thisrow.sigma_down_2025		= std::stod(rowVec[7]);
				thisrow.sigma_down_2022		= std::stod(rowVec[8]);
				thisrow.tau_2025		= std::stod(rowVec[9]);
				thisrow.tau_2022		= std::stod(rowVec[10]);
				gradientmap[rowVec[0]] = thisrow;
			} catch (const std::exception& e) {
				std::cerr << "Error parsing line: " << line << "\nReason: " << e.what() << std::endl;
			}
		}
	}
	return gradientmap;
}

struct Bvector{
	double rho = 0;
	double phi = 0;
	double z = 0;
};

Bvector HarmonicPolynomials_Cylindrical(int l, int m, double rho, double phi, double z){
	double Pi_rho = 0;
	double Pi_phi = 0;
	double Pi_z = 0;
	double Dl = 1;

	// Precompute some powers for convenience
	double rho7 = std::pow(rho,7);
	double rho6 = std::pow(rho,6);
	double rho5 = std::pow(rho,5);
	double rho4 = std::pow(rho,4);
	double rho3 = std::pow(rho,3);
	double rho2 = std::pow(rho,2);
	double z6 = std::pow(z,6);
	double z5 = std::pow(z,5);
	double z4 = std::pow(z,4);
	double z3 = std::pow(z,3);
	double z2 = std::pow(z,2);

	// --- Degree l = 0 ---
	if (l == 0) {
		if (m == -1) {
			Pi_rho = std::sin(phi);
			Pi_phi = std::cos(phi);
			Pi_z   = 0.0;
		} else if (m == 0) {
			Pi_rho = 0.0;
			Pi_phi = 0.0;
			Pi_z   = 1.0;
		} else if (m == 1) {
			Pi_rho = std::cos(phi);
			Pi_phi = -std::sin(phi);
			Pi_z   = 0.0;
		}
		Dl = 1;
	}
	// --- Degree l = 1 ---
	else if (l == 1) {
		if (m == -2) {
			Pi_rho = rho * std::sin(2 * phi);
			Pi_phi = rho * std::cos(2 * phi);
			Pi_z   = 0.0;
		} else if (m == -1) {
			Pi_rho = z * std::sin(phi);
			Pi_phi = z * std::cos(phi);
			Pi_z   = rho * std::sin(phi);
		} else if (m == 0) {
			Pi_rho = -0.5 * rho;
			Pi_phi = 0.0;
			Pi_z   = z;
		} else if (m == 1) {
			Pi_rho = z * std::cos(phi);
			Pi_phi = -z * std::sin(phi);
			Pi_z   = rho * std::cos(phi);
		} else if (m == 2) {
			Pi_rho = rho * std::cos(2 * phi);
			Pi_phi = -rho * std::sin(2 * phi);
			Pi_z   = 0.0;
		}
		Dl = 1;
	}
	// --- Degree l = 2 ---
	else if (l == 2) {
		if (m == -3) {
			Pi_rho = rho2 * std::sin(3 * phi);
			Pi_phi = rho2 * std::cos(3 * phi);
			Pi_z   = 0.0;
		} else if (m == -2) {
			Pi_rho = 2 * rho * z * std::sin(2 * phi);
			Pi_phi = 2 * rho * z * std::cos(2 * phi);
			Pi_z   = rho2 * std::sin(2 * phi);
		} else if (m == -1) {
			Pi_rho = 0.25 * (4 * z2 - 3 * rho2) * std::sin(phi);
			Pi_phi = 0.25 * (4 * z2 - rho2) * std::cos(phi);
			Pi_z   = 2 * rho * z * std::sin(phi);
		} else if (m == 0) {
			Pi_rho = -rho * z;
			Pi_phi = 0.0;
			Pi_z   = -0.5 * rho2 + z2;
		} else if (m == 1) {
			Pi_rho = 0.25 * (4 * z2 - 3 * rho2) * std::cos(phi);
			Pi_phi = 0.25 * (rho2 - 4 * z2) * std::sin(phi);
			Pi_z   = 2 * rho * z * std::cos(phi);
		} else if (m == 2) {
			Pi_rho = 2 * rho * z * std::cos(2 * phi);
			Pi_phi = -2 * rho * z * std::sin(2 * phi);
			Pi_z   = rho2 * std::cos(2 * phi);
		} else if (m == 3) {
			Pi_rho = rho2 * std::cos(3 * phi);
			Pi_phi = -rho2 * std::sin(3 * phi);
			Pi_z   = 0.0;
		}
		Dl = 18.; // cm
	}
	// --- Degree l = 3 ---
	else if (l == 3) {
		if (m == -4) {
			Pi_rho = rho3 * std::sin(4 * phi);
			Pi_phi = rho3 * std::cos(4 * phi);
			Pi_z   = 0.0;
		} else if (m == -3) {
			Pi_rho = 3 * rho2 * z * std::sin(3 * phi);
			Pi_phi = 3 * rho2 * z * std::cos(3 * phi);
			Pi_z   = rho3 * std::sin(3 * phi);
		} else if (m == -2) {
			Pi_rho = rho * (3 * z2 - rho2) * std::sin(2 * phi);
			Pi_phi = 0.5 * rho * (6 * z2 - rho2) * std::cos(2 * phi);
			Pi_z   = 3 * rho2 * z * std::sin(2 * phi);
		} else if (m == -1) {
			Pi_rho = 0.25 * z * (4 * z2 - 9 * rho2) * std::sin(phi);
			Pi_phi = 0.25 * z * (4 * z2 - 3 * rho2) * std::cos(phi);
			Pi_z   = rho * (3 * z2 - 0.75 * rho2) * std::sin(phi);
		} else if (m == 0) {
			Pi_rho = 0.375 * rho * (rho2 - 4 * z2); // 3/8 = 0.375
			Pi_phi = 0.0;
			Pi_z   = 0.5 * z * (2 * z2 - 3 * rho2);
		} else if (m == 1) {
			Pi_rho = 0.25 * z * (4 * z2 - 9 * rho2) * std::cos(phi);
			Pi_phi = 0.25 * z * (3 * rho2 - 4 * z2) * std::sin(phi);
			Pi_z   = rho * (3 * z2 - 0.75 * rho2) * std::cos(phi);
		} else if (m == 2) {
			Pi_rho = rho * (3 * z2 - rho2) * std::cos(2 * phi);
			Pi_phi = 0.5 * rho * (rho2 - 6 * z2) * std::sin(2 * phi);
			Pi_z   = 3 * rho2 * z * std::cos(2 * phi);
		} else if (m == 3) {
			Pi_rho = 3 * rho2 * z * std::cos(3 * phi);
			Pi_phi = -3 * rho2 * z * std::sin(3 * phi);
			Pi_z   = rho3 * std::cos(3 * phi);
		} else if (m == 4) {
			Pi_rho = rho3 * std::cos(4 * phi);
			Pi_phi = -rho3 * std::sin(4 * phi);
			Pi_z   = 0.0;
		}
		Dl = 23.7*23.7; //cm^2
	}
	// --- Degree l = 4 ---
	else if (l == 4) {
		if (m == -5) {
			Pi_rho = rho4 * std::sin(5 * phi);
			Pi_phi = rho4 * std::cos(5 * phi);
			Pi_z   = 0.0;
		} 
		else if (m == -4) {
			Pi_rho = 4.0 * rho3 * z * std::sin(4 * phi);
			Pi_phi = 4.0 * rho3 * z * std::cos(4 * phi);
			Pi_z   = rho4 * std::sin(4 * phi);
		} 
		else if (m == -3) {
			Pi_rho = 0.25 * (24.0 * rho2 * z2 - 5.0 * rho4) * std::sin(3 * phi);
			Pi_phi = 0.75 * (8.0 * rho2 * z2 - rho4) * std::cos(3 * phi);
			Pi_z   = 4.0 * rho3 * z * std::sin(3 * phi);
		} 
		else if (m == -2) {
			Pi_rho = 4.0 * (rho * z3 - rho3 * z) * std::sin(2 * phi);
			Pi_phi = 2.0 * (2.0 * rho * z3 - rho3 * z) * std::cos(2 * phi);
			Pi_z   = (6.0 * rho2 * z2 - rho4) * std::sin(2 * phi);
		} 
		else if (m == -1) {
			Pi_rho = 0.125 * (8.0 * z4 - 36.0 * rho2 * z2 + 5.0 * rho4) * std::sin(phi);
			Pi_phi = 0.125 * (8.0 * z4 - 12.0 * rho2 * z2 + rho4) * std::cos(phi);
			Pi_z   = (4.0 * rho * z3 - 3.0 * rho3 * z) * std::sin(phi);
		} 
		else if (m == 0) {
			Pi_rho = 0.5 * (3.0 * rho3 * z - 4.0 * rho * z3);
			Pi_phi = 0.0;
			Pi_z   = 0.125 * (8.0 * z4 - 24.0 * rho2 * z2 + 3.0 * rho4);
		} 
		else if (m == 1) {
			Pi_rho = 0.125 * (8.0 * z4 - 36.0 * rho2 * z2 + 5.0 * rho4) * std::cos(phi);
			Pi_phi = -0.125 * (8.0 * z4 - 12.0 * rho2 * z2 + rho4) * std::sin(phi);
			Pi_z   = (4.0 * rho * z3 - 3.0 * rho3 * z) * std::cos(phi);
		} 
		else if (m == 2) {
			Pi_rho = 4.0 * (rho * z3 - rho3 * z) * std::cos(2 * phi);
			Pi_phi = -2.0 * (2.0 * rho * z3 - rho3 * z) * std::sin(2 * phi);
			Pi_z   = (6.0 * rho2 * z2 - rho4) * std::cos(2 * phi);
		} 
		else if (m == 3) {
			Pi_rho = 0.25 * (24.0 * rho2 * z2 - 5.0 * rho4) * std::cos(3 * phi);
			Pi_phi = -0.75 * (8.0 * rho2 * z2 - rho4) * std::sin(3 * phi);
			Pi_z   = 4.0 * rho3 * z * std::cos(3 * phi);
		} 
		else if (m == 4) {
			Pi_rho = 4.0 * rho3 * z * std::cos(4 * phi);
			Pi_phi = -4.0 * rho3 * z * std::sin(4 * phi);
			Pi_z   = rho4 * std::cos(4 * phi);
		} 
		else if (m == 5) {
			Pi_rho = rho4 * std::cos(5 * phi);
			Pi_phi = -rho4 * std::sin(5 * phi);
			Pi_z   = 0.0;
		}
		Dl = std::pow(-29.1,3); //cm^3
	}
	// --- Degree l = 5 ---
	else if (l == 5) {
		if (m == -6) {
			Pi_rho = rho5 * std::sin(6 * phi);
			Pi_phi = rho5 * std::cos(6 * phi);
			Pi_z   = 0.0;
		} 
		else if (m == -5) {
			Pi_rho = 5.0 * rho4 * z * std::sin(5 * phi);
			Pi_phi = 5.0 * rho4 * z * std::cos(5 * phi);
			Pi_z   = rho5 * std::sin(5 * phi);
		} 
		else if (m == -4) {
			Pi_rho = 0.5 * (20.0 * rho3 * z2 - 3.0 * rho5) * std::sin(4 * phi);
			Pi_phi = rho3 * (10.0 * z2 - rho2) * std::cos(4 * phi);
			Pi_z   = 5.0 * rho4 * z * std::sin(4 * phi);
		} 
		else if (m == -3) {
			Pi_rho = 1.25 * (8.0 * rho2 * z3 - 5.0 * rho4 * z) * std::sin(3 * phi);
			Pi_phi = 1.25 * (8.0 * rho2 * z3 - 3.0 * rho4 * z) * std::cos(3 * phi);
			Pi_z   = 1.25 * (8.0 * rho3 * z2 - rho5) * std::sin(3 * phi);
		} 
		else if (m == -2) {
			Pi_rho = 0.3125 * (16.0 * rho * z4 - 32.0 * rho3 * z2 + 3.0 * rho5) * std::sin(2 * phi);
			Pi_phi = 0.3125 * (16.0 * rho * z4 - 16.0 * rho3 * z2 + rho5) * std::cos(2 * phi);
			Pi_z   = 5.0 * (2.0 * rho2 * z3 - rho4 * z) * std::sin(2 * phi);
		} 
		else if (m == -1) {
			Pi_rho = 0.125 * (8.0 * z5 - 60.0 * rho2 * z3 + 25.0 * rho4 * z) * std::sin(phi);
			Pi_phi = 0.125 * (8.0 * z5 - 20.0 * rho2 * z3 + 5.0 * rho4 * z) * std::cos(phi);
			Pi_z   = 0.625 * (8.0 * rho * z4 - 12.0 * rho3 * z2 + rho5) * std::sin(phi);
		} 
		else if (m == 0) {
			Pi_rho = 0.3125 * (-8.0 * rho * z4 + 12.0 * rho3 * z2 - rho5);
			Pi_phi = 0.0;
			Pi_z   = 0.125 * (8.0 * z5 - 40.0 * rho2 * z3 + 15.0 * rho4 * z);
		} 
		else if (m == 1) {
			Pi_rho = 0.125 * (8.0 * z5 - 60.0 * rho2 * z3 + 25.0 * rho4 * z) * std::cos(phi);
			Pi_phi = -0.125 * (8.0 * z5 - 20.0 * rho2 * z3 + 5.0 * rho4 * z) * std::sin(phi);
			Pi_z   = 0.625 * (8.0 * rho * z4 - 12.0 * rho3 * z2 + rho5) * std::cos(phi);
		} 
		else if (m == 2) {
			Pi_rho = 0.3125 * (16.0 * rho * z4 - 32.0 * rho3 * z2 + 3.0 * rho5) * std::cos(2 * phi);
			Pi_phi = -0.3125 * (16.0 * rho * z4 - 16.0 * rho3 * z2 + rho5) * std::sin(2 * phi);
			Pi_z   = 5.0 * (2.0 * rho2 * z3 - rho4 * z) * std::cos(2 * phi);
		} 
		else if (m == 3) {
			Pi_rho = 1.25 * (8.0 * rho2 * z3 - 5.0 * rho4 * z) * std::cos(3 * phi);
			Pi_phi = -1.25 * (8.0 * rho2 * z3 - 3.0 * rho4 * z) * std::sin(3 * phi);
			Pi_z   = 1.25 * (8.0 * rho3 * z2 - rho5) * std::cos(3 * phi);
		} 
		else if (m == 4) {
			Pi_rho = 0.5 * (20.0 * rho3 * z2 - 3.0 * rho5) * std::cos(4 * phi);
			Pi_phi = -rho3 * (10.0 * z2 - rho2) * std::sin(4 * phi);
			Pi_z   = 5.0 * rho4 * z * std::cos(4 * phi);
		} 
		else if (m == 5) {
			Pi_rho = 5.0 * rho4 * z * std::cos(5 * phi);
			Pi_phi = -5.0 * rho4 * z * std::sin(5 * phi);
			Pi_z   = rho5 * std::cos(5 * phi);
		} 
		else if (m == 6) {
			Pi_rho = rho5 * std::cos(6 * phi);
			Pi_phi = -rho5 * std::sin(6 * phi);
			Pi_z   = 0.0;
		}
		Dl = std::pow(31.8,4); //cm^4
	}
	// --- Degree l = 6 ---
	else if (l == 6) {
		if (m == -7) {
			Pi_rho = rho6 * std::sin(7 * phi);
			Pi_phi = rho6 * std::cos(7 * phi);
			Pi_z   = 0.0;
		} 
		else if (m == -6) {
			Pi_rho = 6.0 * rho5 * z * std::sin(6 * phi);
			Pi_phi = 6.0 * rho5 * z * std::cos(6 * phi);
			Pi_z   = rho6 * std::sin(6 * phi);
		} 
		else if (m == -5) {
			Pi_rho = 0.25 * rho4 * (60.0 * z2 - 7.0 * rho2) * std::sin(5 * phi);
			Pi_phi = 1.25 * rho4 * (12.0 * z2 - rho2) * std::cos(5 * phi);
			Pi_z   = 6.0 * rho5 * z * std::sin(5 * phi);
		} 
		else if (m == -4) {
			Pi_rho = rho3 * z * (20.0 * z2 - 9.0 * rho2) * std::cos(4 * phi);
			Pi_phi = 2.0 * rho3 * z * (10.0 * z2 - 3.0 * rho2) * std::cos(4 * phi);
			Pi_z   = 1.5 * rho4 * (10.0 * z2 - rho2) * std::sin(4 * phi);
		} 
		else if (m == -3) {
			Pi_rho = 0.1875 * rho2 * (80.0 * z4 - 100.0 * rho2 * z2 + 7.0 * rho4) * std::cos(3 * phi);
			Pi_phi = 0.1875 * rho2 * (80.0 * z4 - 60.0 * rho2 * z2 + 3.0 * rho4) * std::cos(3 * phi);
			Pi_z   = 2.5 * rho3 * z * (8.0 * z2 - 3.0 * rho2) * std::sin(3 * phi);
		} 
		else if (m == -2) {
			Pi_rho = 0.125 * rho * z * (48.0 * z4 - 160.0 * rho2 * z2 + 45.0 * rho4) * std::cos(2 * phi);
			Pi_phi = 0.125 * rho * z * (48.0 * z4 - 80.0 * rho2 * z2 + 15.0 * rho4) * std::cos(2 * phi);
			Pi_z   = 0.9375 * rho2 * (16.0 * z4 - 16.0 * rho2 * z2 + rho4) * std::sin(2 * phi);
		} 
		else if (m == -1) {
			Pi_rho = (1.0/64.0) * (64.0 * z6 - 720.0 * rho2 * z4 + 600.0 * rho4 * z2 - 35.0 * rho6) * std::cos(phi);
			Pi_phi = (1.0/64.0) * (64.0 * z6 - 240.0 * rho2 * z4 + 120.0 * rho4 * z2 - 5.0 * rho6) * std::cos(phi);
			Pi_z   = 0.75 * rho * z * (8.0 * z4 - 20.0 * rho2 * z2 + 5.0 * rho4) * std::sin(phi);
		} 
		else if (m == 0) {
			Pi_rho = 0.375 * rho * (-8.0 * z5 + 20.0 * rho2 * z3 - 5.0 * rho4 * z);
			Pi_phi = 0.0;
			Pi_z   = 0.0625 * (16.0 * z6 - 120.0 * rho2 * z4 + 90.0 * rho4 * z2 - 5.0 * rho6);
		} 
		else if (m == 1) {
			Pi_rho = (1.0/64.0) * (64.0 * z6 - 720.0 * rho2 * z4 + 600.0 * rho4 * z2 - 35.0 * rho6) * std::sin(phi);
			Pi_phi = -(1.0/64.0) * (64.0 * z6 - 240.0 * rho2 * z4 + 120.0 * rho4 * z2 - 5.0 * rho6) * std::sin(phi);
			Pi_z   = 0.75 * rho * z * (8.0 * z4 - 20.0 * rho2 * z2 + 5.0 * rho4) * std::cos(phi);
		} 
		else if (m == 2) {
			Pi_rho = 0.125 * rho * z * (48.0 * z4 - 160.0 * rho2 * z2 + 45.0 * rho4) * std::sin(2 * phi);
			Pi_phi = -0.125 * rho * z * (48.0 * z4 - 80.0 * rho2 * z2 + 15.0 * rho4) * std::sin(2 * phi);
			Pi_z   = 0.9375 * rho2 * (16.0 * z4 - 16.0 * rho2 * z2 + rho4) * std::cos(2 * phi);
		} 
		else if (m == 3) {
			Pi_rho = 0.1875 * rho2 * (80.0 * z4 - 100.0 * rho2 * z2 + 7.0 * rho4) * std::sin(3 * phi);
			Pi_phi = -0.1875 * rho2 * (80.0 * z4 - 60.0 * rho2 * z2 + 3.0 * rho4) * std::sin(3 * phi);
			Pi_z   = 2.5 * rho3 * z * (8.0 * z2 - 3.0 * rho2) * std::cos(3 * phi);
		} 
		else if (m == 4) {
			Pi_rho = rho3 * z * (20.0 * z2 - 9.0 * rho2) * std::sin(4 * phi);
			Pi_phi = -2.0 * rho3 * z * (10.0 * z2 - 3.0 * rho2) * std::sin(4 * phi);
			Pi_z   = 1.5 * rho4 * (10.0 * z2 - rho2) * std::cos(4 * phi);
		} 
		else if (m == 5) {
			Pi_rho = 0.25 * rho4 * (60.0 * z2 - 7.0 * rho2) * std::cos(5 * phi);
			Pi_phi = -1.25 * rho4 * (12.0 * z2 - rho2) * std::sin(5 * phi);
			Pi_z   = 6.0 * rho5 * z * std::cos(5 * phi);
		} 
		else if (m == 6) {
			Pi_rho = 6.0 * rho5 * z * std::cos(6 * phi);
			Pi_phi = -6.0 * rho5 * z * std::sin(6 * phi);
			Pi_z   = rho6 * std::cos(6 * phi);
		} 
		else if (m == 7) {
			Pi_rho = rho6 * std::cos(7 * phi);
			Pi_phi = -rho6 * std::sin(7 * phi);
			Pi_z   = 0.0;
		}
		Dl = std::pow(39.7,5); //cm^5
	}
	// --- Degree l = 7 ---
	else if (l == 7) {
		if (m == -8) {
			Pi_rho = rho7 * std::sin(8 * phi);
			Pi_phi = rho7 * std::cos(8 * phi);
			Pi_z   = 0.0;
		} 
		else if (m == -7) {
			Pi_rho = 7.0 * rho6 * z * std::sin(7 * phi);
			Pi_phi = 7.0 * rho6 * z * std::cos(7 * phi);
			Pi_z   = rho7 * std::sin(7 * phi);
		} 
		else if (m == -6) {
			Pi_rho = rho5 * (21.0 * z2 - 2.0 * rho2) * std::sin(6 * phi);
			Pi_phi = 1.5 * rho5 * (14.0 * z2 - rho2) * std::cos(6 * phi);
			Pi_z   = 7.0 * rho6 * z * std::sin(6 * phi);
		} 
		else if (m == -5) {
			Pi_rho = 1.75 * rho4 * z * (20.0 * z2 - 7.0 * rho2) * std::sin(5 * phi);
			Pi_phi = 8.75 * rho4 * z * (4.0 * z2 - rho2) * std::cos(5 * phi);
			Pi_z   = 1.75 * rho5 * (12.0 * z2 - rho2) * std::sin(5 * phi);
		} 
		else if (m == -4) {
			Pi_rho = 1.75 * rho3 * (20.0 * z4 - 18.0 * rho2 * z2 + rho4) * std::sin(4 * phi);
			Pi_phi = 0.875 * rho3 * (40.0 * z4 - 24.0 * rho2 * z2 + 3.0 * rho4) * std::cos(4 * phi);
			Pi_z   = 3.5 * rho4 * z * (10.0 * z2 - 3.0 * rho2) * std::sin(4 * phi);
		} 
		else if (m == -3) {
			Pi_rho = (7.0/16.0) * rho2 * z * (48.0 * z4 - 100.0 * rho2 * z2 + 21.0 * rho4) * std::sin(3 * phi);
			Pi_phi = (21.0/16.0) * rho2 * z * (16.0 * z4 - 20.0 * rho2 * z2 + 3.0 * rho4) * std::cos(3 * phi);
			Pi_z   = (7.0/16.0) * rho3 * (80.0 * z4 - 60.0 * rho2 * z2 + 3.0 * rho4) * std::sin(3 * phi);
		} 
		else if (m == -2) {
			Pi_rho = (7.0/16.0) * rho * (16.0 * z6 - 80.0 * rho2 * z4 + 45.0 * rho4 * z2 - 2.0 * rho6) * std::sin(2 * phi);
			Pi_phi = (7.0/32.0) * rho * (32.0 * z6 - 80.0 * rho2 * z4 + 30.0 * rho4 * z2 - rho6) * std::cos(2 * phi);
			Pi_z   = (7.0/16.0) * rho2 * z * (48.0 * z4 - 80.0 * rho2 * z2 + 15.0 * rho4) * std::sin(2 * phi);
		} 
		else if (m == -1) {
			Pi_rho = (1.0/64.0) * z * (64.0 * z6 - 1008.0 * rho2 * z4 + 1400.0 * rho4 * z2 - 245.0 * rho6) * std::sin(phi);
			Pi_phi = (1.0/64.0) * z * (64.0 * z6 - 336.0 * rho2 * z4 + 280.0 * rho4 * z2 - 35.0 * rho6) * std::cos(phi);
			Pi_z   = (7.0/64.0) * rho * (64.0 * z6 - 240.0 * rho2 * z4 + 120.0 * rho4 * z2 - 5.0 * rho6) * std::sin(phi);
		} 
		else if (m == 0) {
			Pi_rho = (7.0/128.0) * rho * (-64.0 * z6 + 240.0 * rho2 * z4 - 120.0 * rho4 * z2 + 5.0 * rho6);
			Pi_phi = 0.0;
			Pi_z   = (1.0/16.0) * z * (16.0 * z6 - 168.0 * rho2 * z4 + 210.0 * rho4 * z2 - 35.0 * rho6);
		} 
		else if (m == 1) {
			Pi_rho = (1.0/64.0) * z * (64.0 * z6 - 1008.0 * rho2 * z4 + 1400.0 * rho4 * z2 - 245.0 * rho6) * std::cos(phi);
			Pi_phi = -(1.0/64.0) * z * (64.0 * z6 - 336.0 * rho2 * z4 + 280.0 * rho4 * z2 - 35.0 * rho6) * std::sin(phi);
			Pi_z   = (7.0/64.0) * rho * (64.0 * z6 - 240.0 * rho2 * z4 + 120.0 * rho4 * z2 - 5.0 * rho6) * std::cos(phi);
		} 
		else if (m == 2) {
			Pi_rho = (7.0/16.0) * rho * (16.0 * z6 - 80.0 * rho2 * z4 + 45.0 * rho4 * z2 - 2.0 * rho6) * std::cos(2 * phi);
			Pi_phi = -(7.0/32.0) * rho * (32.0 * z6 - 80.0 * rho2 * z4 + 30.0 * rho4 * z2 - rho6) * std::sin(2 * phi);
			Pi_z   = (7.0/16.0) * rho2 * z * (48.0 * z4 - 80.0 * rho2 * z2 + 15.0 * rho4) * std::cos(2 * phi);
		} 
		else if (m == 3) {
			Pi_rho = (7.0/16.0) * rho2 * z * (48.0 * z4 - 100.0 * rho2 * z2 + 21.0 * rho4) * std::cos(3 * phi);
			Pi_phi = -(21.0/16.0) * rho2 * z * (16.0 * z4 - 20.0 * rho2 * z2 + 3.0 * rho4) * std::sin(3 * phi);
			Pi_z   = (7.0/16.0) * rho3 * (80.0 * z4 - 60.0 * rho2 * z2 + 3.0 * rho4) * std::cos(3 * phi);
		} 
		else if (m == 4) {
			Pi_rho = 1.75 * rho3 * (20.0 * z4 - 18.0 * rho2 * z2 + rho4) * std::cos(4 * phi);
			Pi_phi = -0.875 * rho3 * (40.0 * z4 - 24.0 * rho2 * z2 + 3.0 * rho4) * std::sin(4 * phi);
			Pi_z   = 3.5 * rho4 * z * (10.0 * z2 - 3.0 * rho2) * std::cos(4 * phi);
		} 
		else if (m == 5) {
			Pi_rho = 1.75 * rho4 * z * (20.0 * z2 - 7.0 * rho2) * std::cos(5 * phi);
			Pi_phi = -8.75 * rho4 * z * (4.0 * z2 - rho2) * std::sin(5 * phi);
			Pi_z   = 1.75 * rho5 * (12.0 * z2 - rho2) * std::cos(5 * phi);
		} 
		else if (m == 6) {
			Pi_rho = rho5 * (21.0 * z2 - 2.0 * rho2) * std::cos(6 * phi);
			Pi_phi = -1.5 * rho5 * (14.0 * z2 - rho2) * std::sin(6 * phi);
			Pi_z   = 7.0 * rho6 * z * std::cos(6 * phi);
		} 
		else if (m == 7) {
			Pi_rho = 7.0 * rho6 * z * std::cos(7 * phi);
			Pi_phi = -7.0 * rho6 * z * std::sin(7 * phi);
			Pi_z   = rho7 * std::cos(7 * phi);
		} 
		else if (m == 8) {
			Pi_rho = rho7 * std::cos(8 * phi);
			Pi_phi = -rho7 * std::sin(8 * phi);
			Pi_z   = 0.0;
		}
		Dl = std::pow(33.8,6); //cm^6
	}
	Bvector Pis;
	Pis.rho = Pi_rho / Dl;
	Pis.phi = Pi_phi / Dl;
	Pis.z 	= Pi_z   / Dl;
	return Pis;
}

Bvector CalculateField( map<string,GradientInfo> GradientMap, double r, double phi, double z, string polarity){

	// Gradients from the table have been normalized, Glm'
	// Units of gradient table are pT for l=0, fT/cm for all higher orders.
	// The polynomial harmonics take into account the D_(l-1) factor so:
	// 	Glm' * Pi(l,m) = fT/cm * 1/cm^(l-1) 
	//
	// Input should be [r]=cm, [phi]=deg, [z]=cm

	Bvector Btotal;

	// Loop over the gradient table to do our l,m iteration
	for( const auto& gradients: GradientMap){

		// Grab l,m from the name flag in the given row
		int l,m;
		sscanf( gradients.first.c_str(), "G_{%d,%d}", &l, &m);

		// Now given the l,m, add up the field for our given position
		phi *= M_PI/180.;
		Bvector Pis = HarmonicPolynomials_Cylindrical( l, m, r, phi, z);

		double Glm = 0;
		if( polarity == "up" ){
			Glm = gradients.second.mu_up_2025;
		}
		else if(polarity == "down"){
			Glm = gradients.second.mu_down_2025;
		}
		else{
			cerr << "undefined polarity for field calculation\n";
			return Btotal;
		}

		if( l == 0 ) Glm *= 1000; // pT -> fT to be consistent with rest of table

		//cout << "l,m " << l << "," << m << ", Glm: " << Glm << ", Pis at position: " << Pis.rho << " " << Pis.phi << " " << Pis.z << "\n";


		Btotal.rho += Glm*Pis.rho / 1e3;
		Btotal.phi += Glm*Pis.phi / 1e3;
		Btotal.z   += Glm*Pis.z   / 1e3;
	}

	return Btotal;
}

Bvector CsCellPosition( int daq_ch ){
	Bvector pos;

	// Top1_1 - Top1_4
	if( daq_ch == 1 ){
		pos.rho = 17.5;
		pos.phi = 135.;
		pos.z 	= 40.;
	}
	else if( daq_ch == 2 ){
		pos.rho = 39.0;
		pos.phi = 135.;
		pos.z 	= 23.;
	}
	else if( daq_ch == 3 ){
		pos.rho = 49.5;
		pos.phi = 135.;
		pos.z 	= 28.;
	}
	else if( daq_ch == 4 ){
		pos.rho = 55.0;
		pos.phi = 135.;
		pos.z 	= 40.;
	}
	// Top6_1 - Top6_4
	else if( daq_ch == 5 ){
		pos.rho = 17.5;
		pos.phi = 270.;
		pos.z 	= 40.;
	}
	else if( daq_ch == 6 ){
		pos.rho = 39.0;
		pos.phi = 270.;
		pos.z 	= 23.;
	}
	else if( daq_ch == 7 ){
		pos.rho = 49.5;
		pos.phi = 270.;
		pos.z 	= 28.;
	}
	else if( daq_ch == 8 ){
		pos.rho = 55.0;
		pos.phi = 270.;
		pos.z 	= 40.;
	}
	// Bot4_1 - Bot4_4
	else if( daq_ch == 9 ){
		pos.rho = 17.5;
		pos.phi = 0.;
		pos.z 	= -40.;
	}
	else if( daq_ch == 10 ){
		pos.rho = 39.0;
		pos.phi = 0.;
		pos.z 	= -23.;
	}
	else if( daq_ch == 11 ){
		pos.rho = 49.5;
		pos.phi = 0.;
		pos.z 	= -28.;
	}
	else if( daq_ch == 12 ){
		pos.rho = 55.0;
		pos.phi = 0.;
		pos.z 	= -40.;
	}
	// Bot6_1 - Bot6_4
	else if( daq_ch == 13 ){
		pos.rho = 17.5;
		pos.phi = 270.;
		pos.z 	= -40.;
	}
	else if( daq_ch == 14 ){
		pos.rho = 39.0;
		pos.phi = 270.;
		pos.z 	= -23.;
	}
	else if( daq_ch == 15 ){
		pos.rho = 49.5;
		pos.phi = 270.;
		pos.z 	= -28.;
	}
	else if( daq_ch == 16 ){
		pos.rho = 55.0;
		pos.phi = 270.;
		pos.z 	= -40.;
	}


	return pos;
}




#endif
