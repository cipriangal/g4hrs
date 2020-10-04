#include "g4hrsDatabase.hh"
#include "G4SystemOfUnits.hh"

#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
using std::cout;
using std::endl;

g4hrsDatabase::g4hrsDatabase(G4String targmat) {

         Table = true;  

	if(targmat == "Pb208") {
		E_min = 0.55*GeV;
		E_step = 0.05*GeV;
		n_E = 14;
		n_Th = 66;
		LoadTable("horpb.dat",0);
		LoadTable("horpb1.dat",1);
	}
	else if(targmat == "C12") {
		E_min = 0.9*GeV;
		E_step = 0.05*GeV;
		n_E = 3;
		n_Th = 71;
		LoadTable("c12_fsu.dat",0);
		// There is no "strech" table for 12C
		// Loading the same table twice will just cause sensitivity = 0
		LoadTable("c12_fsu.dat",1);
 	} else if(targmat == "Diamond") {
                E_min = 0.9*GeV;
                E_step = 0.05*GeV;
                n_E = 3;
                n_Th = 71;
                LoadTable("c12_fsu.dat",0);
                LoadTable("c12_fsu.dat",1);
        } 
	else if(targmat == "Ca48") {
		E_min = 0.5*GeV;
		E_step = 0.05*GeV;
		n_E = 62;
		n_Th = 141;
		LoadTable("ca48_fsu.dat",0);
		LoadTable("ca48_fsu_stretched.dat",1);
	}

	else if(targmat == "Ca40") {
		E_min = 0.5*GeV;
		E_step = 0.05*GeV;
		n_E = 62;
		n_Th = 141;
		LoadTable("ca40_fsu.dat",0);
		LoadTable("ca40_fsu_stretched.dat",1);
	}

        else { Table = false; cout << "No Table for this material" << endl; }
 


// end constructer 
} 

g4hrsDatabase::~g4hrsDatabase()	{
}

void g4hrsDatabase::LoadTable(string filename, int stretch) {
	// stretch = {0, 1} = {R_n not stretched, R_n stretched 1%)
	
	ifstream datafile(filename, std::ifstream::in);

        if( !datafile ){
            fprintf(stderr, "ERROR %s line %d:  %s opening %s failed\n", __FILE__, __LINE__, __FUNCTION__, filename.c_str());
            exit(1);
        }
	
	double thisEnergy, thisAngle, thisXS, thisAsym, ignore;
	string dummy;

	for(int i=0; i<n_E; i++) {

		datafile >> dummy;
		thisEnergy = E_min + double(i)*E_step;
		energy.push_back(thisEnergy);

		vector<double> row_XS;
		vector<double> row_asym;

		for(int j=0; j<n_Th; j++) {
			datafile >> thisAngle >> thisXS >> thisAsym;
			angle.push_back(thisAngle*deg);			
			row_XS.push_back(thisXS);
			row_asym.push_back(thisAsym);
		} // end j (angle) for loop
	
		if(stretch == 0) {
			XS.push_back(row_XS);
			asym.push_back(row_asym);
		} 
		
		if(stretch == 1) {
			XS_str.push_back(row_XS);
			asym_str.push_back(row_asym);
		}

	} // end i (energy) for loop 

	cout << "Table " << filename << " loaded" << endl;

// end LoadTable
}

double g4hrsDatabase::Interpolate(double thisE, double thisTh, int stretch, int value) { 
	// stretch = {0, 1} = {R_n not stretched, R_n stretched 1%)
	// value = {0, 1} = {form factor squared, asymmetry}	
	double th0, th1, e0, e1;
	int i0, i1, j0, j1;

        if(!Table) return 0;

        // Check ranges
        if( thisTh < angle[0] ) return 0;
        if( thisTh > angle[angle.size()-1] ) return 0;
        if( thisE < energy[0] ) return 0;
        if( thisE > energy[energy.size()-1] ) return 0;


	for(int i=1; i<n_E; i++) {
		e0 = energy[i-1];
		e1 = energy[i];
		i0 = i-1;
		i1 = i;
		if(energy[i-1] < thisE && energy[i] >= thisE) {
			break; 
		}  
	} 	
	for(int j=1; j<n_Th; j++) {
		th0 = angle[j-1];
		th1 = angle[j];
		j0 = j-1;
		j1 = j;
		if(angle[j-1] < thisTh && angle[j] >= thisTh) {
			break;
		}
	}

	double XS_e0_th0; 
	double XS_e0_th1; 
	double XS_e1_th0;
	double XS_e1_th1;
	double asym_e0_th0;
	double asym_e0_th1;
	double asym_e1_th0;
	double asym_e1_th1;

	if(stretch == 0) {
		XS_e0_th0 = XS[i0][j0];
		XS_e0_th1 = XS[i0][j1];
		XS_e1_th0 = XS[i1][j0];
		XS_e1_th1 = XS[i1][j1];
		asym_e0_th0 = asym[i0][j0];
		asym_e0_th1 = asym[i0][j1];
		asym_e1_th0 = asym[i1][j0];
		asym_e1_th1 = asym[i1][j1];
	}
	
	if(stretch == 1) {
		XS_e0_th0 = XS_str[i0][j0];
		XS_e0_th1 = XS_str[i0][j1];
		XS_e1_th0 = XS_str[i1][j0];
		XS_e1_th1 = XS_str[i1][j1];
		asym_e0_th0 = asym_str[i0][j0];
		asym_e0_th1 = asym_str[i0][j1];
		asym_e1_th0 = asym_str[i1][j0];
		asym_e1_th1 = asym_str[i1][j1];
	}
	
	double answer;

	if(value == 0) {
		double m_e0 = (XS_e0_th1 - XS_e0_th0)/(th1 - th0);
		double m_e1 = (XS_e1_th1 - XS_e1_th0)/(th1 - th0);
		double XS0 = XS_e0_th0 + m_e0*(thisTh - th0); 
		double XS1 = XS_e1_th0 + m_e1*(thisTh - th0);
		double m_th = (XS1 - XS0)/(e1 - e0);
		answer = XS0 + m_th*(thisE - e0); 
	}	
	if(value == 1) {
		double m_e0 = (asym_e0_th1 - asym_e0_th0)/(th1 - th0);
		double m_e1 = (asym_e1_th1 - asym_e1_th0)/(th1 - th0);
		double asym0 = asym_e0_th0 + m_e0*(thisTh - th0); 
		double asym1 = asym_e1_th0 + m_e1*(thisTh - th0);
		double m_th = (asym1 - asym0)/(e1 - e0);
		answer = asym0 + m_th*(thisE - e0); 
	}	
	
	return answer;

}












