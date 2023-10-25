#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>
// #include <math.h>
#include <algorithm>
#include <string>
#include <chrono>
#include <TF1.h>
#include "TStopwatch.h"

using namespace std::chrono;
#include "/w/halla-scshelf2102/sbs/jboyd/include/include_files.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/GEM_lookups.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/beam_variables.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/MC_lookups.h"

void list_files(TString directory, vector<TString> &filenames, TString ext){
	const char *dirname = directory.Data();

	TSystemDirectory dir(dirname, dirname);
	TList *files = dir.GetListOfFiles();

	if(files){
		TSystemFile *file;
		TString fname;
		TIter next(files);
		while( (file = (TSystemFile*)next() ) ){
			fname = file->GetName();
			if( !file->IsDirectory() && fname.EndsWith( ext.Data() )) {
				filenames.push_back( Form("%s/%s", directory.Data(), fname.Data()) );
				// cout << fname.Data() << endl;
			}
		}
	}
}

Double_t fit_gaus(Double_t * x, Double_t *par){

	Double_t g = 0.0;

	g = par[0]*exp((-0.5)*pow(((x[0] -  par[1])/par[2]),2));
	
	return g;
}

template<typename T>
double VectorMean(std::vector<T> const& v){
	if(v.empty()){
		return 0;
	}
	return std::accumulate(v.begin(), v.end(), 0.0)/v.size();
}

bool single_run = true;
bool multi_run = false;

TString hadron_select = "both";
bool parsed_run_files = false;

bool calc_W = true;
bool plot_dxdy = true;
bool dxdy_only = false;

bool use_heavy_cut = false;


//Run info and lookups
int runnum = 11989;
// vector<int> runnum_vec = {13585, 13586, 13587, 13581, 13582, 13583, 13584};
// vector<int> runnum_vec = {13566, 13558, 13559, 13560, 13561, 13562, 13563, 13564, 13565, 13567, 13568, 13569, 13570, 13571};
vector<int> runnum_vec = {11449, 11451, 11452, 11456, 11493, 11494, 11495, 11496, 11551, 11554, 11562, 11563, 11564, 11565, 11568, 11570, 11571};
// vector<int> runnum_vec;
TString runs_string;
TString experiment = "gmn";
int pass = 0;

//Experimental Lookup Parameters
int kine = 8; //Kinematic setting: 4 for SBS-4, 9 for SBS-9, etc.
TString run_target = "LD2"; //Target --> "LH2" or "LD2"

double E_beam = lookup_beam_energy_from_kine( kine ); //Electron beam energy (electron energy) in GeV
int sbsfieldscale = 70;
double SBS_field = 1.0*sbsfieldscale/100.0;

TString rootfile_dir;
// TString input_rootfile;

TFile *outfile;
TChain *TC = new TChain("T");
vector<TString> master_cut_vec;
TString master_cut_string;

TString elastic_yield_str = "";
TCut master_cut = "";


//Experimental Constants, Thresholds, cuts, etc DEFINITIONS
const double pi = TMath::Pi();
const double Mp = 0.938272; //Mass of proton [GeV]
const double Mn = 0.939565; //Mass of neutron [GeV]
double W2_mean; //Invariant Mass-squared (mean val) {With perfect optics W2 = Mp. Can be calculated run-by-run}
double W2_sigma; //Invariant Mass-squared sigma {Reasonable default/guess. Can be calculated run-by-run from W plot}

//HCal constants and stuff
double tdiff = 510;		//Time difference between BBCal and HCal signals
double tdiff_max = 10;	//Maximum time difference from coincidences through tdctrig cut
double HCal_dist; 	//Distace from HCal face to target chamber
double HCal_theta;		//Theta angle for HCal from downstream beamline

const Double_t sampfrac = 0.077; 	//Estimate of the sampling fraction from MC
const Int_t kNcell = 288; // Total number of HCal modules
const Int_t kNrows = 24; // Total number of HCal rows
const Int_t kNcols = 12; // Total number of HCal columns
const Int_t kNtrack = 100; // Reasonable max number of tracks per event
const Int_t kNtdc = 1000; // Reasonable max number of tdc signals per event
const Double_t Xi = -2.20; // Distance from beam center to top of HCal in m
const Double_t Xf = 1.47; // Distance from beam center to bottom of HCal in m
const Double_t Yi = -0.853; // Distance from beam center to opposite-beam side of HCal in m
const Double_t Yf = 0.853; // Distance from beam center to beam side of HCal in m

//Declare vars
Double_t atime[kNcell], row[kNcell], col[kNcell], tdctime[kNcell], cblkid[kNcell], cblke[kNcell];
Double_t nblk, nclus, SH_nclus, PS_nclus, hcal_x, hcal_y, hcal_e;

Double_t par[3];

Double_t bb_tr_p[kNtrack], bb_tr_px[kNtrack], bb_tr_py[kNtrack], bb_tr_pz[kNtrack];
Double_t bb_tr_vz[kNtrack];
Double_t bb_tr_n, bb_ps_x, bb_ps_y, bb_ps_e, bb_sh_x, bb_sh_y, bb_sh_e;

Double_t TDCT_id[kNtdc], TDCT_tdc[kNtdc], hodo_tmean[kNtdc]; 
Int_t TDCTndata;

Long64_t Nevents;

//INITIALIZE ALL HISTOGRAMS:
TH1D *h_simc_Ep, *h_simc_PS, *h_simc_HCal_e, *h_simc_SHPS, *h_simc_W2, *h_simc_W2recon, *h_simc_vert, *h_simc_W, *h_simc_Wrecon;
TH1D *hin_Ep, *hin_PS, *hin_HCal_e, *hin_SHPS, *hin_W2, *hin_W2recon, *hin_W, *hin_Wrecon;
TH1D *h_simc_dx_cut, *h_simc_dx_fcut, *h_simc_dx_wcut_fcut, *h_simc_Q2;
TH2D *hin_dxdy_wcut, *hin_dxdy_all;

TH1D *h_simc_dx, *h_simc_dx_wcut, *h_simc_dy, *h_simc_Y;
TH1D *h_simc_dy_wcut, *h_simc_dx_p_wcut, *h_simc_dx_n_wcut;
TH1D *h_simc_dx_p, *h_simc_dx_n, *h_simc_dx_p_fcut, *h_simc_dx_n_fcut, *h_simc_dx_p_wcut_fcut, *h_simc_dx_n_wcut_fcut;
TH2D *h_simc_dxdy_all, *h_simc_dxdy_wcut, *h_simc_dxdy, *h_simc_dxdy_fcut;

double Ep_center, Ep_sigma, PS_center, PS_sigma, PS_min, SHPS_center, SHPS_sigma, HCal_e;
double HCal_e_min, W2_fit_center, W2_fit_sigma, W_fit_center, W_fit_sigma;
vector<double> cuts_from_fits, dxdy_vector;
double mc_omega, mc_sigma, mc_fnucl, mc_Weight, luminosity;
double pn_weight = 0.0;
bool is_n = false, is_p = false;
double dx, dy;

// TH1D *h_simc_atime, *h_simc_E, *h_simc_E_cut, *h_simc_E_exp, *h_simc_vert, *h_simc_W2, *h_simc_W2recon, *h_simc_dx, *h_simc_dy, *h_simc_KE_p, *h_simc_KElow, *h_simc_diff, *h_simc_X, *h_simc_Y;

double target_int = -1;

TString proton_infile, neutron_infile;

void simc_calibrate(){
	auto total_time_start = high_resolution_clock::now();

	TStopwatch *StopWatch = new TStopwatch();

	gStyle->SetPalette(55);
	cout << "--------------------------------------" << endl;
	cout << "Analysis started. " << endl;
	cout << "--------------------------------------" << endl;
	cout << "Cut booleans: " << endl;
	cout << "Calc W: " << calc_W << endl;
	cout << "Heavy cut: " << use_heavy_cut << endl;

	if( strcmp(run_target.Data(), "LH2") ){
		target_int = 0;
	}
	if( strcmp(run_target.Data(), "LD2" ) ){
		target_int = 1;
	}

	if( parsed_run_files ){
		cout << "Building runnum vector. " << endl;
		for(int i = 0; i < lookup_parsed_runs_cnt(run_target.Data(), kine, sbsfieldscale); i++){
			runnum_vec.push_back(lookup_parsed_runnums(run_target.Data(), kine, sbsfieldscale, i));
		}
	}


	if( !calc_W ){
		outfile = new TFile(Form("rootfiles/MC_SBS%i_%s_mag%i_prime_calibration_15_08_2023.root", kine, run_target.Data(), sbsfieldscale), "RECREATE");
	}
	if( calc_W ){
		outfile = new TFile(Form("rootfiles/MC_SBS%i_%s_mag%i_full_calibration_15_08_2023.root", kine, run_target.Data(), sbsfieldscale), "RECREATE");
	}		


	//Define Histograms:
	h_simc_Ep = new TH1D("h_simc_Ep", Form("E/p (E/p > 0)- SBS%i = %i%%, %s, %s", kine, sbsfieldscale, run_target.Data(), runs_string.Data()), 200, 0, 2);
	h_simc_PS = new TH1D("h_simc_PS", Form("Pre-Shower Clus. E (PS > 0) - SBS%i = %i%%, %s, %s; GeV", kine, sbsfieldscale, run_target.Data(), runs_string.Data()), 300, 0, 3);
	h_simc_HCal_e = new TH1D("h_simc_HCal_e", Form("HCal Clus. E (HCal_E > 0) - SBS%i = %i%%, %s, %s; GeV", kine, sbsfieldscale, run_target.Data(), runs_string.Data()), 200, 0, 0.4);
	h_simc_SHPS = new TH1D("h_simc_SHPS", Form("SH + PS Clus. E (SHPS > 0)- SBS%i = %i%%, %s, %s; GeV", kine, sbsfieldscale, run_target.Data(), runs_string.Data()), 500, 0, 5);
	h_simc_W = new TH1D( "h_simc_W", Form("Invariant Mass, W (No Cuts) - SBS%i = %i%%, %s, %s; GeV", kine, sbsfieldscale, run_target.Data(), runs_string.Data()), 200, 0.0, 2.0 );
	h_simc_Wrecon = new TH1D( "h_simc_Wrecon", Form("Recon. Invariant Mass, W  - SBS%i = %i%%, %s, %s; GeV", kine, sbsfieldscale, run_target.Data(), runs_string.Data()), 300, 0.0, 3.0 );
	h_simc_W2 = new TH1D( "h_simc_W2", Form("Invariant Mass, W^{2} - SBS%i = %i%%, %s, %s; GeV", kine, sbsfieldscale, run_target.Data(), runs_string.Data()), 400, 0.0, 4.0 );
	h_simc_W2recon = new TH1D( "h_simc_W2recon", Form("Reconstructed W^{2}s - SBS%i = %i%%, %s, %s; GeV", kine, sbsfieldscale, run_target.Data(), runs_string.Data()), 400, 0.0, 4.0 );
	h_simc_vert = new TH1D( "h_simc_vert", "Vertex Position; m", 200, -.1, .1 );	
	h_simc_Q2 = new TH1D("h_Q2", Form("Momentum Transfer Q^2 - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 750, 0.5, 9.0);

	h_simc_dxdy_wcut = new TH2D("h_simc_dxdy_wcut",Form("HCal dxdy (W cuts) - SBS%i = %i%%, %s, %s;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data(), runs_string.Data()), 400, -2, 2, 400, -2.0, 2.0 );
	h_simc_dxdy_all = new TH2D("h_simc_dxdy_all",Form("HCal dxdy (NO CUTS) - SBS%i = %i%%, %s, %s;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data(), runs_string.Data()),400,-2,2,400,-2,2);
	h_simc_dx = new TH1D( "h_simc_dx", Form("HCal dx - SBS%i = %i%%, %s, %s; m", kine, sbsfieldscale, run_target.Data(), runs_string.Data()), 400, -2.0, 2.0 );
	h_simc_dy = new TH1D( "h_simc_dy", Form("HCal dy - SBS%i = %i%%, %s, %s; m", kine, sbsfieldscale, run_target.Data(), runs_string.Data()), 400, -2, 2 );
	h_simc_Y = new TH1D( "h_simc_Y", Form("HCal Y - SBS%i = %i%%, %s, %s; m", kine, sbsfieldscale, run_target.Data(), runs_string.Data()), 400, -2, 2 );
	h_simc_dx_wcut = new TH1D("h_simc_dx_wcut",Form("dx (W cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);

	h_simc_dx = new TH1D("h_simc_dx",Form("dx (NO CUTS) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);
	h_simc_dx_p = new TH1D("h_simc_dx_p",Form("dx - proton only -(NO CUTS) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);
	h_simc_dx_n = new TH1D("h_simc_dx_n",Form("dx - neutron only -(NO CUTS) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);
	h_simc_dx_p_wcut = new TH1D("h_simc_dx_p_wcut",Form("dx - proton only - (W cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);
	h_simc_dx_n_wcut = new TH1D("h_simc_dx_n_wcut",Form("dx - neutron only - (W cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);
	h_simc_dx_p_fcut = new TH1D("h_simc_dx_p_fcut",Form("dx - proton only - (F cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);
	h_simc_dx_n_fcut = new TH1D("h_simc_dx_n_fcut",Form("dx - neutron only - (F cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);
	h_simc_dx_p_wcut_fcut = new TH1D("h_simc_dx_p_wcut_fcut",Form("dx - proton only - (W & F cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);
	h_simc_dx_n_wcut_fcut = new TH1D("h_simc_dx_n_wcut_fcut",Form("dx - neutron only - (W & F cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);

	h_simc_dx_fcut = new TH1D("h_simc_dx_fcut",Form("dx (f cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);
	h_simc_dx_wcut_fcut = new TH1D("h_simc_dx_wcut_fcut",Form("dx (W & f cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);
	h_simc_dy = new TH1D("h_simc_dy",Form("dy (NO CUTS) - SBS%i = %i%%, %s; y_{HCal} - y_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);
	h_simc_dy_wcut = new TH1D("h_simc_dy_wcut",Form("dy (W Cuts) - SBS%i = %i%%, %s; y_{HCal} - y_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 250, -1.25, 1.25);  

	h_simc_dxdy = new TH2D("h_simc_dxdy", Form("Hadron Spot(s) on HCal (NO CUTS) - SBS%i = %i%%, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 300, -1.5, 1.5, 500, -2.5, 2.5 );
	h_simc_dxdy_wcut = new TH2D("h_simc_dxdy_wcut", Form("Hadron Spot(s) on HCal (W cut) - SBS%i = %i%%, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 300, -1.5, 1.5, 500, -2.5, 2.5 );
	h_simc_dxdy_fcut = new TH2D("h_simc_dxdy_fcut", Form("Hadron Spot(s) on HCal (f cut) - SBS%i = %i%%, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 250, -1.25, 1.25, 500, -2.5, 2.5 );

	HCal_dist = lookup_HCal_dist_by_kine(kine); 	//Distace from HCal face to target chamber
	HCal_theta = (pi/180.0)*lookup_HCal_angle_by_kine(kine, "deg");		//Theta angle for HCal from downstream beamline
	cout << endl << "-----------------------------------" << endl;
	cout << "Run parameters: " << endl;
	cout << "Run: " << runnum << endl;
	cout << "Kinematic: SBS" << kine << endl;
	cout << "Target: " << run_target.Data() << endl;
	cout << "Beam Energy: " << E_beam << endl;
	cout << "SBS Field: " << sbsfieldscale << "%" << endl;
	cout << "-----------------------------------" << endl;
	cout << "BB angle: " << lookup_BB_angle_by_kine(kine, "deg") << endl;
	cout << "BB distance: " << lookup_BB_dist(kine) << endl;
	cout << "SBS angle: " << lookup_SBS_angle_by_kine(kine, "deg") << endl;
	cout << "HCal angle: " << HCal_theta*TMath::RadToDeg() << endl;
	cout << "HCal distance: " << HCal_dist << endl;
	cout << "-----------------------------------" << endl << endl;;


	// rootfile_dir = Form("/lustre19/expphy/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass%i/SBS%i/%s/rootfiles", pass, kine, run_target.Data());
	rootfile_dir = "/volatile/halla/sbs/jboyd/simulation/out_dir/MC_REPLAY_OUT_DIR/simc";
	// rootfile_dir = "/volatile/halla/sbs/jboyd/simulation/Rootfiles/";


	if( single_run ){
		cout << "Running in single run mode: " << run_target.Data() <<" SBS" << kine << "mag" << sbsfieldscale << endl;
		cout << "--------------------------------------" << endl;
		cout << "Adding files to TChain from: " << rootfile_dir.Data() << endl;

		// TC->Add(Form("%s/*%i*.root", rootfile_dir.Data(), runnum));	
		// TC->Add("/volatile/halla/sbs/jboyd/simulation/out_dir/MC_REPLAY_OUT_DIR/replayed_jb_gmn_SBS8_Mag70_500k.root");
		// TC->Add(Form("%s/replayed_gmn_sbs4_ld2_30p_job_*", rootfile_dir.Data()));
		// TC->Add(Form("%s/replayed_jb_gmn_SBS4_LD2_job*", rootfile_dir.Data()));	
		// TC->Add(Form("%s/replayed_jb_gmn_SBS7_LH2_mag70_1000uA_elas_100k_job*", rootfile_dir.Data()));	
		// TC->Add(Form("%s/replayed_jb_gmn_SBS8_LD2_mag70mod3_500uA_elas_100k_19_07_2023_job*", rootfile_dir.Data()));
		
		//Proton file:
		 
		if( hadron_select == "proton" || hadron_select == "both" ){
			if( kine == 9 ){
				proton_infile = "replayed_digitized_jb_FARM_SIMC_gmn_SBS9_LD2_proton_mag70port0301_1200uA_elas_250k_19_09_2023_job*.root";
			}

			if( kine == 8 ){
				proton_infile = "replayed_digitized_jb_FARM_SIMC_gmn_SBS8_LD2_proton_mag70port0301_500uA_elas_250k_24_09_2023_job*.root";
			}
			//proton_infile = "replayed_digitized_jb_FARM_SIMC_gmn_SBS8_LD2_proton_mag70port0308_500uA_elas_250k_20_08_2023_job*.root";

			cout << "Adding proton infile: " << proton_infile.Data() << endl;
			TC->Add(Form("%s/%s", rootfile_dir.Data(), proton_infile.Data()));
		}

		//Neutron file:
		if( hadron_select == "neutron" || hadron_select == "both" ){
			if( kine == 9 ){
				neutron_infile = "replayed_digitized_jb_FARM_SIMC_gmn_SBS9_LD2_neutron_mag70port0314_1200uA_elas_250k_19_09_2023_job*.root";
			}

			if(kine == 8 ){
				neutron_infile = "replayed_digitized_jb_FARM_SIMC_gmn_SBS8_LD2_neutron_mag70port0301_500uA_elas_250k_24_09_2023_job*.root";
			}
			// neutron_infile = "replayed_digitized_jb_FARM_SIMC_gmn_SBS8_LD2_neutron_mag70port0308_500uA_elas_250k_20_08_2023_job*.root";

			cout << "Adding neutron infile: " << neutron_infile.Data() << endl;

			TC->Add(Form("%s/%s", rootfile_dir.Data(),neutron_infile.Data()));
		}

		runs_string = Form("Run: %i", runnum);
	}
	cout << "--------------------------------------" << endl;
	if( multi_run ){
		runnum = runnum_vec[0];
		cout << "Running in multi-run mode for runs: "  << endl;
		for(size_t run = 0; run < runnum_vec.size(); run++){
			cout << runnum_vec[run] << " ";
		}
		cout << endl;
		cout << "--------------------------------------" << endl;
		cout << "Adding files to TChain from: " << rootfile_dir.Data() << endl;
		for(size_t run = 0; run < runnum_vec.size(); run++){

			TC->Add(Form("%s/*%i*.root", rootfile_dir.Data(), runnum_vec[run]));

		}
		runs_string = Form("Runs: %i thru %i", *min_element(runnum_vec.begin(), runnum_vec.end()), *max_element(runnum_vec.begin(), runnum_vec.end()));
	}
	cout << "Finished adding files to TChain. " << endl;
	cout << "--------------------------------------" << endl;

	// if( calc_W || plot_dxdy ){
	W2_mean = lookup_simc_cut(kine, sbsfieldscale, run_target, "W2");
	W2_sigma = lookup_simc_cut(kine, sbsfieldscale, run_target, "W2_sigma");
	// }
	// else{
	// 	W2_mean = 0.92;
	// 	W2_sigma = 0.3;
	// }

	//DEFINE CUTS

	if( !calc_W ){
		master_cut_vec = {
			"sbs.hcal.nclus>0",
			"bb.ps.nclus>0",
			"bb.sh.nclus>0",
			"abs(bb.tr.vz[0])<=0.075",
			"bb.gem.track.nhits[0]>3",
			"bb.tr.n==1",
		};
	}
	if( calc_W ){
		master_cut_vec = {
			"sbs.hcal.nclus>0",
			"bb.ps.nclus>0",
			"bb.sh.nclus>0",
			"abs(bb.tr.vz)<0.08",
			"bb.gem.track.nhits[0]>3",
			"bb.tr.n==1",
			// Form("bb.ps.e>%f", lookup_cut(runnum, "PS_clus_e_cut")),
			// Form("((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))>(%f)&&((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))<(%f)", lookup_cut(runnum, "Ep") - (lookup_cut(runnum, "Ep_sigma")), lookup_cut(runnum, "Ep") + (lookup_cut(runnum, "Ep_sigma"))),
			Form("((abs(((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))))-%f)<%f", lookup_simc_cut(kine, sbsfieldscale, run_target, "Ep"), 2.0*lookup_simc_cut(kine, sbsfieldscale, run_target, "Ep_sigma")),
			// Form("sbs.hcal.e>%f", lookup_cut(runnum, "HCal_clus_e_cut")),
			Form("(bb.sh.e+bb.ps.e)>%f", lookup_simc_cut(kine, sbsfieldscale, run_target, "SH_PS_mean") - 3.0*lookup_simc_cut(kine, sbsfieldscale, run_target, "SH_PS_sigma")),
			//11449:
			Form("bb.ps.e>%f", 0.15),
			// Form("((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))>(%f)&&((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))<(%f)", 0.75, 1.2),
			Form("sbs.hcal.e>%f", 0.01)
			// Form("(bb.sh.e+bb.ps.e)>%f", 1.25)
		};
	}
	if( use_heavy_cut ){
		master_cut_vec = {
			"sbs.hcal.nclus>0",
			"bb.ps.nclus>0",
			"bb.sh.nclus>0",
			"abs(bb.tr.vz[0])<=0.075",
			"bb.gem.track.nhits[0]>4",
			"bb.tr.n==1",
			// Form("bb.ps.e>%f", 0.3),
			// Form("((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))>(%f)&&((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))<(%f)", lookup_cut(runnum, "Ep") - (lookup_cut(runnum, "Ep_sigma")/4), lookup_cut(runnum, "Ep") + (lookup_cut(runnum, "Ep_sigma")/4)),
			// Form("sbs.hcal.e>%f", 0.03),
			// Form("(bb.sh.e+bb.ps.e)>%f", lookup_cut(runnum, "SH_PS_clus_e_cut") - lookup_cut(runnum, "SH_PS_sigma")/6)
	

		};
	}
	for(size_t cut = 0; cut < master_cut_vec.size(); cut++){
		if(cut == master_cut_vec.size() - 1){
			master_cut_string.Append(Form("%s", master_cut_vec[cut].Data()));
		}
		else{
			master_cut_string.Append(Form("%s%s", master_cut_vec[cut].Data(), "&&"));
		}
	}

	master_cut = Form("%s", master_cut_string.Data());

	cout << "--------------------------------------" << endl;
	cout << "Setting up branches... ";
	//Setup Branches
	TC->SetMakeClass(1); //Allows for viewing of general detector params from Event_Branch
	TC->SetBranchStatus( "*", 0 );

	//MC values for normalization
	TC->SetBranchStatus( "MC.mc_omega", 1 );
	TC->SetBranchStatus( "MC.simc_sigma", 1 );
	TC->SetBranchStatus( "MC.simc_fnucl", 1);
	TC->SetBranchStatus( "MC.simc_Weight", 1);

	TC->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );  
	TC->SetBranchStatus( "sbs.hcal.clus_blk.row", 1 );
	TC->SetBranchStatus( "sbs.hcal.clus_blk.col", 1 );
	TC->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
	TC->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
	TC->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
	TC->SetBranchStatus( "sbs.hcal.x", 1 );
	TC->SetBranchStatus( "sbs.hcal.y", 1 );
	TC->SetBranchStatus( "sbs.hcal.e", 1 );
	TC->SetBranchStatus( "sbs.hcal.nblk", 1 );
	TC->SetBranchStatus( "sbs.hcal.nclus", 1 );
	TC->SetBranchStatus( "bb.sh.nclus", 1 );
	TC->SetBranchStatus( "bb.ps.nclus", 1 );
	TC->SetBranchStatus( "bb.tr.n", 1 );
	TC->SetBranchStatus( "bb.tr.px", 1 );
	TC->SetBranchStatus( "bb.tr.py", 1 );
	TC->SetBranchStatus( "bb.tr.pz", 1 );
	TC->SetBranchStatus( "bb.tr.p", 1 );
	TC->SetBranchStatus( "bb.tr.vz", 1 );
	TC->SetBranchStatus( "bb.ps.e", 1 );
	TC->SetBranchStatus( "bb.ps.x", 1 );
	TC->SetBranchStatus( "bb.ps.y", 1 );
	TC->SetBranchStatus( "bb.sh.e", 1 );
	TC->SetBranchStatus( "bb.sh.x", 1 );
	TC->SetBranchStatus( "bb.sh.y", 1 );
	// TC->SetBranchStatus( "bb.tdctrig.tdc", 1 );
	// TC->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
	// TC->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );
	TC->SetBranchStatus( "bb.gem.track.nhits", 1);

	//MC normalization variables
	TC->SetBranchAddress( "MC.mc_omega", &mc_omega );
	TC->SetBranchAddress( "MC.simc_sigma", &mc_sigma );
	TC->SetBranchAddress( "MC.simc_fnucl", &mc_fnucl );
	TC->SetBranchAddress( "MC.simc_Weight", &mc_Weight );

	TC->SetBranchAddress( "sbs.hcal.clus_blk.atime", atime );
	TC->SetBranchAddress( "sbs.hcal.clus_blk.row", row );
	TC->SetBranchAddress( "sbs.hcal.clus_blk.col", col );
	TC->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", tdctime );
	TC->SetBranchAddress( "sbs.hcal.clus_blk.id", cblkid );
	TC->SetBranchAddress( "sbs.hcal.clus_blk.e", cblke );
	TC->SetBranchAddress( "sbs.hcal.x", &hcal_x );
	TC->SetBranchAddress( "sbs.hcal.y", &hcal_y );
	TC->SetBranchAddress( "sbs.hcal.e", &hcal_e );
	TC->SetBranchAddress( "sbs.hcal.nblk", &nblk );
	TC->SetBranchAddress( "sbs.hcal.nclus", &nclus );
	TC->SetBranchAddress( "bb.sh.nclus", &SH_nclus );
	TC->SetBranchAddress( "bb.ps.nclus", &PS_nclus );
	TC->SetBranchAddress( "bb.tr.n", &bb_tr_n );
	TC->SetBranchAddress( "bb.tr.px", bb_tr_px );
	TC->SetBranchAddress( "bb.tr.py", bb_tr_py );
	TC->SetBranchAddress( "bb.tr.pz", bb_tr_pz );
	TC->SetBranchAddress( "bb.tr.vz", bb_tr_vz );
	TC->SetBranchAddress( "bb.tr.p", bb_tr_p );
	TC->SetBranchAddress( "bb.ps.e", &bb_ps_e );
	TC->SetBranchAddress( "bb.ps.x", &bb_ps_x );
	TC->SetBranchAddress( "bb.ps.y", &bb_ps_y );
	TC->SetBranchAddress( "bb.sh.e", &bb_sh_e );
	TC->SetBranchAddress( "bb.sh.x", &bb_sh_x );
	TC->SetBranchAddress( "bb.sh.y", &bb_sh_y );
	// TC->SetBranchAddress( "bb.tdctrig.tdcelemID", TDCT_id );
	// TC->SetBranchAddress( "bb.tdctrig.tdc", TDCT_tdc );
	// TC->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &TDCTndata );
	cout << " done. " << endl;
	cout << "--------------------------------------" << endl;
	cout << "Applying master cut: " << endl;
	cout << master_cut_string.Data() << endl;
	cout << "--------------------------------------" << endl;

	TEventList *ev_list = new TEventList("ev_list", "Elastic Event List");

//Event list with cuts:
	TC->Draw(">>ev_list");
	cout << "--------------------------------------" << endl;
	cout << "--------------------------------------" << endl;
	
	Nevents = ev_list->GetN();

	cout << "Number of events to analyze: " << Nevents << endl;
	cout << "--------------------------------------" << endl;
	cout << "--------------------------------------" << endl;
	cout << "Starting analysis loop on events..... " << endl;

	//Beam and detector vars including definition of HCal coordinate axis
	TLorentzVector Pbeam( 0, 0, E_beam, E_beam );
	TLorentzVector Ptarg( 0, 0, 0, 0.5*(Mp+Mn) ); //Average of proton and neutron rest mass
	TVector3 hcal_origin( -HCal_dist*sin(HCal_theta), 0, HCal_dist*cos(HCal_theta));
	TVector3 hcal_zaxis = hcal_origin.Unit();
	TVector3 hcal_xaxis(0,-1,0);
	TVector3 hcal_yaxis = hcal_zaxis.Cross( hcal_xaxis ).Unit();

	Int_t elastic_yield = 0;

	int watch_cnt = 0;
	int five_percent = int(0.05*Nevents);
	vector<double> time_for_five;
	double average_time = 0.0, time_remaining;
	StopWatch->Start();

	for(Long64_t nevent = 0; nevent < Nevents; nevent++){
		TC->GetEntry( ev_list->GetEntry( nevent ));
		if( plot_dxdy ){
			elastic_yield_str = Form("Elastic yield = %i.", elastic_yield);
		}

		if( nevent%five_percent == 0){
			StopWatch->Stop();

			if( watch_cnt == 0){
				cout << "Evt: " << nevent <<"/" << Nevents << Form("(%.0f/100%%)", 100.0*double(1.0*nevent/Nevents)) << ". " << elastic_yield_str.Data() << endl;
			}

			if( watch_cnt > 0 ){
				time_for_five.push_back(StopWatch->RealTime());	
				average_time = VectorMean(time_for_five)/60;
				time_remaining = average_time*100*( 1.0 - double(nevent)/double(Nevents));
				cout << "Evt: " << nevent <<"/" << Nevents << Form("(%.0f/100%%)", 100.0*double(1.0*nevent/Nevents)) << ". " << elastic_yield_str.Data() << " Time left: " << time_remaining << " minutes." << endl;
			}
			StopWatch->Reset();
			StopWatch->Continue();
			watch_cnt++;
		}
		cout.flush();

		if( int(mc_fnucl) == 0 ){ 
			is_n = true; 
			is_p = false;
		}
		if( int(mc_fnucl) == 1 ){ 
			is_p = true; 
			is_n = false;
		}

		pn_weight = 0.0;
		pn_weight = mc_Weight;
	
		double bbcal_time=0.0, hcal_time=0.0;

		for(int ihit=0; ihit<TDCTndata; ihit++){
			if(TDCT_id[ihit]==5){
				bbcal_time=TDCT_tdc[ihit];
			}
			if(TDCT_id[ihit]==0){
				hcal_time=TDCT_tdc[ihit];
			}
		}

		double diff = hcal_time - bbcal_time; 
		// if( fabs(diff - tdiff)>tdiff_max ){
		// 	continue;
		// }

		//Sanity check
      	// if( (int)bb_tr_n!=1 ){
      	// 	cout << "**************************************************************" << endl;
      	// 	cout << "--------------------------------------------------------------" << endl;
      	// 	cout << endl << endl << "WARNING: Total tracks not as expected from global cut. Check globalcut for errors." << endl << endl;
      	// 	cout << "--------------------------------------------------------------" << endl;
      	// 	cout << "**************************************************************" << endl;
      	// } 
	
		Double_t Ep = (bb_ps_e + bb_sh_e)/(bb_tr_p[0]);
		Double_t PS = bb_ps_e;
		Double_t SHPS = bb_ps_e + bb_sh_e;
		Double_t HCal_e = hcal_e;

		if( Ep > 0.0 ){
			h_simc_Ep->Fill(Ep);			
		}
		if( PS > 0.0 ){
			h_simc_PS->Fill(PS);		
		}
		if( SHPS > 0.0 ){
			h_simc_SHPS->Fill(SHPS);			
		}
		if( hcal_e > 0.0 ){
			h_simc_HCal_e->Fill(hcal_e);		
		}


//----------------------------------------------------
//---------------INVARIANT MASS -----------------------
		if( calc_W ){
			TLorentzVector kprime( bb_tr_px[0], bb_tr_py[0], bb_tr_pz[0], bb_tr_p[0] );
			TLorentzVector q = Pbeam - kprime; //Standard q-vector
			TVector3 vertex( 0, 0, bb_tr_vz[0] ); //Location of scattering event in target
			TVector3 qunit = q.Vect().Unit(); //q-vector direction

			Double_t sintersect = (hcal_origin-vertex).Dot( hcal_zaxis )/qunit.Dot( hcal_zaxis );
			TVector3 hcal_intersect = vertex + sintersect*qunit;
			Double_t x_expect = hcal_intersect.Dot( hcal_xaxis );
			Double_t y_expect = hcal_intersect.Dot( hcal_yaxis );

			Double_t W2recon = (Ptarg + q).M2();
			Double_t Wrecon = (Ptarg + q).M();
			Double_t E_ep = bb_tr_p[0]; // Obtain the scattered electron energy, neglect mass e
			Double_t p_ep = bb_tr_p[0]; // Obtain the magnitude of scattered electron momentum
			Double_t Q2 = 2*E_beam*E_ep*( 1-(bb_tr_pz[0]/p_ep) );
			Double_t nu = E_beam-E_ep; // Obtain energy transfer
			Double_t W2 = pow( Mp,2 )+2*Mp*nu-Q2; // Obtain W2 from Q2 and nu
			Double_t W = sqrt(W2);

		    h_simc_vert->Fill(bb_tr_vz[0]);
		    h_simc_W2->Fill(W2, pn_weight);
		    h_simc_W2recon->Fill(W2recon, pn_weight);

		    h_simc_W->Fill(W, pn_weight);
		    h_simc_Wrecon->Fill(Wrecon, pn_weight);

		    dx = hcal_x - x_expect;
		    dy = hcal_y - y_expect;

			if( plot_dxdy ){
				h_simc_dxdy_all->Fill( dy, dx, pn_weight );
			}

			if( plot_dxdy ){

			// Preliminary HCal projections with single cut on W
				if( fabs(W2 - W2_mean) < W2_sigma ){
					h_simc_dx_wcut->Fill( dx, pn_weight );
					h_simc_dy_wcut->Fill ( dy, pn_weight );
					h_simc_dxdy_wcut->Fill( dy, dx, pn_weight );

					//Histogram by hadron
					if( is_p ){
						h_simc_dx_p_wcut->Fill(dx, pn_weight);
					}
					if( is_n ){
						h_simc_dx_n_wcut->Fill(dx, pn_weight);
					}
				}

				h_simc_dx->Fill( dx, pn_weight );
				h_simc_dy->Fill( dy, pn_weight );
				h_simc_Y->Fill( hcal_y, pn_weight );
			}
			
			elastic_yield++;

		}

	}

	cout << "---------------------------------------" << endl;
	cout << "-----Finished going through events-----" << endl;
	cout << "---------------------------------------" << endl;
	cout << "--- Writing to output file and saving ---" << endl;
	cout << "---------------------------------------" << endl;

	outfile->Write();
	outfile->Close();

	cout << "---------------------------------------" << endl;
	cout << "       Opening written file and plotting" << endl;
	cout << "---------------------------------------" << endl;

	TFile *infile;

	if( !calc_W ){
		infile = new TFile(outfile->GetName(), "READ");
	}

	if( calc_W ){
		infile = new TFile(outfile->GetName(), "READ");
		hin_W2recon = static_cast<TH1D*>(infile->Get("h_simc_W2recon"));
		hin_W2 = static_cast<TH1D*>(infile->Get("h_simc_W2"));
		hin_Wrecon = static_cast<TH1D*>(infile->Get("h_simc_Wrecon"));
		hin_W = static_cast<TH1D*>(infile->Get("h_simc_W"));
	}
	if( plot_dxdy){
		hin_dxdy_wcut = static_cast<TH2D*>(infile->Get("h_simc_dxdy_wcut"));
		hin_dxdy_all = static_cast<TH2D*>(infile->Get("h_simc_dxdy_all"));
	}


	hin_Ep = static_cast<TH1D*>(infile->Get("h_simc_Ep"));
	hin_PS = static_cast<TH1D*>(infile->Get("h_simc_PS"));
	hin_SHPS = static_cast<TH1D*>(infile->Get("h_simc_SHPS"));
	hin_HCal_e = static_cast<TH1D*>(infile->Get("h_simc_HCal_e"));

//---------------------------------------
//---------------- Ep -----------------------

	TCanvas *c_Ep = new TCanvas("c_Ep", "E/p", 600, 500);
	hin_Ep->Draw();

	if( !calc_W || calc_W){
		TF1 *fit_Ep = new TF1("fit_Ep", fit_gaus, 0.5, 1.5, 3);
		fit_Ep->SetParName(0, "Ep Norm");
		fit_Ep->SetParName(1, "Ep Center");
		fit_Ep->SetParName(2, "Ep Sigma");

		fit_Ep->SetParLimits(0, 0, hin_Ep->GetMaximum());
		fit_Ep->SetParLimits(1, 0.95, 1.);
		fit_Ep->SetParLimits(2, 0.01, 0.4);

		hin_Ep->Fit("fit_Ep", "R");
		fit_Ep->Draw("same");
		Ep_center = fit_Ep->GetParameter(1);
		Ep_sigma = fit_Ep->GetParameter(2);

		par[0] = 0.0; par[1] = 0.0, par[2] = 0.0;
	}


//---------------------------------------
//---------------- PS -----------------------

	TCanvas *c_PS = new TCanvas("c_PS", "PS", 600, 500);
	hin_PS->Draw();

	if( !calc_W || calc_W){
		TF1 *fit_PS = new TF1("fit_PS", fit_gaus, 0.15, 3, 3);
		fit_PS->SetParName(0, "PS Norm");
		fit_PS->SetParName(1, "PS Center");
		fit_PS->SetParName(2, "PS Sigma");

		fit_PS->SetParLimits(0, 0, hin_PS->GetMaximum());
		fit_PS->SetParLimits(1, 0.6, 1.1);
		fit_PS->SetParLimits(2, 0.01, 0.8);

		hin_PS->Fit("fit_PS", "R");
		// fit_PS->Draw("same");
		PS_center = fit_PS->GetParameter(1);
		PS_sigma = fit_PS->GetParameter(2);

		par[0] = 0.0; par[1] = 0.0, par[2] = 0.0;
	

		cout << "Setting range of PS to : " << (3.0/300.0)*hin_PS->FindFirstBinAbove((0.75)*(hin_PS->GetMaximum())) << " - " << PS_center << endl;
		hin_PS->GetXaxis()->SetRangeUser((3.0/300.0)*hin_PS->FindFirstBinAbove((0.75)*(hin_PS->GetMaximum())), PS_center);
		cout << "---------------------------------------" << endl;
		PS_min = (3.0/300.0)*hin_PS->GetMinimumBin();
		cout << "PS_min before rounding: " << PS_min << endl;
		// PS_min = ceil(100*PS_min)/100.0;
		cout << "PS_min set to: " << PS_min << endl;
		cout << "---------------------------------------" << endl;

		hin_PS->GetXaxis()->SetRangeUser(0, 3);
		TLine *tl_PS = new TLine(PS_min, 0, PS_min, hin_PS->GetMaximum());
		tl_PS->SetLineColor(6);
		tl_PS->Draw("same");
	}


//---------------------------------------
//---------------- SHPS -----------------------
	TCanvas *c_SHPS = new TCanvas("c_SHPS", "SHPS", 600, 500);
	hin_SHPS->Draw();
	if( !calc_W || calc_W){
		TF1 *fit_SHPS = new TF1("fit_SHPS", fit_gaus, 0.5, 5, 3);
		fit_SHPS->SetParName(0, "SHPS Norm");
		fit_SHPS->SetParName(1, "SHPS Center");
		fit_SHPS->SetParName(2, "SHPS Sigma");

		fit_SHPS->SetParLimits(0, 0, hin_SHPS->GetMaximum());
		fit_SHPS->SetParLimits(1, 1.93, 4.0);
		fit_SHPS->SetParLimits(2, 0.01, 0.23);

		hin_SHPS->Fit("fit_SHPS", "R");
		fit_SHPS->Draw("same");
		SHPS_center = fit_SHPS->GetParameter(1);
		SHPS_sigma = fit_SHPS->GetParameter(2);

		par[0] = 0.0; par[1] = 0.0, par[2] = 0.0;
	}

//---------------------------------------
//---------------- HCal_e -----------------------
	TCanvas *c_HCal_e = new TCanvas("c_HCal_e", "HCal_e", 600, 500);
	hin_HCal_e->Draw();
	if( !calc_W || calc_W){
		cout << "---------------------------------------" << endl;
		cout << "HCal_e_min before rounding: " << HCal_e_min << endl;
		HCal_e_min = (0.4/200.0)*hin_HCal_e->FindFirstBinAbove((0.15)*(hin_HCal_e->GetMaximum()));
		HCal_e_min = ceil(100*HCal_e_min)/100.0;
		cout << "HCal_e_min set to: " << HCal_e_min << endl;
		cout << "---------------------------------------" << endl;

		TLine *tl_HCal_e = new TLine(HCal_e_min, 0, HCal_e_min, hin_HCal_e->GetMaximum());
		tl_HCal_e->SetLineColor(6);
		tl_HCal_e->Draw("same");

		par[0] = 0.0; par[1] = 0.0, par[2] = 0.0;
	}
//---------------------------------------
//---------------- W2 Recon -----------------------	
	if( calc_W ){
		TCanvas *c_W2recon = new TCanvas("c_W2recon", "W2recon", 600, 500);
		hin_W2recon->Draw();
		TF1 *fit_W2recon = new TF1("fit_W2recon", fit_gaus, 0.5, 1.25, 3);
		fit_W2recon->SetParName(0, "W2recon Norm");
		fit_W2recon->SetParName(1, "W2recon Center");
		fit_W2recon->SetParName(2, "W2recon Sigma");

		fit_W2recon->SetParLimits(0, 0, hin_W2recon->GetMaximum());
		fit_W2recon->SetParLimits(1, 0.85, 1.08);
		fit_W2recon->SetParLimits(2, 0.01, 0.4);

		hin_W2recon->Fit("fit_W2recon", "R");
		fit_W2recon->Draw("same");
		W2_fit_center = fit_W2recon->GetParameter(1);
		W2_fit_sigma = fit_W2recon->GetParameter(2);

		hin_W2->SetLineColor(6);
		hin_W2->SetLineStyle(6);
		hin_W2->Draw("same");
		par[0] = 0.0; par[1] = 0.0, par[2] = 0.0;
	}

//---------------------------------------
//---------------- W Recon -----------------------	
	if( calc_W ){
		TCanvas *c_Wrecon = new TCanvas("c_Wrecon", "Wrecon", 600, 500);
		hin_Wrecon->Draw();
		TF1 *fit_Wrecon = new TF1("fit_Wrecon", fit_gaus, 0.5, 1.25, 3);
		fit_Wrecon->SetParName(0, "Wrecon Norm");
		fit_Wrecon->SetParName(1, "Wrecon Center");
		fit_Wrecon->SetParName(2, "Wrecon Sigma");

		fit_Wrecon->SetParLimits(0, 0, hin_Wrecon->GetMaximum());
		fit_Wrecon->SetParLimits(1, 0.94, 1.005);
		fit_Wrecon->SetParLimits(2, 0.01, 0.4);

		hin_Wrecon->Fit("fit_Wrecon", "R");
		fit_Wrecon->Draw("same");
		W_fit_center = fit_Wrecon->GetParameter(1);
		W_fit_sigma = fit_Wrecon->GetParameter(2);

		hin_W->SetLineColor(6);
		hin_W->SetLineStyle(6);
		hin_W->Draw("same");
		par[0] = 0.0; par[1] = 0.0, par[2] = 0.0;
	}

//---------------------------------------
//---------------- dxdy -----------------------	
	if( plot_dxdy ){
		TCanvas *c_dxdy_wcut = new TCanvas("c_dxdy_wcut", "c_dxdy_wcut", 600, 500);
		hin_dxdy_wcut->Draw("colz");
	
		TCanvas *c_dx_wcut = new TCanvas("c_dx_wcut", "c_dx_wcut", 600, 500);
		
		TH1D *h_simc_dx_wcut = hin_dxdy_wcut->ProjectionY();
		TH1D *h_simc_dy_wcut = hin_dxdy_wcut->ProjectionX();
		// hin_dxdy_wcut->ProjectionY()->Draw();
		h_simc_dx_wcut->Draw();

		double dx_p_min = -1.5;
		double dx_p_max = -0.8;
		double dx_p_sigma_max = 0.18;

		TF1 *tf_dx_wcut_p = new TF1("tf_dx_wcut_p", fit_gaus, dx_p_min, dx_p_max, 3);
		tf_dx_wcut_p->SetLineColor(6);
		tf_dx_wcut_p->SetParLimits(0, 0.9*h_simc_dx_wcut->GetMaximum(), h_simc_dx_wcut->GetMaximum());
		tf_dx_wcut_p->SetParLimits(1, dx_p_min, dx_p_max);
		tf_dx_wcut_p->SetParLimits(2, 0.01, dx_p_sigma_max);
		h_simc_dx_wcut->Fit("tf_dx_wcut_p", "R+");

		dxdy_vector.push_back(tf_dx_wcut_p->GetParameter(1));
		dxdy_vector.push_back(tf_dx_wcut_p->GetParameter(2));

		double dx_n_min = -0.7;
		double dx_n_max = 0.0;
		double dx_n_sigma_max = 0.17;

		TF1 *tf_dx_wcut_n = new TF1("tf_dx_wcut_n", fit_gaus, dx_n_min, dx_n_max, 3);
		tf_dx_wcut_n->SetLineColor(3);
		tf_dx_wcut_n->SetParLimits(0, 0.9*h_simc_dx_wcut->GetMaximum(), h_simc_dx_wcut->GetMaximum());
		tf_dx_wcut_n->SetParLimits(1, dx_n_min, dx_n_max);
		tf_dx_wcut_n->SetParLimits(2, 0.01, dx_n_sigma_max);
		h_simc_dx_wcut->Fit("tf_dx_wcut_n", "R+");

		dxdy_vector.push_back(tf_dx_wcut_n->GetParameter(1));
		dxdy_vector.push_back(tf_dx_wcut_n->GetParameter(2));

		TCanvas *c_dxdy_all = new TCanvas("c_dxdy_all", "c_dxdy_all", 600, 500);
		hin_dxdy_all->Draw("colz");
	
		TCanvas *c_dx_all = new TCanvas("c_dx_all", "c_dx_all", 600, 500);
		hin_dxdy_all->ProjectionY()->Draw();

		TCanvas *c_dy_wcut = new TCanvas("c_dy_wcut", "c_dy_wcut", 600, 500);
		h_simc_dy_wcut->Draw();

		double dy_min = -1.0;
		double dy_max = 1.0;
		double dy_sigma_max = 0.26;

		TF1 *tf_dy_wcut = new TF1("tf_dy_wcut", fit_gaus, dy_min, dy_max, 3);
		tf_dy_wcut->SetLineColor(9);
		tf_dy_wcut->SetParLimits(0, 0.9*h_simc_dy_wcut->GetMaximum(), h_simc_dy_wcut->GetMaximum());
		tf_dy_wcut->SetParLimits(1, dy_min, dy_max);
		tf_dy_wcut->SetParLimits(2, 0.01, dy_sigma_max);
		h_simc_dy_wcut->Fit("tf_dy_wcut", "R+");

		dxdy_vector.push_back(tf_dy_wcut->GetParameter(1));
		dxdy_vector.push_back(tf_dy_wcut->GetParameter(2));

	}	

//---------------------------------------
//---------------- CUT VECTOR -----------------------
// 	[0]target_int, [1]kine, [2]sbs_field [3]PS_clus_e_cut, [4]SH_PS_mean, [5]SH_PS_sigma, [6]HCal_clus_e_cut, [7]Ep, [8]Ep_sigma, [9]W2, [10]W2_sigma, [11]W, [12]W_sigma, 

	cuts_from_fits.push_back(target_int);
	cuts_from_fits.push_back(kine);	
	cuts_from_fits.push_back(sbsfieldscale);
	cuts_from_fits.push_back(PS_min);
	cuts_from_fits.push_back(SHPS_center);
	cuts_from_fits.push_back(SHPS_sigma);
	cuts_from_fits.push_back(HCal_e_min);
	cuts_from_fits.push_back(Ep_center);
	cuts_from_fits.push_back(Ep_sigma);
	if( calc_W ){
		cuts_from_fits.push_back(W2_fit_center);
		cuts_from_fits.push_back(W2_fit_sigma);
		cuts_from_fits.push_back(W_fit_center);
		cuts_from_fits.push_back(W_fit_sigma);
	}


	cout << "---------------------------------------" << endl;
	cout << "---------------------------------------" << endl << endl;	
	cout << "Vector with cuts: " << endl;
	cout << "{";
	for(size_t cut = 0; cut < cuts_from_fits.size(); cut++){
		if( cut < cuts_from_fits.size() -1 ){
			cout << cuts_from_fits[cut] << ", ";
		}
		if( cut == cuts_from_fits.size() -1 ){
			cout << cuts_from_fits[cut] << " }" << endl;
		}
	}


	// TFile *inputfile = new TFile(Form("rootfiles/%i_dxdy", runnum), "READ");



	// TH1D *h_simc_atime = new TH1D( "atime", Form("HCal ADC Time, All Channels - SBS = %i%%, %s, Run %i; ns", sbsfieldscale, run_target.Data(), runnum), 160, 0, 160 );
	// TH1D *h_simc_E = new TH1D( "E", Form("HCal Cluster E, All Channels - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 500, 0, 0.5 );
	// TH1D *h_simc_E_cut = new TH1D( "E_cut", Form("HCal Cluster E All Cuts, All Channels - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 100, 0, 0.5 );
	// TH1D *h_simc_E_exp = new TH1D( "E_exp", Form("Expected Energy Dep in HCal - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 100, 0, 0.2 );
	// TH1D *h_simc_vert = new TH1D( "vert", Form("Vertex Position - SBS = %i%%, %s, Run %i; m", sbsfieldscale, run_target.Data(), runnum), 200, -1.0, 1.0 );
	// TH1D *hdx = new TH1D( "dx", Form("HCal dx - SBS = %i%%, %s, Run %i; m", sbsfieldscale, run_target.Data(), runnum), 200, -4.0, 2.0 );
	// TH1D *hdy = new TH1D( "dy", Form("HCal dy - SBS = %i%%, %s, Run %i; m", sbsfieldscale, run_target.Data(), runnum), 100, -1.2, 1.2 );
	// TH1D *hKE_p = new TH1D( "KE_p", Form("Scattered Proton Kinetic Energy - SBS = %i%%, %s, Run %i", sbsfieldscale, run_target.Data(), runnum), 500, 0.0, 5.0 );
	// TH1D *hKElow = new TH1D( "KElow", Form("Lowest Elastic E Sampled in HCal - SBS = %i%%, %s, Run %i; GeV", sbsfieldscale, run_target.Data(), runnum), 500, 0.0, 0.2 );
	// TH1D *hDiff = new TH1D( "hDiff",Form("HCal time - BBCal time - SBS = %i%%, %s, Run %i; ns", sbsfieldscale, run_target.Data(), runnum), 1300, -500, 800 );
	// TH1D *hX = new TH1D( "X", Form("HCal X - SBS = %i%%, %s, Run %i; m", sbsfieldscale, run_target.Data(), runnum), 100, -4.0, 2.0 );
	// TH1D *hY = new TH1D( "Y", Form("HCal Y - SBS = %i%%, %s, Run %i; m", sbsfieldscale, run_target.Data(), runnum), 100, -1.2, 1.2 );
	// TH2D *h2_CvCh = new TH2D( "CvCh", Form("HCal Coeff Single Block Clusters - SBS = %i%%, %s, Run %i; channel, GeV", sbsfieldscale, run_target.Data(), runnum), kNcell, 0, kNcell, 200, 0, 1.0 );
	// TH2D *h2_EvCh = new TH2D( "EvCh", Form("HCal Cluster E Single Block Clusters - SBS = %i%%, %s, Run %i; channel, GeV", sbsfieldscale, run_target.Data(), runnum), kNcell, 0, kNcell, 50, 0, 0.5 );
	// TH2D *h2_dxdy = new TH2D("dxdy",Form("Delta Plot (q-proj. dxdy) - SBS = %i%%, %s, Run %i;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", sbsfieldscale, run_target.Data(), runnum), 50, -2, 2, 100, -4.0, 2.0 );
	// TH2D *h2_dxdy_all = new TH2D("dxdy_all",Form("Delta Plot (NO CUT) - SBS = %i%%, %s, Run %i;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", sbsfieldscale, run_target.Data(), runnum),125,-2,2,125,-4,6);
	// TH2D *h2_dxVE = new TH2D("dxVE",Form("dxVE - SBS = %i%%, %s, Run %i;x_{HCAL}-x_{expect} (m); GeV", sbsfieldscale, run_target.Data(), runnum), 100, -4.0, 2.0, 100, 0, 0.5 );
	// TH2D *h2_rowcol = new TH2D( "hrowcol", Form("HCal Block Position Elastics, HCal - SBS = %i%%, %s, Run %i; Col; Row", sbsfieldscale, run_target.Data(), runnum), kNcols, 0, kNcols, kNrows, -kNrows, 0 );
cout << "------------------------------------------------------------------"<< endl;
cout << "                       ANALYSIS FINISHED" << endl;
cout << "------------------------------------------------------------------"<< endl;
	cout << "Run parameters: " << endl;
	cout << "Kinematic: SBS" << kine << endl;
	cout << "Target: " << run_target.Data() << endl;
	cout << "Beam Energy: " << E_beam << endl;
	cout << "SBS Field: " << sbsfieldscale << "%" << endl;
	cout << "-----------------------------------" << endl;
	cout << "BB angle: " << lookup_BB_angle_by_kine(kine, "deg") << endl;
	cout << "SBS angle: " << lookup_SBS_angle_by_kine(kine, "deg") << endl;
	cout << "HCal angle: " << HCal_theta << endl;
	cout << "HCal distance: " << HCal_dist << endl;
	cout << "-----------------------------------" << endl << endl;
	cout << "Elastic yield: " << elastic_yield << endl << endl;
	cout << "---------------------------------------" << endl << endl;	
	cout << "Vector with cuts: " << endl;

	if( !calc_W ){
		cout << "{ kine, sbsfieldscale, PS, Sh_simc_PS_center, SH_PS_sigma, HCal, Ep, Ep_sigma, W2, W2_sigma, W, W_sigma }" << endl << endl;
		cout << "{ ";
		for(size_t cut = 0; cut < cuts_from_fits.size(); cut++){
			if( cut < cuts_from_fits.size() -1 ){
				cout << cuts_from_fits[cut] << ", ";
			}
			if( cut == cuts_from_fits.size() -1 ){
				cout << cuts_from_fits[cut] << " }" << endl;
			}
		}
	}
	if( calc_W ){
		cout << "{ kine, sbsfieldscale, PS, SH_PS_center, SH_PS_sigma, HCal, Ep, Ep_sigma, W2, W2_sigma, W, W_sigma }" << endl << endl;
		cout << "{";
		// cout << "SBS" << kine << ", " << lookup_simc_cut(kine, sbsfieldscale, run_target, "PS_clus_e_cut") << ", ";
		// cout << lookup_simc_cut(kine, sbsfieldscale, run_target, "SH_PSmean") << ", " << lookup_simc_cut(kine, sbsfieldscale, run_target, "SH_PS_sigma") << ", ";
		// cout << lookup_simc_cut(kine, sbsfieldscale, run_target, "HCal_clus_e_cut") << ", " << lookup_simc_cut(kine, sbsfieldscale, run_target, "Ep") << ", ";
		// cout << lookup_simc_cut(kine, sbsfieldscale, run_target, "Ep_sigma") << ", " << W2_fit_center <<", " << W2_fit_sigma << ", ";
		// cout << W_fit_center << ", " << W_fit_sigma << " }" << endl;
		for(size_t cut = 0; cut < cuts_from_fits.size(); cut++){
			if( cut < cuts_from_fits.size() -1 ){
				cout << cuts_from_fits[cut] << ", ";
			}
			if( cut == cuts_from_fits.size() -1 ){
				cout << cuts_from_fits[cut] << " }" << endl;
			}
		}	
	}

	if( plot_dxdy ){
		cout << "{ dx_p_mean, dx_p_sigma, dx_n_mean, dx_n_sigma, dy_mean, dy_sigma }" << endl;
		cout << "{";
		for( size_t dxdy_val = 0; dxdy_val < dxdy_vector.size(); dxdy_val++ ){
			if( dxdy_val < dxdy_vector.size() -1 ){
				cout << dxdy_vector[dxdy_val] << ", ";
			}
			if( dxdy_val == dxdy_vector.size() -1 ){
				cout << dxdy_vector[dxdy_val] << " }" << endl;
			}
		}
	}



	cout << "Output file created: " << endl;
	cout << outfile->GetName() << endl;

	auto total_time_end = high_resolution_clock::now();
	auto total_time_duration = duration_cast<minutes>(total_time_end - total_time_start);
	cout << "Total time for analysis: " << total_time_duration.count() << " minutes. " << endl;

}

