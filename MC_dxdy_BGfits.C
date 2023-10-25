//Jboyd 27 July 2023
//Script to help find the bg function (pol4) for the MC data

#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <string>
#include <chrono>
#include <TF1.h>
#include "TStopwatch.h"

using namespace std::chrono;
#include "/w/halla-scshelf2102/sbs/jboyd/include/include_files.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/GEM_lookups.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/beam_variables.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/calc_functions.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/MC_lookups.h"


Double_t par[3], par_full[11];
Bool_t reject;

Double_t background_with_reject(Double_t *x, Double_t *par_full){
	double x_min = -1.6;
	double x_max = 1.75;

	if( reject && (x[0] >= x_min || x[0] <= x_max) ){
		TF1::RejectPoint();
		return 0;
	}

	return par_full[0] + par_full[1]*x[0] + par_full[2]*x[0]*x[0] + par_full[3]*x[0]*x[0]*x[0] + par_full[4]*x[0]*x[0]*x[0]*x[0];

}

Double_t background(Double_t *x, Double_t *par_full){
	return par_full[0] + par_full[1]*x[0] + par_full[2]*x[0]*x[0] + par_full[3]*x[0]*x[0]*x[0] + par_full[4]*x[0]*x[0]*x[0]*x[0];
}

Double_t yieldPeaks(Double_t *x, Double_t *par){
	return par_full[0]*exp((-0.5)*pow(((x[0] -  par_full[1])/par_full[2]),2)) + par_full[3]*exp((-0.5)*pow(((x[0] -  par_full[4])/par_full[5]),2));
}


Double_t fullFitFunction(Double_t *x, Double_t *par_full){
	return yieldPeaks(x,par_full) + background(x, &par_full[6]);
}

Double_t fit_gaus(Double_t * x, Double_t *par){

	Double_t g = 0.0;

	g = par[0]*exp((-0.5)*pow(((x[0] -  par[1])/par[2]),2));
	
	return g;
}

bool fiducial_cut = true;
bool apply_fcut;
bool calc_pn_weight = true;

bool correct_beam_energy = false;
bool calibrate = true;
bool fit_full_function = true;
bool fit_peaks_only;

//RUN Info/Parameters
int kine = 8;
int sbsfieldscale = 70;
TString run_target = "LD2";
// double I_beam = 5.0;
// SBS8
double I_beam_uA = 5.00;

//SBS9
// double I_beam_uA = 12.00;

double I_beam;
TString I_beam_str;
double E_beam = lookup_beam_energy_from_kine( kine );

double SBS_field = 1.0*sbsfieldscale/100.0;
double ngen_total = 500000.0;
// double ngen_total = 4000000.0;

TString rootfile_dir;
TFile *outfile, *fits_outfile;
TString fits_outfile_name;
TChain *TC = new TChain("T");
TString infile;

vector<TString> master_cut_vec;
TString master_cut_string;
TCut master_cut = "";

//Experimental Parameters
	//BigBite, BBCal
	double BB_dist = lookup_BB_dist_by_kine( kine );  //Distance to BigBite magnet
	double BB_theta = lookup_BB_angle_by_kine( kine, "deg" ); //degrees, BB arm angle
	double BB_theta_rad = lookup_BB_angle_by_kine( kine, "rad" ); //radians, BB arm angle

	//SBS, HCal
	double HCal_dist = lookup_HCal_dist_by_kine( kine ); //Distance to Hcal face form target chamber
	double HCal_theta = lookup_HCal_angle_by_kine (kine, "deg" ); //degrees, Angle for downsream arm to HCal
	double HCal_theta_rad = lookup_HCal_angle_by_kine (kine, "rad" ); //radians, Angle for downsream arm to HCal

	double SBS_theta = lookup_SBS_angle_by_kine( kine, "deg" ); //degrees
	double SBS_theta_rad = lookup_SBS_angle_by_kine( kine, "rad" ); //radians
	double scint_intersect, x_expected_HCal, y_expected_HCal;

	//HCal Fiducial cut parameters:
	double hcal_y_fmin = -0.75;
	double hcal_y_fmax = 0.75;
	double hcal_x_fmin = -2.015;
	double hcal_x_fmax = 1.285;

//Scattered kinematics
double e_prime_theta; //Scattered electron angle, theta
double e_prime_phi; //Scattered electrong angle, phi
double p_el, nu, pp, p_ep, dx, dy;
double nucleon_theta, nucleon_phi;
double Q2, W, W2, E_nucleon, KE_p, E_pp, E_ep;
Double_t Wrecon, W2recon, x_recon_expect, y_recon_expect;
double W_mean, W_sigma, W2_mean, W2_sigma;

double pn_weight = 0;

//BEAM Variables, etc.
double p_Beam, E_loss_outgoing;
double Eloss, E_beam_final, p_corr;

///////    CONSTANTS  ////////
//Physical Constants
const double pi = TMath::Pi();
const double Mp = 0.938272; //Mass of proton [GeV]
const double Mn = 0.939565; //Mass of neutron [GeV]
const double Me = 0.00051; //Mass of electron [GeV]

//Static Detector Parameters
const int maxTracks = 1000; // Reasonable limit on tracks to be stored per event
const int maxTdcChan = 10; // Set to accomodate original 5 TDCTrig channels with buffer
// const double HCal_height = -0.2897; // Height of HCal above beamline
// const double HCal_height = 0.135; // Height of HCal above beamline
double HCal_height;

//Static Target Parameters
const double l_tgt = 0.15; // Length of the target (m)
const double rho_tgt = 0.0723; // Density of target (g/cc)
const double rho_Al = 2.7; // Density of aluminum windows (g/cc)
const double cell_diameter = 1.6*2.54; //cm, right now this is a guess
const double Ztgt = 1.0;
const double Atgt = 1.0;
const double Mmol_tgt = 1.008; //g/mol

//For energy-loss correction to beam energy:
const double dEdx_tgt=0.00574; //According to NIST ESTAR, the collisional stopping power of hydrogen is about 5.74 MeV*cm2/g at 2 GeV energy
const double dEdx_Al = 0.0021; //According to NIST ESTAR, the collisional stopping power of Aluminum is about 2.1 MeV*cm2/g between 1-4 GeV
const double uwallthick_LH2 = 0.0145; //cm
const double dwallthick_LH2 = 0.015; //cm
const double cellthick_LH2 = 0.02; //cm, this is a guess;
const double Alshieldthick = 2.54/8.0; //= 1/8 inch * 2.54 cm/inch
int useAlshield = 0;

//EXTRACTED VARIABLES
double n_integral, p_integral, n_center, n_sigma, p_center, p_sigma;
int n_counts, p_counts, elastic_yield;

//dxdy variables
double dx_p, dx_p_sigma, dy_p, dy_p_sigma, dx_n, dx_n_sigma, dy_n, dy_n_sigma, dx_pn_max;

//HISTOGRAMS
TH1D *h_Ep, *h_PS, *h_HCal_e, *h_SHPS; 

TH1D *h_W, *h_W2, *h_W_cut, *h_W_fcut, *h_vz_cut;
TH1D *h_Wrecon, *h_W2recon;
TH1D *h_KE_p, *h_KE_low, *h_X, *h_Y, *h_E_eloss;
TH1D *h_Q2, *h_E_ep, *h_E_pp;

TH1D *h_dx, *h_dx_cut, *h_dx_wcut, *h_dx_fcut, *h_dx_wcut_fcut;
TH1D *h_dy, *h_dy_cut, *h_dy_wcut;
TH1D *hin_dx_wcut, *hin_dy_wcut;
TH2D *h_dxdy, *h_dxdy_cut, *h_dxdy_wcut, *h_dxdy_ncut, *h_dxdy_pcut, *h_dxdy_fcut;
 
TH2D *h_E_ecorr_vs_vert;

TH2D *h_xy, *h_xy_cut, *h_xy_fcut, *h_xy_cut_p, *h_xy_cut_n, *h_PAngleCorr_theta, *h_PAngleCorr_phi;
TH1D *h_dx_BG;
TH2D *hin_dxdy_Main, *hin_dxdy_BeyondDyCut,*hin_dxdy_midCut;
double dy_min, dy_max;
int dy_min_bin, dy_max_bin;

//BRANCH VARIABLES

double bb_tr_p[maxTracks], bb_tr_px[maxTracks], bb_tr_py[maxTracks], bb_tr_pz[maxTracks];
double bb_tr_vx[maxTracks], bb_tr_vy[maxTracks], bb_tr_vz[maxTracks], bb_tr_chi2[maxTracks];
double bb_fp_x[maxTracks], bb_fp_y[maxTracks], bb_fp_th[maxTracks], bb_fp_ph[maxTracks];
double bb_tgt_x[maxTracks], bb_tgt_y[maxTracks], bb_tgt_th[maxTracks], bb_tgt_ph[maxTracks];
double bb_tr_n, bb_ps_x, bb_ps_y, bb_ps_e, bb_sh_x, bb_sh_y, bb_sh_e;
Double_t nblk, nclus, SH_nclus, PS_nclus, hcal_x, hcal_y, hcal_e;
double mc_omega, mc_sigma, mc_fnucl, luminosity;

Long64_t Nevents;

double dx_p_scale = 1.0;
double dx_n_scale = 1.0;

TCanvas *c_dx, *c_dy;
TF1 *tf_dx_bg, *fit_dx_dyBG, *fit_dx_full;
double BG_scale_factor;
double fit_dx_par1min, fit_dx_par1max, fit_dx_par2min, fit_dx_par2max, fit_dx_par4min, fit_dx_par4max;
double fit_dx_par5min, fit_dx_par5max, fit_dx_par6min, fit_dx_par6max, fit_dx_par7min, fit_dx_par7max;
double fit_dx_par8min, fit_dx_par8max, fit_dx_par9min, fit_dx_par9max, fit_dx_par10min, fit_dx_par10max;

void MC_dxdy_BGfits(){

	auto total_time_start = high_resolution_clock::now();
	TStopwatch *StopWatch = new TStopwatch();

	gStyle->SetPalette(55);
	cout << "--------------------------------------" << endl;
	cout << "Analysis started. " << endl;
	cout << "--------------------------------------" << endl;

	W_mean = lookup_MC_cut(kine, sbsfieldscale, run_target, "W");
	W_sigma = lookup_MC_cut(kine, sbsfieldscale, run_target, "W_sigma");
	W2_mean = lookup_MC_cut(kine, sbsfieldscale, run_target, "W2");
	W2_sigma = lookup_MC_cut(kine, sbsfieldscale, run_target, "W2_sigma");

	//Set defaults???
	if( W_mean == -1 ){ W_mean = Mp; }
	if( W_sigma == -1 ){ W_sigma = 0.15; }
	if( W2_mean == -1){ W2_mean = pow(Mp, 2); }
	if( W2_sigma == -1){ W2_sigma = 0.25; }

	// HCal_height = -0.2897; // Height of HCal above beamline     ---- ORIGINAL ----
	if( kine == 4 ){
		HCal_height = 0.0; // Modified to calibrate peak centers for MC	
	}
	if( kine == 8 ){
		HCal_height = 0.0; // Modified to calibrate peak centers for MC	
	}


	// rootfile_dir = Form("/lustre19/expphy/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass%i/SBS%i/%s/rootfiles", pass, kine, run_target.Data());
	// rootfile_dir = "/volatile/halla/sbs/jboyd/simulation/Rootfiles/";
	rootfile_dir = "/lustre19/expphy/volatile/halla/sbs/jboyd/simulation/out_dir/MC_REPLAY_OUT_DIR";
	// rootfile_dir = "/lustre19/expphy/volatile/halla/sbs/pdbforce/g4sbs_output/sdr/sbs4-sbs50p/simc";

	// rootfile_dir = "/lustre19/expphy/volatile/halla/sbs/adr/sim_results/MC_replay_out";
	I_beam = I_beam_uA*(1.0e-6);
	I_beam_str = Form("%0.2f", I_beam_uA);
	I_beam_str.ReplaceAll(".", "");
	outfile = new TFile(Form("rootfiles/MC_SBS%i_%s_mag%i_%suA_dxdy.root", kine, run_target.Data(), sbsfieldscale, I_beam_str.Data()), "RECREATE");	


//INITIALIZE HISTOGRAMS
	cout << "Initiliazing histograms...";

	h_Ep = new TH1D("h_Ep", Form("E/p - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 200, 0, 2);
	h_PS = new TH1D("h_PS", Form("Pre-Shower Clus. E - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 300, 0, 3);
	h_HCal_e = new TH1D("h_HCal_e", Form("HCal Clus. E - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 200, 0, 0.4);
	h_SHPS = new TH1D("h_SHPS", Form("SH + PS Clus. E - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 500, 0, 5);

	h_E_eloss = new TH1D("E_eloss", Form("Scattered Electron Energy Loss in Target - SBS%i = %i%%, %s", kine, sbsfieldscale, run_target.Data()), 500, 0.0, (0.1)*E_beam);
	h_E_ecorr_vs_vert = new TH2D("h_E_ecorr_vs_vert", Form("Corrected Beam Energy vs Vertex - SBS%i = %i%%, %s; E_{e} (GeV); Z_{vertex} (m)", kine, sbsfieldscale, run_target.Data()), 250, -0.125, 0.125, 500, 0, 0.001);
	h_Q2 = new TH1D("h_Q2", Form("Momentum Transfer Q^2 - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 750, 0.5, 9.0);
	h_E_ep = new TH1D("h_E_ep", Form("Scattered Electron Energy - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 500, 0.0, 1.5*E_beam);
	h_E_pp = new TH1D("h_E_pp", Form("Scattered Proton Energy - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 500, 0.0, 1.5*E_beam);
	h_KE_p = new TH1D("h_KE_p", Form("Scattered Proton Kinetic Energy - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 500, 0.0, 1.5*E_beam);

	h_W = new TH1D("h_W", Form("Invariant Mass W - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 300, 0.0, 3.0);
	h_W2 = new TH1D("h_W2", Form("Invariant Mass W^2 - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 300, 0.0, 3.0);
	h_W_cut = new TH1D("h_W_cut", Form("Invariant Mass W (Coin & Vert Cuts) - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 300, 0.0, 3.0);
	h_W_fcut = new TH1D("h_W_fcut", Form("Invariant Mass W (Fiduc. Cuts) - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 300, 0.0, 3.0);

	h_Wrecon = new TH1D("h_Wrecon", Form("Invariant Mass W recon - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 300, 0.0, 3.0);
	h_W2recon = new TH1D("h_W2recon", Form("Invariant Mass W^2 recon - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 300, 0.0, 3.0);

	h_dx = new TH1D("h_dx",Form("dx (NO CUTS) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 300, -3.5, 2.5);
	h_dx_cut = new TH1D("h_dx_cut",Form("dx (Basic CUTS) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 300, -3.5, 2.5);
	h_dx_wcut = new TH1D("h_dx_wcut",Form("dx (W cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 300, -3.5, 2.5);
	h_dx_fcut = new TH1D("h_dx_fcut",Form("dx (f cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 300, -3.5, 2.5);
	h_dx_wcut_fcut = new TH1D("h_dx_wcut_fcut",Form("dx (W & f cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 300, -3.5, 2.5);
	h_dy = new TH1D("h_dy",Form("dy (NO CUTS) - SBS%i = %i%%, %s; y_{HCal} - y_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);
	h_dy_cut = new TH1D("h_dy_cut",Form("dy (Basic Cuts) - SBS%i = %i%%, %s; y_{HCal} - y_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 250, -1.25, 1.25);  
	h_dy_wcut = new TH1D("h_dy_wcut",Form("dy (W Cuts) - SBS%i = %i%%, %s; y_{HCal} - y_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 250, -1.25, 1.25);  

	h_dxdy = new TH2D("h_dxdy", Form("Hadron Spot(s) on HCal (NO CUTS) - SBS%i = %i%%, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 300, -1.5, 1.5, 500, -2.5, 2.5 );
	h_dxdy_wcut = new TH2D("h_dxdy_wcut", Form("Hadron Spot(s) on HCal (W cut) - SBS%i = %i%%, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 300, -1.5, 1.5, 500, -2.5, 2.5 );
	h_dxdy_cut = new TH2D("h_dxdy_cut", Form("Hadron Spot(s) on HCal (Basic cuts) - SBS%i = %i%%, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 300, -1.5, 1.5, 500, -2.5, 2.5 );
	h_dxdy_ncut = new TH2D("h_dxdy_ncut", Form("Hadron Spot(s) on HCal (n cut) - SBS%i = %i%%, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 250, -1.25, 1.25, 250, -2.5, 2.5 );
	h_dxdy_pcut = new TH2D("h_dxdy_pcut", Form("Hadron Spot(s) on HCal (p cut) - SBS%i = %i%%, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 250, -1.25, 1.25, 250, -2.5, 2.5 );
	h_dxdy_fcut = new TH2D("h_dxdy_fcut", Form("Hadron Spot(s) on HCal (f cut) - SBS%i = %i%%, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 250, -1.25, 1.25, 250, -2.5, 2.5 );

	h_xy = new TH2D("h_xy",Form("HCal Hadron Spots (x, y) (NO CUTS) - SBS%i = %i%%, %s;y_{HCal} (m); x_{HCal} (m)", kine, sbsfieldscale, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);
	h_xy_cut = new TH2D("h_xy_cut", Form("HCal Hadron Spots (x, y) (BASIC CUTS) - SBS%i = %i%%, %s;y_{HCal} (m); x_{HCal} (m)", kine, sbsfieldscale, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);
	h_xy_fcut = new TH2D("h_xy_fcut", Form("HCal Hadron Spots (x, y) (Fiduc. CUTS) - SBS%i = %i%%, %s;y_{HCal} (m); x_{HCal} (m)", kine, sbsfieldscale, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);
	h_xy_cut_p = new TH2D("h_xy_cut_p", Form("HCal Hadron Spots (x, y) (p CUT) - SBS%i = %i%%, %s;y_{HCal} (m); x_{HCal} (m)", kine, sbsfieldscale, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);
	h_xy_cut_n = new TH2D("h_xy_cut_n", Form("HCal Hadron Spots (x, y) (n CUT) - SBS%i = %i%%, %s;y_{HCal} (m); x_{HCal} (m)", kine, sbsfieldscale, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);

	h_PAngleCorr_theta = new TH2D( "h_PAngCorr_theta",Form("BB theta vs HCal theta - SBS%i = %i%%, %s", kine, sbsfieldscale, run_target.Data()), 200, 0.55, 0.75, 300, 0.35, 0.65 );
	h_PAngleCorr_phi = new TH2D( "h_PAngCorr_phi",Form("BB phi vs HCal phi - SBS%i = %i%%, %s", kine, sbsfieldscale, run_target.Data()), 500, -0.4, 0.1, 500, 2.7, 3.2 );
	h_vz_cut = new TH1D("h_vz_cut",Form("BB phi vs HCal phi - SBS%i = %i%%, %s; vertex z (m);", kine, sbsfieldscale, run_target.Data()), 250,-0.125,0.125);	

	cout << "finished. " << endl << endl;

	cout << "--------------------------------------" << endl;
	cout << "Adding files to TChain from: " << rootfile_dir.Data() << endl;

	if( kine == 4 ){
		infile = Form("%s/replayed_jb_gmn_SBS%i_%s_mag%i_175uA_elas_100k_job*", rootfile_dir.Data(), kine, run_target.Data(), sbsfieldscale);
		cout << "Adding file: " << endl;
		cout << infile.Data() << endl << endl;
		TC->Add(infile.Data());
		// TC->Add(Form("%s/replayed_adr_gmn_SBS4_LD2_mag30_elas_1M_job0_1.root", rootfile_dir.Data()));	
	}
	if( kine == 8 ){
		// infile = Form("%s/replayed_jb_gmn_SBS%i_%s_mag%i_%suA_elas_100k_job*", rootfile_dir.Data(), kine, run_target.Data(), sbsfieldscale, I_beam_str.Data());	
		infile = Form("%s/replayed_jb_gmn_SBS8_LD2_mag70mod3_500uA_elas_100k_27_07_2023_job*", rootfile_dir.Data());
		cout << "Adding file: " << endl;
		cout << infile.Data() << endl << endl;
		TC->Add(infile.Data());
		// TC->Add(Form("%s/replayed_adr_gmn_SBS4_LD2_mag30_elas_1M_job0_1.root", rootfile_dir.Data()));	
	}
	if( kine == 9 ){
		infile = Form("%s/replayed_jb_gmn_SBS%i_%s_mag%i_%suA_elas_100k_job*", rootfile_dir.Data(), kine, run_target.Data(), sbsfieldscale, I_beam_str.Data());
		cout << "Adding file: " << endl;
		cout << infile.Data() << endl << endl;
		TC->Add(infile.Data());
		// TC->Add(Form("%s/replayed_adr_gmn_SBS4_LD2_mag30_elas_1M_job0_1.root", rootfile_dir.Data()));	
	}
	// else{
	// 	infile = Form("%s/replayed_jb_gmn_SBS%i_%s_mag%i_%suA_elas_100k_job*", rootfile_dir.Data(), kine, run_target.Data(), sbsfieldscale, I_beam_str.Data());
	// 	cout << "Adding file: " << endl;
	// 	cout << infile.Data() << endl << endl;
	// 	TC->Add(infile.Data());		
	// }

	// TC->Add(Form("%s/replayed_simc_sbs4_sbs50p_89T_elas_job*", rootfile_dir.Data()));
	// TC->Add(Form("%s/replayed_gmn_sbs4_ld2_30p_job_*", rootfile_dir.Data()));
	// TC->Add(Form("%s/replayed_jb_gmn_SBS8_Mag70_500k.root", rootfile_dir.Data()));

	cout << "--------------------------------------" << endl;
	cout << "Setting up branches... ";
	TC->SetBranchStatus( "*", 0 );

	//MC values for normalization
	TC->SetBranchStatus( "MC.mc_omega", 1 );
	TC->SetBranchStatus( "MC.mc_sigma", 1 );
	TC->SetBranchStatus( "MC.mc_nucl", 1);

	// HCal
	TC->SetBranchStatus( "sbs.hcal.x", 1 );
	TC->SetBranchStatus( "sbs.hcal.y", 1 );
	TC->SetBranchStatus( "sbs.hcal.e", 1 );
	TC->SetBranchStatus( "sbs.hcal.nclus", 1);

	// BB track
	TC->SetBranchStatus( "bb.tr.chi2", 1 );
	TC->SetBranchStatus( "bb.tr.n", 1 );
	TC->SetBranchStatus( "bb.tr.px", 1 );
	TC->SetBranchStatus( "bb.tr.py", 1 );
	TC->SetBranchStatus( "bb.tr.pz", 1 );    
	TC->SetBranchStatus( "bb.tr.p", 1 );
	TC->SetBranchStatus( "bb.tr.vx", 1 );
	TC->SetBranchStatus( "bb.tr.vy", 1 );
	TC->SetBranchStatus( "bb.tr.vz", 1 );
	TC->SetBranchStatus( "bb.tr.r_x", 1 );
	TC->SetBranchStatus( "bb.tr.r_y", 1 );
	TC->SetBranchStatus( "bb.tr.r_th", 1 );
	TC->SetBranchStatus( "bb.tr.r_ph", 1 );
	TC->SetBranchStatus( "bb.tr.tg_x", 1 );
	TC->SetBranchStatus( "bb.tr.tg_y", 1 );
	TC->SetBranchStatus( "bb.tr.tg_th", 1 );
	TC->SetBranchStatus( "bb.tr.tg_ph", 1 );
	TC->SetBranchStatus( "bb.gem.track.nhits", 1);

	// BBCal shower preshower
	TC->SetBranchStatus( "bb.ps.e", 1 );
	TC->SetBranchStatus( "bb.ps.x", 1 );
	TC->SetBranchStatus( "bb.ps.y", 1 );
	TC->SetBranchStatus( "bb.sh.e", 1 );
	TC->SetBranchStatus( "bb.sh.x", 1 );
	TC->SetBranchStatus( "bb.sh.y", 1 );
	TC->SetBranchStatus( "bb.sh.nclus", 1 );
	TC->SetBranchStatus( "bb.ps.nclus", 1 );

// Set BRANCH ADDRESSES
	//MC normalization variables
	TC->SetBranchAddress( "MC.mc_omega", &mc_omega );
	TC->SetBranchAddress( "MC.mc_sigma", &mc_sigma );
	TC->SetBranchAddress( "MC.mc_fnucl", &mc_fnucl );

	// HCal
	TC->SetBranchAddress( "sbs.hcal.x", &hcal_x );
	TC->SetBranchAddress( "sbs.hcal.y", &hcal_y );
	TC->SetBranchAddress( "sbs.hcal.e", &hcal_e );
	TC->SetBranchAddress( "sbs.hcal.nclus", &nclus );

	// BB track
	TC->SetBranchAddress( "bb.tr.chi2", bb_tr_chi2 );
	TC->SetBranchAddress( "bb.tr.n", &bb_tr_n );
	TC->SetBranchAddress( "bb.tr.px", bb_tr_px );
	TC->SetBranchAddress( "bb.tr.py", bb_tr_py );
	TC->SetBranchAddress( "bb.tr.pz", bb_tr_pz );
	TC->SetBranchAddress( "bb.tr.p", bb_tr_p );
	TC->SetBranchAddress( "bb.tr.vx", bb_tr_vx );
	TC->SetBranchAddress( "bb.tr.vy", bb_tr_vy );
	TC->SetBranchAddress( "bb.tr.vz", bb_tr_vz );
	TC->SetBranchAddress( "bb.tr.r_x", bb_fp_x );
	TC->SetBranchAddress( "bb.tr.r_y", bb_fp_y );
	TC->SetBranchAddress( "bb.tr.r_th", bb_fp_th );
	TC->SetBranchAddress( "bb.tr.r_ph", bb_fp_ph );
	TC->SetBranchAddress( "bb.tr.tg_x", bb_tgt_x );
	TC->SetBranchAddress( "bb.tr.tg_y", bb_tgt_y );
	TC->SetBranchAddress( "bb.tr.tg_th", bb_tgt_th );
	TC->SetBranchAddress( "bb.tr.tg_ph", bb_tgt_ph );

	// BBCal shower preshower
	TC->SetBranchAddress( "bb.ps.e", &bb_ps_e );
	TC->SetBranchAddress( "bb.ps.x", &bb_ps_x );
	TC->SetBranchAddress( "bb.ps.y", &bb_ps_y );
	TC->SetBranchAddress( "bb.sh.e", &bb_sh_e );
	TC->SetBranchAddress( "bb.sh.x", &bb_sh_x );
	TC->SetBranchAddress( "bb.sh.y", &bb_sh_y );
	TC->SetBranchAddress( "bb.sh.nclus", &SH_nclus );
	TC->SetBranchAddress( "bb.ps.nclus", &PS_nclus );

	cout << "finished with branches. " << endl << endl;

	cout << "-------- Starting scan on TChain --------" << endl;

	master_cut_vec = {
		"sbs.hcal.nclus>0",
		"bb.ps.nclus>0",
		"bb.sh.nclus>0",
		"abs(bb.tr.vz[0])<=0.075",
		"bb.gem.track.nhits[0]>3",
		"bb.tr.n==1",
	};

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
	cout << "Applying Master Cut: " << endl;
	cout << master_cut << endl;

	TEventList *ev_list = new TEventList("ev_list", "Elastic Events List");
	TC->Draw(">>ev_list", master_cut);
	Nevents = ev_list->GetN();

	cout << "--------------------------------------" << endl;
	cout << "Number of events to analyze: " << Nevents << endl;
	cout << "--------------------------------------" << endl;
	cout << "--------------------------------------" << endl;
	cout << "Starting analysis loop on events..... " << endl;

	elastic_yield = 0;

	int watch_cnt = 0;	
	int five_percent = int(0.05*Nevents);
	vector<double> time_for_five;
	double average_time = 0.0, time_remaining;
	StopWatch->Start();

	for(Long64_t nevent = 0; nevent < Nevents; nevent++){
		TC->GetEntry( ev_list->GetEntry( nevent ));

		if( nevent%five_percent == 0){		
			StopWatch->Stop();

			if( watch_cnt == 0){
				cout << "Evt: " << nevent <<"/" << Nevents << Form("(%.0f/100%%)", 100.0*double(1.0*nevent/Nevents)) << ". Elastic yield = " << elastic_yield << ". "<< endl;
			}

			if( watch_cnt > 0 ){
				time_for_five.push_back(StopWatch->RealTime());	
				average_time = VectorMean(time_for_five);
				cout << "average time for 5 = " << average_time << endl;
				time_remaining = average_time*( 1.0 - double(nevent)/double(Nevents));
				cout << "Evt: " << nevent <<"/" << Nevents << Form("(%.0f/100%%)", 100.0*double(1.0*nevent/Nevents)) << ". Elastic yield = " << elastic_yield << ". Time left: " << time_remaining << endl;
			}
			StopWatch->Reset();
			StopWatch->Continue();
		}

      	if( !correct_beam_energy){
      		Eloss = (bb_tr_vz[0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
      		h_E_eloss->Fill( Eloss );

      		E_beam_final = E_beam;
      		h_E_ecorr_vs_vert->Fill( bb_tr_vz[0], E_beam_final);
      	}

		if( correct_beam_energy ){
      		Eloss = (bb_tr_vz[0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
      		h_E_eloss->Fill( Eloss );

      		E_beam_final = E_beam - Eloss;
      		h_E_ecorr_vs_vert->Fill( bb_tr_vz[0], E_beam_final);    	
      	}

      	p_corr = bb_tr_p[0] - E_loss_outgoing; //Neglecting mass of e'

      	//Proceed only if at least one track exists in BB arm - lowest chi2 track always first element
      	if( false ){ //Null-ing to check
	      	if( bb_tr_n > 1){
	      		continue;
	      	}
	    }

//ENERGY CALCULATIONS -- CALIBRATIONS
	    Double_t Ep = (bb_ps_e + bb_sh_e)/(bb_tr_p[0]);
		Double_t PS = bb_ps_e;
		Double_t SHPS = bb_ps_e + bb_sh_e;
		Double_t HCal_e = hcal_e;

		// temporary fix is to set luminosity to 1
		// luminosity = 1.0;
		luminosity = calc_luminosity(I_beam, run_target );
		if( calc_pn_weight ){
			pn_weight = ((mc_omega*mc_sigma)*luminosity)/ngen_total;			
		}
		else{ pn_weight = 1.0; }


		h_Ep->Fill(Ep);
		h_PS->Fill(PS);
		h_SHPS->Fill(SHPS);
		h_HCal_e->Fill(hcal_e);

// ---------------------------------------------------
//  KINEMATIC VARIABLE CALCULATIONS
// ---------------------------------------------------

       	TVector3 vertex( 0, 0, bb_tr_vz[0] ); // z location of vertex in hall coordinates
		TLorentzVector P_beam( 0, 0, E_beam_final, E_beam_final ); //Mass of e negligable
		TLorentzVector k_prime( bb_tr_px[0], bb_tr_py[0], bb_tr_pz[0], bb_tr_p[0] );
		TLorentzVector P_targ( 0, 0, 0, Mp );     	

		TLorentzVector q = P_beam - k_prime; //Standard q-vector
		TVector3 qunit = q.Vect().Unit(); //q-vector direction

		//Define HCal coordinate system
		TVector3 HCal_zaxis( sin(-HCal_theta_rad ), 0, cos(-HCal_theta_rad) );
		TVector3 HCal_xaxis( 0, -1, 0 );
		TVector3 HCal_yaxis = HCal_zaxis.Cross(HCal_xaxis).Unit();

		TVector3 HCal_origin = HCal_dist*HCal_zaxis + (HCal_height)*HCal_xaxis;

		//Define intersection points for hadron vector
		Double_t sintersect = (HCal_origin-vertex).Dot( HCal_zaxis )/qunit.Dot( HCal_zaxis ); //Scintillator face
		TVector3 HCal_intersect = vertex + sintersect*qunit; //HCal Face


//---------"dxdy" method

		x_expected_HCal = (HCal_intersect - HCal_origin).Dot( HCal_xaxis );
		y_expected_HCal = (HCal_intersect - HCal_origin).Dot( HCal_yaxis );

		E_ep = sqrt( pow(Me,2) + pow(bb_tr_p[0],2) ); // Obtain the scattered electron energy
		h_E_ep->Fill( E_ep );

		p_ep = bb_tr_p[0];

		Q2 = 2*E_beam_final*E_ep*( 1-(bb_tr_pz[0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta
		h_Q2->Fill( Q2 );

		//Get invariant mass transfer W from the four-momentum of the scattered nucleon
		TLorentzVector P_gammaN = P_targ + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)
		W = P_gammaN.M(); 
		h_W->Fill( W );
		h_W2->Fill(pow(W, 2));

		//Use the electron kinematics to predict the proton momedntum assuming elastic scattering on free proton at rest (will need to correct for fermi motion):
		E_pp = nu + Mp; // Get energy of the proton
		E_nucleon = sqrt(pow(pp,2)+pow(Mp,2)); // Check on E_pp, same
		h_E_pp->Fill( E_pp ); // Fill histogram

		KE_p = nu; // For elastics
		h_KE_p->Fill( KE_p );

		dx = hcal_x - x_expected_HCal;
		dy = hcal_y - y_expected_HCal;

	//Resolve the hadron spots without cuts
		h_dx->Fill( dx, pn_weight);
		h_dy->Fill( dy, pn_weight);
		h_dxdy->Fill( dy, dx, pn_weight );
		h_xy->Fill( hcal_y, hcal_x );

	// Preliminary HCal projections with single cut on W
		if( fabs(W - W_mean) < W_sigma ){
			h_dx_wcut->Fill( dx, pn_weight );
			h_dy_wcut->Fill ( dy, pn_weight );
			h_dxdy_wcut->Fill( dy, dx, pn_weight );
		}

	//Populate BB/HCal correlation histograms from elastics
		h_PAngleCorr_phi->Fill( e_prime_phi, nucleon_phi );
		h_PAngleCorr_theta->Fill( e_prime_theta, nucleon_theta );

	//Fill vertex position histogram for cut on tracks
    	h_vz_cut->Fill( bb_tr_vz[0] );

	//FIDUCIAL Cut
		//Check "elastic" events on center HCal for id with spot checks
		bool HCal_on = false;
		bool is_p = false;
		bool is_n = false;

		if( fiducial_cut ){
			if( kine == 4 && sbsfieldscale == 30 ){
				dx_p_scale = 0.6;
				dx_n_scale = 0.6;
			}
			if( kine == 8 && sbsfieldscale == 70 ){
				dx_p_scale = 1.0;
				dx_n_scale = 1.0;
			}

			dx_p = lookup_MC_dxdy(kine, sbsfieldscale, run_target, "dx_p");
			dx_p_sigma = dx_p_scale*lookup_MC_dxdy(kine, sbsfieldscale, run_target, "dx_p_sigma");
			dy_p = lookup_MC_dxdy(kine, sbsfieldscale, run_target, "dy");
			dy_p_sigma = lookup_MC_dxdy(kine, sbsfieldscale, run_target, "dy_sigma");
			dx_n = lookup_MC_dxdy(kine, sbsfieldscale, run_target, "dx_n");
			dx_n_sigma = dx_n_scale*lookup_MC_dxdy(kine, sbsfieldscale, run_target, "dx_n_sigma");
			dy_n = lookup_MC_dxdy(kine, sbsfieldscale, run_target, "dy");
			dy_n_sigma = lookup_MC_dxdy(kine, sbsfieldscale, run_target, "dy_sigma");
			dx_pn_max = abs( dx_p ) + abs( dx_n );
			
		
			if( hcal_y > hcal_y_fmin && hcal_y < hcal_y_fmax && hcal_x >hcal_x_fmin && hcal_x < hcal_x_fmax ){
				HCal_on = true;
			}

			apply_fcut = ((y_expected_HCal - dy_p_sigma) > hcal_y_fmin) && ((y_expected_HCal + dy_p_sigma) < hcal_y_fmax) && ((x_expected_HCal - dx_pn_max - dx_p_sigma) > hcal_x_fmin) && ((x_expected_HCal + dx_p_sigma) < hcal_x_fmax);

			// if( pow( (hcal_x - x_expected_HCal - dx_p)/dx_p_sigma, 2) + pow( (hcal_y - y_expected_HCal - dy_p)/dy_p_sigma,2) <= pow(2.5,2) ){
			// 	is_p = true;
			// 	if( !calc_pn_weight ){
			// 		pn_weight = 1.0;					
			// 	}

			// }
			// if( pow( (hcal_x - x_expected_HCal - dx_n)/dx_n_sigma,2) + pow( (hcal_y - y_expected_HCal - dy_n)/dy_n_sigma,2) <= pow(2.5,2) ){
			// 	is_n = true;
			// 	if( !calc_pn_weight ){
			// 		pn_weight = 0.33333333;					
			// 	}
			// }
			if( int(mc_fnucl) == 0 ){ 
				is_n = true; 
				is_p = false;
			}
			if( int(mc_fnucl) == 1 ){ 
				is_p = true; 
				is_n = false;
			}

		//Fill respective histograms for these checks.
			if( HCal_on && is_n && apply_fcut ) h_dxdy_ncut->Fill( dy, dx );
			if( HCal_on && is_p && apply_fcut ) h_dxdy_pcut->Fill( dy, dx );

		//----------neutron
			if( HCal_on && is_n && apply_fcut ){
				if( !calc_pn_weight ){
					pn_weight = (1.0/3.0);					
				}

				if( (hcal_x - dx_pn_max )>hcal_x_fmin ){
					h_dxdy_fcut->Fill( dy, dx, pn_weight );
					h_dx_fcut->Fill( dx, pn_weight );
					h_W_fcut->Fill( W, pn_weight );
					h_xy_fcut->Fill( hcal_y, hcal_x, pn_weight );
					h_xy_cut_n->Fill( hcal_y, hcal_x );
					if( fabs(W - W_mean) < W_sigma ){
						h_dx_wcut_fcut->Fill(dx, pn_weight);
					}
					elastic_yield++;
				}
			}		
		//----------proton
			else if( HCal_on && is_p && apply_fcut ){
				if( !calc_pn_weight ){
					pn_weight = 1.0;					
				}

				if( (hcal_x + dx_pn_max)<hcal_x_fmax ){
					h_dxdy_fcut->Fill( dy, dx, pn_weight );
					h_dx_fcut->Fill( dx, pn_weight );
					h_W_fcut->Fill( W, pn_weight );
					h_xy_fcut->Fill( hcal_y, hcal_x, pn_weight );
					h_xy_cut_p->Fill( hcal_y, hcal_x );
					if( fabs(W - W_mean) < W_sigma ){
						h_dx_wcut_fcut->Fill(dx, pn_weight);
					}
					elastic_yield++;
				}
			}
		//END OF FIDUCIAL CUT
		}

		//Still should count elastic yields if we got this far.....
		if( !fiducial_cut ){
			elastic_yield++;
		}


//------(Calibrate method) RECONSTRUCTED
	    e_prime_theta = acos( bb_tr_pz[0]/bb_tr_p[0] ); //track momenutm to reconstruct e' theta
      	e_prime_phi = atan2( bb_tr_py[0], bb_tr_px[0]);
      	TLorentzVector Ptarg_recon( 0, 0, 0, 0.5*(Mp+Mn) ); //Average of proton and neutron rest mass
		
		//FIDUCIAL STUFF
		//Define the expected position of hadron on HCal from BB track 
		x_recon_expect = HCal_intersect.Dot( HCal_xaxis );
		y_recon_expect = HCal_intersect.Dot( HCal_yaxis );
//-------------
		
		Double_t E_ep = bb_tr_p[0]; // Obtain the scattered electron energy, neglect mass e
		Double_t p_ep = bb_tr_p[0]; // Obtain the magnitude of scattered electron momentum
		Double_t Q2 = 2*E_beam_final*E_ep*( 1-(bb_tr_pz[0]/p_ep) );
		Double_t nu = E_beam_final-E_ep; // Obtain energy transfer
		W2recon = pow( Mp,2 )+2*Mp*nu-Q2; // Obtain W2 from Q2 and nu
		Wrecon = sqrt(W2);		

		h_Wrecon->Fill(Wrecon);
		h_W2recon->Fill(W2recon);

	//end of events loop  
    }


	cout << "---------------------------------------" << endl;
	cout << "-----Finished going through events-----" << endl;
	cout << "---------------------------------------" << endl;

	outfile->Write();
	outfile->Close();

	TFile *infile = new TFile(outfile->GetName(), "UPDATE");
	TH1D *hin_dx_fcut;
	TH2D *hin_dxdy, *hin_dxdy_fcut, *hin_dxdy_wcut;


	hin_dx_fcut = static_cast<TH1D*>(infile->Get("h_dx_fcut"));
	hin_dxdy = static_cast<TH2D*>(infile->Get("h_dxdy"));
	hin_dxdy_fcut = static_cast<TH2D*>(infile->Get("h_dxdy_fcut"));
	hin_dxdy_wcut = static_cast<TH2D*>(infile->Get("h_dxdy_wcut"));

    TCanvas *c_dxdy = new TCanvas("c_dxdy", "c_dxdy", 600, 500);
    hin_dxdy->Draw("colz");

    TCanvas *c_dxdy_fcut = new TCanvas("c_dxdy_fcut", "c_dxdy_fcut", 600, 500);
    hin_dxdy_fcut->Draw("colz");

    TCanvas *c_dxdy_wcut = new TCanvas("c_dxdy_wcut", "c_dxdy_wcut", 600, 500);
    hin_dxdy_wcut->Draw("colz");

    TCanvas *c_dx_fcut = new TCanvas("c_dx_fcut", "c_dx_fcut", 600, 500);
    hin_dx_fcut->Draw("hist");

	if( calibrate ){
    	TFile *infile = new TFile(outfile->GetName(), "READ");

    	fits_outfile_name = Form("%s", outfile->GetName());
		fits_outfile_name.ReplaceAll(".root", "_fits.root");
		// fits_outfile = new TFile(Form("%s", fits_outfile_name.Data()), "RECREATE");

    	hin_dx_wcut = static_cast<TH1D*>(infile->Get("h_dx_wcut"));
    	hin_dy_wcut = static_cast<TH1D*>(infile->Get("h_dy_wcut"));

    	c_dx = new TCanvas("c_dx", "c_dx", 600, 500);
    	hin_dx_wcut->GetYaxis()->SetRangeUser(0, 1.01*hin_dx_wcut->GetMaximum());
    	hin_dx_wcut->Draw();
  	
  		if( !fit_full_function ){
  			fit_peaks_only = true;
  		}

  		if( fit_peaks_only ){



	  	//------ p -------
	    	TF1 *fit_dx_p = new TF1("fit_dx_p", fit_gaus, -1.5, -0.2, 3);
	  
	    	fit_dx_p->SetParName(0, "dx_p Norm");
			fit_dx_p->SetParName(1, "dx_p Center");
			fit_dx_p->SetParName(2, "dx_p Sigma");
			fit_dx_p->SetLineColor(2);

			if( kine == 4 && sbsfieldscale == 30){
				fit_dx_p->SetParLimits(0, 0, hin_dx_wcut->GetMaximum());
				fit_dx_p->SetParLimits(1, -0.7, -0.45);
				fit_dx_p->SetParLimits(2, 0.1, 0.178);
			}

			if( kine == 4 && sbsfieldscale == 50){
				fit_dx_p->SetParLimits(0, 0, hin_dx_wcut->GetMaximum());
				fit_dx_p->SetParLimits(1, -0.9, -0.8);
				fit_dx_p->SetParLimits(2, 0.1, 0.19);
			}
			//for custom field setting to get to correct 70%:
			if( kine == 8 && sbsfieldscale == 70){
				fit_dx_p->SetParLimits(0, 0, hin_dx_wcut->GetMaximum());
				fit_dx_p->SetParLimits(1, -1.05, -0.95);
				fit_dx_p->SetParLimits(2, 0.22, 0.24);
			}

			// if( kine == 8 && sbsfieldscale == 70){
			// 	fit_dx_p->SetParLimits(0, 0, hin_dx_wcut->GetMaximum());
			// 	fit_dx_p->SetParLimits(1, -0.9, -0.7);
			// 	fit_dx_p->SetParLimits(2, 0.1, 0.16);			
			// }
		
		
			hin_dx_wcut->Fit("fit_dx_p", "R+");
			dx_p = fit_dx_p->GetParameter(1);
			dx_p_sigma = fit_dx_p->GetParameter(2);	

		//------ n -------
	    	TF1 *fit_dx_n = new TF1("fit_dx_n", fit_gaus, -0.6, 0.5, 3);
	  
	    	fit_dx_n->SetParName(0, "dx_n Norm");
			fit_dx_n->SetParName(1, "dx_n Center");
			fit_dx_n->SetParName(2, "dx_n Sigma");
			fit_dx_n->SetLineColor(3);

			if( kine == 4 && sbsfieldscale == 30){
				fit_dx_n->SetParLimits(0, 0, (0.45)*hin_dx_wcut->GetMaximum());
				fit_dx_n->SetParLimits(1, -0.05, 0.05);
				fit_dx_n->SetParLimits(2, 0.1, 0.16);
			}	

			if( kine == 4 && sbsfieldscale == 50){
				fit_dx_n->SetParLimits(0, 0, (0.45)*hin_dx_wcut->GetMaximum());
				fit_dx_n->SetParLimits(1, 0.0, 0.13);
				fit_dx_n->SetParLimits(2, 0.1, 0.18);
			}	

			if( kine == 8 && sbsfieldscale == 70){
				fit_dx_n->SetParLimits(0, 0, (0.35)*hin_dx_wcut->GetMaximum());
				fit_dx_n->SetParLimits(1, -0.05, 0.05);
				fit_dx_n->SetParLimits(2, 0.1, 0.16);
			}	
		
			hin_dx_wcut->Fit("fit_dx_n", "R+");
			dx_n = fit_dx_n->GetParameter(1);
			dx_n_sigma = fit_dx_n->GetParameter(2);	
		}

		if( fit_full_function ){

			cout << "Plotting full dx functions... " << endl;

			fit_dx_full = new TF1("fit_dx_full", fullFitFunction, -2.0, 1.0, 11);

			fit_dx_full->SetParName(0, "dx - p Norm");
			fit_dx_full->SetParName(1, "dx - p Center");
			fit_dx_full->SetParName(2, "dx - p Sigma");
			fit_dx_full->SetParName(3, "dx - n Norm");
			fit_dx_full->SetParName(4, "dx - n Center");
			fit_dx_full->SetParName(5, "dx - n Sigma");
			fit_dx_full->SetParName(6, "dx - BG a");
			fit_dx_full->SetParName(7, "dx - BG b");
			fit_dx_full->SetParName(8, "dx - BG c");
			fit_dx_full->SetParName(9, "dx - BG d");	
			fit_dx_full->SetParName(10, "dx - BG e");

			if( kine == 8 && sbsfieldscale == 70 ){
				fit_dx_full->SetParLimits(0, 0.9*hin_dx_wcut->GetMaximum(), hin_dx_wcut->GetMaximum());
				fit_dx_full->SetParLimits(1, -1.05, -0.95);
				fit_dx_full->SetParLimits(2, 0.22, 0.24);
				fit_dx_full->SetParLimits(3, (0.3)*hin_dx_wcut->GetMaximum(), (0.33)*hin_dx_wcut->GetMaximum());
				fit_dx_full->SetParLimits(4, -0.05, 0.05);
				fit_dx_full->SetParLimits(5, 0.1, 0.16);

			}

			if( kine == 4 ){
				fit_dx_full->SetParLimits(0, 0.9*hin_dx_wcut->GetMaximum(), hin_dx_wcut->GetMaximum());
				fit_dx_full->SetParLimits(1, -0.95, -0.85);
				fit_dx_full->SetParLimits(2, 0.1, 0.16);
				fit_dx_full->SetParLimits(3, (0.3)*hin_dx_wcut->GetMaximum(), (0.33)*hin_dx_wcut->GetMaximum());
				fit_dx_full->SetParLimits(4, -0.060, -0.055);
				fit_dx_full->SetParLimits(5, 0.1, 0.16);				
			}

			// fit_dx_full->SetParLimits(6, 30, 100);
			// fit_dx_full->SetParLimits(7, -80, -30);
			// fit_dx_full->SetParLimits(8, -50, -20);
			// fit_dx_full->SetParLimits(9, 5, 20);
			// fit_dx_full->SetParLimits(10, 0, 20);
			// fit_dx_full->SetParLimits(6, 1.4e-08, 3.5e-08);
			// fit_dx_full->SetParLimits(7, -2.32e-06, 2.32e-08);
			// fit_dx_full->SetParLimits(8, -1.23e-06, -1.23e-10);
			// fit_dx_full->SetParLimits(9, 1.54e-08, 1.54e-06);
			// fit_dx_full->SetParLimits(10, 6.985e-09, 6.985e-07);	
			//so-so
			// fit_dx_full->SetParLimits(6, 0.00000001, 0.000000059);
			// fit_dx_full->SetParLimits(7, -0.000000055, -0.000000015);
			// fit_dx_full->SetParLimits(8, -0.000000030, -0.000000010);
			// fit_dx_full->SetParLimits(9, 0.0000000060, 0.000000040);
			// fit_dx_full->SetParLimits(10, 0.0000000020, 0.00000001);

			// fit_dx_full->SetParLimits(6, 0.00000001, 0.000000059);
			// fit_dx_full->SetParLimits(7, -0.000000105, -0.000000015);
			// fit_dx_full->SetParLimits(8, -0.000000040, -0.000000010);
			// fit_dx_full->SetParLimits(9, 0.0000000060, 0.000000070);
			// fit_dx_full->SetParLimits(10, 0.0000000020, 0.00000002);

			// fit_dx_full->SetParLimits(6, 2.0e-04, 3.0e-04);
			// fit_dx_full->SetParLimits(7, -4.0e-04, -3.0e-04);
			// fit_dx_full->SetParLimits(8, -3.0e-04, -2.0e-04);
			// fit_dx_full->SetParLimits(9, -10.0e-04, -9.0e-04);
			// fit_dx_full->SetParLimits(10, 1.0e-05, 1.9e-05);

			// a: 74.5009
			// b: -34.5542
			// c: -21.1193
			// d: 6.13723
			// e: 2.5271
			fit_dx_full->SetParLimits(6, 74.25, 74.75);
			fit_dx_full->SetParLimits(7, -34.75, -34.25);
			fit_dx_full->SetParLimits(8, -21.2, -21.0);
			fit_dx_full->SetParLimits(9, 6.0, 6.2);
			fit_dx_full->SetParLimits(10, 2.4, 2.6);

	
   // 1  p0           2.40000e-04     fixed    
   // 2  p1          -3.30000e-04     fixed    
   // 3  p2          -2.64080e-04     fixed    
   // 4  p3          -9.40000e-06     fixed    
   // 5  p4           1.14000e-05     fixed  			

			double bg_par6 = 74.5009;
			double bg_scale = (bg_par6)/3.06383e-03;

			// fit_dx_full->FixParameter(6, bg_par6);
			// fit_dx_full->SetParameter(7, bg_scale*(-3.98150e-03));
			// fit_dx_full->SetParameter(8, bg_scale*(-1.81918e-03));
			// fit_dx_full->SetParameter(9, bg_scale*(2.32083e-03));
			// fit_dx_full->SetParameter(10, bg_scale*(9.43538e-04));

			// fit_dx_full->FixParameter(6, bg_par6);
			// fit_dx_full->FixParameter(7, bg_scale*(-3.98150e-03));
			// fit_dx_full->FixParameter(8, bg_scale*(-1.81918e-03));
			// fit_dx_full->FixParameter(9, bg_scale*(2.32083e-03));
			// fit_dx_full->FixParameter(10, bg_scale*(9.43538e-04));


			hin_dx_wcut->Fit("fit_dx_full", "R");
			fit_dx_full->GetParameters(&par_full[0]);
			fit_dx_full->SetParameters(par_full);
			// fit_dx_full->Draw

			dx_p = fit_dx_full->GetParameter(1);
			dx_p_sigma = fit_dx_full->GetParameter(2);
			dx_n = fit_dx_full->GetParameter(4);
			dx_n_sigma = fit_dx_full->GetParameter(5);
		}

	//------- dy -------
    	c_dy = new TCanvas("c_dy", "c_dy", 600, 500);
    	hin_dy_wcut->Draw();

    	TF1 *fit_dy = new TF1("fit_dy", fit_gaus, -1.5, 1.5, 3);
  
    	fit_dy->SetParName(0, "dy Norm");
		fit_dy->SetParName(1, "dy Center");
		fit_dy->SetParName(2, "dy Sigma");
		fit_dy->SetLineColor(3);

		if( kine == 8 && sbsfieldscale == 70){
			fit_dy->SetParLimits(0, 0, hin_dy_wcut->GetMaximum());
			fit_dy->SetParLimits(1, -0.15, 0.15);
			fit_dy->SetParLimits(2, 0.1, 0.21);	
		}
		else{
			fit_dy->SetParLimits(0, 0, hin_dy_wcut->GetMaximum());
			fit_dy->SetParLimits(1, -0.15, 0.15);
			fit_dy->SetParLimits(2, 0.1, 0.23);				
		}

	
		hin_dy_wcut->Fit("fit_dy", "R+");
		dy_p = fit_dy->GetParameter(1);
		dy_p_sigma = fit_dy->GetParameter(2);	

//----------------------------------------------

	//BACKGROUND FITTING
	    TString BG_fit_type = "wcut";
	    double dy_sigma = lookup_MC_dxdy(kine, sbsfieldscale, run_target, "dy_sigma");
	    int num_par = 11;

	    //Choose some appropriate "scales" and min/max values:
		double dy_mult = 3.5;
		dy_min = -dy_mult*dy_sigma;
		dy_max = dy_mult*dy_sigma;
		dy_min_bin = int(100*(1.5 + dy_min));
		dy_max_bin = int(100*(1.5 + dy_max));

		c_dy->cd();

		TLine *tl_dy_min = new TLine(dy_min, 0, dy_min, hin_dy_wcut->GetMaximum());
		tl_dy_min->SetLineStyle(8);
		tl_dy_min->Draw("same");

		TLine *tl_dy_max = new TLine(dy_max, 0, dy_max, hin_dy_wcut->GetMaximum());
		tl_dy_max->SetLineStyle(8);
		tl_dy_max->Draw("same");

		cout << "---------------------------" << endl;
		cout << "dy_sigma: " << dy_sigma << endl;
		cout << "dy_min: " << dy_min << endl;
		cout << "dy_min_bin: " << dy_min_bin << endl;
		cout << "dy_max: " << dy_max << endl;
		cout << "dy_max_bin: " << dy_max_bin << endl;
		cout << "---------------------------" << endl;

	    //Make a copy of the dxdy histogram to use for beyond central dxdy region out
		if( BG_fit_type == "wcut"){
			hin_dxdy_Main = (TH2D*)hin_dxdy_wcut->Clone();
			hin_dxdy_midCut = (TH2D*)hin_dxdy_wcut->Clone();
		}

		hin_dxdy_BeyondDyCut = new TH2D("hin_dxdy_BeyondDyCut", "dxdy excluding regions beyond dy limits", hin_dxdy_Main->GetXaxis()->GetNbins(), -1.5, 1.5, hin_dxdy_Main->GetYaxis()->GetNbins(), -2.5, 2.5);

	//Lets cut the central primary region of the p and n peaks --> Left with only background from outer regions
		for( int x = dy_min_bin; x < dy_max_bin; x++ ){
			for( int y = 0; y < hin_dxdy_midCut->GetYaxis()->GetNbins(); y++ ){
				hin_dxdy_midCut->SetBinContent(x, y, 0);
			}
		}
		TCanvas *c_dxdy_midCut_BG = new TCanvas("c_dxdy_midCut_BG", "c_dxdy_midcut_BG", 600, 500);
		h_dx_BG = hin_dxdy_midCut->ProjectionY();
		h_dx_BG->Scale(1/(h_dx_BG->Integral()));
		h_dx_BG->Draw("hist");

	// FIT POLYNOMIAL TO BG
	tf_dx_bg = new TF1("tf_dx_bg", background, -2.0, 1.0, 5);

	tf_dx_bg->SetParName(0, "dxBG dyCut - a");
	tf_dx_bg->SetParName(1, "dxBG dyCut - b");
	tf_dx_bg->SetParName(2, "dxBG dyCut - c");
	tf_dx_bg->SetParName(3, "dxBG dyCut - d");
	tf_dx_bg->SetParName(4, "dxBG dyCut - e");

	if( dy_mult == 3.0 ){
		tf_dx_bg->SetParLimits(0, 27, 33);
		tf_dx_bg->SetParLimits(1, -25, -15);
		tf_dx_bg->SetParLimits(2, -18, -12);
		tf_dx_bg->SetParLimits(3, 3, 9);
		tf_dx_bg->SetParLimits(4, 0, 6);		
	}
	if( dy_mult == 3.5 && kine ==8 ){
		tf_dx_bg->SetParLimits(0, 0.00045, 0.0005);
		tf_dx_bg->SetParLimits(1, -0.0008, 0.0007);
		tf_dx_bg->SetParLimits(2, -0.0003, -0.0002);
		tf_dx_bg->SetParLimits(3, 0.00030, 0.0004);
		tf_dx_bg->SetParLimits(4, 0.00005, 0.00015);
	}

	// if( dy_mult == 3.5 && kine ==8 ){
	// 	tf_dx_bg->SetParLimits(0, 0.0001, 0.0010);
	// 	tf_dx_bg->SetParLimits(1, -0.002, 0.000);
	// 	tf_dx_bg->SetParLimits(2, -0.001, -0.000);
	// 	tf_dx_bg->SetParLimits(3, 0.0, 0.00075);
	// 	tf_dx_bg->SetParLimits(4, 0.0, 0.0002);
	// }

	// if( dy_mult == 3.5 && kine ==8 ){
	// 	tf_dx_bg->SetParLimits(0, 0.001, 0.0050);
	// 	tf_dx_bg->SetParLimits(1, -0.0055, -0.0015);
	// 	tf_dx_bg->SetParLimits(2, -0.0030, -0.0010);
	// 	tf_dx_bg->SetParLimits(3, 0.00060, 0.0040);
	// 	tf_dx_bg->SetParLimits(4, 0.00020, 0.001);
	// }

	h_dx_BG->Fit("tf_dx_bg", "R");
	tf_dx_bg->Draw("same");

	cout << "Fitting BG, p, and n peaks..." << endl;

	cout << "sbsfield scale = " << sbsfieldscale << ". Num params: " << num_par << endl;

	fit_dx_dyBG = new TF1("fit_dx_dyBG", fullFitFunction, -2.5, 1.5, 11);
	
	fit_dx_dyBG->SetParName(0, "dx_dyBG - p Norm");
	fit_dx_dyBG->SetParName(1, "dx_dyBG - p Center");
	fit_dx_dyBG->SetParName(2, "dx_dyBG - p Sigma");
	fit_dx_dyBG->SetParName(3, "dx_dyBG - n Norm");
	fit_dx_dyBG->SetParName(4, "dx_dyBG - n Center");
	fit_dx_dyBG->SetParName(5, "dx_dyBG - n Sigma");
	fit_dx_dyBG->SetParName(6, "dx_dyBG - BG a");
	fit_dx_dyBG->SetParName(7, "dx_dyBG - BG b");
	fit_dx_dyBG->SetParName(8, "dx_dyBG - BG c");
	fit_dx_dyBG->SetParName(9, "dx_dyBG - BG d");	
	fit_dx_dyBG->SetParName(10, "dx_dyBG - BG e");		

   // 1  p0           2.03285e-03
   // 2  p1          -1.34790e-03
   // 3  p2          -1.16418e-03
   // 4  p3           4.56436e-04
   // 5  p4           2.43961e-04

   // 1  p0           2.75772e-03
   // 2  p1          -2.00000e-03
   // 3  p2          -1.63567e-03
   // 4  p3           7.01072e-04
   // 5  p4           3.59674e-04

   // 1  p0           2.03001e-03
   // 2  p1          -1.34422e-03
   // 3  p2          -1.16196e-03
   // 4  p3           4.54717e-04
   // 5  p4           2.43194e-04

   //    1  p0           3.50000e-03   9.57108e-05   2.37160e-03** at limit **
   // 2  p1          -3.00000e-03   1.64378e-05   4.59861e-04** at limit **
   // 3  p2          -1.84727e-03   6.50597e-05   2.80599e-04   5.00377e-03
   // 4  p3           1.00000e-03   1.05445e-05   7.79765e-04** at limit **
   // 5  p4           4.13884e-04   1.47928e-05   1.24276e-05  -2.43631e-02

   // 1  p0           2.02986e-03   5.60479e-04   1.25392e-04   5.33703e-03
   // 2  p1          -1.34414e-03   5.76365e-04   7.22534e-04  -2.57406e-03
   // 3  p2          -1.16185e-03   3.83789e-04  -7.48149e-04   1.54623e-03
   // 4  p3           4.54684e-04   1.83786e-04   2.91212e-05   7.19153e-02
   // 5  p4           2.43169e-04   1.00585e-04   1.00585e-04  -2.85330e-02


		// fits_outfile->Write();
		// fits_outfile->Close();
    //END OF CALIBRATION
    }
	cout << "------------------------------------------------------------------"<< endl;
	cout << "                       ANALYSIS FINISHED" << endl;
	cout << "------------------------------------------------------------------"<< endl;
	cout << "Run parameters: " << endl;
	cout << "Kinematic: SBS" << kine << endl;
	cout << "Target: " << run_target.Data() << endl;
	cout << "Beam Energy: " << E_beam << endl;
	cout << "SBS Field: " << sbsfieldscale << "%" << endl;
	cout << "-----------------------------------" << endl;
	cout << "BB angle: " << BB_theta << " degrees; " << BB_theta_rad << " radians " << endl;
	cout << "SBS angle: " << SBS_theta << " degrees; " << SBS_theta_rad << " radians "  << endl;
	cout << "HCal angle: " << HCal_theta << " degrees; " << HCal_theta_rad << " radians " << endl;
	cout << "HCal distance: " << HCal_dist << " m " << endl;
	cout << "-----------------------------------" << endl << endl;
	cout << "Elastic yield: " << elastic_yield << endl << endl;
	cout << "---------------------------------------" << endl << endl;	
	cout << "------------------------------------------------------------------"<< endl;
	cout << "Cut info:" << endl;
	cout << "W_mean: " << W_mean << "; W_sigma: " << W_sigma << endl;
	cout << "W2_mean: " << W2_mean << "; W2_sigma: " << W2_sigma << endl;
	cout << "------------------------------------------------------------------"<< endl << endl;

	if( calibrate ){
		cout << "dx_p = : " << dx_p << "; dx_p_sigma = " << dx_p_sigma << endl;
		cout << "dx_n = : " << dx_n << "; dx_n_sigma = " << dx_n_sigma << endl;
		cout << "dy_p = " << dy_p << "; dy_p_sigma = " << dy_p_sigma << endl;
	}
	cout << "------------------------------------------------------------------"<< endl;
	cout << "Output file: " << outfile->GetName() << endl;
	cout << "------------------------------------------------------------------"<< endl;
	// cout << "Fits Output file: " << fits_outfile->GetName() << endl << endl;
	cout << "------------------------------------------------------------------"<< endl;
	
	auto total_time_end = high_resolution_clock::now();
	auto total_time_duration = duration_cast<minutes>(total_time_end - total_time_start);
	cout << "Total time for analysis: " << total_time_duration.count() << " minutes. " << endl;

}