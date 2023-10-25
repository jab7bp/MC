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
#include "/work/halla/sbs/jboyd/include/MC_lookups.h"

double par[3];

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


//RUN Info/Parameters
int kine = 8;
int sbsfieldscale = 70;
TString run_target = "LD2";
bool use_particle_gun = false;
double I_beam_uA;

double I_beam;
TString I_beam_str;
double E_beam = lookup_beam_energy_from_kine( kine );

double SBS_field = 1.0*sbsfieldscale/100.0;
double ngen_total = 500000.0;
// double ngen_total = 4000000.0;

TString rootfile_dir;
TFile *outfile;
TChain *TC = new TChain("T");
vector<TString> infile;

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

	double SBS_theta = 35.0; //lookup_HCal_angle_by_kine( kine, "deg" ); //degrees
	double SBS_theta_rad = 35.0*TMath::DegToRad(); //lookup_HCal_angle_by_kine( kine, "rad" ); //radians
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
TH1D *h_Ep, *h_PS, *h_HCal_e, *h_SHPS, *h_weighted_W; 

TH1D *h_W, *h_W2, *h_W_cut, *h_W_fcut, *h_vz_cut;
TH1D *h_Wrecon, *h_W2recon;
TH1D *h_KE_p, *h_KE_low, *h_X, *h_Y, *h_E_eloss;
TH1D *h_Q2, *h_E_ep, *h_E_pp;

TH1D *h_BG_dx, *h_BG_dx_cut, *h_BG_dx_wcut, *h_BG_dx_fcut, *h_BG_dx_wcut_fcut;
TH1D *h_BG_dx_p, *h_BG_dx_n, *h_BG_dx_p_wcut, *h_BG_dx_n_wcut, *h_BG_dx_p_fcut, *h_BG_dx_n_fcut;
TH1D *h_BG_dx_p_wcut_fcut, *h_BG_dx_n_wcut_fcut;
TH1D *h_BG_dy, *h_BG_dy_cut, *h_BG_dy_wcut;
TH2D *h_BG_dxdy, *h_BG_dxdy_cut, *h_BG_dxdy_wcut, *h_BG_dxdy_ncut, *h_BG_dxdy_pcut, *h_BG_dxdy_fcut;
 
TH2D *h_E_ecorr_vs_vert;

TH2D *h_xy, *h_xy_cut, *h_xy_fcut, *h_xy_cut_p, *h_xy_cut_n, *h_PAngleCorr_theta, *h_PAngleCorr_phi;

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

bool is_p = false;
bool is_n = false;

TString magmod = "";

void MC_background(){

	auto total_time_start = high_resolution_clock::now();
	TStopwatch *StopWatch = new TStopwatch();

	gStyle->SetPalette(55);
	cout << "--------------------------------------" << endl;
	cout << "Analysis started. " << endl;
	cout << "--------------------------------------" << endl;

//SBS4 LD2
// double I_beam_uA = 1.75;

//SBS4 LH2
// double I_beam_uA = 3.7393;

// SBS8
//LD2 && LH2
// double I_beam_uA = 5.0;

//SBS9 
//LH2
// double I_beam_uA = 12.0;

//LD2
// double I_beam_uA = 12.0;

	if( kine == 4 ){
		if( run_target == "LH2" ){
			I_beam_uA = 3.5;
		}
		if( run_target == "LD2" ){
			I_beam_uA = 1.75;
		}
	}

	if( kine == 8 ){
		I_beam_uA = 5.0;
	}

	if( kine == 9 ){
		if( run_target == "LH2" ){
			I_beam_uA = 15.0;
		}
		if( run_target == "LD2" ){
			I_beam_uA = 12.0;
		}
	}

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
		HCal_height = 0.0; // Modified to calibrate peak centers for MC: 0.168695, 0.25235289
	}
	if( kine == 9 ){
		HCal_height = 0.0;
	}


	// rootfile_dir = Form("/lustre19/expphy/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass%i/SBS%i/%s/rootfiles", pass, kine, run_target.Data());
	// rootfile_dir = "/volatile/halla/sbs/jboyd/simulation/Rootfiles/";
	rootfile_dir = "/lustre19/expphy/volatile/halla/sbs/jboyd/simulation/out_dir/MC_REPLAY_OUT_DIR";
	// rootfile_dir = "/lustre19/expphy/volatile/halla/sbs/pdbforce/g4sbs_output/sdr/sbs4-sbs50p/simc";

	// rootfile_dir = "/lustre19/expphy/volatile/halla/sbs/adr/sim_results/MC_replay_out";
	I_beam = I_beam_uA*(1.0e-6);
	I_beam_str = Form("%0.2f", I_beam_uA);
	I_beam_str.ReplaceAll(".", "");
	magmod = "0308";
	outfile = new TFile(Form("rootfiles/MC_BACKGROUND_SBS%i_%s_mag%imod%s_%suA_dxdy_23_08_2023.root", kine, run_target.Data(), sbsfieldscale, magmod.Data(), I_beam_str.Data()), "RECREATE");	

//INITIALIZE HISTOGRAMS
	cout << "Initiliazing histograms...";

	h_Ep = new TH1D("h_Ep", Form("E/p (E/p > 0) - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 200, 0, 2);
	h_PS = new TH1D("h_PS", Form("Pre-Shower Clus. E (PS > 0) - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 300, 0, 3);
	h_HCal_e = new TH1D("h_HCal_e", Form("HCal Clus. E (HCal_E > 0) - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 500, 0, 1.0);
	h_SHPS = new TH1D("h_SHPS", Form("SH + PS Clus. E (SH+PS > 0) - SBS%i = %i%%, %s; GeV", kine, sbsfieldscale, run_target.Data()), 500, 0, 5);

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

	h_BG_dx = new TH1D("h_BG_dx",Form("dx (NO CUTS) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);
	h_BG_dx_p = new TH1D("h_BG_dx_p",Form("dx - proton only - (NO CUTS) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);
	h_BG_dx_n = new TH1D("h_BG_dx_n",Form("dx - neutron only - (NO CUTS) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);	
	h_BG_dx_p_wcut = new TH1D("h_BG_dx_p_wcut",Form("dx - proton only - (W Cuts) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);
	h_BG_dx_n_wcut = new TH1D("h_BG_dx_n_wcut",Form("dx - neutron only - (W Cuts) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);	
	h_BG_dx_p_fcut = new TH1D("h_BG_dx_p_fcut",Form("dx - proton only - (F Cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);
	h_BG_dx_n_fcut = new TH1D("h_BG_dx_n_fcut",Form("dx - neutron only - (F Cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);	
	h_BG_dx_p_wcut_fcut = new TH1D("h_BG_dx_p_wcut_fcut",Form("dx - proton only - (W & F Cuts) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);
	h_BG_dx_n_wcut_fcut = new TH1D("h_BG_dx_n_wcut_fcut",Form("dx - neutron only - (W & F Cuts) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);	

	h_BG_dx_cut = new TH1D("h_BG_dx_cut",Form("dx (Basic CUTS) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);
	h_BG_dx_wcut = new TH1D("h_BG_dx_wcut",Form("dx (W cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);
	h_BG_dx_fcut = new TH1D("h_BG_dx_fcut",Form("dx (f cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);
	h_BG_dx_wcut_fcut = new TH1D("h_BG_dx_wcut_fcut",Form("dx (W & f cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);
	h_BG_dy = new TH1D("h_BG_dy",Form("dy (NO CUTS) - SBS%i = %i%%, %s; y_{HCal} - y_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 500, -2.5, 2.5);
	h_BG_dy_cut = new TH1D("h_BG_dy_cut",Form("dy (Basic Cuts) - SBS%i = %i%%, %s; y_{HCal} - y_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 250, -1.25, 1.25);  
	h_BG_dy_wcut = new TH1D("h_BG_dy_wcut",Form("dy (W Cuts) - SBS%i = %i%%, %s; y_{HCal} - y_{exp} (m);", kine, sbsfieldscale, run_target.Data()), 250, -1.25, 1.25);  

	h_BG_dxdy = new TH2D("h_BG_dxdy", Form("Hadron Spot(s) on HCal (NO CUTS) - SBS%i = %i%%, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 300, -1.5, 1.5, 500, -2.5, 2.5 );
	h_BG_dxdy_wcut = new TH2D("h_BG_dxdy_wcut", Form("Hadron Spot(s) on HCal (W cut) - SBS%i = %i%%, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 300, -1.5, 1.5, 500, -2.5, 2.5 );
	h_BG_dxdy_cut = new TH2D("h_BG_dxdy_cut", Form("Hadron Spot(s) on HCal (Basic cuts) - SBS%i = %i%%, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 300, -1.5, 1.5, 500, -2.5, 2.5 );
	h_BG_dxdy_ncut = new TH2D("h_BG_dxdy_ncut", Form("Hadron Spot(s) on HCal (n cut) - SBS%i = %i%%, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 250, -1.25, 1.25, 500, -2.5, 2.5 );
	h_BG_dxdy_pcut = new TH2D("h_BG_dxdy_pcut", Form("Hadron Spot(s) on HCal (p cut) - SBS%i = %i%%, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 250, -1.25, 1.25, 500, -2.5, 2.5 );
	h_BG_dxdy_fcut = new TH2D("h_BG_dxdy_fcut", Form("Hadron Spot(s) on HCal (f cut) - SBS%i = %i%%, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine, sbsfieldscale, run_target.Data()), 250, -1.25, 1.25, 500, -2.5, 2.5 );

	h_xy = new TH2D("h_xy",Form("HCal Hadron Spots (x, y) (NO CUTS) - SBS%i = %i%%, %s;y_{HCal} (m); x_{HCal} (m)", kine, sbsfieldscale, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);
	h_xy_cut = new TH2D("h_xy_cut", Form("HCal Hadron Spots (x, y) (BASIC CUTS) - SBS%i = %i%%, %s;y_{HCal} (m); x_{HCal} (m)", kine, sbsfieldscale, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);
	h_xy_fcut = new TH2D("h_xy_fcut", Form("HCal Hadron Spots (x, y) (Fiduc. CUTS) - SBS%i = %i%%, %s;y_{HCal} (m); x_{HCal} (m)", kine, sbsfieldscale, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);
	h_xy_cut_p = new TH2D("h_xy_cut_p", Form("HCal Hadron Spots (x, y) (p CUT) - SBS%i = %i%%, %s;y_{HCal} (m); x_{HCal} (m)", kine, sbsfieldscale, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);
	h_xy_cut_n = new TH2D("h_xy_cut_n", Form("HCal Hadron Spots (x, y) (n CUT) - SBS%i = %i%%, %s;y_{HCal} (m); x_{HCal} (m)", kine, sbsfieldscale, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);

	h_PAngleCorr_theta = new TH2D( "h_PAngCorr_theta",Form("BB theta vs HCal theta - SBS%i = %i%%, %s", kine, sbsfieldscale, run_target.Data()), 200, 0.55, 0.75, 300, 0.35, 0.65 );
	h_PAngleCorr_phi = new TH2D( "h_PAngCorr_phi",Form("BB phi vs HCal phi - SBS%i = %i%%, %s", kine, sbsfieldscale, run_target.Data()), 500, -0.4, 0.1, 500, 2.7, 3.2 );
	h_vz_cut = new TH1D("h_vz_cut",Form("BB phi vs HCal phi - SBS%i = %i%%, %s; vertex z (m);", kine, sbsfieldscale, run_target.Data()), 250,-0.125,0.125);	

	h_weighted_W = new TH1D("h_weighted_W", Form("MC Weight - SBS%i = %i%%, %s; vertex z (m);", kine, sbsfieldscale, run_target.Data()), 200, 0.0, 2.0);

	cout << "finished. " << endl << endl;

	cout << "--------------------------------------" << endl;
	cout << "Adding files to TChain from: " << rootfile_dir.Data() << endl;

	if( kine == 4 ){
		if( run_target == "LD2"){
			infile = {Form("%s/replayed_jb_gmn_SBS4_LD2_mag70_175uA_elas_100k_job*", rootfile_dir.Data())};			
		}
		else{
			infile = {Form("%s/replayed_jb_gmn_SBS4_LH2_mag0_350uA_pGun_100k_job*", rootfile_dir.Data())};
		}
		// infile = Form("%s/replayed_jb_gmn_SBS%i_%s_mag%imod1_175uA_elas_100k_job*", rootfile_dir.Data(), kine, run_target.Data(), sbsfieldscale);
		// TC->Add(infile.Data());
		
		for( size_t file = 0; file < infile.size(); file++ ){
			cout << "Adding file: " << endl;
			cout << infile[file].Data() << endl << endl;
			TC->Add( infile[file].Data() );	
		}
		
	}

	if( kine == 8 ){
		// if( run_target == "LD2" ){
		// 	infile = Form("%s/replayed_jb_gmn_SBS%i_%s_mag%imod3_%suA_elas_100k_job*", rootfile_dir.Data(), kine, run_target.Data(), sbsfieldscale, I_beam_str.Data());
		// }
		// if( run_target == "LH2" ){
		// 	// infile = Form("%s/replayed_jb_gmn_SBS%i_%s_mag%i_%suA_elas_100k_job*", rootfile_dir.Data(), kine, run_target.Data(), sbsfieldscale, I_beam_str.Data());
		// 	infile = "/lustre19/expphy/volatile/halla/sbs/jboyd/simulation/out_dir/MC_REPLAY_OUT_DIR/replayed_jb_gmn_SBS8_LH2_mag70_elas_500uA_50k_19_07_2023_job*";
		// }
		// infile = Form("/lustre19/expphy/volatile/halla/sbs/jboyd/simulation/out_dir/MC_REPLAY_OUT_DIR/replayed_jb_gmn_SBS8_LD2_mag70mod3_500uA_elas_100k_27_07_2023_job*");
		// infile = Form("/lustre19/expphy/volatile/halla/sbs/jboyd/simulation/out_dir/MC_REPLAY_OUT_DIR/replayed_jb_gmn_SBS8_LH2_mag70mod0815_500uA_elas_100k_07_08_2023_job*");
		
		//field a litle off but god
		// infile = Form("/lustre19/expphy/volatile/halla/sbs/jboyd/simulation/out_dir/MC_REPLAY_OUT_DIR/replayed_jb_gmn_SBS8_LD2_mag70mod0815_500uA_elas_100k_08_08_2023_job*.root");
		// infile = Form("/lustre19/expphy/volatile/halla/sbs/jboyd/simulation/out_dir/MC_REPLAY_OUT_DIR/replayed_jb_gmn_SBS8_LD2_mag70mod%s_500uA_elas_100k_09_08_2023_job*", magmod.Data());
		infile = {"/lustre19/expphy/volatile/halla/sbs/jboyd/simulation/out_dir/MC_REPLAY_OUT_DIR/replayed_jb_gmn_SBS8_LD2_inelastic_mag70port0308_500uA_elas_voffsetFix_100k_22_08_2023_job*.root",
			"/lustre19/expphy/volatile/halla/sbs/jboyd/simulation/out_dir/MC_REPLAY_OUT_DIR/replayed_jb_gmn_SBS8_LD2_inelastic_mag70port0308_500uA_elas_voffsetFix_500k_22_08_2023_job*.root"
		};
		// infile = "/lustre19/expphy/volatile/halla/sbs/seeds/sbs8_inelastic/replayed_inelastic_sbs8_0p_g4sbs_deeN_job_*.root";

		// TC->Add(infile.Data());
		// TC->Add(Form("%s/replayed_adr_gmn_SBS4_LD2_mag30_elas_1M_job0_1.root", rootfile_dir.Data()));

		for( size_t file = 0; file < infile.size(); file++ ){
			cout << "Adding file: " << endl;
			cout << infile[file].Data() << endl << endl;
			TC->Add( infile[file].Data() );	
		}
	
	}

	if( kine == 9 ){
		//John's
		infile = {Form("%s/replayed_jb_gmn_SBS9_LH2_mag70mod0815_sbsAng225_1500uA_elas_100k_09_08_2023_job_*.root", rootfile_dir.Data())};
		//Anu's
		// infile = Form("%s/replayed_gmn_SBS9_job_*", rootfile_dir.Data());
		// infile = Form("%s/replayed_jb_gmn_SBS%i_%s_mag%imod1_175uA_elas_100k_job*", rootfile_dir.Data(), kine, run_target.Data(), sbsfieldscale);
		for( size_t file = 0; file < infile.size(); file++ ){
			cout << "Adding file: " << endl;
			cout << infile[file].Data() << endl << endl;
			TC->Add( infile[file].Data() );
	}
		// TC->Add(Form("%s/replayed_adr_gmn_SBS4_LD2_mag30_elas_1M_job0_1.root", rootfile_dir.Data()));	
	}
	// else{
	// 	infile = Form("%s/replayed_jb_gmn_SBS%i_%s_mag%imod3_%suA_elas_100k_job*", rootfile_dir.Data(), kine, run_target.Data(), sbsfieldscale, I_beam_str.Data());
	// 	cout << "Adding file: " << endl;
	// 	cout << infile.Data() << endl << endl;
	// 	TC->Add(infile.Data());		
	// }
	if( use_particle_gun ){
		cout << endl << "Using particle gun.... " << endl;
		cout << "Adding file: " << endl;
		infile = {Form("%s/replayed_jb_gmn_SBS4_LD2_mag0_100uA_pGun_Px_100k_v2_job0.root", rootfile_dir.Data())};
		cout << infile[0].Data() << endl << endl;
	}
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

	if( kine == 8 ){
		master_cut_vec = {
			"sbs.hcal.nclus>0",
			"bb.ps.nclus>0",
			"bb.sh.nclus>0",
			"abs(bb.tr.vz[0])<=0.075",
			"bb.gem.track.nhits[0]>3",
			"bb.tr.n==1",
			Form("((abs(((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))))-%f)<%f", lookup_MC_cut(kine, sbsfieldscale, run_target, "Ep"), 1.0*lookup_MC_cut(kine, sbsfieldscale, run_target, "Ep_sigma")),
			// Form("(bb.sh.e+bb.ps.e)>%f", lookup_MC_cut(kine, sbsfieldscale, run_target, "SH_PS_mean") - 3.0*lookup_MC_cut(kine, sbsfieldscale, run_target, "SH_PS_sigma")),
			Form("bb.ps.e>%f", 0.15),
			Form("sbs.hcal.e>%f", 0.01)
		};		
	}
	if( kine == 9 ){
		master_cut_vec = {
			"sbs.hcal.nclus>0",
			"bb.ps.nclus>0",
			"bb.sh.nclus>0",
			"abs(bb.tr.vz[0])<=0.075",
			"bb.gem.track.nhits[0]>3",
			"bb.tr.n==1",
			// Form("((abs(((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))))-%f)<%f", lookup_MC_cut(kine, sbsfieldscale, run_target, "Ep"), 1.0*lookup_MC_cut(kine, sbsfieldscale, run_target, "Ep_sigma")),
			// Form("(bb.sh.e+bb.ps.e)>%f", lookup_MC_cut(kine, sbsfieldscale, run_target, "SH_PS_mean") - 3.0*lookup_MC_cut(kine, sbsfieldscale, run_target, "SH_PS_sigma")),
			// Form("bb.ps.e>%f", 0.15),
			// Form("sbs.hcal.e>%f", 0.01)
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
	cout << "Applying Master Cut: " << endl;
	cout << master_cut << endl;

	TEventList *ev_list = new TEventList("ev_list", "Elastic Events List");
	TC->Draw(">>ev_list", master_cut);
	Nevents = ev_list->GetN();

	cout << "---- Raw events in TC: " << TC->GetEntries() << endl;
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
		
		if( int(mc_fnucl) == 0 ){ 
			is_n = true; 
			is_p = false;
		}
		if( int(mc_fnucl) == 1 ){ 
			is_p = true; 
			is_n = false;
		}

		if( calc_pn_weight ){
			pn_weight = ((mc_omega*mc_sigma)*luminosity)/ngen_total;		

			// cout << "pn weight: " << pn_weight << endl;	
		}
		else{ pn_weight = 1.0; }

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

		luminosity = calc_luminosity(I_beam, run_target );
		// luminosity = 1.0;

		if( Ep > 0.0 ){
			h_Ep->Fill(Ep);
		}

		if( PS > 0.0 ){
			h_PS->Fill(PS);			
		}

		if( SHPS > 0.0 ){
			h_SHPS->Fill(SHPS);			
		}

		if( hcal_e > 0.0 ){
			h_HCal_e->Fill(hcal_e);		
		}

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
		if( W > 0 ){
			h_weighted_W->Fill(W*pn_weight);			
		}

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
		h_BG_dx->Fill( dx, pn_weight);
		//Selectively fill for proton or neutron:
		if( is_p ){
			h_BG_dx_p->Fill( dx, pn_weight );
		}
		if( is_n ){
			h_BG_dx_n->Fill( dx, pn_weight );
		}

		h_BG_dy->Fill( dy, pn_weight);
		h_BG_dxdy->Fill( dy, dx, pn_weight );
		h_xy->Fill( hcal_y, hcal_x );

	// Preliminary HCal projections with single cut on W
		if( fabs(W - W_mean) < W_sigma ){
			h_BG_dx_wcut->Fill( dx, pn_weight );
			h_BG_dy_wcut->Fill ( dy, pn_weight );
			h_BG_dxdy_wcut->Fill( dy, dx, pn_weight );

			//Histogram by hadron
			if( is_p ){
				h_BG_dx_p_wcut->Fill(dx, pn_weight);
			}
			if( is_n ){
				h_BG_dx_n_wcut->Fill(dx, pn_weight);
			}
		}

	//Populate BB/HCal correlation histograms from elastics
		h_PAngleCorr_phi->Fill( e_prime_phi, nucleon_phi );
		h_PAngleCorr_theta->Fill( e_prime_theta, nucleon_theta );

	//Fill vertex position histogram for cut on tracks
    	h_vz_cut->Fill( bb_tr_vz[0] );

	//FIDUCIAL Cut
		//Check "elastic" events on center HCal for id with spot checks
		bool HCal_on = false;

		if( fiducial_cut ){
			if( kine == 4 && sbsfieldscale == 30 ){
				dx_p_scale = 1.0;
				dx_n_scale = 1.0;
			}
			if( kine == 8 && sbsfieldscale == 70 ){
				dx_p_scale = 1.0;
				dx_n_scale = 1.0;
			}

			if( sbsfieldscale != 0 ){
				dx_p = lookup_MC_dxdy(kine, sbsfieldscale, run_target, "dx_p");
				dx_p_sigma = dx_p_scale*lookup_MC_dxdy(kine, sbsfieldscale, run_target, "dx_p_sigma");
				dy_p = lookup_MC_dxdy(kine, sbsfieldscale, run_target, "dy");
				dy_p_sigma = lookup_MC_dxdy(kine, sbsfieldscale, run_target, "dy_sigma");
				dx_n = lookup_MC_dxdy(kine, sbsfieldscale, run_target, "dx_n");
				dx_n_sigma = dx_n_scale*lookup_MC_dxdy(kine, sbsfieldscale, run_target, "dx_n_sigma");
				dy_n = lookup_MC_dxdy(kine, sbsfieldscale, run_target, "dy");
				dy_n_sigma = lookup_MC_dxdy(kine, sbsfieldscale, run_target, "dy_sigma");
				dx_pn_max = abs( dx_p ) + abs( dx_n );
			}	
			else{
				dx_p = 0;
				dx_p_sigma = 0.16;
				dy_p = 0;
				dy_p_sigma = 0.21;
				dx_n = 0;
				dx_n_sigma = 0.16;
				dy_n = 0;
				dy_n_sigma = 0.21;
				dx_pn_max = 0;
			}
		
			if( hcal_y > hcal_y_fmin && hcal_y < hcal_y_fmax && hcal_x >hcal_x_fmin && hcal_x < hcal_x_fmax ){
				HCal_on = true;
			}

			// apply_fcut = ((y_expected_HCal - dy_p_sigma) > hcal_y_fmin) && ((y_expected_HCal + dy_p_sigma) < hcal_y_fmax) && ((x_expected_HCal - dx_pn_max - dx_p_sigma) > hcal_x_fmin) && ((x_expected_HCal + dx_p_sigma) < hcal_x_fmax);
		//Modified Aug 10, 2023:
			apply_fcut = ((y_expected_HCal - dy_p_sigma) > hcal_y_fmin) && ((y_expected_HCal + dy_p_sigma) < hcal_y_fmax) && ((x_expected_HCal - dx_pn_max - dx_p_sigma) > hcal_x_fmin) && ((x_expected_HCal + dx_pn_max + dx_p_sigma) < hcal_x_fmax);

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


		//Fill respective histograms for these checks.
			if( HCal_on && is_n && apply_fcut ) h_BG_dxdy_ncut->Fill( dy, dx );
			if( HCal_on && is_p && apply_fcut ) h_BG_dxdy_pcut->Fill( dy, dx );

		//----------neutron
			if( HCal_on && is_n && apply_fcut ){
				if( !calc_pn_weight ){
					pn_weight = (1.0/3.0);					
				}
				// if( (hcal_x - dx_pn_max )>hcal_x_fmin ){
					h_BG_dxdy_fcut->Fill( dy, dx, pn_weight );
					h_BG_dx_fcut->Fill( dx, pn_weight );
					h_W_fcut->Fill( W, pn_weight );
					h_xy_fcut->Fill( hcal_y, hcal_x, pn_weight );
					h_xy_cut_n->Fill( hcal_y, hcal_x );
					h_BG_dx_n_fcut->Fill(dx, pn_weight);
					if( fabs(W - W_mean) < W_sigma ){
						h_BG_dx_wcut_fcut->Fill(dx, pn_weight);
						h_BG_dx_n_wcut_fcut->Fill(dx, pn_weight);
					}
					elastic_yield++;
				// }
			}		
		//----------proton
			else if( HCal_on && is_p && apply_fcut ){
				if( !calc_pn_weight ){
					pn_weight = 1.0;					
				}

				// if( (hcal_x + dx_pn_max)<hcal_x_fmax ){
					h_BG_dxdy_fcut->Fill( dy, dx, pn_weight );
					h_BG_dx_fcut->Fill( dx, pn_weight );
					h_W_fcut->Fill( W, pn_weight );
					h_xy_fcut->Fill( hcal_y, hcal_x, pn_weight );
					h_xy_cut_p->Fill( hcal_y, hcal_x );
					h_BG_dx_p_fcut->Fill(dx, pn_weight);
					if( fabs(W - W_mean) < W_sigma ){
						h_BG_dx_wcut_fcut->Fill(dx, pn_weight);
						h_BG_dx_p_wcut_fcut->Fill(dx, pn_weight);
					}
					elastic_yield++;
				// }
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

	TFile *infile = new TFile(outfile->GetName(), "READ");
	TH1D *hin_dx, *hin_dx_fcut, *hin_dx_wcut_fcut;
	TH2D *hin_dxdy, *hin_dxdy_fcut, *hin_dxdy_wcut;

	hin_dx = static_cast<TH1D*>(infile->Get("h_BG_dx"));
	hin_dx_fcut = static_cast<TH1D*>(infile->Get("h_BG_dx_fcut"));
	hin_dx_wcut_fcut = static_cast<TH1D*>(infile->Get("h_BG_dx_wcut_fcut"));
	hin_dxdy = static_cast<TH2D*>(infile->Get("h_BG_dxdy"));
	hin_dxdy_fcut = static_cast<TH2D*>(infile->Get("h_BG_dxdy_fcut"));
	hin_dxdy_wcut = static_cast<TH2D*>(infile->Get("h_BG_dxdy_wcut"));

    TCanvas *c_dxdy = new TCanvas("c_dxdy", "c_dxdy", 600, 500);
    hin_dxdy->Draw("colz");

    TCanvas *c_dxdy_fcut = new TCanvas("c_dxdy_fcut", "c_dxdy_fcut", 600, 500);
    hin_dxdy_fcut->Draw("colz");

    TCanvas *c_dxdy_wcut = new TCanvas("c_dxdy_wcut", "c_dxdy_wcut", 600, 500);
    hin_dxdy_wcut->Draw("colz");

    TCanvas *c_dx_fcut = new TCanvas("c_dx_fcut", "c_dx_fcut", 600, 500);
    hin_dx_fcut->Draw("hist");

    TCanvas *c_dx_wcut_fcut = new TCanvas("c_dx_wcut_fcut", "c_dx_wcut_fcut", 600, 500);
    hin_dx_wcut_fcut->Draw("hist");

	if( calibrate ){
    	TFile *infile = new TFile(outfile->GetName(), "READ");
    	TH1D *hin_dx_wcut = static_cast<TH1D*>(infile->Get("h_BG_dx_wcut"));
    	TH1D *hin_dy_wcut = static_cast<TH1D*>(infile->Get("h_BG_dy_wcut"));

    	TCanvas *c_dx = new TCanvas("c_dx", "c_dx", 600, 500);
    	hin_dx->Draw("hist");

       	TCanvas *c_dx_wcut = new TCanvas("c_dx_wcut", "c_dx_wcut", 600, 500);
    	hin_dx_wcut->Draw("hist");
  	
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
			fit_dx_p->SetParLimits(1, -1.2, -0.8);
			fit_dx_p->SetParLimits(2, 0.1, 0.14);			
		}

		if( kine == 9 && sbsfieldscale == 70 && run_target == "LH2"){
			fit_dx_p->SetParLimits(0, 0.8*hin_dx_wcut->GetMaximum(), hin_dx_wcut->GetMaximum());
			fit_dx_p->SetParLimits(1, -0.8, -0.7);
			fit_dx_p->SetParLimits(2, 0.1, 0.16);			
		}
		if( kine == 9 && sbsfieldscale == 70 && run_target == "LD2"){
			fit_dx_p->SetParLimits(0, 0.8*hin_dx_wcut->GetMaximum(), hin_dx_wcut->GetMaximum());
			fit_dx_p->SetParLimits(1, -0.85, -0.7);
			fit_dx_p->SetParLimits(2, 0.1, 0.18);			
		}
	
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
			fit_dx_n->SetParLimits(1, -0.20, 0.05);
			fit_dx_n->SetParLimits(2, 0.1, 0.14);
		}
		if( kine == 9 && sbsfieldscale == 70){
			fit_dx_n->SetParLimits(0, 0, (0.35)*hin_dx_wcut->GetMaximum());
			fit_dx_n->SetParLimits(1, -0.05, 0.05);
			fit_dx_n->SetParLimits(2, 0.1, 0.16);
		}	
		
		if( run_target != "LH2" ){
			hin_dx_wcut->Fit("fit_dx_n", "R+");
			dx_n = fit_dx_n->GetParameter(1);
			dx_n_sigma = fit_dx_n->GetParameter(2);				
		}


	//------- dy -------
    	TCanvas *c_dy = new TCanvas("c_dy", "c_dy", 600, 500);
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
		if( kine == 9 && sbsfieldscale == 70){
			fit_dy->SetParLimits(0, 0.8*hin_dy_wcut->GetMaximum(), hin_dy_wcut->GetMaximum());
			fit_dy->SetParLimits(1, -0.15, 0.15);
			fit_dy->SetParLimits(2, 0.0, 0.15);	
		}
		else{
			fit_dy->SetParLimits(0, 0, hin_dy_wcut->GetMaximum());
			fit_dy->SetParLimits(1, -0.15, 0.15);
			fit_dy->SetParLimits(2, 0.1, 0.23);				
		}

	
		hin_dy_wcut->Fit("fit_dy", "R+");
		dy_p = fit_dy->GetParameter(1);
		dy_p_sigma = fit_dy->GetParameter(2);	
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
	cout << "Output file: " << outfile->GetName() << endl << endl;
	cout << "------------------------------------------------------------------"<< endl;
	
	auto total_time_end = high_resolution_clock::now();
	auto total_time_duration = duration_cast<minutes>(total_time_end - total_time_start);
	cout << "Total time for analysis: " << total_time_duration.count() << " minutes. " << endl;

}