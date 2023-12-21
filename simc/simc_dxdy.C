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
#include "/w/halla-scshelf2102/sbs/jboyd/include/utility_functions.h"
#include "/work/halla/sbs/jboyd/include/MC_lookups.h"

double par[3], par2gaus[6];

Double_t fit_gaus(Double_t * x, Double_t *par){

	Double_t g = 0.0;

	g = par[0]*exp((-0.5)*pow(((x[0] -  par[1])/par[2]),2));
	
	return g;
}

Double_t fit_two_gaus(Double_t *x, Double_t *par2gaus){
	Double_t g1 = 0.0;
	Double_t g2 = 0.0;

	g1 = par2gaus[0]*exp((-0.5)*pow(((x[0] -  par2gaus[1])/par2gaus[2]),2));
	g2 = par2gaus[3]*exp((-0.5)*pow(((x[0] -  par2gaus[4])/par2gaus[5]),2));

	return g1 + g2;
}

bool match_file_cnts = true;

//---debugging:
bool manually_select_file_cnt = false;
int selected_file_cnt = 25;
TString manual_select_file_cnt_str = "";
//------

bool use_hist_file = true;
bool plot_only = false;

bool fiducial_cut = true;
bool apply_fcut;
bool calc_pn_weight = true;	

bool correct_beam_energy = false; ///DONT USE FOR SIMULATION
bool calibrate = true;

TString hadron_select = "both";

//RUN Info/Parameters
int kine = 9;
int sbsfieldscale = 70;
TString run_target = "LD2";
bool use_particle_gun = false;

double I_beam_uA;

double I_beam;
TString I_beam_str;
double E_beam;


double ngen_total = 100000.0;
// double ngen_total = 4000000.0;

TString rootfile_dir;
TString dateInfile = "";

TFile *outfile;
TChain *TC = new TChain("T");
TString infile, histfile;
TFile *currentChainFile;
TString currentChainFileName = "";
int current_file_index = -1, current_nfile_index = -1, current_pfile_index = -1;
int current_file_jobNum_index = -1, current_file_jobNum = -1, current_nfile_jobNum = -1, current_pfile_jobNum = -1;

vector<TString> master_cut_vec;
TString master_cut_string;
TCut master_cut = "";

TString histfile_dir = "/lustre19/expphy/volatile/halla/sbs/jboyd/simulation/out_dir/MC_REPLAY_OUT_DIR/hist";


//Experimental Parameters
//BigBite, BBCal
double BB_dist, BB_theta, BB_theta_rad;

//SBS, HCal
double HCal_dist, HCal_theta, HCal_theta_rad;

double SBS_theta, SBS_theta_rad, scint_intersect, x_expected_HCal, y_expected_HCal;

	//HCal Fiducial cut parameters:
	// double hcal_y_fmin = -0.75;
	// double hcal_y_fmax = 0.75;
	// double hcal_x_fmin = -2.015;
	// double hcal_x_fmax = 1.285;

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

TH1D *h_simc_dx, *h_simc_dx_cut, *h_simc_dx_wcut, *h_simc_dx_fcut, *h_simc_dx_wcut_fcut;
TH1D *h_simc_dy, *h_simc_dy_cut, *h_simc_dy_wcut;
TH2D *h_simc_W_dx, *h_simc_W2_dx;
TH2D *h_simc_dxdy, *h_simc_dxdy_cut, *h_simc_dxdy_wcut, *h_simc_dxdy_ncut, *h_simc_dxdy_pcut, *h_simc_dxdy_fcut;

//SEPRATED HISTOGRAMS BASED ON HADRON
TH1D *h_simc_dx_p, *h_simc_dx_n, *h_simc_dx_p_wcut, *h_simc_dx_n_wcut;
TH1D *h_simc_dx_p_fcut, *h_simc_dx_n_fcut, *h_simc_dx_p_wcut_fcut, *h_simc_dx_n_wcut_fcut;
 
TH2D *h_E_ecorr_vs_vert;

TH2D *h_xy, *h_xy_cut, *h_xy_fcut, *h_xy_cut_p, *h_xy_cut_n, *h_PAngleCorr_theta, *h_PAngleCorr_phi;

TH1D *h_mc_p_weight, *h_mc_n_weight, *h_mc_p_final_weight, *h_mc_n_final_weight, *h_mc_n_p_weight_Ratio, *h_mc_n_p_final_weight_Ratio;
TH1D *h_mc_p_sigma, *h_mc_n_sigma, *h_mc_np_sigma_ratio;
TH1D *h_mc_p_sigma_weighted, *h_mc_n_sigma_weighted, *h_mc_np_sigma_ratio_weighted;
TH1D *h_mc_simc_Q2, *h_mc_simc_nu, *h_mc_simc_epsilon, *h_mc_simc_Ebeam, *h_mc_simc_p_n, *h_mc_simc_theta_n, *h_mc_simc_p_e, *h_mc_simc_theta_e;
TH1D *h_mc_simc_Q2_wcut, *h_mc_simc_nu_wcut, *h_mc_simc_epsilon_wcut, *h_mc_simc_Ebeam_wcut;
TH1D *h_mc_simc_p_n_wcut, *h_mc_simc_theta_n_wcut, *h_mc_simc_p_e_wcut, *h_mc_simc_theta_e_wcut;


//BRANCH VARIABLES

double bb_tr_p[maxTracks], bb_tr_px[maxTracks], bb_tr_py[maxTracks], bb_tr_pz[maxTracks];
double bb_tr_vx[maxTracks], bb_tr_vy[maxTracks], bb_tr_vz[maxTracks], bb_tr_chi2[maxTracks];
double bb_fp_x[maxTracks], bb_fp_y[maxTracks], bb_fp_th[maxTracks], bb_fp_ph[maxTracks];
double bb_tgt_x[maxTracks], bb_tgt_y[maxTracks], bb_tgt_th[maxTracks], bb_tgt_ph[maxTracks];
double bb_tr_n, bb_ps_x, bb_ps_y, bb_ps_e, bb_sh_x, bb_sh_y, bb_sh_e;
Double_t nblk, nclus, SH_nclus, PS_nclus, hcal_x, hcal_y, hcal_e;
double mc_omega, mc_sigma, mc_fnucl, mc_Weight, luminosity, mc_luminosity, mc_genvol, mc_Final_Weight;
double mc_simc_Q2, mc_simc_nu, mc_simc_epsilon, mc_simc_Ebeam, mc_simc_p_n, mc_simc_theta_n, mc_simc_p_e, mc_simc_theta_e;
//variables for SIMC weights
	//protons
double mc_p_luminosity, mc_p_genvol, mc_p_Final_Weight;
vector<double> mc_p_Final_Weight_vec = {}, mc_p_luminosity_vec = {}, mc_p_genvol_vec = {};
int mc_p_Ntried;
vector<int> mc_p_Ntried_vec = {};
vector<int> mc_p_rootfile_jobnum = {};
vector<int> mc_p_histfile_jobnum = {};

	//neutrons
double mc_n_luminosity, mc_n_genvol, mc_n_Final_Weight;
vector<double> mc_n_Final_Weight_vec = {}, mc_n_luminosity_vec = {}, mc_n_genvol_vec = {};
int mc_n_Ntried;
vector<int> mc_n_Ntried_vec = {};
vector<int> mc_n_rootfile_jobnum = {};
vector<int> mc_n_histfile_jobnum = {};

Long64_t Nevents;

double dx_p_scale = 1.0;
double dx_n_scale = 1.0;

bool is_p = false;
bool is_n = false;

TString SBS4_infile_basename, SBS8_infile_basename, SBS9_infile_basename;
TString proton_infile, proton_infile_basename, neutron_infile, neutron_infile_basename;
vector<TString> proton_infile_vec = {};
vector<TString> proton_histfile_vec = {};
vector<Int_t> proton_infile_jobNum_vec = {};

vector<TString> neutron_infile_vec = {};
vector<TString> neutron_histfile_vec = {};
vector<Int_t> neutron_infile_jobNum_vec = {};

int nucleon_with_min_file_cnt; //indicator of which nucleon (p or n) has less simulation files. 0: nucleon, 1: proton
int proton_infile_cnt = 0, neutron_infile_cnt = 0;
vector<int> nucleon_infile_cnt = {0, 0};
TString nucleon_with_min_file_cnt_string;

TH1D *hin_dx_wcut, *hin_dy_wcut;

TString outfilename = "";
TString portField = "";
TString numEvents_simFile_str = "";

void simc_dxdy(int kine_select = -1, int sbsfieldscale_select = -1, TString portField_select = "" ){
	
	auto total_time_start = high_resolution_clock::now();
	TStopwatch *StopWatch = new TStopwatch();

	gStyle->SetPalette(55);
	cout << "--------------------------------------" << endl;
	cout << "Analysis started. " << endl;
	cout << "--------------------------------------" << endl;

	if( kine_select != -1 ){
		kine_select = kine;
		cout << "Terminal-line-defined kinematic: " << kine_select << endl;
	}
	else{
		kine_select = kine;
	}
	if( sbsfieldscale_select != -1 ){
		sbsfieldscale_select = sbsfieldscale;
		cout << "Terminal-line-defined sbsfieldscale: " << sbsfieldscale_select << endl;
	}
	else{
		sbsfieldscale_select = sbsfieldscale;
	}

	I_beam_uA =  lookup_beam_current( kine_select, run_target );
	E_beam = lookup_beam_energy_from_kine( kine_select );

	BB_dist = lookup_BB_dist_by_kine( kine );  //Distance to BigBite magnet
	BB_theta = lookup_BB_angle_by_kine( kine, "deg" ); //degrees, BB arm angle
	BB_theta_rad = lookup_BB_angle_by_kine( kine, "rad" ); //radians, BB arm angle

	//SBS, HCal
	HCal_dist = lookup_HCal_dist_by_kine( kine ); //Distance to Hcal face form target chamber
	HCal_theta = lookup_SBS_angle_by_kine (kine, "deg" ); //degrees, Angle for downsream arm to HCal
	HCal_theta_rad = lookup_SBS_angle_by_kine (kine, "rad" ); //radians, Angle for downsream arm to HCal

	SBS_theta = lookup_SBS_angle_by_kine( kine, "deg" ); //degrees
	SBS_theta_rad = lookup_SBS_angle_by_kine( kine, "rad" ); //radians
	scint_intersect, x_expected_HCal, y_expected_HCal;



	W_mean = lookup_simc_cut(kine_select, sbsfieldscale_select, run_target, "W");
	W_sigma = lookup_simc_cut(kine_select, sbsfieldscale_select, run_target, "W_sigma");
	W2_mean = lookup_simc_cut(kine_select, sbsfieldscale_select, run_target, "W2");
	W2_sigma = lookup_simc_cut(kine_select, sbsfieldscale_select, run_target, "W2_sigma");

	//Set defaults???
	if( W_mean == -1 ){ W_mean = Mp; }
	if( W_sigma == -1 ){ W_sigma = 0.15; }
	if( W2_mean == -1){ W2_mean = pow(Mp, 2); }
	if( W2_sigma == -1){ W2_sigma = 0.25; }

	// HCal_height = -0.2897; // Height of HCal above beamline     ---- ORIGINAL ----
	if( kine_select == 4 ){
		HCal_height = 0.0; // Modified to calibrate peak centers for MC	
	}
	if( kine_select == 8 ){
		HCal_height = 0.0; // Modified to calibrate peak centers for MC: 0.168695, 0.25235289
	}
	if( kine_select == 9 ){
		HCal_height = 0.0;
	}


	// rootfile_dir = Form("/lustre19/expphy/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass%i/SBS%i/%s/rootfiles", pass, kine_select, run_target.Data());
	// rootfile_dir = "/volatile/halla/sbs/jboyd/simulation/Rootfiles/";
	rootfile_dir = Form( "/lustre19/expphy/volatile/halla/sbs/jboyd/simulation/out_dir/MC_REPLAY_OUT_DIR/simc/SBS%i", kine_select );
	// rootfile_dir = "/lustre19/expphy/volatile/halla/sbs/pdbforce/g4sbs_output/sdr/sbs4-sbs50p/simc";

	// rootfile_dir = "/lustre19/expphy/volatile/halla/sbs/adr/sim_results/MC_replay_out";
	I_beam = I_beam_uA*(1.0e-6);
	I_beam_str = Form("%0.2f", I_beam_uA);
	I_beam_str.ReplaceAll(".", "");

	if( manually_select_file_cnt ){
		manual_select_file_cnt_str =  Form("_manualFilCnt%i", selected_file_cnt);
	}

	if( kine_select == 4 ){
		portField = "0387";
		dateInfile = "10_12_2023";
	}
	if( kine_select == 8 ){
		portField = "0653";
		dateInfile = "04_11_2023";		
	}
	if( kine_select == 9 ){
		portField = "0666";		
		dateInfile = "12_19_2023";	
	}

	if( portField_select != "" ){
		portField = portField_select;
		cout << "Terminal-line-defined portField: " << portField_select.Data() << endl;
		cout << "---------------------------------------" << endl;
		cout << endl;
	}
	else{
		portField_select = portField;
	}

	outfilename = Form("rootfiles/simc_SBS%i_%s_mag%i_%suA_port%s_%s%s.root", kine_select, run_target.Data(), sbsfieldscale_select, I_beam_str.Data(), portField_select.Data(), dateInfile.Data(), manual_select_file_cnt_str.Data() );

	if( !plot_only ){

		outfile = new TFile(outfilename.Data(), "RECREATE");	

	//INITIALIZE HISTOGRAMS
		cout << "Initiliazing histograms...";

		h_Ep = new TH1D("h_Ep", Form("E/p - SBS%i = %i%%, %s; GeV", kine_select, sbsfieldscale_select, run_target.Data()), 200, 0, 2);
		h_PS = new TH1D("h_PS", Form("Pre-Shower Clus. E - SBS%i = %i%%, %s; GeV", kine_select, sbsfieldscale_select, run_target.Data()), 300, 0, 3);
		h_HCal_e = new TH1D("h_HCal_e", Form("HCal Clus. E - SBS%i = %i%%, %s; GeV", kine_select, sbsfieldscale_select, run_target.Data()), 200, 0, 0.4);
		h_SHPS = new TH1D("h_SHPS", Form("SH + PS Clus. E - SBS%i = %i%%, %s; GeV", kine_select, sbsfieldscale_select, run_target.Data()), 500, 0, 5);

		h_E_eloss = new TH1D("E_eloss", Form("Scattered Electron Energy Loss in Target - SBS%i = %i%%, %s", kine_select, sbsfieldscale_select, run_target.Data()), 500, 0.0, (0.1)*E_beam);
		h_E_ecorr_vs_vert = new TH2D("h_E_ecorr_vs_vert", Form("Corrected Beam Energy vs Vertex - SBS%i = %i%%, %s; E_{e} (GeV); Z_{vertex} (m)", kine_select, sbsfieldscale_select, run_target.Data()), 250, -0.125, 0.125, 500, 0, 0.001);
		h_Q2 = new TH1D("h_Q2", Form("Momentum Transfer Q^2 - SBS%i = %i%%, %s; GeV", kine_select, sbsfieldscale_select, run_target.Data()), 750, 0.5, 9.0);
		h_E_ep = new TH1D("h_E_ep", Form("Scattered Electron Energy - SBS%i = %i%%, %s; GeV", kine_select, sbsfieldscale_select, run_target.Data()), 500, 0.0, 1.5*E_beam);
		h_E_pp = new TH1D("h_E_pp", Form("Scattered Proton Energy - SBS%i = %i%%, %s; GeV", kine_select, sbsfieldscale_select, run_target.Data()), 500, 0.0, 1.5*E_beam);
		h_KE_p = new TH1D("h_KE_p", Form("Scattered Proton kine_selecttic Energy - SBS%i = %i%%, %s; GeV", kine_select, sbsfieldscale_select, run_target.Data()), 500, 0.0, 1.5*E_beam);

		h_W = new TH1D("h_W", Form("Invariant Mass W - SBS%i = %i%%, %s; GeV", kine_select, sbsfieldscale_select, run_target.Data()), 300, 0.0, 3.0);
		h_W2 = new TH1D("h_W2", Form("Invariant Mass W^2 - SBS%i = %i%%, %s; GeV", kine_select, sbsfieldscale_select, run_target.Data()), 300, 0.0, 3.0);
		h_W_cut = new TH1D("h_W_cut", Form("Invariant Mass W (Coin & Vert Cuts) - SBS%i = %i%%, %s; GeV", kine_select, sbsfieldscale_select, run_target.Data()), 300, 0.0, 3.0);
		h_W_fcut = new TH1D("h_W_fcut", Form("Invariant Mass W (Fiduc. Cuts) - SBS%i = %i%%, %s; GeV", kine_select, sbsfieldscale_select, run_target.Data()), 300, 0.0, 3.0);

		h_simc_W_dx = new TH2D("h_simc_W_dx", Form("Invariant Mass W vs dx - SBS%i, mag%i, %s; x_{HCal} - x_{exp} (m); GeV", kine_select, sbsfieldscale_select, run_target.Data()), 500, -2.5, 2.5, 300, 0.0, 3.0);
		h_simc_W2_dx = new TH2D("h_simc_W2_dx", Form("Invariant Mass W^{2} vs dx - SBS%i, mag%i, %s; x_{HCal} - x_{exp} (m); GeV", kine_select, sbsfieldscale_select, run_target.Data()), 500, -2.5, 2.5, 300, 0.0, 3.0);
		h_Wrecon = new TH1D("h_Wrecon", Form("Invariant Mass W recon - SBS%i = %i%%, %s; GeV", kine_select, sbsfieldscale_select, run_target.Data()), 300, 0.0, 3.0);
		h_W2recon = new TH1D("h_W2recon", Form("Invariant Mass W^2 recon - SBS%i = %i%%, %s; GeV", kine_select, sbsfieldscale_select, run_target.Data()), 300, 0.0, 3.0);

		h_simc_dx = new TH1D("h_simc_dx",Form("dx (NO CUTS) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine_select, sbsfieldscale_select, run_target.Data()), 500, -2.5, 2.5);
		h_simc_dx_p = new TH1D("h_simc_dx_p",Form("dx - proton only -(NO CUTS) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine_select, sbsfieldscale_select, run_target.Data()), 500, -2.5, 2.5);
		h_simc_dx_n = new TH1D("h_simc_dx_n",Form("dx - neutron only -(NO CUTS) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine_select, sbsfieldscale_select, run_target.Data()), 500, -2.5, 2.5);
		h_simc_dx_cut = new TH1D("h_simc_dx_cut",Form("dx (Basic CUTS) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine_select, sbsfieldscale_select, run_target.Data()), 500, -2.5, 2.5);
		h_simc_dx_wcut = new TH1D("h_simc_dx_wcut",Form("dx (W cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine_select, sbsfieldscale_select, run_target.Data()), 500, -2.5, 2.5);
		h_simc_dx_p_wcut = new TH1D("h_simc_dx_p_wcut",Form("dx - proton only - (W cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine_select, sbsfieldscale_select, run_target.Data()), 500, -2.5, 2.5);
		h_simc_dx_n_wcut = new TH1D("h_simc_dx_n_wcut",Form("dx - neutron only - (W cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine_select, sbsfieldscale_select, run_target.Data()), 500, -2.5, 2.5);
		h_simc_dx_p_fcut = new TH1D("h_simc_dx_p_fcut",Form("dx - proton only - (F cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine_select, sbsfieldscale_select, run_target.Data()), 500, -2.5, 2.5);
		h_simc_dx_n_fcut = new TH1D("h_simc_dx_n_fcut",Form("dx - neutron only - (F cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine_select, sbsfieldscale_select, run_target.Data()), 500, -2.5, 2.5);
		h_simc_dx_p_wcut_fcut = new TH1D("h_simc_dx_p_wcut_fcut",Form("dx - proton only - (W & F cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine_select, sbsfieldscale_select, run_target.Data()), 500, -2.5, 2.5);
		h_simc_dx_n_wcut_fcut = new TH1D("h_simc_dx_n_wcut_fcut",Form("dx - neutron only - (W & F cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine_select, sbsfieldscale_select, run_target.Data()), 500, -2.5, 2.5);

		h_simc_dx_fcut = new TH1D("h_simc_dx_fcut",Form("dx (f cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine_select, sbsfieldscale_select, run_target.Data()), 500, -2.5, 2.5);
		h_simc_dx_wcut_fcut = new TH1D("h_simc_dx_wcut_fcut",Form("dx (W & f cut) - SBS%i = %i%%, %s; x_{HCal} - x_{exp} (m);", kine_select, sbsfieldscale_select, run_target.Data()), 500, -2.5, 2.5);
		h_simc_dy = new TH1D("h_simc_dy",Form("dy (NO CUTS) - SBS%i = %i%%, %s; y_{HCal} - y_{exp} (m);", kine_select, sbsfieldscale_select, run_target.Data()), 500, -2.5, 2.5);
		h_simc_dy_cut = new TH1D("h_simc_dy_cut",Form("dy (Basic Cuts) - SBS%i = %i%%, %s; y_{HCal} - y_{exp} (m);", kine_select, sbsfieldscale_select, run_target.Data()), 250, -1.25, 1.25);  
		h_simc_dy_wcut = new TH1D("h_simc_dy_wcut",Form("dy (W Cuts) - SBS%i = %i%%, %s; y_{HCal} - y_{exp} (m);", kine_select, sbsfieldscale_select, run_target.Data()), 250, -1.25, 1.25);  

		h_simc_dxdy = new TH2D("h_simc_dxdy", Form("Hadron Spot(s) on HCal (NO CUTS) - SBS%i = %i%%, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine_select, sbsfieldscale_select, run_target.Data()), 300, -1.5, 1.5, 500, -2.5, 2.5 );
		h_simc_dxdy_wcut = new TH2D("h_simc_dxdy_wcut", Form("Hadron Spot(s) on HCal (W cut) - SBS%i = %i%%, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine_select, sbsfieldscale_select, run_target.Data()), 300, -1.5, 1.5, 500, -2.5, 2.5 );
		h_simc_dxdy_cut = new TH2D("h_simc_dxdy_cut", Form("Hadron Spot(s) on HCal (Basic cuts) - SBS%i = %i%%, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine_select, sbsfieldscale_select, run_target.Data()), 300, -1.5, 1.5, 500, -2.5, 2.5 );
		h_simc_dxdy_ncut = new TH2D("h_simc_dxdy_ncut", Form("Hadron Spot(s) on HCal (n cut) - SBS%i = %i%%, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine_select, sbsfieldscale_select, run_target.Data()), 250, -1.25, 1.25, 500, -2.5, 2.5 );
		h_simc_dxdy_pcut = new TH2D("h_simc_dxdy_pcut", Form("Hadron Spot(s) on HCal (p cut) - SBS%i = %i%%, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine_select, sbsfieldscale_select, run_target.Data()), 250, -1.25, 1.25, 500, -2.5, 2.5 );
		h_simc_dxdy_fcut = new TH2D("h_simc_dxdy_fcut", Form("Hadron Spot(s) on HCal (f cut) - SBS%i = %i%%, %s;y_{HCal}-y_{expect} (m); x_{HCal}-x_{expect} (m)", kine_select, sbsfieldscale_select, run_target.Data()), 250, -1.25, 1.25, 500, -2.5, 2.5 );

		h_xy = new TH2D("h_xy",Form("HCal Hadron Spots (x, y) (NO CUTS) - SBS%i = %i%%, %s;y_{HCal} (m); x_{HCal} (m)", kine_select, sbsfieldscale_select, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);
		h_xy_cut = new TH2D("h_xy_cut", Form("HCal Hadron Spots (x, y) (BASIC CUTS) - SBS%i = %i%%, %s;y_{HCal} (m); x_{HCal} (m)", kine_select, sbsfieldscale_select, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);
		h_xy_fcut = new TH2D("h_xy_fcut", Form("HCal Hadron Spots (x, y) (Fiduc. CUTS) - SBS%i = %i%%, %s;y_{HCal} (m); x_{HCal} (m)", kine_select, sbsfieldscale_select, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);
		h_xy_cut_p = new TH2D("h_xy_cut_p", Form("HCal Hadron Spots (x, y) (p CUT) - SBS%i = %i%%, %s;y_{HCal} (m); x_{HCal} (m)", kine_select, sbsfieldscale_select, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);
		h_xy_cut_n = new TH2D("h_xy_cut_n", Form("HCal Hadron Spots (x, y) (n CUT) - SBS%i = %i%%, %s;y_{HCal} (m); x_{HCal} (m)", kine_select, sbsfieldscale_select, run_target.Data()),12,-0.9,0.9,24,-2.165,1.435);

		h_PAngleCorr_theta = new TH2D( "h_PAngCorr_theta",Form("BB theta vs HCal theta - SBS%i = %i%%, %s", kine_select, sbsfieldscale_select, run_target.Data()), 200, 0.55, 0.75, 300, 0.35, 0.65 );
		h_PAngleCorr_phi = new TH2D( "h_PAngCorr_phi",Form("BB phi vs HCal phi - SBS%i = %i%%, %s", kine_select, sbsfieldscale_select, run_target.Data()), 500, -0.4, 0.1, 500, 2.7, 3.2 );
		h_vz_cut = new TH1D("h_vz_cut",Form("BB phi vs HCal phi - SBS%i = %i%%, %s; vertex z (m);", kine_select, sbsfieldscale_select, run_target.Data()), 250,-0.125,0.125);	

		h_mc_p_weight = new TH1D("h_mc_p_weight", "mc_weight from tree for proton", 2000, 0.0, 0.000000002);
		h_mc_n_weight = new TH1D("h_mc_n_weight", "mc_weight from tree for neutron", 2000, 0.0, 0.000000001);
		h_mc_n_p_weight_Ratio = new TH1D("h_mc_n_p_weight_Ratio", "ratio of mc weights: neutron/proton", 100, 0, 1);
		h_mc_p_sigma = new TH1D("h_mc_p_sigma", "mc sigma from tree for proton", 1000, 0.0, 10.0);
		h_mc_n_sigma = new TH1D("h_mc_n_sigma", "mc sigma from tree for neutron", 1000, 0.0, 10.0);
		h_mc_np_sigma_ratio = new TH1D("h_mc_np_sigma_ratio", "Ratio of sigmas for n/p from mc tree", 100, 0.0, 1.0);
		h_mc_p_sigma_weighted = new TH1D("h_mc_p_sigma_weighted", "mc sigma from tree for proton - weighted", 1000, 0.0, 10.0);
		h_mc_n_sigma_weighted = new TH1D("h_mc_n_sigma_weighted", "mc sigma from tree for neutron- weighted", 1000, 0.0, 10.0);
		h_mc_np_sigma_ratio_weighted = new TH1D("h_mc_np_sigma_ratio_weighted", "Ratio of sigmas for n/p from mc tree- weighted", 100, 0.0, 1.0);

		if( use_hist_file ){
			h_mc_p_final_weight = new TH1D("h_mc_p_final_weight", "mc_Final_weight calculated for proton", 100000, 0.0, 0.1);
			h_mc_n_final_weight = new TH1D("h_mc_n_final_weight", "mc_Final_weight calculated for neutron", 100000, 0.0, 0.1);
			h_mc_n_p_final_weight_Ratio = new TH1D("h_mc_n_p_final_weight_Ratio", "ratio of mc final weights: neutron/proton", 250, 0, 5);			
		}
		//Various SIMC variables
		h_mc_simc_Q2 = new TH1D("h_mc_simc_Q2", Form("simc calculated Q^{2} momentum transfer - SBS%i Mag%i %s; Momentum-Transfer Q^{2} (GeV); Entries (N)", kine_select, sbsfieldscale_select, run_target.Data()), 1000, 0, 10);
		h_mc_simc_epsilon = new TH1D("h_mc_simc_epsilon", Form("simc calcualted virtual photon polarization, #epsilon - SBS%i Mag%i %s; Virtual Photon Polarization, #epsilon (-); Entries (N)", kine_select, sbsfieldscale_select, run_target.Data()), 150, 0.0, 1.5);
		h_mc_simc_Ebeam = new TH1D("h_mc_simc_Ebeam", Form("simc calcualted beam energy, E_{Beam} - SBS%i Mag%i %s; Beam Energy (GeV); Entries (N)", kine_select, sbsfieldscale_select, run_target.Data()), 1000, 0, 10.0);
		h_mc_simc_p_n = new TH1D("h_mc_simc_p_n", Form("simc calculated scattered nucleon momentum, p_{N} - SBS%i Mag%i %s; Nucleon Momentum (GeV); Entries (N)", kine_select, sbsfieldscale_select, run_target.Data()), 1000, 0., 10.0);
		h_mc_simc_theta_n = new TH1D("h_mc_simc_theta_n", Form("simc calculated scattered nucleon polar angle (DEG), #theta_{N} - SBS%i Mag%i %s; Nucleon Polar Angle (Deg); Entries (N)", kine_select, sbsfieldscale_select, run_target.Data()), 900, 0.0, 90);
		h_mc_simc_p_e = new TH1D("h_mc_simc_p_e", Form("simc calcualted scattered electron momentum - SBS%i Mag%i %s; Electron Momentum (GeV); Entries (N)", kine_select, sbsfieldscale_select, run_target.Data()), 1000, 0, 10.0);
		h_mc_simc_theta_e = new TH1D("h_MC_simc_theta_e", Form("simc calcualted scattered electron polar angle (DEG) - SBS%i Mag%i %s; Electron Polar Angle (Deg); Entries (N)", kine_select, sbsfieldscale_select, run_target.Data()), 900, 0.0, 90);

		//with cut on W
		h_mc_simc_Q2_wcut = new TH1D("h_mc_simc_Q2_wcut", Form("simc calculated Q^{2} momentum transfer (W cut) - SBS%i Mag%i %s; Momentum-Transfer Q^{2} (GeV); Entries (N)", kine_select, sbsfieldscale_select, run_target.Data()), 1000, 0, 10);
		h_mc_simc_epsilon_wcut = new TH1D("h_mc_simc_epsilon_wcut", Form("simc calcualted virtual photon polarization, #epsilon (W cut) - SBS%i Mag%i %s; Virtual Photon Polarization, #epsilon (-); Entries (N)", kine_select, sbsfieldscale_select, run_target.Data()), 150, 0.0, 1.5);
		h_mc_simc_Ebeam_wcut = new TH1D("h_mc_simc_Ebeam_wcut", Form("simc calcualted beam energy, E_{Beam} (W cut) - SBS%i Mag%i %s; Beam Energy (GeV); Entries (N)", kine_select, sbsfieldscale_select, run_target.Data()), 1000, 0, 10.0);
		h_mc_simc_p_n_wcut = new TH1D("h_mc_simc_p_n_wcut", Form("simc calculated scattered nucleon momentum, p_{N} (W cut) - SBS%i Mag%i %s; Nucleon Momentum (GeV); Entries (N)", kine_select, sbsfieldscale_select, run_target.Data()), 1000, 0., 10.0);
		h_mc_simc_theta_n_wcut = new TH1D("h_mc_simc_theta_n_wcut", Form("simc calcualted scattered nucleon polar angle (DEG), #theta_{N} (W cut) - SBS%i Mag%i %s; Nucleon Polar Angle (Deg); Entries (N)", kine_select, sbsfieldscale_select, run_target.Data()), 900, 0.0, 90);
		h_mc_simc_p_e_wcut = new TH1D("h_mc_simc_p_e_wcut", Form("simc calcualted scattered electron momentum (W cut) - SBS%i Mag%i %s; Electron Momentum (GeV); Entries (N)", kine_select, sbsfieldscale_select, run_target.Data()), 1000, 0, 10.0);
		h_mc_simc_theta_e_wcut = new TH1D("h_MC_simc_theta_e_wcut", Form("simc calcualted scattered electron polar angle (DEG) (W cut) - SBS%i Mag%i %s; Electron Polar Angle (Deg); Entries (N)", kine_select, sbsfieldscale_select, run_target.Data()), 900, 0.0, 90);


		cout << "finished. " << endl << endl;

		cout << "--------------------------------------" << endl;
		cout << "Adding files to TChain from: " << rootfile_dir.Data() << endl;

		if( kine_select == 4 ){

			SBS4_infile_basename = Form("replayed_digitized_jb_FARM_SIMC_gmn_SBS4_LD2_NUCLEON_mag30port%s_175uA_elas_250k_%s_job*.root", portField_select.Data(), dateInfile.Data() );

			// if( run_target == "LD2"){
			// 	infile = Form("%s/replayed_jb_gmn_SBS4_LD2_mag70_175uA_elas_100k_job*", rootfile_dir.Data());			
			// }
			// else{
			// 	infile = Form("%s/replayed_jb_gmn_SBS4_LH2_mag0_350uA_pGun_100k_job*", rootfile_dir.Data());
			// }

			// infile = Form("%s/replayed_jb_gmn_SBS%i_%s_mag%imod1_175uA_elas_100k_job*", rootfile_dir.Data(), kine_select, run_target.Data(), sbsfieldscale_select);
			//Proton file:
			 
			if( hadron_select == "proton" || hadron_select == "both" ){
				// proton_infile = "replayed_jb_SIMC_gmn_SBS8_LD2_proton_mag70mod0312_500uA_elas_100k_31_07_2023_job*";
				// proton_infile = "replayed_jb_SIMC_gmn_SBS8_LD2_proton_mag70mod0312_500uA_elas_500k_02_08_2023_job*";

				// proton_infile = Form("replayed_digitized_jb_FARM_SIMC_gmn_SBS4_LD2_proton_mag30port%s_175uA_elas_%s_10_12_2023_job*.root", portField_select.Data(), numEvents_simFile_str.Data());
				proton_infile = SBS4_infile_basename;
				proton_infile.ReplaceAll("_NUCLEON_", "_proton_");

				// proton_infile = "replayed_jb_SIMC_gmn_SBS8_LD2_proton_mag70mod0312_500uA_elas_500k_02_08_2023_job10.root*";
				proton_infile_basename = Form("%s", proton_infile.Data() );
				proton_infile_basename.ReplaceAll("*.root", "");
				cout << "Adding proton infile names to vector. " << endl;
				cout << "Adding proton infiles with basename: "<< proton_infile_basename.Data() << endl;

				FillVectorByFilesAndDirWithPattern( rootfile_dir.Data(), proton_infile_basename.Data(), proton_infile_vec );

				cout << "Sorting proton infile vector... " << endl;
				sort( proton_infile_vec.begin(), proton_infile_vec.end(), CompareStringsWithNumbers );

				proton_infile_cnt = proton_infile_vec.size();
				nucleon_infile_cnt[1] = proton_infile_cnt;

				cout << "Number of files found for proton simulation: " << proton_infile_cnt << endl;
			}

			//Neutron file:
			if( hadron_select == "neutron" || hadron_select == "both" ){
				// neutron_infile = "replayed_jb_SIMC_gmn_SBS8_LD2_neutron_mag70mod0312_500uA_elas_100k_31_07_2023_job*";
				// neutron_infile = "replayed_jb_SIMC_gmn_SBS8_LD2_neutron_mag70mod0312_500uA_elas_500k_02_08_2023_job*";
				
				// neutron_infile = Form("replayed_digitized_jb_FARM_SIMC_gmn_SBS4_LD2_neutron_mag30port%s_175uA_elas_%s_10_12_2023_job*.root", portField_select.Data(), numEvents_simFile_str.Data());
				neutron_infile = SBS4_infile_basename;
				neutron_infile.ReplaceAll("_NUCLEON_", "_neutron_");

				// neutron_infile = "replayed_jb_SIMC_gmn_SBS8_LD2_neutron_mag70mod0312_500uA_elas_500k_02_08_2023_job10.root";
				neutron_infile_basename = Form("%s", neutron_infile.Data() );
				neutron_infile_basename.ReplaceAll("*.root", "");
				cout << "Adding neutron infile names to vector. " << endl;
				cout << "Adding neutron infiles with basename: "<< neutron_infile_basename.Data() << endl;

				FillVectorByFilesAndDirWithPattern( rootfile_dir.Data(), neutron_infile_basename.Data(), neutron_infile_vec );

				cout << "Sorting neutron infile vector... " << endl;
				sort( neutron_infile_vec.begin(), neutron_infile_vec.end(), CompareStringsWithNumbers );

				neutron_infile_cnt = neutron_infile_vec.size();
				nucleon_infile_cnt[0] = neutron_infile_cnt;

				cout << "Number of files found for neutron simulation: " << neutron_infile_cnt << endl;

			}
		}

		if( kine_select == 8 ){
			
			SBS8_infile_basename = Form("replayed_digitized_jb_FARM_SIMC_gmn_SBS8_LD2_NUCLEON_mag70port%s_500uA_elas_250k_%s_job*.root", portField_select.Data(), dateInfile.Data());


			if( hadron_select == "proton" || hadron_select == "both" ){
				// proton_infile = "replayed_jb_SIMC_gmn_SBS8_LD2_proton_mag70mod0312_500uA_elas_100k_31_07_2023_job*";
				// proton_infile = "replayed_digitized_jb_FARM_SIMC_gmn_SBS8_LD2_proton_mag70port0301_500uA_elas_250k_24_09_2023_job*.root";
				// proton_infile = "replayed_digitized_jb_FARM_SIMC_gmn_SBS8_LD2_proton_mag70port0301_500uA_elas_250k_08_10_2023_job*.root";
				// proton_infile = "replayed_digitized_jb_FARM_SIMC_gmn_SBS8_LD2_proton_mag70port0310_500uA_elas_250k_16_10_2023_job*.root";
				//***
				// proton_infile = "replayed_digitized_jb_FARM_SIMC_gmn_SBS8_LD2_proton_mag70port0310_500uA_elas_250k_03_12_2023_job*.root";
				// proton_infile = Form("replayed_digitized_jb_FARM_SIMC_gmn_SBS8_LD2_proton_mag70port%s_500uA_elas_250k_04_11_2023_job*.root", portField_select.Data());

				proton_infile = SBS8_infile_basename;
				proton_infile.ReplaceAll("_NUCLEON_", "_proton_");

				// proton_infile = "replayed_digitized_jb_FARM_SIMC_gmn_SBS8_LD2_proton_mag70port0800_500uA_elas_100k_04_11_2023_job*.root";
				proton_infile_basename = Form("%s", proton_infile.Data() );
				proton_infile_basename.ReplaceAll("*.root", "");
				cout << "Adding proton infile names to vector. " << endl;
				cout << "Adding proton infiles with basename: "<< proton_infile_basename.Data() << endl;

				FillVectorByFilesAndDirWithPattern( rootfile_dir.Data(), proton_infile_basename.Data(), proton_infile_vec );

				cout << "Sorting proton infile vector... " << endl;
				sort( proton_infile_vec.begin(), proton_infile_vec.end(), CompareStringsWithNumbers );

				proton_infile_cnt = proton_infile_vec.size();
				nucleon_infile_cnt[1] = proton_infile_cnt;

				cout << "Number of files found for proton simulation: " << proton_infile_cnt << endl;

			}

			//Neutron file:
			if( hadron_select == "neutron" || hadron_select == "both" ){
				// neutron_infile = "replayed_jb_SIMC_gmn_SBS8_LD2_neutron_mag70mod0312_500uA_elas_100k_31_07_2023_job*";
				// neutron_infile = "replayed_digitized_jb_FARM_SIMC_gmn_SBS8_LD2_neutron_mag70port0301_500uA_elas_250k_24_09_2023_job*.root";
				// neutron_infile = "replayed_digitized_jb_FARM_SIMC_gmn_SBS8_LD2_neutron_mag70port0301_500uA_elas_250k_08_10_2023_job*.root";
				// neutron_infile = "replayed_digitized_jb_FARM_SIMC_gmn_SBS8_LD2_neutron_mag70port0310_500uA_elas_250k_16_10_2023_job*.root";
				//***
				// neutron_infile = "replayed_digitized_jb_FARM_SIMC_gmn_SBS8_LD2_neutron_mag70port0310_500uA_elas_250k_03_12_2023_job*.root";
				// neutron_infile = Form("replayed_digitized_jb_FARM_SIMC_gmn_SBS8_LD2_neutron_mag70port%s_500uA_elas_250k_04_11_2023_job*.root", portField_select.Data());

				neutron_infile = SBS8_infile_basename;
				neutron_infile.ReplaceAll("_NUCLEON_", "_neutron_");

				neutron_infile_basename = Form("%s", neutron_infile.Data() );
				neutron_infile_basename.ReplaceAll("*.root", "");
				cout << "Adding neutron infile names to vector. " << endl;
				cout << "Adding neutron infiles with basename: "<< neutron_infile_basename.Data() << endl;

				FillVectorByFilesAndDirWithPattern( rootfile_dir.Data(), neutron_infile_basename.Data(), neutron_infile_vec );

				cout << "Sorting neutron infile vector... " << endl;
				sort( neutron_infile_vec.begin(), neutron_infile_vec.end(), CompareStringsWithNumbers );

				neutron_infile_cnt = neutron_infile_vec.size();
				nucleon_infile_cnt[0] = neutron_infile_cnt;

				cout << "Number of files found for neutron simulation: " << neutron_infile_cnt << endl;

			}
		
		}

		if( kine_select == 9 ){
			
			SBS9_infile_basename = Form("replayed_digitized_jb_FARM_SIMC_gmn_SBS9_LD2_NUCLEON_mag70port%s_1200uA_elas_250k_%s_job*.root", portField_select.Data(), dateInfile.Data());

			if( hadron_select == "proton" || hadron_select == "both" ){
				cout << "------------------------------------------" << endl;
				// proton_infile = "replayed_jb_SIMC_gmn_SBS8_LD2_proton_mag70mod0312_500uA_elas_100k_31_07_2023_job*";
				// proton_infile = "replayed_digitized_jb_FARM_SIMC_gmn_SBS9_LD2_proton_mag70port0314_1200uA_elas_250k_19_09_2023_job*.root";
				//
				// proton_infile = "replayed_digitized_jb_FARM_SIMC_gmn_SBS9_LD2_proton_mag70port0314_1200uA_elas_250k_08_10_2023_job*.root";
				// proton_infile = Form("replayed_digitized_jb_FARM_SIMC_gmn_SBS9_LD2_proton_mag70port%s_1200uA_elas_250k_08_10_2023_job*.root", portField_select.Data());

				proton_infile = SBS9_infile_basename;
				proton_infile.ReplaceAll("_NUCLEON_", "_proton_");

				proton_infile_basename = Form("%s", proton_infile.Data() );
				proton_infile_basename.ReplaceAll("*.root", "");
				cout << "Adding proton infile names to vector. " << endl;
				cout << "Adding proton infiles with basename: "<< proton_infile_basename.Data() << endl;

				FillVectorByFilesAndDirWithPattern( rootfile_dir.Data(), proton_infile_basename.Data(), proton_infile_vec );

				cout << "Sorting proton infile vector... " << endl;
				sort( proton_infile_vec.begin(), proton_infile_vec.end(), CompareStringsWithNumbers );

				proton_infile_cnt = proton_infile_vec.size();
				nucleon_infile_cnt[1] = proton_infile_cnt;

				cout << "Number of files found for proton simulation: " << proton_infile_cnt << endl;
				cout << "------------------------------------------" << endl;
			}
			//Neutron file:
			if( hadron_select == "neutron" || hadron_select == "both" ){
				cout << "------------------------------------------" << endl;
				// neutron_infile = "replayed_jb_SIMC_gmn_SBS8_LD2_neutron_mag70mod0312_500uA_elas_100k_31_07_2023_job*";
				// neutron_infile = "replayed_digitized_jb_FARM_SIMC_gmn_SBS9_LD2_neutron_mag70port0314_1200uA_elas_250k_19_09_2023_job*.root";
				///
				// neutron_infile = "replayed_digitized_jb_FARM_SIMC_gmn_SBS9_LD2_neutron_mag70port0314_1200uA_elas_250k_08_10_2023_job*.root";
				// neutron_infile = Form("replayed_digitized_jb_FARM_SIMC_gmn_SBS9_LD2_neutron_mag70port%s_1200uA_elas_250k_08_10_2023_job*.root", portField_select.Data());

				neutron_infile = SBS9_infile_basename;
				neutron_infile.ReplaceAll("_NUCLEON_", "_neutron_");

				neutron_infile_basename = Form("%s", neutron_infile.Data() );
				neutron_infile_basename.ReplaceAll("*.root", "");
				cout << "Adding neutron infile names to vector. " << endl;
				cout << "Adding neutron infiles with basename: "<< neutron_infile_basename.Data() << endl;

				FillVectorByFilesAndDirWithPattern( rootfile_dir.Data(), neutron_infile_basename.Data(), neutron_infile_vec );

				cout << "Sorting neutron infile vector... " << endl;
				sort( neutron_infile_vec.begin(), neutron_infile_vec.end(), CompareStringsWithNumbers );

				neutron_infile_cnt = neutron_infile_vec.size();
				nucleon_infile_cnt[0] = neutron_infile_cnt;

				cout << "Number of files found for proton simulation: " << neutron_infile_cnt << endl;
				cout << "------------------------------------------" << endl;
			}
		}

		if( use_hist_file ){
			cout << "------------------------------------------" << endl;
			cout << "Gathering .hist files corresponding to the proton infiles...." << endl;
			if( proton_infile_vec.size() == 0 || neutron_infile_vec.size() == 0 ){
				cout << "Empty infile vector...." << endl;
				cout << "Sleeping to allow for a ctrl-c catch..." << endl;
				sleep(20);
				cout << "Exiting...." << endl;
				exit(1);				
			}
			cout << "Number of proton files: " << proton_infile_vec.size() << ".....";
			vector<int> missing_pfile_index = {};
			vector<int> missing_nfile_index = {};

			//We will store the job number of each input file to check against the .hist file
			//We will store the index/position of the string _job in each nucleon filename:
			int p_file_job_index = proton_infile_vec[0].Index("_job");
			int n_file_job_index = neutron_infile_vec[0].Index("_job");

			//Now let us a build a corresponding vector of SIMC hist files for each input:

			for( size_t pFile = 0; pFile < proton_infile_vec.size(); pFile++ ){
				cout << pFile << " ";
				TString temp_filename = proton_infile_vec[pFile];

				temp_filename.ReplaceAll("replayed_digitized_", "");
				temp_filename.ReplaceAll(Form("/simc/SBS%i/", kine_select), "/hist/");
				temp_filename.ReplaceAll(".root", ".hist");
				// temp_filename.Prepend( Form( "%s/", histfile_dir.Data() ));
				proton_histfile_vec.push_back(temp_filename);

				//Let us also store the job number from this file so that we can match things up later...
				//Extract the substring after "_job"
				TString p_jobNumberStr = "";
				p_jobNumberStr = proton_infile_vec[pFile]( p_file_job_index + 4, proton_infile_vec[pFile].Length() - p_file_job_index - 9); //9 accounts for lenght of "_job" and ".root"
				
				int p_jobNumber = -1;
				p_jobNumber = p_jobNumberStr.Atoi();

				proton_infile_jobNum_vec.push_back( p_jobNumber );

			//SET UP SOME PROVISIONS IN CASE THE HIST FILE IS MISSING. 
				double hist_test_var = searchSimcHistFile( "luminosity", proton_histfile_vec[pFile] );
				//If searchSimcHistFile can't access the file it returns the value -99
				if( hist_test_var == -99 ){
					missing_pfile_index.push_back(pFile);
					cout << "ERRRRRROOOOOORRRRRR ------ .hist file not found for pFile: " << pFile << endl;
					cout << "Sleeping to allow for a ctrl-c catch..." << endl;
					sleep(20);
					cout << "Exiting...." << endl;
					exit(1);
				}

				mc_p_Ntried_vec.push_back( int(searchSimcHistFile( "Ntried", proton_histfile_vec[pFile] ) ) );
				mc_p_luminosity_vec.push_back( searchSimcHistFile( "luminosity", proton_histfile_vec[pFile] ) );
				mc_p_genvol_vec.push_back( searchSimcHistFile( "genvol", proton_histfile_vec[pFile] ) );
				mc_p_Final_Weight_vec.push_back( mc_p_luminosity_vec[pFile]*mc_p_genvol_vec[pFile]*(1.0/(1.0*mc_p_Ntried_vec[pFile])) );
				temp_filename = "";
			}


			cout << endl;
			cout << "done." << endl;
			cout << "- - - - - - - - - - - - - - - - - - - - -" << endl;
			cout << "Gathering .hist files corresponding to the neutron infiles...." << endl;
			cout << "Number of neutron files: " << neutron_infile_vec.size() << endl;
			for( size_t nFile = 0; nFile < neutron_infile_vec.size(); nFile++ ){
				cout << nFile << " ";
				TString temp_filename = neutron_infile_vec[nFile];
				temp_filename.ReplaceAll("replayed_digitized_", "");
				temp_filename.ReplaceAll(Form("/simc/SBS%i/", kine_select), "/hist/");
				temp_filename.ReplaceAll(".root", ".hist");
				// temp_filename.Prepend( Form( "%s/", histfile_dir.Data() ));
				neutron_histfile_vec.push_back(temp_filename);

				//Let us also store the job number from this file so that we can match things up later...
				//Extract the substring after "_job"
				TString n_jobNumberStr = "";
				n_jobNumberStr = neutron_infile_vec[nFile]( n_file_job_index + 4, neutron_infile_vec[nFile].Length() - n_file_job_index - 9); //9 accounts for lenght of "_job" and ".root"
				
				int n_jobNumber = -1;
				n_jobNumber = n_jobNumberStr.Atoi();

				neutron_infile_jobNum_vec.push_back( n_jobNumber );

			//SET UP SOME PROVISIONS IN CASE THE HIST FILE IS MISSING. 
				double hist_test_var = searchSimcHistFile( "luminosity", neutron_histfile_vec[nFile] );
				//If searchSimcHistFile can't access the file it returns the value -99
				if( hist_test_var == -99 ){
					missing_nfile_index.push_back(nFile);
					cout << "ERRRRRROOOOOORRRRRR ------ .hist file not found for nFile: " << nFile << endl;
					cout << "Sleeping to allow for a ctrl-c catch..." << endl;
					sleep(20);
					cout << "Exiting...." << endl;
					exit(1);
				}

				mc_n_Ntried_vec.push_back( int(searchSimcHistFile( "Ntried", neutron_histfile_vec[nFile] ) ) );
				mc_n_luminosity_vec.push_back( searchSimcHistFile( "luminosity", neutron_histfile_vec[nFile] ) );
				mc_n_genvol_vec.push_back( searchSimcHistFile( "genvol", neutron_histfile_vec[nFile] ) );
				mc_n_Final_Weight_vec.push_back( mc_n_luminosity_vec[nFile]*mc_n_genvol_vec[nFile]*(1.0/(1.0*mc_n_Ntried_vec[nFile])) );
				temp_filename = "";
			}
			cout << endl;
			cout << "------------------------------------------" << endl;	
		}		

		if( !match_file_cnts ){
			cout << "-------------------NOT MATCHING FILE COUNTS---------------------" << endl;
			int p_max_files = 0;
			int n_max_files = 0;

			if( manually_select_file_cnt ){
				p_max_files = selected_file_cnt;
				n_max_files = selected_file_cnt;
				cout << "Max file cnt selected as: " << selected_file_cnt << endl;
			}

			if( !manually_select_file_cnt ){
				cout << "Adding all proton and neutron simulation files to TChain... " << endl;
				p_max_files = proton_infile_vec.size();
				n_max_files = neutron_infile_vec.size();
				cout << "Max file cnt matches infile vector sizes: p = " << proton_infile_vec.size() << ", n = " << neutron_infile_vec.size() << endl;
			}

			for( int p_file = 0; p_file < p_max_files; p_file++ ){
				TC->Add(proton_infile_vec[p_file].Data());
			}
			for( int n_file = 0; n_file < n_max_files; n_file++ ){
				TC->Add(neutron_infile_vec[n_file].Data());
			}
			// TC->Add(Form("%s/%s", rootfile_dir.Data(), neutron_infile.Data()));
			// TC->Add(Form("%s/%s", rootfile_dir.Data(), proton_infile.Data()));

			if( kine_select == 9 ){
				for( int p = 0; p < 200; p++ ){
					TC->Add(proton_infile_vec[p].Data());
				}
				for( int n = 0; n < 233; n++ ){
					TC->Add(neutron_infile_vec[n].Data());
				}				
			}
			// if( true ){
			// 	for( size_t p = 0; p < proton_infile_vec.size(); p++ ){
			// 		TC->Add(proton_infile_vec[p].Data());
			// 		TC->Add(neutron_infile_vec[p].Data());
			// 	}
			// }

		}

		if( match_file_cnts ){

			if( !manually_select_file_cnt ){
				cout << "Adding matching number of proton and neutron simulation files to TChain... " << endl;
				if( neutron_infile_cnt < proton_infile_cnt ){ 
					nucleon_with_min_file_cnt = 0; 
					nucleon_with_min_file_cnt_string = Form("Neutron: %i (Proton has %i)", neutron_infile_cnt, proton_infile_cnt);
				}
				if( neutron_infile_cnt > proton_infile_cnt ){ 
					nucleon_with_min_file_cnt = 1; 
					nucleon_with_min_file_cnt_string = Form("Proton: %i (Neutron has %i)", proton_infile_cnt, neutron_infile_cnt);
				}
				if( neutron_infile_cnt == proton_infile_cnt ){ 
					nucleon_with_min_file_cnt = 0; 
					nucleon_with_min_file_cnt_string = Form("Input file counts are equal: %i and %i", neutron_infile_cnt, proton_infile_cnt );
				} //Just deafult to neutron
	//Now that we know which nucleon has the lower number input files we can only add as many files as that count:
				cout << "Nucleon with less file cnts: " << nucleon_with_min_file_cnt_string.Data() << endl;
				cout << "Number of files for each nucleon to add: " << nucleon_infile_cnt[nucleon_with_min_file_cnt] << endl;


				cout << "Adding files to TChain...." << endl;
				for( size_t infile = 0; infile < nucleon_infile_cnt[nucleon_with_min_file_cnt]; infile++ ){
					TC->Add( neutron_infile_vec[infile].Data() );
					TC->Add( proton_infile_vec[infile].Data() );
				}
			}

			if( manually_select_file_cnt ){
				cout << "Adding a fixed/manually selected file count: " << selected_file_cnt << endl;
				cout << "Existing file counts: p: " << proton_infile_cnt << ", n: " << neutron_infile_cnt << endl;
				cout << endl;
				cout << "Infile: ";
				for( int infile = 0; infile < selected_file_cnt; infile++ ){
					cout << infile << " ";
					TC->Add( neutron_infile_vec[infile].Data() );
					TC->Add( proton_infile_vec[infile].Data() );				
				}
				cout << endl;
			}

			cout << "Finishing adding files to TChain..." << endl;
		}


		// else{
		// 	infile = Form("%s/replayed_jb_gmn_SBS%i_%s_mag%imod3_%suA_elas_100k_job*", rootfile_dir.Data(), kine_select, run_target.Data(), sbsfieldscale_select, I_beam_str.Data());
		// 	cout << "Adding file: " << endl;
		// 	cout << infile.Data() << endl << endl;
		// 	TC->Add(infile.Data());		
		// }
		if( use_particle_gun ){
			cout << endl << "Using particle gun.... " << endl;
			cout << "Adding file: " << endl;
			infile = Form("%s/replayed_jb_gmn_SBS4_LD2_mag0_100uA_pGun_Px_100k_v2_job0.root", rootfile_dir.Data());
			cout << infile.Data() << endl << endl;
		}
		// TC->Add(Form("%s/replayed_simc_sbs4_sbs50p_89T_elas_job*", rootfile_dir.Data()));
		// TC->Add(Form("%s/replayed_gmn_sbs4_ld2_30p_job_*", rootfile_dir.Data()));
		// TC->Add(Form("%s/replayed_jb_gmn_SBS8_Mag70_500k.root", rootfile_dir.Data()));

		cout << "--------------------------------------" << endl;
		cout << "Setting up branches... ";
		TC->SetBranchStatus( "*", 0 );

		//MC values for normalization
		TC->SetBranchStatus( "MC.mc_omega", 1 );
		TC->SetBranchStatus( "MC.simc_sigma", 1 );
		TC->SetBranchStatus( "MC.simc_fnucl", 1);
		TC->SetBranchStatus( "MC.simc_Weight", 1);
		TC->SetBranchStatus( "MC.simc_Q2", 1);
		TC->SetBranchStatus( "MC.simc_nu", 1);
		TC->SetBranchStatus( "MC.simc_epsilon", 1);
		TC->SetBranchStatus( "MC.simc_Ebeam", 1);
		TC->SetBranchStatus( "MC.simc_p_n", 1);
		TC->SetBranchStatus( "MC.simc_theta_n", 1);
		TC->SetBranchStatus( "MC.simc_p_e", 1);
		TC->SetBranchStatus( "MC.simc_theta_e", 1);


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
		TC->SetBranchAddress( "MC.simc_sigma", &mc_sigma );
		TC->SetBranchAddress( "MC.simc_fnucl", &mc_fnucl );
		TC->SetBranchAddress( "MC.simc_Weight", &mc_Weight );
		TC->SetBranchAddress( "MC.simc_Q2", &mc_simc_Q2);
		TC->SetBranchAddress( "MC.simc_nu", &mc_simc_nu);
		TC->SetBranchAddress( "MC.simc_epsilon", &mc_simc_epsilon);
		TC->SetBranchAddress( "MC.simc_Ebeam", &mc_simc_Ebeam);
		TC->SetBranchAddress( "MC.simc_p_n", &mc_simc_p_n);
		TC->SetBranchAddress( "MC.simc_theta_n", &mc_simc_theta_n);
		TC->SetBranchAddress( "MC.simc_p_e", &mc_simc_p_e);
		TC->SetBranchAddress( "MC.simc_theta_e", &mc_simc_theta_e);

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
			Form("bb.tr.p[0]>%f", 1.00*lookup_parsed_cut(run_target.Data(), kine_select, sbsfieldscale_select, "SH_PS_mean") ),
			// Form("((abs(((bb.sh.e+bb.ps.e)/(bb.tr.p[0]))))-%f)<%f", lookup_simc_cut(kine_select, sbsfieldscale_select, run_target, "Ep"), 3.0*lookup_simc_cut(kine_select, sbsfieldscale_select, run_target, "Ep_sigma")),
			// Form("(bb.sh.e+bb.ps.e)>%f", lookup_simc_cut(kine_select, sbsfieldscale_select, run_target, "SH_PS_mean") - 3.0*lookup_simc_cut(kine_select, sbsfieldscale_select, run_target, "SH_PS_sigma")),
			Form("bb.ps.e>%f", 0.20),
			// Form("sbs.hcal.e>%f", 0.01)
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

		///TO Calculate the appropriate MC weight using the values from the .hist file we need to know which file the current event comes from. 
			//After we know what file that event comes from we can look up the corresponding hist file values. 

			if( use_hist_file ){
				currentChainFile = TC->GetCurrentFile();
				currentChainFileName = currentChainFile->GetName();

				//We should compare job numbers now with the index of the infile_vec so that we pull the right value
				current_file_jobNum_index = -1;
				current_file_jobNum_index = currentChainFileName.Index("_job");

				TString current_file_jobNumStr = "";
				current_file_jobNumStr = currentChainFileName(current_file_jobNum_index + 4, currentChainFileName.Length() - current_file_jobNum_index - 9);

				int current_file_jobNum = current_file_jobNumStr.Atoi();

				//From here we can either use mc_fnucl to determine if this is from a proton or neutron
				//or we can just search through each proton/neutron infile vector for a metching value. 
				if( is_n ){
					for( size_t nFileIndex = 0; nFileIndex < neutron_infile_vec.size(); nFileIndex++ ){
						if( neutron_infile_vec[nFileIndex] == currentChainFileName ){
							// cout << "Checking job numbers between neutron hist and input files: " << endl;
							// cout << "neutron infile: " << current_file_jobNum << ", hist: " << neutron_infile_jobNum_vec[nFileIndex] << endl;
							
							if( neutron_infile_jobNum_vec[nFileIndex] == current_file_jobNum ){
								current_file_index = nFileIndex;
								current_nfile_index = nFileIndex;								
							}
							else{
								cout << "ERROR: MISMATCH in current neutron file jobNum and current histFile jobNum. nFile index: " << nFileIndex << ", current filename: " << endl;
								cout << currentChainFileName.Data() << endl;
							}
							break;
						}
					}
					if(current_nfile_index == -1){
						cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
						cout << "Current neutron file index is -1..... no matching hist file was found!!! Check Script!!" << endl;
						cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
						exit(1);
					}
				}
				if( is_p ){
					for( size_t pFileIndex = 0; pFileIndex < proton_infile_vec.size(); pFileIndex++ ){
						if( proton_infile_vec[pFileIndex] == currentChainFileName ){
							// cout << "Checking job numbers between proton hist and input files: " << endl;
							// cout << "proton infile: " << current_file_jobNum << ", hist: " << proton_infile_jobNum_vec[pFileIndex] << endl;

							if( proton_infile_jobNum_vec[pFileIndex] == current_file_jobNum ){
								current_file_index = pFileIndex;
								current_pfile_index = pFileIndex;
							}
							else{
								cout << "ERROR: MISMATCH in current proton file jobNum and current histFile jobNum. pFile index: " << pFileIndex << ", current filename: " << endl;
								cout << currentChainFileName.Data() << endl;
							}

							break;
						}
					}
					if(current_pfile_index == -1){
						cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
						cout << "Current proton file index is -1..... no matching hist file was found!!! Check Script!!" << endl;
						cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
						exit(1);
					}
				}				
			}
			if( calc_pn_weight ){
				// pn_weight = ((mc_omega*mc_sigma)*luminosity)/ngen_total;
				pn_weight = mc_Weight;

				if( is_n ){				
					h_mc_n_weight->Fill( mc_Weight );
					// if( mc_Weight>0 ){
					// 	cout << "mc_n_weight = " << mc_Weight << endl;
					// }
				}
				if( is_p ){
					h_mc_p_weight->Fill( mc_Weight );
					// if( mc_Weight>0 ){
					// 	cout << "mc_p_weight = " << mc_Weight << endl;
					// }
				}
				// cout << "pn weight: " << pn_weight << endl;

		//IF using the .hist file then we need to references the values pulled from each .hist file
				if( use_hist_file ){
					//mc_Final_Weight = mc_Weight*mc_luminosity*mc_genvol*(1/mc_Ntried)
					//we already calculate most of the final weight in the mc_X_final_weight_vec vector, everything but mc_Weight
					if( is_n ){
						pn_weight = mc_Weight*mc_n_Final_Weight_vec[current_nfile_index];
						// cout << "n final weight: " << pn_weight << endl;
						h_mc_n_final_weight->Fill( pn_weight ); 
					}
					if( is_p ){
						pn_weight = mc_Weight*mc_p_Final_Weight_vec[current_pfile_index];
						// cout <<"pfile index: " << current_pfile_index << endl;
						// cout << "mcWeight: " << mc_Weight << endl;
						// cout << "mc_p_final_weightvec-val: " << mc_p_Final_Weight_vec[current_pfile_index]<< endl;
						// cout << "p final weight: " << pn_weight << endl;
						h_mc_p_final_weight->Fill( pn_weight ); 
					}
				}
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
					// cout << "average time for 5 = " << average_time << endl;
					time_remaining = average_time*( 1.0 - double(nevent)/double(Nevents));
					cout << "Evt: " << nevent <<"/" << Nevents << Form("(%.0f/100%%)", 100.0*double(1.0*nevent/Nevents)) << ". Elastic yield = " << elastic_yield << ". Time left: " << time_remaining << endl;
				}
				StopWatch->Reset();
				StopWatch->Continue();
				watch_cnt++;
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
		//SIMC branch variables
			h_mc_simc_Q2->Fill( mc_simc_Q2 );
			h_mc_simc_epsilon->Fill( mc_simc_epsilon );
			h_mc_simc_Ebeam->Fill( mc_simc_Ebeam );
			h_mc_simc_p_n->Fill( mc_simc_p_n );
			h_mc_simc_theta_n->Fill( TMath::RadToDeg()*mc_simc_theta_n );
			h_mc_simc_p_e->Fill( mc_simc_p_e );
			h_mc_simc_theta_e->Fill( TMath::RadToDeg()*mc_simc_theta_e );


	//ENERGY CALCULATIONS -- CALIBRATIONS
		    Double_t Ep = (bb_ps_e + bb_sh_e)/(bb_tr_p[0]);
			Double_t PS = bb_ps_e;
			Double_t SHPS = bb_ps_e + bb_sh_e;
			Double_t HCal_e = hcal_e;

			luminosity = calc_luminosity(I_beam, run_target );
			// luminosity = 1.0;


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
			h_W->Fill( W, pn_weight );
			h_W2->Fill(pow(W, 2), pn_weight);

			//Use the electron kinematics to predict the proton momedntum assuming elastic scattering on free proton at rest (will need to correct for fermi motion):
			E_pp = nu + Mp; // Get energy of the proton
			E_nucleon = sqrt(pow(pp,2)+pow(Mp,2)); // Check on E_pp, same
			h_E_pp->Fill( E_pp ); // Fill histogram

			KE_p = nu; // For elastics
			h_KE_p->Fill( KE_p );

			dx = hcal_x - x_expected_HCal;
			dy = hcal_y - y_expected_HCal;

		//Resolve the hadron spots without cuts
			h_simc_dx->Fill( dx, pn_weight);
			h_simc_dy->Fill( dy, pn_weight);
			h_simc_dxdy->Fill( dy, dx, pn_weight );
			h_xy->Fill( hcal_y, hcal_x );

			h_simc_W_dx->Fill(dx, W, pn_weight);
			h_simc_W2_dx->Fill(dx, pow(W, 2), pn_weight);

		//Histograms by hadron
			if( is_p ){
				h_simc_dx_p->Fill(dx, pn_weight);
				h_mc_p_sigma->Fill( mc_sigma);
				h_mc_p_sigma_weighted->Fill( mc_sigma, pn_weight );
			}
			if( is_n ){
				h_simc_dx_n->Fill(dx, pn_weight);
				h_mc_n_sigma->Fill( mc_sigma);
				h_mc_n_sigma_weighted->Fill( mc_sigma, pn_weight );
			}

		// Preliminary HCal projections with single cut on W
			if( fabs(W - W_mean) < W_sigma ){
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

			//SIMC branch variables with Wcut
				h_mc_simc_Q2_wcut->Fill( mc_simc_Q2 );
				h_mc_simc_epsilon_wcut->Fill( mc_simc_epsilon );
				h_mc_simc_Ebeam_wcut->Fill( mc_simc_Ebeam );
				h_mc_simc_p_n_wcut->Fill( mc_simc_p_n );
				h_mc_simc_theta_n_wcut->Fill( TMath::RadToDeg()*mc_simc_theta_n );
				h_mc_simc_p_e_wcut->Fill( mc_simc_p_e );
				h_mc_simc_theta_e_wcut->Fill( TMath::RadToDeg()*mc_simc_theta_e );

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
				if( kine_select == 4 && sbsfieldscale_select == 30 ){
					dx_p_scale = 1.0;
					dx_n_scale = 1.0;
				}
				if( kine_select == 8 && sbsfieldscale_select == 70 ){
					dx_p_scale = 1.0;
					dx_n_scale = 1.0;
				}

				if( sbsfieldscale_select != 0 ){
					dx_p = lookup_simc_dxdy(kine_select, sbsfieldscale_select, run_target, "dx_p");
					dx_p_sigma = dx_p_scale*lookup_simc_dxdy(kine_select, sbsfieldscale_select, run_target, "dx_p_sigma");
					dy_p = lookup_simc_dxdy(kine_select, sbsfieldscale_select, run_target, "dy");
					dy_p_sigma = lookup_simc_dxdy(kine_select, sbsfieldscale_select, run_target, "dy_sigma");
					dx_n = lookup_simc_dxdy(kine_select, sbsfieldscale_select, run_target, "dx_n");
					dx_n_sigma = dx_n_scale*lookup_simc_dxdy(kine_select, sbsfieldscale_select, run_target, "dx_n_sigma");
					dy_n = lookup_simc_dxdy(kine_select, sbsfieldscale_select, run_target, "dy");
					dy_n_sigma = lookup_simc_dxdy(kine_select, sbsfieldscale_select, run_target, "dy_sigma");
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
				if( HCal_on && is_n && apply_fcut ) h_simc_dxdy_ncut->Fill( dy, dx );
				if( HCal_on && is_p && apply_fcut ) h_simc_dxdy_pcut->Fill( dy, dx );

			//----------neutron
				if( HCal_on && is_n && apply_fcut ){
					if( !calc_pn_weight ){
						pn_weight = (1.0/3.0);					
					}

					// if( (hcal_x - dx_pn_max )>hcal_x_fmin ){
						h_simc_dxdy_fcut->Fill( dy, dx, pn_weight );
						h_simc_dx_fcut->Fill( dx, pn_weight );
						h_W_fcut->Fill( W, pn_weight );
						h_xy_fcut->Fill( hcal_y, hcal_x, pn_weight );
						h_xy_cut_n->Fill( hcal_y, hcal_x );
						h_simc_dx_n_fcut->Fill(dx, pn_weight);
						if( fabs(W - W_mean) < W_sigma ){
							h_simc_dx_wcut_fcut->Fill(dx, pn_weight);
							h_simc_dx_n_wcut_fcut->Fill(dx, pn_weight);
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
						h_simc_dxdy_fcut->Fill( dy, dx, pn_weight );
						h_simc_dx_fcut->Fill( dx, pn_weight );
						h_W_fcut->Fill( W, pn_weight );
						h_xy_fcut->Fill( hcal_y, hcal_x, pn_weight );
						h_xy_cut_p->Fill( hcal_y, hcal_x );
						h_simc_dx_p_fcut->Fill(dx, pn_weight);
						if( fabs(W - W_mean) < W_sigma ){
							h_simc_dx_wcut_fcut->Fill(dx, pn_weight);
							h_simc_dx_p_wcut_fcut->Fill(dx, pn_weight);
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

			h_Wrecon->Fill(Wrecon, pn_weight);
			h_W2recon->Fill(W2recon, pn_weight);

		//end of events loop  
	    }
		cout << "---------------------------------------" << endl;
		cout << "-----Finished going through events-----" << endl;
		cout << "---------------------------------------" << endl;
		for( int bin = 1; bin < h_mc_n_weight->GetXaxis()->GetNbins(); bin++ ){
			double n_weight = h_mc_n_weight->GetBinContent(bin);
			double p_weight = h_mc_p_weight->GetBinContent(bin);

			h_mc_n_p_weight_Ratio->Fill( p_weight/n_weight );

			if( use_hist_file ){
				double n_final_weight = h_mc_n_final_weight->GetBinContent(bin);
				double p_final_weight = h_mc_p_final_weight->GetBinContent(bin);
				h_mc_n_p_final_weight_Ratio->Fill( p_final_weight/n_final_weight );
			}
		}

		//get the MC sigma n/p ratio:
		for( int bin = 1; bin < h_mc_n_sigma->GetNbinsX(); bin++ ){
			if( h_mc_p_sigma->GetBinContent(bin) != 0.0 ){
				double np_ratio = h_mc_n_sigma->GetBinContent(bin)/h_mc_p_sigma->GetBinContent(bin);
				double np_ratio_weighted = h_mc_n_sigma_weighted->GetBinContent(bin)/h_mc_p_sigma_weighted->GetBinContent(bin);

				h_mc_np_sigma_ratio->SetBinContent(bin, np_ratio);
				h_mc_np_sigma_ratio_weighted->SetBinContent(bin, np_ratio_weighted);
			}
		}

		outfile->Write();
		outfile->Close();
	}


	TFile *infile = new TFile(outfilename.Data(), "READ");
	TH1D *hin_dx_fcut, *hin_dx, *hin_dx_wcut_fcut;

	TH2D *hin_dxdy, *hin_dxdy_fcut, *hin_dxdy_wcut;

	hin_dx = static_cast<TH1D*>(infile->Get("h_simc_dx"));
	hin_dx_fcut = static_cast<TH1D*>(infile->Get("h_simc_dx_fcut"));
	hin_dx_wcut_fcut = static_cast<TH1D*>(infile->Get("h_simc_dx_wcut_fcut"));
	hin_dxdy = static_cast<TH2D*>(infile->Get("h_simc_dxdy"));
	hin_dxdy_fcut = static_cast<TH2D*>(infile->Get("h_simc_dxdy_fcut"));
	hin_dxdy_wcut = static_cast<TH2D*>(infile->Get("h_simc_dxdy_wcut"));

    TCanvas *c_dxdy = new TCanvas("c_dxdy", "c_dxdy", 600, 500);
    hin_dxdy->Draw("colz");

    TCanvas *c_dxdy_fcut = new TCanvas("c_dxdy_fcut", "c_dxdy_fcut", 600, 500);
    hin_dxdy_fcut->Draw("colz");

    TCanvas *c_dxdy_wcut = new TCanvas("c_dxdy_wcut", "c_dxdy_wcut", 600, 500);
    hin_dxdy_wcut->Draw("colz");

    TCanvas *c_dx = new TCanvas("c_dx", "c_dx", 600, 500);
    hin_dx->Draw("hist");

    TCanvas *c_dx_fcut = new TCanvas("c_dx_fcut", "c_dx_fcut", 600, 500);
    hin_dx_fcut->Draw("hist");

    TCanvas *c_dx_wcut_fcut = new TCanvas("c_dx_wcut_fcut", "c_dx_wcut_fcut", 600, 500);

	if( calibrate ){
    	TFile *infile = new TFile(outfilename.Data(), "READ");
    	hin_dx_wcut = static_cast<TH1D*>(infile->Get("h_simc_dx_wcut"));
    	hin_dy_wcut = static_cast<TH1D*>(infile->Get("h_simc_dy_wcut"));

    	TCanvas *c_dx_wcut = new TCanvas("c_dx_wcut", "c_dx_wcut", 600, 500);
    	hin_dx_wcut->Draw("hist");
  	
  	//------ p -------
    	TF1 *fit_dx_p = new TF1("fit_dx_p", fit_gaus, -1.5, -0.2, 3);
  
    	fit_dx_p->SetParName(0, "dx_p Norm");
		fit_dx_p->SetParName(1, "dx_p Center");
		fit_dx_p->SetParName(2, "dx_p Sigma");
		fit_dx_p->SetLineColor(2);

		if( kine_select == 4 && sbsfieldscale_select == 30){
			fit_dx_p->SetParLimits(0, 0, hin_dx_wcut->GetMaximum());
			fit_dx_p->SetParLimits(1, -0.7, -0.45);
			fit_dx_p->SetParLimits(2, 0.1, 0.178);
		}

		if( kine_select == 4 && sbsfieldscale_select == 50){
			fit_dx_p->SetParLimits(0, 0, hin_dx_wcut->GetMaximum());
			fit_dx_p->SetParLimits(1, -0.9, -0.8);
			fit_dx_p->SetParLimits(2, 0.1, 0.19);
		}
		//for custom field setting to get to correct 70%:
		if( kine_select == 8 && sbsfieldscale_select == 70){
			fit_dx_p->SetParLimits(0, 0, hin_dx_wcut->GetMaximum());
			fit_dx_p->SetParLimits(1, -1.2, -0.8);
			fit_dx_p->SetParLimits(2, 0.1, 0.14);			
		}

		if( kine_select == 9 && sbsfieldscale_select == 70 && run_target == "LH2"){
			fit_dx_p->SetParLimits(0, 0.8*hin_dx_wcut->GetMaximum(), hin_dx_wcut->GetMaximum());
			fit_dx_p->SetParLimits(1, -0.8, -0.7);
			fit_dx_p->SetParLimits(2, 0.1, 0.16);			
		}
		if( kine_select == 9 && sbsfieldscale_select == 70 && run_target == "LD2"){
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

		if( kine_select == 4 && sbsfieldscale_select == 30){
			fit_dx_n->SetParLimits(0, 0, (0.45)*hin_dx_wcut->GetMaximum());
			fit_dx_n->SetParLimits(1, -0.05, 0.05);
			fit_dx_n->SetParLimits(2, 0.1, 0.16);
		}	

		if( kine_select == 4 && sbsfieldscale_select == 50){
			fit_dx_n->SetParLimits(0, 0, (0.45)*hin_dx_wcut->GetMaximum());
			fit_dx_n->SetParLimits(1, 0.0, 0.13);
			fit_dx_n->SetParLimits(2, 0.1, 0.18);
		}	

		if( kine_select == 8 && sbsfieldscale_select == 70){
			fit_dx_n->SetParLimits(0, 0, (0.35)*hin_dx_wcut->GetMaximum());
			fit_dx_n->SetParLimits(1, -0.20, 0.05);
			fit_dx_n->SetParLimits(2, 0.1, 0.14);
		}
		if( kine_select == 9 && sbsfieldscale_select == 70){
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

		if( kine_select == 8 && sbsfieldscale_select == 70){
			fit_dy->SetParLimits(0, 0, hin_dy_wcut->GetMaximum());
			fit_dy->SetParLimits(1, -0.15, 0.15);
			fit_dy->SetParLimits(2, 0.1, 0.21);	
		}
		if( kine_select == 9 && sbsfieldscale_select == 70){
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
	cout << "Kinematic: SBS" << kine_select << endl;
	cout << "Target: " << run_target.Data() << endl;
	cout << "Beam Energy: " << E_beam << endl;
	cout << "SBS Field: " << sbsfieldscale_select << "%" << endl;
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
	cout << "Output file: " << outfilename.Data() << endl << endl;
	cout << "------------------------------------------------------------------"<< endl;
	
	auto total_time_end = high_resolution_clock::now();
	auto total_time_duration = duration_cast<minutes>(total_time_end - total_time_start);
	cout << "Total time for analysis: " << total_time_duration.count() << " minutes. " << endl;

}