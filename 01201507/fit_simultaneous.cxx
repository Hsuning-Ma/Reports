//******************************************************************************
//*********************************************************
//
//  Copyright (c) 2015, Nankai Univ. and IHEP
//  All rights reserved
//  Redistribution and use in source and binary forms,
//  with or without modification, are permitted
//
//  Name:   fit_simultaneous.cxx
//  Description
//      This is a ROOT micro
//      Which Fit the signal of exclusive and inclusive
//      at sqrt(s) = 4.23GeV, 4.26GeV, 4.36GeV 4.42GeV
//      simultaneously
//
//  Version:            10.1
//  Author:             Ma Hsuning
//  Email:              maxn@ihep.ac.cn
//  Modification Date:  2015.06.02
//
//  Version:            10.0
//  Author:             Ma Hsuning
//  Email:              maxn@ihep.ac.cn
//  Modification Date:  2015.05.19
//
//  Version:            9.0
//  Author:             Ma Hsuning
//  Email:              maxn@ihep.ac.cn
//  Modification Date:  2015.05.19
//
//  Version:            8.0
//  Author:             Ma Hsuning
//  Email:              maxn@ihep.ac.cn
//  Modification Date:  2015.05.19
//
//  Version:            7.0
//  Author:             Ma Hsuning
//  Email:              maxn@ihep.ac.cn
//  Modification Date:  2015.05.18
//
//  Version:            6.0
//  Author:             Ma Hsuning
//  Email:              maxn@ihep.ac.cn
//  Modification Date:  2015.05.18
//
//  Version:            5.0
//  Author:             Ma Hsuning
//  Email:              maxn@ihep.ac.cn
//  Modification Date:  2015.05.18
//
//  Version:            4.0
//  Author:             Ma Hsuning
//  Email:              maxn@ihep.ac.cn
//  Modification Date:  2015.05.16
//
//  Version:            3.0
//  Author:             Ma Hsuning
//  Email:              maxn@ihep.ac.cn
//  Modification Date:  2015.05.02
//
//  Version Based on:   archive/fit_simultaneous_v2.cxx
//  Author:             Ma Hsuning
//  Email:              maxn@ihep.ac.cn
//  Modification Date:  2015.04.30
//
//  Version Based on:   $SCRIPTS/Reference/guoaq/simulfit_etac_16ch_hc_effline_forsyse.cxx
//  Author:             Guo Aiqiang
//  Email:              guoaq@ihep.ac.cn
//  Modification Date:  unknown
//
//  Change Log:
//      v10.1
//          Add code to calculate the branching fraction
//      v10.0
//          Use the results of resolution_v2.cxx of Exclusive process
//          As the resolutions of both inclusive and exclusive processes
//          at each energy point
//          Convolve a unfixed common Gaussian to describe the difference
//          between data and Monto Carlo at each energy point
//      
//      Use n-order polynomial to describe the background of the inclusive process
//      Use common peak value
// 		Add sideband as the background pdf
//
//**********************************************************
//******************************************************************************

#include <map>
#include <math.h>
#include <vector>
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include <iostream>
#include <fstream>
#include "TStyle.h"
#include "TDirectory.h"
#include "TCut.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TPostScript.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "RooRealVar.h"
#include "RooNovosibirsk.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooArgusBG.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooMCStudy.h"
//#include "Get_FWHM.C"
//#include "PlotDataMC.h"
#include "RooGaussian.h"
#include "RooKeysPdf.h"
#include "RooGaussian.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "/besfs/users/maxn/scripts/Project_Charlie/RooBWMod/RooBWMod_cxx.so"

using namespace RooFit;

int fit_simultaneous()
{
    Int_t Nbin = 33;
    Double_t left = 2.54;
    Double_t right = 3.20;
    //Double_t right = 3.20;
    Int_t bg_pdf = 0;
    Bool_t b_draw_dataset = 0;
	Bool_t b_BW_Mod = 1;
    Bool_t b_exc_4230 = 1;
    Bool_t b_exc_4260 = 1;
    Bool_t b_exc_4360 = 1;
    Bool_t b_exc_4420 = 1;
    Bool_t b_inc_4230 = 1;
    Bool_t b_inc_4260 = 1;
    Bool_t b_inc_4360 = 1;
    Bool_t b_inc_4420 = 1;
    Bool_t b_sideband = 0;
    Bool_t b_common_m = 1;
    Bool_t b_inc_log_Y = 0;
    Int_t type_bg_inc = 9;
    // type_bg_inc = 0  -->   0-order Chebychev polynomial
    // type_bg_inc = 1  -->   1-order Chebychev polynomial
    // type_bg_inc = 2  -->   2-order Chebychev polynomial
    // type_bg_inc = 3  -->   3-order Chebychev polynomial
    // type_bg_inc = 4  -->   4-order Chebychev polynomial
    // type_bg_inc = 5  -->   5-order Chebychev polynomial
    // type_bg_inc = 6  -->   exponential
    // type_bg_inc = 7  -->   exponential + 3-order Chebychev polynomial
    // type_bg_inc = 8  -->   exponential + 4-order Chebychev polynomial
    // type_bg_inc = 9  -->   sideband

    TString TreeName( "kskp" );
    TString BranchName( "rec_pipigam_kskp" );

    TCut exc_4230_cut_chi2 = "chisq4c<40";
    TCut exc_4260_cut_chi2 = "chisq4c<30";
    TCut exc_4360_cut_chi2 = "chisq4c<30";
    TCut exc_4420_cut_chi2 = "chisq4c<35";
    TCut cut_signal = "rec_pipi_kskp<3.535&&rec_pipi_kskp>3.515";
    TCut cut_sideband = "(rec_pipi_kskp>3.485&&rec_pipi_kskp<3.505)||(rec_pipi_kskp>3.545&&rec_pipi_kskp<3.565)";

    RooMsgService::instance().deleteStream(0);
    RooMsgService::instance().deleteStream(1);
    RooMsgService::instance().Print();
    gErrorIgnoreLevel = kWarning;
    //gROOT->ProcessLine( ".L $SCRIPTS/Project_Charlie/RooBWMod/RooBWMod.cxx+" );

    RooRealVar mes( BranchName, "M_{#pi^{+}#pi^{-}#gamma}^{recoil}",
                    left, right );
    //RooRealVar mes( "mes", "M_{#pi^{+}#pi^{-}#gamma}^{recoil}", 2.98, left, right );

    //**************************************************************************
    // Parameters definition
    //**************************************************************************
    RooRealVar etacwidth( "_eta_c__width", "",           3.22000e-02 );
    RooRealVar etac_peak( "_eta_c_peak__", "",           2.98360e+00 );
    //RooRealVar etacwidth( "_eta_c__width", "",          3.22000e-02, 2.90000e-02, 3.40000e-02 );
    //RooRealVar etac_peak( "_eta_c_peak__", "",          2.98360e+00, 2.80000e+00, 3.20000e+00 );

    RooRealVar dif_4230_mean( "dif_4230_mean", "",      0.000000e-03, -0.01, 0.01 );
    RooRealVar dif_4260_mean( "dif_4260_mean", "",      0.000000e-03, -0.01, 0.01 );
    RooRealVar dif_4360_mean( "dif_4360_mean", "",      0.000000e-03, -0.01, 0.01 );
    RooRealVar dif_4420_mean( "dif_4420_mean", "",      0.000000e-03, -0.01, 0.01 );
    RooRealVar dif_4230_sigm( "dif_4230_sigm", "",      1.000000e-03, 0.0, 0.1 );
    RooRealVar dif_4260_sigm( "dif_4260_sigm", "",      1.000000e-03, 0.0, 0.1 );
    RooRealVar dif_4360_sigm( "dif_4360_sigm", "",      1.000000e-03, 0.0, 0.1 );
    RooRealVar dif_4420_sigm( "dif_4420_sigm", "",      1.000000e-03, 0.0, 0.1 );
    RooRealVar exc_4230_Nbkg( "exc_4230_Nbkg", "",      3.50379e+01, 0, 50 );
    RooRealVar exc_4230_Nsig( "exc_4230_Nsig", "",      5.79712e+01, 0, 80 );
    RooRealVar exc_4260_Nbkg( "exc_4260_Nbkg", "",      9.45969e+00, 0, 50 );
    RooRealVar exc_4260_Nsig( "exc_4260_Nsig", "",      4.75472e+01, 0, 60 );
    RooRealVar exc_4360_Nbkg( "exc_4360_Nbkg", "",      6.21590e+00, 0, 50 );
    RooRealVar exc_4360_Nsig( "exc_4360_Nsig", "",      4.77905e+01, 0, 60 );
    RooRealVar exc_4420_Nbkg( "exc_4420_Nbkg", "",      1.35743e+01, 0, 50 );
    RooRealVar exc_4420_Nsig( "exc_4420_Nsig", "",      6.24344e+01, 0, 80 );
    RooRealVar inc_4230_Nbkg( "inc_4230_Nbkg", "",      1.40181e+06, 0, 1500000 );
    RooRealVar inc_4230_Nsig( "inc_4230_Nsig", "",      1.19227e+04, 0, 25000 );
    RooRealVar inc_4260_Nbkg( "inc_4260_Nbkg", "",      1.02740e+06, 0, 1500000 );
    RooRealVar inc_4260_Nsig( "inc_4260_Nsig", "",      8.02902e+03, 0, 25000 );
    RooRealVar inc_4360_Nbkg( "inc_4360_Nbkg", "",      7.14899e+05, 0, 1500000 );
    RooRealVar inc_4360_Nsig( "inc_4360_Nsig", "",      7.17645e+03, 0, 15000 );
    RooRealVar inc_4420_Nbkg( "inc_4420_Nbkg", "",      1.40353e+06, 0, 1500000 );
    RooRealVar inc_4420_Nsig( "inc_4420_Nsig", "",      1.24777e+04, 0, 25000 );

    //**************************************************************************
    //  parameters
    //**************************************************************************
    //  4230
    //**************************************************************************
    RooRealVar exc_4230_frac_DG( "exc_4230_frac_DG", "",      6.44614e-01 );
    RooRealVar exc_4230_mean1( "exc_4230_mean1", "",          2.60966e-03 );
    RooRealVar exc_4230_mean2( "exc_4230_mean2", "",          2.36083e-02 );
    RooRealVar exc_4230_sig_sigma1( "exc_4230_massRes1", "",  1.12871e-02 );
    RooRealVar exc_4230_sig_sigma2( "exc_4230_massRes2", "",  2.63651e-02 );
    RooRealVar exc_4230_eff_c0( "exc_4230_eff_c0", "",        3.00342e-02 );
    RooRealVar exc_4230_eff_c1( "exc_4230_eff_c1", "",       -1.51907e-02 );
    RooRealVar exc_4230_eff_c2( "exc_4230_eff_c2", "",        2.12448e-03 );
    RooRealVar exc_4230_eff_c3( "exc_4230_eff_c3", "",       -9.06519e-03 );
    //**************************************************************************
    RooRealVar inc_4230_frac_DG( "inc_4230_frac_DG", "",      6.44614e-01 );
    RooRealVar inc_4230_mean1( "inc_4230_mean1", "",          2.60966e-03 );
    RooRealVar inc_4230_mean2( "inc_4230_mean2", "",          2.36083e-02 );
    RooRealVar inc_4230_sig_sigma1( "inc_4230_massRes1", "",  1.12871e-02 );
    RooRealVar inc_4230_sig_sigma2( "inc_4230_massRes2", "",  2.63651e-02 );
    RooRealVar inc_4230_eff_c0( "inc_4230_eff_c0", "",       -9.21247e-05 );
    RooRealVar inc_4230_eff_c1( "inc_4230_eff_c1", "",       -8.77962e-03 );
    RooRealVar inc_4230_eff_c2( "inc_4230_eff_c2", "",       -5.29260e-03 );
    RooRealVar inc_4230_eff_c3( "inc_4230_eff_c3", "",        1.56562e-03 );
    //**************************************************************************
    //  4260
    //**************************************************************************
    RooRealVar exc_4260_frac_DG( "exc_4260_frac_DG", "",      6.04471e-01 );
    RooRealVar exc_4260_mean1( "exc_4260_mean1", "",          1.73044e-03 );
    RooRealVar exc_4260_mean2( "exc_4260_mean2", "",          2.01322e-02 );
    RooRealVar exc_4260_sig_sigma1( "exc_4260_massRes1", "",  1.06856e-02 );
    RooRealVar exc_4260_sig_sigma2( "exc_4260_massRes2", "",  2.37019e-02 );
    RooRealVar exc_4260_eff_c0( "exc_4260_eff_c0", "",        8.49044e-03 );
    RooRealVar exc_4260_eff_c1( "exc_4260_eff_c1", "",       -3.13517e-02 );
    RooRealVar exc_4260_eff_c2( "exc_4260_eff_c2", "",        3.89154e-03 );
    RooRealVar exc_4260_eff_c3( "exc_4260_eff_c3", "",       -9.74477e-03 );
    //**************************************************************************
    RooRealVar inc_4260_frac_DG( "inc_4260_frac_DG", "",      6.04471e-01 );
    RooRealVar inc_4260_mean1( "inc_4260_mean1", "",          1.73044e-03 );
    RooRealVar inc_4260_mean2( "inc_4260_mean2", "",          2.01322e-02 );
    RooRealVar inc_4260_sig_sigma1( "inc_4260_massRes1", "",  1.06856e-02 );
    RooRealVar inc_4260_sig_sigma2( "inc_4260_massRes2", "",  2.37019e-02 );
    RooRealVar inc_4260_eff_c0( "inc_4260_eff_c0", "",        1.05395e-03 );
    RooRealVar inc_4260_eff_c1( "inc_4260_eff_c1", "",       -1.00394e-02 );
    RooRealVar inc_4260_eff_c2( "inc_4260_eff_c2", "",       -9.97755e-03 );
    RooRealVar inc_4260_eff_c3( "inc_4260_eff_c3", "",       -1.94744e-03 );
    //**************************************************************************
    //  4360
    //**************************************************************************
    RooRealVar exc_4360_frac_DG( "exc_4360_frac_DG", "",      6.01291e-01 );
    RooRealVar exc_4360_mean1( "exc_4360_mean1", "",          1.64410e-03 );
    RooRealVar exc_4360_mean2( "exc_4360_mean2", "",          2.05414e-02 );
    RooRealVar exc_4360_sig_sigma1( "exc_4360_massRes1", "",  1.07271e-02 );
    RooRealVar exc_4360_sig_sigma2( "exc_4360_massRes2", "",  2.35260e-02 );
    RooRealVar exc_4360_eff_c0( "exc_4360_eff_c0", "",        1.07313e-02 );
    RooRealVar exc_4360_eff_c1( "exc_4360_eff_c1", "",       -5.37816e-02 );
    RooRealVar exc_4360_eff_c2( "exc_4360_eff_c2", "",        1.63318e-02 );
    RooRealVar exc_4360_eff_c3( "exc_4360_eff_c3", "",       -2.21815e-02 );
    //**************************************************************************
    RooRealVar inc_4360_frac_DG( "inc_4360_frac_DG", "",      6.01291e-01 );
    RooRealVar inc_4360_mean1( "inc_4360_mean1", "",          1.64410e-03 );
    RooRealVar inc_4360_mean2( "inc_4360_mean2", "",          2.05414e-02 );
    RooRealVar inc_4360_sig_sigma1( "inc_4360_massRes1", "",  1.07271e-02 );
    RooRealVar inc_4360_sig_sigma2( "inc_4360_massRes2", "",  2.35260e-02 );
    RooRealVar inc_4360_eff_c0( "inc_4360_eff_c0", "",       -4.21721e-03 );
    RooRealVar inc_4360_eff_c1( "inc_4360_eff_c1", "",       -1.96160e-03 );
    RooRealVar inc_4360_eff_c2( "inc_4360_eff_c2", "",       -8.05649e-03 );
    RooRealVar inc_4360_eff_c3( "inc_4360_eff_c3", "",       -4.08510e-05 );
    //**************************************************************************
    //  4420
    //**************************************************************************
    RooRealVar exc_4420_frac_DG( "exc_4420_frac_DG", "",      6.34061e-01 );
    RooRealVar exc_4420_mean1( "exc_4420_mean1", "",          2.45160e-03 );
    RooRealVar exc_4420_mean2( "exc_4420_mean2", "",          2.20984e-02 );
    RooRealVar exc_4420_sig_sigma1( "exc_4420_massRes1", "",  1.12767e-02 );
    RooRealVar exc_4420_sig_sigma2( "exc_4420_massRes2", "",  2.57629e-02 );
    RooRealVar exc_4420_eff_c0( "exc_4420_eff_c0", "",        1.26035e-02 );
    RooRealVar exc_4420_eff_c1( "exc_4420_eff_c1", "",       -3.53234e-02 );
    RooRealVar exc_4420_eff_c2( "exc_4420_eff_c2", "",        1.13393e-02 );
    RooRealVar exc_4420_eff_c3( "exc_4420_eff_c3", "",       -6.15503e-03 );
    //**************************************************************************
    RooRealVar inc_4420_frac_DG( "inc_4420_frac_DG", "",      6.34061e-01 );
    RooRealVar inc_4420_mean1( "inc_4420_mean1", "",          2.45160e-03 );
    RooRealVar inc_4420_mean2( "inc_4420_mean2", "",          2.20984e-02 );
    RooRealVar inc_4420_sig_sigma1( "inc_4420_massRes1", "",  1.12767e-02 );
    RooRealVar inc_4420_sig_sigma2( "inc_4420_massRes2", "",  2.57629e-02 );
    RooRealVar inc_4420_eff_c0( "inc_4420_eff_c0", "",       -2.32829e-03 );
    RooRealVar inc_4420_eff_c1( "inc_4420_eff_c1", "",       -7.74328e-03 );
    RooRealVar inc_4420_eff_c2( "inc_4420_eff_c2", "",       -2.82231e-03 );
    RooRealVar inc_4420_eff_c3( "inc_4420_eff_c3", "",       -5.87680e-03 );
    //**************************************************************************
    //**************************************************************************
	// 	Prepare for fit
    //**************************************************************************
    //**************************************************************************
    //  4230
    //**************************************************************************
    //**************************************************************************
    //  Data Sets
    //**************************************************************************
    //  exclusive data set
    //**************************************************************************
    TFile *exc_4230_file = new TFile( "$PCHARLIE/merge/PipihcExclusive_4230_data_without_cut.root" );
    TTree *exc_4230_tree = ( TTree* )exc_4230_file->Get( TreeName );
    TH1F *h_exc_4230_sig = new TH1F( "h_exc_4230_sig", "", Nbin, left, right );
    exc_4230_tree->Project( "h_exc_4230_sig", BranchName,
                           exc_4230_cut_chi2 + cut_signal );
    RooDataHist exc_4230_sig_dataset( "exc_4230_sig_dataset", "dataset",
                                      mes, h_exc_4230_sig );
    //**************************************************************************
    // inclusive dataset
    //**************************************************************************
    TFile *inc_4230_file = new TFile( "$PCHARLIE/merge/PipihcInclusive_4230_data_without_cut.root" );
    TTree *inc_4230_tree = ( TTree* )inc_4230_file->Get( TreeName );
    TH1F *h_inc_4230_sig = new TH1F( "h_inc_4230_sig", "", Nbin, left, right );
    TH1F *h_inc_4230_sb = new TH1F( "h_inc_4230_sb", "", Nbin, left, right );
    inc_4230_tree->Project( "h_inc_4230_sig", BranchName,
                            cut_signal );
    inc_4230_tree->Project( "h_inc_4230_sb",
                            BranchName, cut_sideband );
    RooDataHist inc_4230_sig_dataset( "inc_4230_sig_dataset", "dataset", 
                                      mes, h_inc_4230_sig );
    RooDataHist inc_4230_sb_dataset( "inc_4230_sb_dataset", "dataset", 
                                     mes, h_inc_4230_sb );
    //**************************************************************************
    //  Construct exclusive signal pdf
    //**************************************************************************
    RooGaussian exc_4230_gaus1( "exc_4230_gaus1", "gaus1",
                                mes, exc_4230_mean1, exc_4230_sig_sigma1 );
    RooGaussian exc_4230_gaus2( "exc_4230_gaus2", "gaus2",
                                mes, exc_4230_mean2, exc_4230_sig_sigma2 );
    RooAddPdf exc_4230_DG( "exc_4230_DG", "DG",
                           RooArgList( exc_4230_gaus1, exc_4230_gaus2 ),
                           exc_4230_frac_DG );
    RooBWMod exc_4230_BW( "exc_4230_BW", "BW",
                          mes, etac_peak, etacwidth,
                          RooArgList( exc_4230_eff_c0, exc_4230_eff_c1,
                                      exc_4230_eff_c2, exc_4230_eff_c3 ) );
    RooFFTConvPdf exc_4230_sig( "exc_4230_sig", "exc_4230_sig",
                                mes, exc_4230_BW, exc_4230_DG );

    //**************************************************************************
    //  Construct exclusive background p.d.f.
    //**************************************************************************
    RooChebychev exc_4230_bg( "exc_4230_bg", "background pdf", mes,
                              RooArgList( ) );

    //**************************************************************************
    //  Construct exclusive composite pdf
    //**************************************************************************
    RooAddPdf exc_4230_model( "exc_4230_model", "model",
                              RooArgList( exc_4230_sig, exc_4230_bg ),
                              RooArgList( exc_4230_Nsig, exc_4230_Nbkg ) );

    //**************************************************************************
    //  Construct inclusive signal pdf
    //**************************************************************************
    RooGaussian inc_4230_gaus1( "inc_4230_gaus1", "gaus1",
                                mes, inc_4230_mean1, inc_4230_sig_sigma1 );
    RooGaussian inc_4230_gaus2( "inc_4230_gaus2", "gaus2",
                                mes, inc_4230_mean2, inc_4230_sig_sigma2 );
    RooAddPdf inc_4230_DG( "inc_4230_DG", "DG",
                           RooArgList( inc_4230_gaus1, inc_4230_gaus2 ),
                           inc_4230_frac_DG );
    RooBWMod inc_4230_BW( "inc_4230_BW", "BW",
                          mes, etac_peak, etacwidth,
                          RooArgList( inc_4230_eff_c0, inc_4230_eff_c1,
                                      inc_4230_eff_c2, inc_4230_eff_c3 ) );
    RooFFTConvPdf inc_4230_sig( "inc_4230_sig", "inc_4230_sig",
                                mes, inc_4230_BW, inc_4230_DG );

    //**************************************************************************
    //  Construct inclusive background pdf
    //**************************************************************************
    RooAbsPdf *inc_4230_bg = new RooHistPdf( "inc_4230_bg", "background pdf",
                                            mes, inc_4230_sb_dataset, 3 );

    //**************************************************************************
    //  Construct inclusive composite pdf
    //**************************************************************************
    RooAddPdf inc_4230_model( "inc_4230_model", "model",
                              RooArgList( inc_4230_sig, *inc_4230_bg ),
                              RooArgList( inc_4230_Nsig, inc_4230_Nbkg ) );

    //**************************************************************************
    //**************************************************************************
    //  4260
    //**************************************************************************
    //**************************************************************************
    //  Data Sets
    //**************************************************************************
    //  exclusive data set
    //**************************************************************************
    TFile *exc_4260_file = new TFile( "$PCHARLIE/merge/PipihcExclusive_4260_data_without_cut.root" );
    TTree *exc_4260_tree = ( TTree* )exc_4260_file->Get( TreeName );
    TH1F *h_exc_4260_sig = new TH1F( "h_exc_4260_sig", "", Nbin, left, right );
    exc_4260_tree->Project( "h_exc_4260_sig", BranchName,
                           exc_4260_cut_chi2 + cut_signal );
    RooDataHist exc_4260_sig_dataset( "exc_4260_sig_dataset", "dataset",
                                      mes, h_exc_4260_sig );
    //**************************************************************************
    // inclusive dataset
    //**************************************************************************
    TFile *inc_4260_file = new TFile( "$PCHARLIE/merge/PipihcInclusive_4260_data_without_cut.root" );
    TTree *inc_4260_tree = ( TTree* )inc_4260_file->Get( TreeName );
    TH1F *h_inc_4260_sig = new TH1F( "h_inc_4260_sig", "", Nbin, left, right );
    TH1F *h_inc_4260_sb = new TH1F( "h_inc_4260_sb", "", Nbin, left, right );
    inc_4260_tree->Project( "h_inc_4260_sig", BranchName,
                            cut_signal );
    inc_4260_tree->Project( "h_inc_4260_sb",
                            BranchName, cut_sideband );
    RooDataHist inc_4260_sig_dataset( "inc_4260_sig_dataset", "dataset", 
                                      mes, h_inc_4260_sig );
    RooDataHist inc_4260_sb_dataset( "inc_4260_sb_dataset", "dataset", 
                                     mes, h_inc_4260_sb );
    //**************************************************************************
    //  Construct exclusive signal pdf
    //**************************************************************************
    RooGaussian exc_4260_gaus1( "exc_4260_gaus1", "gaus1",
                                mes, exc_4260_mean1, exc_4260_sig_sigma1 );
    RooGaussian exc_4260_gaus2( "exc_4260_gaus2", "gaus2",
                                mes, exc_4260_mean2, exc_4260_sig_sigma2 );
    RooAddPdf exc_4260_DG( "exc_4260_DG", "DG",
                           RooArgList( exc_4260_gaus1, exc_4260_gaus2 ),
                           exc_4260_frac_DG );
    RooBWMod exc_4260_BW( "exc_4260_BW", "BW",
                          mes, etac_peak, etacwidth,
                          RooArgList( exc_4260_eff_c0, exc_4260_eff_c1,
                                      exc_4260_eff_c2, exc_4260_eff_c3 ) );
    RooFFTConvPdf exc_4260_sig( "exc_4260_sig", "exc_4260_sig",
                                mes, exc_4260_BW, exc_4260_DG );

    //**************************************************************************
    //  Construct exclusive background p.d.f.
    //**************************************************************************
    RooChebychev exc_4260_bg( "exc_4260_bg", "background pdf", mes,
                              RooArgList( ) );

    //**************************************************************************
    //  Construct exclusive composite pdf
    //**************************************************************************
    RooAddPdf exc_4260_model( "exc_4260_model", "model",
                              RooArgList( exc_4260_sig, exc_4260_bg ),
                              RooArgList( exc_4260_Nsig, exc_4260_Nbkg ) );

    //**************************************************************************
    //  Construct inclusive signal pdf
    //**************************************************************************
    RooGaussian inc_4260_gaus1( "inc_4260_gaus1", "gaus1",
                                mes, inc_4260_mean1, inc_4260_sig_sigma1 );
    RooGaussian inc_4260_gaus2( "inc_4260_gaus2", "gaus2",
                                mes, inc_4260_mean2, inc_4260_sig_sigma2 );
    RooAddPdf inc_4260_DG( "inc_4260_DG", "DG",
                           RooArgList( inc_4260_gaus1, inc_4260_gaus2 ),
                           inc_4260_frac_DG );
    RooBWMod inc_4260_BW( "inc_4260_BW", "BW",
                          mes, etac_peak, etacwidth,
                          RooArgList( inc_4260_eff_c0, inc_4260_eff_c1,
                                      inc_4260_eff_c2, inc_4260_eff_c3 ) );
    RooFFTConvPdf inc_4260_sig( "inc_4260_sig", "inc_4260_sig",
                                mes, inc_4260_BW, inc_4260_DG );

    //**************************************************************************
    //  Construct inclusive background pdf
    //**************************************************************************
    RooAbsPdf *inc_4260_bg = new RooHistPdf( "inc_4260_bg", "background pdf",
                                            mes, inc_4260_sb_dataset, 3 );

    //**************************************************************************
    //  Construct inclusive composite pdf
    //**************************************************************************
    RooAddPdf inc_4260_model( "inc_4260_model", "model",
                              RooArgList( inc_4260_sig, *inc_4260_bg ),
                              RooArgList( inc_4260_Nsig, inc_4260_Nbkg ) );

    //**************************************************************************
    //**************************************************************************
    //  4360
    //**************************************************************************
    //**************************************************************************
    //  Data Sets
    //**************************************************************************
    //  exclusive data set
    //**************************************************************************
    TFile *exc_4360_file = new TFile( "$PCHARLIE/merge/PipihcExclusive_4360_data_without_cut.root" );
    TTree *exc_4360_tree = ( TTree* )exc_4360_file->Get( TreeName );
    TH1F *h_exc_4360_sig = new TH1F( "h_exc_4360_sig", "", Nbin, left, right );
    exc_4360_tree->Project( "h_exc_4360_sig", BranchName,
                           exc_4360_cut_chi2 + cut_signal );
    RooDataHist exc_4360_sig_dataset( "exc_4360_sig_dataset", "dataset",
                                      mes, h_exc_4360_sig );
    //**************************************************************************
    // inclusive dataset
    //**************************************************************************
    TFile *inc_4360_file = new TFile( "$PCHARLIE/merge/PipihcInclusive_4360_data_without_cut.root" );
    TTree *inc_4360_tree = ( TTree* )inc_4360_file->Get( TreeName );
    TH1F *h_inc_4360_sig = new TH1F( "h_inc_4360_sig", "", Nbin, left, right );
    TH1F *h_inc_4360_sb = new TH1F( "h_inc_4360_sb", "", Nbin, left, right );
    inc_4360_tree->Project( "h_inc_4360_sig", BranchName,
                            cut_signal );
    inc_4360_tree->Project( "h_inc_4360_sb",
                            BranchName, cut_sideband );
    RooDataHist inc_4360_sig_dataset( "inc_4360_sig_dataset", "dataset", 
                                      mes, h_inc_4360_sig );
    RooDataHist inc_4360_sb_dataset( "inc_4360_sb_dataset", "dataset", 
                                     mes, h_inc_4360_sb );
    //**************************************************************************
    //  Construct exclusive signal pdf
    //**************************************************************************
    RooGaussian exc_4360_gaus1( "exc_4360_gaus1", "gaus1",
                                mes, exc_4360_mean1, exc_4360_sig_sigma1 );
    RooGaussian exc_4360_gaus2( "exc_4360_gaus2", "gaus2",
                                mes, exc_4360_mean2, exc_4360_sig_sigma2 );
    RooAddPdf exc_4360_DG( "exc_4360_DG", "DG",
                           RooArgList( exc_4360_gaus1, exc_4360_gaus2 ),
                           exc_4360_frac_DG );
    RooBWMod exc_4360_BW( "exc_4360_BW", "BW",
                          mes, etac_peak, etacwidth,
                          RooArgList( exc_4360_eff_c0, exc_4360_eff_c1,
                                      exc_4360_eff_c2, exc_4360_eff_c3 ) );
    RooFFTConvPdf exc_4360_sig( "exc_4360_sig", "exc_4360_sig",
                                mes, exc_4360_BW, exc_4360_DG );

    //**************************************************************************
    //  Construct exclusive background p.d.f.
    //**************************************************************************
    RooChebychev exc_4360_bg( "exc_4360_bg", "background pdf", mes,
                              RooArgList( ) );

    //**************************************************************************
    //  Construct exclusive composite pdf
    //**************************************************************************
    RooAddPdf exc_4360_model( "exc_4360_model", "model",
                              RooArgList( exc_4360_sig, exc_4360_bg ),
                              RooArgList( exc_4360_Nsig, exc_4360_Nbkg ) );

    //**************************************************************************
    //  Construct inclusive signal pdf
    //**************************************************************************
    RooGaussian inc_4360_gaus1( "inc_4360_gaus1", "gaus1",
                                mes, inc_4360_mean1, inc_4360_sig_sigma1 );
    RooGaussian inc_4360_gaus2( "inc_4360_gaus2", "gaus2",
                                mes, inc_4360_mean2, inc_4360_sig_sigma2 );
    RooAddPdf inc_4360_DG( "inc_4360_DG", "DG",
                           RooArgList( inc_4360_gaus1, inc_4360_gaus2 ),
                           inc_4360_frac_DG );
    RooBWMod inc_4360_BW( "inc_4360_BW", "BW",
                          mes, etac_peak, etacwidth,
                          RooArgList( inc_4360_eff_c0, inc_4360_eff_c1,
                                      inc_4360_eff_c2, inc_4360_eff_c3 ) );
    RooFFTConvPdf inc_4360_sig( "inc_4360_sig", "inc_4360_sig",
                                mes, inc_4360_BW, inc_4360_DG );

    //**************************************************************************
    //  Construct inclusive background pdf
    //**************************************************************************
    RooAbsPdf *inc_4360_bg = new RooHistPdf( "inc_4360_bg", "background pdf",
                                            mes, inc_4360_sb_dataset, 3 );

    //**************************************************************************
    //  Construct inclusive composite pdf
    //**************************************************************************
    RooAddPdf inc_4360_model( "inc_4360_model", "model",
                              RooArgList( inc_4360_sig, *inc_4360_bg ),
                              RooArgList( inc_4360_Nsig, inc_4360_Nbkg ) );

    //**************************************************************************
    //**************************************************************************
    //  4420
    //**************************************************************************
    //**************************************************************************
    //  Data Sets
    //**************************************************************************
    //  exclusive data set
    //**************************************************************************
    TFile *exc_4420_file = new TFile( "$PCHARLIE/merge/PipihcExclusive_4420_data_without_cut.root" );
    TTree *exc_4420_tree = ( TTree* )exc_4420_file->Get( TreeName );
    TH1F *h_exc_4420_sig = new TH1F( "h_exc_4420_sig", "", Nbin, left, right );
    exc_4420_tree->Project( "h_exc_4420_sig", BranchName,
                           exc_4420_cut_chi2 + cut_signal );
    RooDataHist exc_4420_sig_dataset( "exc_4420_sig_dataset", "dataset",
                                      mes, h_exc_4420_sig );
    //**************************************************************************
    // inclusive dataset
    //**************************************************************************
    TFile *inc_4420_file = new TFile( "$PCHARLIE/merge/PipihcInclusive_4420_data_without_cut.root" );
    TTree *inc_4420_tree = ( TTree* )inc_4420_file->Get( TreeName );
    TH1F *h_inc_4420_sig = new TH1F( "h_inc_4420_sig", "", Nbin, left, right );
    TH1F *h_inc_4420_sb = new TH1F( "h_inc_4420_sb", "", Nbin, left, right );
    inc_4420_tree->Project( "h_inc_4420_sig", BranchName,
                            cut_signal );
    inc_4420_tree->Project( "h_inc_4420_sb",
                            BranchName, cut_sideband );
    RooDataHist inc_4420_sig_dataset( "inc_4420_sig_dataset", "dataset", 
                                      mes, h_inc_4420_sig );
    RooDataHist inc_4420_sb_dataset( "inc_4420_sb_dataset", "dataset", 
                                     mes, h_inc_4420_sb );
    //**************************************************************************
    //  Construct exclusive signal pdf
    //**************************************************************************
    RooGaussian exc_4420_gaus1( "exc_4420_gaus1", "gaus1",
                                mes, exc_4420_mean1, exc_4420_sig_sigma1 );
    RooGaussian exc_4420_gaus2( "exc_4420_gaus2", "gaus2",
                                mes, exc_4420_mean2, exc_4420_sig_sigma2 );
    RooAddPdf exc_4420_DG( "exc_4420_DG", "DG",
                           RooArgList( exc_4420_gaus1, exc_4420_gaus2 ),
                           exc_4420_frac_DG );
    RooBWMod exc_4420_BW( "exc_4420_BW", "BW",
                          mes, etac_peak, etacwidth,
                          RooArgList( exc_4420_eff_c0, exc_4420_eff_c1,
                                      exc_4420_eff_c2, exc_4420_eff_c3 ) );
    RooFFTConvPdf exc_4420_sig( "exc_4420_sig", "exc_4420_sig",
                                mes, exc_4420_BW, exc_4420_DG );

    //**************************************************************************
    //  Construct exclusive background p.d.f.
    //**************************************************************************
    RooChebychev exc_4420_bg( "exc_4420_bg", "background pdf", mes,
                              RooArgList( ) );

    //**************************************************************************
    //  Construct exclusive composite pdf
    //**************************************************************************
    RooAddPdf exc_4420_model( "exc_4420_model", "model",
                              RooArgList( exc_4420_sig, exc_4420_bg ),
                              RooArgList( exc_4420_Nsig, exc_4420_Nbkg ) );

    //**************************************************************************
    //  Construct inclusive signal pdf
    //**************************************************************************
    RooGaussian inc_4420_gaus1( "inc_4420_gaus1", "gaus1",
                                mes, inc_4420_mean1, inc_4420_sig_sigma1 );
    RooGaussian inc_4420_gaus2( "inc_4420_gaus2", "gaus2",
                                mes, inc_4420_mean2, inc_4420_sig_sigma2 );
    RooAddPdf inc_4420_DG( "inc_4420_DG", "DG",
                           RooArgList( inc_4420_gaus1, inc_4420_gaus2 ),
                           inc_4420_frac_DG );
    RooBWMod inc_4420_BW( "inc_4420_BW", "BW",
                          mes, etac_peak, etacwidth,
                          RooArgList( inc_4420_eff_c0, inc_4420_eff_c1,
                                      inc_4420_eff_c2, inc_4420_eff_c3 ) );
    RooFFTConvPdf inc_4420_sig( "inc_4420_sig", "inc_4420_sig",
                                mes, inc_4420_BW, inc_4420_DG );

    //**************************************************************************
    //  Construct inclusive background pdf
    //**************************************************************************
    RooAbsPdf *inc_4420_bg = new RooHistPdf( "inc_4420_bg", "background pdf",
                                            mes, inc_4420_sb_dataset, 3 );

    //**************************************************************************
    //  Construct inclusive composite pdf
    //**************************************************************************
    RooAddPdf inc_4420_model( "inc_4420_model", "model",
                              RooArgList( inc_4420_sig, *inc_4420_bg ),
                              RooArgList( inc_4420_Nsig, inc_4420_Nbkg ) );

    //**************************************************************************
    //  Start Simultaneous Fit
    //**************************************************************************
    //  Define category to distinguish inclusive and exclusive events
    //**************************************************************************
    RooCategory sample( "sample", "sample" );
    sample.defineType( "all" );
    if ( b_exc_4230 )   sample.defineType( "exc_4230" );
    if ( b_inc_4230 )   sample.defineType( "inc_4230" );
    if ( b_inc_4230 )   sample.defineType( "inc_4230_sideband" );
    if ( b_exc_4260 )   sample.defineType( "exc_4260" );
    if ( b_inc_4260 )   sample.defineType( "inc_4260" );
    if ( b_inc_4260 )   sample.defineType( "inc_4260_sideband" );
    if ( b_exc_4360 )   sample.defineType( "exc_4360" );
    if ( b_inc_4360 )   sample.defineType( "inc_4360" );
    if ( b_inc_4360 )   sample.defineType( "inc_4360_sideband" );
    if ( b_exc_4420 )   sample.defineType( "exc_4420" );
    if ( b_inc_4420 )   sample.defineType( "inc_4420" );
    if ( b_inc_4420 )   sample.defineType( "inc_4420_sideband" );

    //**************************************************************************
    //  Construct combined dataset
    //**************************************************************************
    TH1F *hist = new TH1F( "hist", "", Nbin, left, right );
    RooDataHist hist_all( "hist_all", "dataset", mes, hist );
    RooDataHist combData( "combData", "combined data",
                          mes, Index(sample),
                          Import( "all", hist_all ) );
    if ( b_exc_4230 )
    {
        RooDataHist exc_4230_sig_combData( "exc_4230_sig_combData",
                                           "combined data",
                                           mes, Index(sample),
                                           Import( "exc_4230",
                                                   exc_4230_sig_dataset ) );
    }
    if ( b_inc_4230 )
    {
        RooDataHist inc_4230_sig_combData( "inc_4230_sig_combData",
                                           "combined data",
                                           mes, Index(sample),
                                           Import( "inc_4230",
                                                   inc_4230_sig_dataset ) );
        RooDataHist inc_4230_sb_combData( "inc_4230_sb_combData",
                                          "combined data",
                                          mes, Index(sample),
                                          Import( "inc_4230_sideband",
                                                  inc_4230_sb_dataset ) );
    }
    if ( b_exc_4260 )
    {
        RooDataHist exc_4260_sig_combData( "exc_4260_sig_combData",
                                           "combined data",
                                           mes, Index(sample),
                                           Import( "exc_4260",
                                                   exc_4260_sig_dataset ) );
    }
    if ( b_inc_4260 )
    {
        RooDataHist inc_4260_sig_combData( "inc_4260_sig_combData",
                                           "combined data",
                                           mes, Index(sample),
                                           Import( "inc_4260",
                                                   inc_4260_sig_dataset ) );
        RooDataHist inc_4260_sb_combData( "inc_4260_sb_combData",
                                          "combined data",
                                          mes, Index(sample),
                                          Import( "inc_4260_sideband",
                                                  inc_4260_sb_dataset ) );
    }
    if ( b_exc_4360 )
    {
        RooDataHist exc_4360_sig_combData( "exc_4360_sig_combData",
                                           "combined data",
                                           mes, Index(sample),
                                           Import( "exc_4360",
                                                   exc_4360_sig_dataset ) );
    }
    if ( b_inc_4360 )
    {
        RooDataHist inc_4360_sig_combData( "inc_4360_sig_combData",
                                           "combined data",
                                           mes, Index(sample),
                                           Import( "inc_4360",
                                                   inc_4360_sig_dataset ) );
        RooDataHist inc_4360_sb_combData( "inc_4360_sb_combData",
                                          "combined data",
                                          mes, Index(sample),
                                          Import( "inc_4360_sideband",
                                                  inc_4360_sb_dataset ) );
    }
    if ( b_exc_4420 )
    {
        RooDataHist exc_4420_sig_combData( "exc_4420_sig_combData",
                                          "combined data",
                                          mes, Index(sample),
                                          Import( "exc_4420",
                                                  exc_4420_sig_dataset ) );
    }
    if ( b_inc_4420 )
    {
        RooDataHist inc_4420_sig_combData( "inc_4420_sig_combData",
                                           "combined data",
                                           mes, Index(sample),
                                           Import( "inc_4420",
                                                  inc_4420_sig_dataset ) );
        RooDataHist inc_4420_sb_combData( "inc_4420_sb_combData",
                                          "combined data",
                                          mes, Index(sample),
                                          Import( "inc_4420_sideband",
                                                  inc_4420_sb_dataset ) );
    }
    if ( b_exc_4230 )   combData.add( exc_4230_sig_combData );
    if ( b_inc_4230 )   combData.add( inc_4230_sig_combData );
    if ( b_inc_4230 && b_sideband )   combData.add( inc_4230_sb_combData );
    if ( b_exc_4260 )   combData.add( exc_4260_sig_combData );
    if ( b_inc_4260 )   combData.add( inc_4260_sig_combData );
    if ( b_inc_4260 && b_sideband )   combData.add( inc_4260_sb_combData );
    if ( b_exc_4360 )   combData.add( exc_4360_sig_combData );
    if ( b_inc_4360 )   combData.add( inc_4360_sig_combData );
    if ( b_inc_4360 && b_sideband )   combData.add( inc_4360_sb_combData );
    if ( b_exc_4420 )   combData.add( exc_4420_sig_combData );
    if ( b_inc_4420 )   combData.add( inc_4420_sig_combData );
    if ( b_inc_4420 && b_sideband )   combData.add( inc_4420_sb_combData );

    //**************************************************************************
    //  Construct a simultaneous pdf
    //**************************************************************************
    RooSimultaneous simPdf( "simPdf", "simultaneous pdf", sample );
    if ( b_exc_4230 )   simPdf.addPdf( exc_4230_model, "exc_4230" );
    if ( b_inc_4230 )   simPdf.addPdf( inc_4230_model, "inc_4230" );
    if ( b_inc_4230 && b_sideband )   simPdf.addPdf( inc_4230_bg, "inc_4230_sideband" );
    if ( b_exc_4260 )   simPdf.addPdf( exc_4260_model, "exc_4260" );
    if ( b_inc_4260 )   simPdf.addPdf( inc_4260_model, "inc_4260" );
    if ( b_inc_4260 && b_sideband )   simPdf.addPdf( inc_4260_bg, "inc_4260_sideband" );
    if ( b_exc_4360 )   simPdf.addPdf( exc_4360_model, "exc_4360" );
    if ( b_inc_4360 )   simPdf.addPdf( inc_4360_model, "inc_4360" );
    if ( b_inc_4360 && b_sideband )   simPdf.addPdf( inc_4360_bg, "inc_4360_sideband" );
    if ( b_exc_4420 )   simPdf.addPdf( exc_4420_model, "exc_4420" );
    if ( b_inc_4420 )   simPdf.addPdf( inc_4420_model, "inc_4420" );
    if ( b_inc_4420 && b_sideband )   simPdf.addPdf( inc_4420_bg, "inc_4420_sideband" );

    //**************************************************************************
    //  Perform a simultaneous fit
    //**************************************************************************
    //simPdf.fitTo( combData, Extended() );
    simPdf.fitTo( combData, Extended(), Minos(kTRUE) );

    //**************************************************************************
    if ( b_exc_4230 )
    {
        RooPlot *exc_4230_frame = mes.frame( Bins( Nbin ),
                                             Title( "Exclusive at #sqrt{s}=4.23GeV" ) );
        combData.plotOn( exc_4230_frame, Binning( Nbin ),
                         //DataError( RooAbsData::SumW2 ),
                         Cut( "sample==sample::exc_4230" ) );
        simPdf.plotOn( exc_4230_frame, Slice( sample, "exc_4230" ),
                       ProjWData( sample, combData ) );
        simPdf.plotOn( exc_4230_frame, Slice( sample, "exc_4230" ),
                       Components( "exc_4230_sig" ),
                       ProjWData( sample, combData ),
                       LineStyle(3), LineColor(3),
                       RooFit::Name( "exc_4230_sig" ) );
        simPdf.plotOn( exc_4230_frame, Slice( sample, "exc_4230" ),
                       Components( "exc_4230_bg" ),
                       ProjWData( sample, combData ),
                       LineStyle(2), LineColor(2),
                       RooFit::Name( "exc_4230_bg" ) );
    }
    //**************************************************************************
    if ( b_exc_4260 )
    {
        RooPlot *exc_4260_frame = mes.frame( Bins( Nbin ),
                                             Title( "Exclusive at #sqrt{s}=4.26GeV" ) );
        combData.plotOn( exc_4260_frame, Binning( Nbin ), 
                         //DataError( RooAbsData::SumW2 ),
                         Cut( "sample==sample::exc_4260" ) );
        simPdf.plotOn( exc_4260_frame, Slice( sample, "exc_4260" ),
                       ProjWData( sample, combData ) );
        simPdf.plotOn( exc_4260_frame, Slice( sample, "exc_4260" ),
                       Components( "exc_4260_sig" ),
                       ProjWData( sample, combData ),
                       LineStyle(3), LineColor(3),
                       RooFit::Name( "exc_4260_sig" ) );
        simPdf.plotOn( exc_4260_frame, Slice( sample, "exc_4260" ),
                       Components( "exc_4260_bg" ),
                       ProjWData( sample, combData ),
                       LineStyle(2), LineColor(2),
                       RooFit::Name( "exc_4260_bg" ) );
    }
    //**************************************************************************
    if ( b_exc_4360 )
    {
        RooPlot *exc_4360_frame = mes.frame( Bins( Nbin ),
                                             Title( "Exclusive at #sqrt{s}=4.36GeV" ) );
        combData.plotOn( exc_4360_frame, Binning( Nbin ),
                         //DataError( RooAbsData::SumW2 ),
                         Cut( "sample==sample::exc_4360" ) );
        simPdf.plotOn( exc_4360_frame, Slice( sample, "exc_4360" ),
                       ProjWData( sample, combData ) );
        simPdf.plotOn( exc_4360_frame, Slice( sample, "exc_4360" ),
                       Components( "exc_4360_sig" ),
                       ProjWData( sample, combData ),
                       LineStyle(3), LineColor(3),
                       RooFit::Name( "exc_4360_sig" ) );
        simPdf.plotOn( exc_4360_frame, Slice( sample, "exc_4360" ),
                       Components( "exc_4360_bg" ),
                       ProjWData( sample, combData ),
                       LineStyle(2), LineColor(2),
                       RooFit::Name( "exc_4360_bg" ) );
    }
    //**************************************************************************
    if ( b_exc_4420 )
    {
        RooPlot *exc_4420_frame = mes.frame( Bins( Nbin ),
                                             Title( "Exclusive at #sqrt{s}=4.42GeV" ) );
        combData.plotOn( exc_4420_frame, Binning( Nbin ), 
                         //DataError( RooAbsData::SumW2 ),
                         Cut( "sample==sample::exc_4420" ) );
        simPdf.plotOn( exc_4420_frame, Slice( sample, "exc_4420" ),
                       ProjWData( sample, combData ) );
        simPdf.plotOn( exc_4420_frame, Slice( sample, "exc_4420" ),
                       Components( "exc_4420_sig" ),
                       ProjWData( sample, combData ),
                       LineStyle(3), LineColor(3),
                       RooFit::Name( "exc_4420_sig" ) );
        simPdf.plotOn( exc_4420_frame, Slice( sample, "exc_4420" ),
                       Components( "exc_4420_bg" ),
                       ProjWData( sample, combData ),
                       LineStyle(2), LineColor(2),
                       RooFit::Name( "exc_4420_bg" ) );
    }
    //**************************************************************************
    if ( b_inc_4230 )
    {
        RooPlot *inc_4230_frame = mes.frame( Bins( Nbin ),
                                             Title( "Inclusive at #sqrt{s}=4.23GeV" ) );
        combData.plotOn( inc_4230_frame, Binning( Nbin ), 
                         //DataError( RooAbsData::SumW2 ),
                         Cut( "sample==sample::inc_4230" ) );
        simPdf.plotOn( inc_4230_frame, Slice( sample, "inc_4230" ),
                       ProjWData( sample, combData ) );
        simPdf.plotOn( inc_4230_frame, Slice( sample, "inc_4230" ),
                       Components( "inc_4230_sig" ),
                       ProjWData( sample, combData ),
                       LineStyle(3), LineColor(3),
                       RooFit::Name( "inc_4230_sig" ) );
        simPdf.plotOn( inc_4230_frame, Slice( sample, "inc_4230" ),
                       Components( "inc_4230_bg" ),
                       ProjWData( sample, combData ),
                       LineStyle(2), LineColor(2),
                       RooFit::Name( "inc_4230_bg" ) );
    }
    //**************************************************************************
    if ( b_inc_4260 )
    {
        RooPlot *inc_4260_frame = mes.frame( Bins( Nbin ),
                                             Title( "Inclusive at #sqrt{s}=4.26GeV" ) );
        combData.plotOn( inc_4260_frame, Binning( Nbin ), 
                         //DataError( RooAbsData::SumW2 ),
                         Cut( "sample==sample::inc_4260" ) );
        simPdf.plotOn( inc_4260_frame, Slice( sample, "inc_4260" ),
                       ProjWData( sample, combData ) );
        simPdf.plotOn( inc_4260_frame, Slice( sample, "inc_4260" ),
                       Components( "inc_4260_sig" ),
                       ProjWData( sample, combData ),
                       LineStyle(3), LineColor(3),
                       RooFit::Name( "inc_4260_sig" ) );
        simPdf.plotOn( inc_4260_frame, Slice( sample, "inc_4260" ),
                       Components( "inc_4260_bg" ),
                       ProjWData( sample, combData ),
                       LineStyle(2), LineColor(2),
                       RooFit::Name( "inc_4260_bg" ) );
    }
    //**************************************************************************
    if ( b_inc_4360 )
    {
        RooPlot *inc_4360_frame = mes.frame( Bins( Nbin ),
                                             Title( "Inclusive at #sqrt{s}=4.36GeV" ) );
        combData.plotOn( inc_4360_frame, Binning( Nbin ), 
                         //DataError( RooAbsData::SumW2 ),
                         Cut( "sample==sample::inc_4360" ) );
        simPdf.plotOn( inc_4360_frame, Slice( sample, "inc_4360" ),
                       ProjWData( sample, combData ) );
        simPdf.plotOn( inc_4360_frame, Slice( sample, "inc_4360" ),
                       Components( "inc_4360_sig" ),
                       ProjWData( sample, combData ),
                       LineStyle(3), LineColor(3),
                       RooFit::Name( "inc_4360_sig" ) );
        simPdf.plotOn( inc_4360_frame, Slice( sample, "inc_4360" ),
                       Components( "inc_4360_bg" ),
                       ProjWData( sample, combData ),
                       LineStyle(2), LineColor(2),
                       RooFit::Name( "inc_4360_bg" ) );
    }
    //**************************************************************************
    if ( b_inc_4420 )
    {
        RooPlot *inc_4420_frame = mes.frame( Bins( Nbin ),
                                             Title( "Incluive at #sqrt{s}=4.42GeV" ) );
        combData.plotOn( inc_4420_frame, Binning( Nbin ), 
                         //DataError( RooAbsData::SumW2 ),
                         Cut( "sample==sample::inc_4420" ) );
        simPdf.plotOn( inc_4420_frame, Slice( sample, "inc_4420" ),
                       ProjWData( sample, combData ) );
        simPdf.plotOn( inc_4420_frame, Slice( sample, "inc_4420" ),
                       Components( "inc_4420_sig" ),
                       ProjWData( sample, combData ),
                       LineStyle(3), LineColor(3),
                       RooFit::Name( "inc_4420_sig" ) );
        simPdf.plotOn( inc_4420_frame, Slice( sample, "inc_4420" ),
                       Components( "inc_4420_bg" ),
                       ProjWData( sample, combData ),
                       LineStyle(2), LineColor(2),
                       RooFit::Name( "inc_4420_bg" ) );
    }

    //**************************************************************************
    //**************************************************************************
    TCanvas *c = new TCanvas( "SimultaneousFit", "SimultaneousFit", 1600, 1200 );
    c->SetFillColor( 10 );
    c->Divide( 3, 4 );

    //**************************************************************************
    //**************************************************************************
    if ( b_exc_4230 )
    {
        c->cd(1);
        gPad->SetLeftMargin( 0.15 );
        exc_4230_frame->GetYaxis()->SetTitleOffset( 1.4 );
        exc_4230_frame->Draw();
        c->Update();
        TLegend *exc_4230_leg = new TLegend( 0.2, 0.6, 0.5, 0.8 );
        exc_4230_leg->SetFillStyle(0);
        exc_4230_leg->SetBorderSize(0);
        exc_4230_leg->AddEntry( exc_4230_frame->findObject( "exc_4230_bg" ),
                                "background", "l" );
        exc_4230_leg->AddEntry( exc_4230_frame->findObject( "exc_4230_sig" ),
                                "signal", "l" );
        //TLatex lat1;
        //lat1.SetNDC();
        //lat1.SetTextSize( 0.04 );
        //lat1.DrawLatex( 0.5, 0.02, Form( "#chi^{2} = ", 1.0*(chi2exc_4230) ) );
        exc_4230_leg->Draw();
    }
    //**************************************************************************
    if ( b_inc_4230 )
    {
        c->cd(2);
        gPad->SetLeftMargin( 0.15 );
        if ( b_inc_log_Y )  gPad->SetLogy();
        inc_4230_frame->GetYaxis()->SetTitleOffset( 1.4 );
        inc_4230_frame->Draw();
        c->Update();
        TLegend *inc_4230_leg = new TLegend( 0.2, 0.6, 0.5, 0.8 );
        inc_4230_leg->SetFillStyle(0);
        inc_4230_leg->SetBorderSize(0);
        inc_4230_leg->AddEntry( inc_4230_frame->findObject( "inc_4230_bg" ),
                                "background", "l" );
        inc_4230_leg->AddEntry( inc_4230_frame->findObject( "inc_4230_sig" ),
                                "signal", "l" );
        //TLatex lat2;
        //lat2.SetNDC();
        //lat2.SetTextSize( 0.04 );
        //lat2.DrawLatex( 0.5, 0.02, Form( "#chi^{2} = ", 1.0*(chi2inc_4230) ) );
        Double_t inc_4230_nsig = inc_4230_Nsig.getVal();
        inc_4230_leg->Draw();
        //**************************************************************************
        RooHist *hresid_4230 = inc_4230_frame->residHist();
        RooHist *pull_4230 = inc_4230_frame->pullHist();
        pull_4230->SetMarkerColor(kRed);
        RooPlot *pull_4230_frame = mes.frame( Title( "Pull Distribution" ) );
        pull_4230_frame->addPlotable( hresid_4230, "P" );
        inc_4230_sig.plotOn( pull_4230_frame, LineWidth(2), LineColor(2),
                             Normalization( inc_4230_nsig, RooAbsReal::NumEvent ) );
        c->cd(3);
        gPad->SetLeftMargin( 0.15 );
        pull_4230_frame->GetYaxis()->SetTitleOffset(1.6);
        pull_4230_frame->Draw();
    }
    
    //**************************************************************************
    //**************************************************************************
    if ( b_exc_4260 )
    {
        c->cd(4);
        gPad->SetLeftMargin( 0.15 );
        exc_4260_frame->GetYaxis()->SetTitleOffset( 1.4 );
        exc_4260_frame->Draw();
        c->Update();
        TLegend *exc_4260_leg = new TLegend( 0.2, 0.6, 0.5, 0.8 );
        exc_4260_leg->SetFillStyle(0);
        exc_4260_leg->SetBorderSize(0);
        exc_4260_leg->AddEntry( exc_4260_frame->findObject( "exc_4260_bg" ),
                                "background", "l" );
        exc_4260_leg->AddEntry( exc_4260_frame->findObject( "exc_4260_sig" ),
                                "signal", "l" );
        //TLatex lat1;
        //lat1.SetNDC();
        //lat1.SetTextSize( 0.04 );
        //lat1.DrawLatex( 0.5, 0.02, Form( "#chi^{2} = ", 1.0*(chi2exc_4260) ) );
        exc_4260_leg->Draw();
    }
    //**************************************************************************
    if ( b_inc_4260 )
    {
        
        c->cd(5);
        gPad->SetLeftMargin( 0.15 );
        if ( b_inc_log_Y )  gPad->SetLogy();
        inc_4260_frame->GetYaxis()->SetTitleOffset( 1.4 );
        inc_4260_frame->Draw();
        c->Update();
        TLegend *inc_4260_leg = new TLegend( 0.2, 0.6, 0.5, 0.8 );
        inc_4260_leg->SetFillStyle(0);
        inc_4260_leg->SetBorderSize(0);
        inc_4260_leg->AddEntry( inc_4260_frame->findObject( "inc_4260_bg" ),
                                "background", "l" );
        inc_4260_leg->AddEntry( inc_4260_frame->findObject( "inc_4260_sig" ),
                                "signal", "l" );
        //TLatex lat2;
        //lat2.SetNDC();
        //lat2.SetTextSize( 0.04 );
        //lat2.DrawLatex( 0.5, 0.02, Form( "#chi^{2} = ", 1.0*(chi2inc_4260) ) );
        Double_t inc_4260_nsig = inc_4260_Nsig.getVal();
        inc_4260_leg->Draw();
        //**************************************************************************
        RooHist *hresid_4260 = inc_4260_frame->residHist();
        RooHist *pull_4260 = inc_4260_frame->pullHist();
        pull_4260->SetMarkerColor(kRed);
        RooPlot *pull_4260_frame = mes.frame( Title( "Pull Distribution" ) );
        pull_4260_frame->addPlotable( hresid_4260, "P" );
        inc_4260_sig.plotOn( pull_4260_frame, LineWidth(2), LineColor(2),
                             Normalization( inc_4260_nsig, RooAbsReal::NumEvent ) );
        c->cd(6);
        gPad->SetLeftMargin( 0.15 );
        pull_4260_frame->GetYaxis()->SetTitleOffset(1.6);
        pull_4260_frame->Draw();
    }
    
    //**************************************************************************
    //**************************************************************************
    if ( b_exc_4360 )
    {
        c->cd(7);
        gPad->SetLeftMargin( 0.15 );
        exc_4360_frame->GetYaxis()->SetTitleOffset( 1.4 );
        exc_4360_frame->Draw();
        c->Update();
        TLegend *exc_4360_leg = new TLegend( 0.2, 0.6, 0.5, 0.8 );
        exc_4360_leg->SetFillStyle(0);
        exc_4360_leg->SetBorderSize(0);
        exc_4360_leg->AddEntry( exc_4360_frame->findObject( "exc_4360_bg" ),
                                "background", "l" );
        exc_4360_leg->AddEntry( exc_4360_frame->findObject( "exc_4360_sig" ),
                                "signal", "l" );
        //TLatex lat1;
        //lat1.SetNDC();
        //lat1.SetTextSize( 0.04 );
        //lat1.DrawLatex( 0.5, 0.02, Form( "#chi^{2} = ", 1.0*(chi2exc_4360) ) );
        exc_4360_leg->Draw();
    }
    //**************************************************************************
    if ( b_inc_4360 )
    {
        c->cd(8);
        gPad->SetLeftMargin( 0.15 );
        if ( b_inc_log_Y )  gPad->SetLogy();
        inc_4360_frame->GetYaxis()->SetTitleOffset( 1.4 );
        inc_4360_frame->Draw();
        c->Update();
        TLegend *inc_4360_leg = new TLegend( 0.2, 0.6, 0.5, 0.8 );
        inc_4360_leg->SetFillStyle(0);
        inc_4360_leg->SetBorderSize(0);
        inc_4360_leg->AddEntry( inc_4360_frame->findObject( "inc_4360_bg" ),
                                "background", "l" );
        inc_4360_leg->AddEntry( inc_4360_frame->findObject( "inc_4360_sig" ),
                                "signal", "l" );
        //TLatex lat2;
        //lat2.SetNDC();
        //lat2.SetTextSize( 0.04 );
        //lat2.DrawLatex( 0.5, 0.02, Form( "#chi^{2} = ", 1.0*(chi2inc_4360) ) );
        Double_t inc_4360_nsig = inc_4360_Nsig.getVal();
        inc_4360_leg->Draw();
        //**************************************************************************
        RooHist *hresid_4360 = inc_4360_frame->residHist();
        RooHist *pull_4360 = inc_4360_frame->pullHist();
        pull_4360->SetMarkerColor(kRed);
        RooPlot *pull_4360_frame = mes.frame( Title( "Pull Distribution" ) );
        pull_4360_frame->addPlotable( hresid_4360, "P" );
        inc_4360_sig.plotOn( pull_4360_frame, LineWidth(2), LineColor(2),
                             Normalization( inc_4360_nsig, RooAbsReal::NumEvent ) );
        c->cd(9);
        gPad->SetLeftMargin( 0.15 );
        pull_4360_frame->GetYaxis()->SetTitleOffset(1.6);
        pull_4360_frame->Draw();
    }
    
    //**************************************************************************
    //**************************************************************************
    if ( b_exc_4420 )
    {
        c->cd(10);
        gPad->SetLeftMargin( 0.15 );
        exc_4420_frame->GetYaxis()->SetTitleOffset( 1.4 );
        exc_4420_frame->Draw();
        c->Update();
        TLegend *exc_4420_leg = new TLegend( 0.2, 0.6, 0.5, 0.8 );
        exc_4420_leg->SetFillStyle(0);
        exc_4420_leg->SetBorderSize(0);
        exc_4420_leg->AddEntry( exc_4420_frame->findObject( "exc_4420_bg" ),
                                "background", "l" );
        exc_4420_leg->AddEntry( exc_4420_frame->findObject( "exc_4420_sig" ),
                                "signal", "l" );
        //TLatex lat1;
        //lat1.SetNDC();
        //lat1.SetTextSize( 0.04 );
        //lat1.DrawLatex( 0.5, 0.02, Form( "#chi^{2} = ", 1.0*(chi2exc_4420) ) );
        exc_4420_leg->Draw();
    }
    //**************************************************************************
    if ( b_inc_4420 )
    {
        c->cd(11);
        gPad->SetLeftMargin( 0.15 );
        if ( b_inc_log_Y )  gPad->SetLogy();
        inc_4420_frame->GetYaxis()->SetTitleOffset( 1.4 );
        inc_4420_frame->Draw();
        c->Update();
        TLegend *inc_4420_leg = new TLegend( 0.2, 0.6, 0.5, 0.8 );
        inc_4420_leg->SetFillStyle(0);
        inc_4420_leg->SetBorderSize(0);
        inc_4420_leg->AddEntry( inc_4420_frame->findObject( "inc_4420_bg" ),
                                "background", "l" );
        inc_4420_leg->AddEntry( inc_4420_frame->findObject( "inc_4420_sig" ),
                                "signal", "l" );
        //TLatex lat2;
        //lat2.SetNDC();
        //lat2.SetTextSize( 0.04 );
        //lat2.DrawLatex( 0.5, 0.02, Form( "#chi^{2} = ", 1.0*(chi2inc_4420) ) );
        Double_t inc_4420_nsig = inc_4420_Nsig.getVal();
        inc_4420_leg->Draw();
        //**************************************************************************
        RooHist *hresid_4420 = inc_4420_frame->residHist();
        RooHist *pull_4420 = inc_4420_frame->pullHist();
        pull_4420->SetMarkerColor(kRed);
        RooPlot *pull_4420_frame = mes.frame( Title( "Pull Distribution" ) );
        pull_4420_frame->addPlotable( hresid_4420, "P" );
        inc_4420_sig.plotOn( pull_4420_frame, LineWidth(2), LineColor(2),
                             Normalization( inc_4420_nsig, RooAbsReal::NumEvent ) );
        c->cd(12);
        gPad->SetLeftMargin( 0.15 );
        pull_4420_frame->GetYaxis()->SetTitleOffset(1.6);
        pull_4420_frame->Draw();
    }

    //**************************************************************************
    //  Add the components together
    //**************************************************************************
    Double_t Nsig_exc = 0;
    Double_t Nbkg_exc = 0;
    Nsig_exc += exc_4230_Nsig.getVal();
    Nsig_exc += exc_4260_Nsig.getVal();
    Nsig_exc += exc_4360_Nsig.getVal();
    Nsig_exc += exc_4420_Nsig.getVal();
    Nbkg_exc += exc_4230_Nbkg.getVal();
    Nbkg_exc += exc_4260_Nbkg.getVal();
    Nbkg_exc += exc_4360_Nbkg.getVal();
    Nbkg_exc += exc_4420_Nbkg.getVal();
    Double_t Ntot_exc = Nsig_exc + Nbkg_exc;

    Double_t Nsig_inc = 0;
    Double_t Nbkg_inc = 0;
    Nsig_inc += inc_4230_Nsig.getVal();
    Nsig_inc += inc_4260_Nsig.getVal();
    Nsig_inc += inc_4360_Nsig.getVal();
    Nsig_inc += inc_4420_Nsig.getVal();
    Nbkg_inc += inc_4230_Nbkg.getVal();
    Nbkg_inc += inc_4260_Nbkg.getVal();
    Nbkg_inc += inc_4360_Nbkg.getVal();
    Nbkg_inc += inc_4420_Nbkg.getVal();
    Double_t Ntot_inc = Nsig_inc + Nbkg_inc;
    Double_t inc_total_nsig = Ntot_inc;
    //Double_t inc_total_nsig = 3*Ntot_inc;
    
    //**************************************************************************
    //**************************************************************************
    RooPlot *exc_total_frame = mes.frame( Bins( Nbin ), Title( "Sum up" ) );
    combData.plotOn( exc_total_frame, Binning( Nbin ),
                    Cut( "sample==sample::exc_4230||sample==sample::exc_4260||sample==sample::exc_4360||sample==sample::exc_4420" ) );
    simPdf.plotOn( exc_total_frame, ProjWData( sample, combData ),
                   Components( "exc_4230_model,exc_4260_model,exc_4360_model,exc_4420_model" ) );
    simPdf.plotOn( exc_total_frame, ProjWData( sample, combData ),
                   LineStyle(2), LineColor(2),
                   Components( "exc_4230_sig,exc_4260_sig,exc_4360_sig,exc_4420_sig" ) );
    simPdf.plotOn( exc_total_frame, ProjWData( sample, combData ),
                   LineStyle(3), LineColor(3),
                   Components( "exc_4230_bg,exc_4260_bg,exc_4360_bg,exc_4420_bg" ) );
    
    RooPlot *inc_total_frame = mes.frame( Bins( Nbin ), Title( "Sum up" ) );
    combData.plotOn( inc_total_frame, Binning( Nbin ),
                     Cut( "sample==sample::inc_4230||sample==sample::inc_4260||sample==sample::inc_4360||sample==sample::inc_4420" ) );
    simPdf.plotOn( inc_total_frame, ProjWData( sample, combData ),
                   Components( "inc_4230_model,inc_4260_model,inc_4360_model,inc_4420_model" ) );
    simPdf.plotOn( inc_total_frame, ProjWData( sample, combData ),
                   LineStyle(2), LineColor(2),
                   Components( "inc_4230_sig,inc_4260_sig,inc_4360_sig,inc_4420_sig" ) );
    simPdf.plotOn( inc_total_frame, ProjWData( sample, combData ),
                   LineStyle(3), LineColor(3),
                   Components( "inc_4230_bg,inc_4260_bg,inc_4360_bg,inc_4420_bg" ) );

    RooHist *hresid_total = inc_total_frame->residHist();
    RooHist *pull_total = inc_total_frame->pullHist();
    pull_total->SetMarkerColor(kRed);
    RooPlot *pull_total_frame = mes.frame( Title( "Pull Distribution" ) );
    pull_total_frame->addPlotable( hresid_total, "P" );
    simPdf.plotOn( pull_total_frame, ProjWData( sample, combData ),
                   LineWidth(2), LineColor(2),
                   Normalization( inc_total_nsig, RooAbsReal::NumEvent ),
                   Components( "inc_4230_sig,inc_4260_sig,inc_4360_sig,inc_4420_sig" ) );
    pull_total_frame->GetYaxis()->SetTitleOffset(1.6);

    TCanvas *canvas = new TCanvas( "Total", "Total", 800, 400 );
    canvas->Divide(3);
    canvas->cd(1);
    exc_total_frame->Draw();
    canvas->cd(2);
    inc_total_frame->Draw();
    canvas->cd(3);
    pull_total_frame->Draw();
    //**************************************************************************
    //**************************************************************************
    c->Print( "$PCHARLIE/pics/Pipihc_data_fit_simultaneous.eps" );
    c->Print( "$PCHARLIE/pics/Pipihc_data_fit_simultaneous.png" );
    canvas->Print( "$PCHARLIE/pics/Pipihc_data_fit_simultaneous_sumup.eps" );
    canvas->Print( "$PCHARLIE/pics/Pipihc_data_fit_simultaneous_sumup.png" );
    //**************************************************************************
    //**************************************************************************
    //  Calculate the Branching fraction
    //**************************************************************************
    //**************************************************************************
    //  4230
    //**************************************************************************
    //**************************************************************************
    Double_t N_exc_4230_eff = 3.13110e+04;
    Double_t E_exc_4230_eff = 1.76719e+02;
    Double_t N_inc_4230_eff = 9.62398e+04;
    Double_t E_inc_4230_eff = 4.73114e+02;
    Double_t N_exc_4230_sig = exc_4230_Nsig.getVal();
    Double_t E_exc_4230_sig = exc_4230_Nsig.getError();
    Double_t N_inc_4230_sig = inc_4230_Nsig.getVal();
    Double_t E_inc_4230_sig = inc_4230_Nsig.getError();
    Double_t N_Br_KS = 0.692;
    Double_t E_Br_KS = 0.005;
    Double_t N_Br_etac_4230 = ( N_exc_4230_sig/N_inc_4230_sig )
                             *( N_inc_4230_eff/N_exc_4230_eff )
                             *( 1/N_Br_KS );
    Double_t E_Br_etac_4230 = sqrt( pow( ( E_exc_4230_sig/N_exc_4230_sig ), 2 )
                                  + pow( ( E_inc_4230_sig/N_inc_4230_sig ), 2 )
                                  + pow( ( E_exc_4230_eff/N_exc_4230_eff ), 2 )
                                  + pow( ( E_inc_4230_eff/N_inc_4230_eff ), 2 )
                                  + pow( ( E_Br_KS/N_Br_KS ), 2 )
                                  )*N_Br_etac_4230;
    //**************************************************************************
    //**************************************************************************
    //  4260
    //**************************************************************************
    //**************************************************************************
    Double_t N_exc_4260_eff = 2.78880e+04;
    Double_t E_exc_4260_eff = 1.66978e+02;
    Double_t N_inc_4260_eff = 8.82858e+04;
    Double_t E_inc_4260_eff = 4.52781e+02;
    Double_t N_exc_4260_sig = exc_4260_Nsig.getVal();
    Double_t E_exc_4260_sig = exc_4260_Nsig.getError();
    Double_t N_inc_4260_sig = inc_4260_Nsig.getVal();
    Double_t E_inc_4260_sig = inc_4260_Nsig.getError();
    Double_t N_Br_etac_4260 = ( N_exc_4260_sig/N_inc_4260_sig )
                             *( N_inc_4260_eff/N_exc_4260_eff )
                             *( 1/N_Br_KS );
    Double_t E_Br_etac_4260 = sqrt( pow( ( E_exc_4260_sig/N_exc_4260_sig ), 2 )
                                  + pow( ( E_inc_4260_sig/N_inc_4260_sig ), 2 )
                                  + pow( ( E_exc_4260_eff/N_exc_4260_eff ), 2 )
                                  + pow( ( E_inc_4260_eff/N_inc_4260_eff ), 2 )
                                  + pow( ( E_Br_KS/N_Br_KS ), 2 )
                                  )*N_Br_etac_4260;
    //**************************************************************************
    //**************************************************************************
    //  4360
    //**************************************************************************
    //**************************************************************************
    Double_t N_exc_4360_eff = 2.78239e+04;
    Double_t E_exc_4360_eff = 1.66785e+02;
    Double_t N_inc_4360_eff = 8.51853e+04;
    Double_t E_inc_4360_eff = 4.44782e+02;
    Double_t N_exc_4360_sig = exc_4360_Nsig.getVal();
    Double_t E_exc_4360_sig = exc_4360_Nsig.getError();
    Double_t N_inc_4360_sig = inc_4360_Nsig.getVal();
    Double_t E_inc_4360_sig = inc_4360_Nsig.getError();
    Double_t N_Br_etac_4360 = ( N_exc_4360_sig/N_inc_4360_sig )
                             *( N_inc_4360_eff/N_exc_4360_eff )
                             *( 1/N_Br_KS );
    Double_t E_Br_etac_4360 = sqrt( pow( ( E_exc_4360_sig/N_exc_4360_sig ), 2 )
                                  + pow( ( E_inc_4360_sig/N_inc_4360_sig ), 2 )
                                  + pow( ( E_exc_4360_eff/N_exc_4360_eff ), 2 )
                                  + pow( ( E_inc_4360_eff/N_inc_4360_eff ), 2 )
                                  + pow( ( E_Br_KS/N_Br_KS ), 2 )
                                  )*N_Br_etac_4360;
    //**************************************************************************
    //**************************************************************************
    //  4420
    //**************************************************************************
    //**************************************************************************
    Double_t N_exc_4420_eff = 3.58073e+04;
    Double_t E_exc_4420_eff = 1.89181e+02;
    Double_t N_inc_4420_eff = 1.02298e+05;
    Double_t E_inc_4420_eff = 3.68825e+02;
    Double_t N_exc_4420_sig = exc_4420_Nsig.getVal();
    Double_t E_exc_4420_sig = exc_4420_Nsig.getError();
    Double_t N_inc_4420_sig = inc_4420_Nsig.getVal();
    Double_t E_inc_4420_sig = inc_4420_Nsig.getError();
    Double_t N_Br_etac_4420 = ( N_exc_4420_sig/N_inc_4420_sig )
                             *( N_inc_4420_eff/N_exc_4420_eff )
                             *( 1/N_Br_KS );
    Double_t E_Br_etac_4420 = sqrt( pow( ( E_exc_4420_sig/N_exc_4420_sig ), 2 )
                                  + pow( ( E_inc_4420_sig/N_inc_4420_sig ), 2 )
                                  + pow( ( E_exc_4420_eff/N_exc_4420_eff ), 2 )
                                  + pow( ( E_inc_4420_eff/N_inc_4420_eff ), 2 )
                                  + pow( ( E_Br_KS/N_Br_KS ), 2 )
                                  )*N_Br_etac_4420;
    //**************************************************************************
    //**************************************************************************
    //  Calculate the average branching fraction using weight-average method
    //**************************************************************************
    //**************************************************************************
    Double_t tmp = 1/pow( E_Br_etac_4230, 2 )
                 + 1/pow( E_Br_etac_4260, 2 )
                 + 1/pow( E_Br_etac_4360, 2 )
                 + 1/pow( E_Br_etac_4420, 2 );
    Double_t N_Br_etac_average, E_Br_etac_average;
    N_Br_etac_average = ( N_Br_etac_4230/pow( E_Br_etac_4230, 2 )
                        + N_Br_etac_4260/pow( E_Br_etac_4260, 2 )
                        + N_Br_etac_4360/pow( E_Br_etac_4360, 2 )
                        + N_Br_etac_4420/pow( E_Br_etac_4420, 2 )
                        )/tmp; 
    E_Br_etac_average = sqrt( 1/tmp );
    //**************************************************************************
    //**************************************************************************
    cout << "######################################" << endl;
    cout << "\t\tchisq_exc\t\tchisq_inc" << endl;
    cout << "4230:\t\t" << exc_4230_frame->chiSquare(2) << "\t\t\t" << inc_4230_frame->chiSquare(2) << endl;
    cout << "4260:\t\t" << exc_4260_frame->chiSquare(2) << "\t\t\t" << inc_4260_frame->chiSquare(2) << endl;
    cout << "4360:\t\t" << exc_4360_frame->chiSquare(2) << "\t\t\t" << inc_4360_frame->chiSquare(2) << endl;
    cout << "4420:\t\t" << exc_4420_frame->chiSquare(2) << "\t\t\t" << inc_4420_frame->chiSquare(2) << endl;
    cout << "######################################" << endl;
    cout << "\t\tBranching fraction\t\tError" << endl;
    cout << "4230:\t\t" << N_Br_etac_4230 << "\t\t\t" << E_Br_etac_4230 << endl;
    cout << "4260:\t\t" << N_Br_etac_4260 << "\t\t\t" << E_Br_etac_4260 << endl;
    cout << "4360:\t\t" << N_Br_etac_4360 << "\t\t\t" << E_Br_etac_4360 << endl;
    cout << "4420:\t\t" << N_Br_etac_4420 << "\t\t\t" << E_Br_etac_4420 << endl;
    cout << "average:\t" << N_Br_etac_average << "\t\t\t" << E_Br_etac_average << endl;
    cout << "######################################" << endl;
}
//######################################
//		chisq_exc		chisq_inc
//4230:	2.06245			15.5738
//4260:	4.10109			9.73613
//4360:	3.55926			10.9923
//4420:	3.6784			16.2379
//######################################
