// A simple fit in energy for signal and a background
#include <string>
#include <vector>
#include <math.h>
#include <Rand.h>
#include <fstream>
#include <iostream>

#include <TCanvas.h>
#include <ROOTNtuple.h>
#include <TRandom3.h>
#include <TH1D.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TLegend.h>

#include <BinnedED.h>
#include <BinnedEDGenerator.h>
#include <SystematicManager.h>
#include <BinnedNLLH.h>
#include <FitResult.h>
#include <Minuit.h>
#include <DistTools.h>
#include <Minuit.h>
#include <Convolution.h>
#include <Scale.h>
#include <BoolCut.h>
#include <BoxCut.h>
#include <Gaussian.h>
#include <ParameterDict.h>
#include <ContainerTools.hpp>
#include <NuOsc.h>
#include <SurvProb.h>

void LHFit(){

	TRandom3 *r1 = new TRandom3();
	r1->SetSeed(0);

	ObsSet dataRep(0);

	// Fill histograms
	AxisCollection axes;
	double Emin = -500;
	double Emax = 500;
	int numbins = 2000;
	axes.AddAxis(BinAxis("ParKE", Emin, Emax, numbins));
	
	const std::string UnOscfile   = "/home/lidgard/antinu_analysis/positronium/processed/htr_positron_modb_12_99_events100_3MeV.root";
	const std::string dataFile = "/home/lidgard/antinu_analysis/positronium/processed/htr_positron_mod_0_events100_3MeV.root";

	TFile *f_in = new TFile(dataFile.c_str());
	BinnedED data_set_pdf("data_set_pdf", axes);
	TTree *data_set_ntp = (TTree*)f_in->Get("ntmc");
	Float_t timeresidual0;
	data_set_ntp->SetBranchAddress("ntpHitTimeResidualsMC", &timeresidual0);
	for(size_t j = 0; j < data_set_ntp->GetEntries(); j++){
		data_set_ntp->GetEntry(j);
		data_set_pdf.Fill(timeresidual0);
	}
	f_in->Close();

	TFile *f_in2 = new TFile(UnOscfile.c_str());
	BinnedED UnOscPdf("UnOscPdf", axes);
	TTree *unosc_ntp = (TTree*)f_in2->Get("ntmc");
	Float_t timeresidual2;
	unosc_ntp->SetBranchAddress("ntpHitTimeResidualsMC", &timeresidual2);
	for(size_t j = 0; j < unosc_ntp->GetEntries(); j++){
		unosc_ntp->GetEntry(j);
		UnOscPdf.Fill(timeresidual2);
	}
	f_in2->Close();

	//UnOscPdf.Normalise();

	std::cout << "Initialised Pdfs" << std::endl;

	// Setting optimisation limits
	ParameterDict minima;
	minima["UnOscPdf_norm"]= 0;

	ParameterDict maxima;
	maxima["UnOscPdf_norm"]= 100000;

	ParameterDict initialval;
	initialval["UnOscPdf_norm"]= 9000;

	ParameterDict initialerr;
	initialerr["UnOscPdf_norm"]= 0.1*initialval["UnOscPdf_norm"];

	BinnedNLLH lhFunction;
	lhFunction.SetDataDist(data_set_pdf); // initialise withe the data set
	lhFunction.SetBufferAsOverflow(true);
	lhFunction.SetBuffer(0,1,1);
	lhFunction.AddDist(UnOscPdf);

	std::cout << "Built LH function " << std::endl;

	// Fit
	Minuit min;
	min.SetMethod("Migrad");
	min.SetMaxCalls(100000);
	min.SetMinima(minima);
	min.SetMaxima(maxima);
	min.SetInitialValues(initialval);
	min.SetInitialErrors(initialerr);

	// FitResult fitResult = min.Optimise(&lhFunction);
	// ParameterDict bestFit = fitResult.GetBestFit();
	// fitResult.Print();

	// plot
	TCanvas* c1 = new TCanvas("c1");
	//gStyle->SetOptStat(0);

	TH1D DataHist = DistTools::ToTH1D(data_set_pdf);
	TH1D FitHist = DistTools::ToTH1D(data_set_pdf);

	DataHist.Sumw2();
	DataHist.SetTitle("Data to Fit");
	DataHist.GetYaxis()->SetTitle("Counts");
	DataHist.Draw();

	FitHist.SetLineColor(kRed);
	FitHist.SetLineWidth(3);
	FitHist.Draw("same e");

	TLegend* leg = new TLegend(0.75,0.8,1.0,1.0);
	leg->AddEntry(&DataHist,"Data","lf");
	leg->AddEntry(&FitHist,"Fit Result","lf");
	leg->Draw();

	// write
	TFile * fitout = new TFile("/home/lidgard/antinu_analysis/positronium/processed/htr_oxsx.root","RECREATE");
	c1->Write();
	FitHist.Write();
	DataHist.Write();
	fitout->Close();

  return;
}

int main(){

  LHFit();

  return 0;
}
