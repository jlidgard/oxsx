// A fit in energy for signal and a background
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

bool test_file_exists(const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    }
    else {
        return false;
    }
}

void readConstraintsInfoFile(const std::string &runConstraintsInfoFileName, std::string reactor_name, Double_t &fit_mean, Double_t &fit_mean_err, Double_t &fit_sigma, Double_t &fit_sigma_err ) {
    //
    // Read constraints info from text file - generated by the 'passes_plot' script which fits number of events per year per reactor
    // reads one at a time, searches for a reactor name if a match is found the fit information is returned
    //

    std::ifstream in;
    in.open(runConstraintsInfoFileName.c_str());
    //std::cout << "opening file: " << runConstraintsInfoFileName.c_str() << std::endl;

    std::string reactor_name2,fit_mean_s,fit_mean_err_s,fit_sigma_s,fit_sigma_err_s;
    ULong64_t line_no = 0, count = 0;

    // read until end of file.
    while(in.good()){
        std::getline(in,reactor_name2,',');
        std::getline(in,fit_mean_s,',');
        std::getline(in,fit_mean_err_s,',');
        std::getline(in,fit_sigma_s,',');
        std::getline(in,fit_sigma_err_s,'\n');

        if (line_no>0){ //skip csv header
            if (strcmp(reactor_name.c_str(),reactor_name2.c_str())==0) {
                fit_mean = atof(fit_mean_s.c_str());
                fit_mean_err = atof(fit_mean_err_s.c_str());
                fit_sigma = atof(fit_sigma_s.c_str());
                fit_sigma_err = atof(fit_sigma_err_s.c_str());
                count++;
                //std::cout << "v: reactor_name: " << reactor_names[line_no-1] << ", distance: " << distances[line_no-1] << ", reactor_type: " << reactor_types[line_no-1] << ", n_core: " << n_cores[line_no-1] << ", power: " << powers[line_no-1] << ", power_err: " << power_errs[line_no-1] << std::endl; //debug check ('-1' for header)
            }
        }
        line_no++;
    }
    in.close();

    // if no reactor is found, or if more than one reactor of the same name is found, show an error
    if (count>1)
        printf("Error: More than 1 of the reactor found, %s (times: %llu)\n", reactor_name.c_str(), count);
    if (count==0)
        printf("Error: Reactor not found in constraint info file, %s\n", reactor_name.c_str());
    //if (count==1)
    //    printf("Reactor found %s\n", reactor_name.c_str());

    // print out read info
    //printf("reactor_name:%s, fit_mean: %.3f, fit_mean_err: %.3f, fit_sigma: %.3f, fit_sigma_err: %.3f\n", reactor_name.c_str(), fit_mean, fit_mean_err, fit_sigma, fit_sigma_err);
}

void readInfoFile(const std::string &runInfoFileName, std::vector<std::string> &reactor_names, std::vector<Double_t> &distances, std::vector<std::string> &reactor_types, std::vector<ULong64_t> &n_cores, std::vector<Double_t> &powers, std::vector<Double_t> &power_errs ) {
    //
    // Read couchDB run-info text file
    //

    std::ifstream in;
    in.open(runInfoFileName.c_str());
    //std::cout << "opening file: " << runInfoFileName.c_str() << std::endl;

    std::fill(reactor_names.begin(), reactor_names.end(), "");
    std::fill(distances.begin(), distances.end(), 0.);
    std::fill(reactor_types.begin(), reactor_types.end(), "");
    std::fill(n_cores.begin(), n_cores.end(), 0);
    std::fill(powers.begin(), powers.end(), 0.);
    std::fill(power_errs.begin(), power_errs.end(), 0.);

    std::string reactor_name,distance,reactor_type,n_core,power,power_err;
    ULong64_t line_no = 0;

    // read until end of file.
    while(in.good()){
        std::getline(in,reactor_name,',');
        std::getline(in,distance,',');
        std::getline(in,reactor_type,',');
        std::getline(in,n_core,',');
        std::getline(in,power,',');
        std::getline(in,power_err,'\n');

        if (line_no>0){ //skip csv header
            if (strcmp(reactor_name.c_str(),"")!=0) {
                reactor_names.push_back(reactor_name);
                distances.push_back(atof(distance.c_str()));
                reactor_types.push_back(reactor_type.c_str());
                n_cores.push_back(atoi(n_core.c_str()));
                powers.push_back(atof(power.c_str()));
                power_errs.push_back(atof(power_err.c_str()));

                //std::cout << "v: reactor_name: " << reactor_names[line_no-1] << ", distance: " << distances[line_no-1] << ", reactor_type: " << reactor_types[line_no-1] << ", n_core: " << n_cores[line_no-1] << ", power: " << powers[line_no-1] << ", power_err: " << power_errs[line_no-1] << std::endl; //debug check ('-1' for header)
            }
        }
        line_no++;
    }
    in.close();

    // print out read info
    //for (size_t i=0; i<(size_t)reactor_names.size(); i++)
    //    printf("i:%llu,reactor_names[i]:%s, distance: %.5f, type: %s, n_cores: %llu, power: %.5f, power_err: %.5f \n", i, reactor_names[i].c_str(), distances[i], reactor_types[i].c_str(), n_cores[i], powers[i], power_errs[i]);
}

void readParameterFile(const std::string &runParameterFileName, std::vector<Double_t> &d_21s, std::vector<Double_t> &s_12s, std::vector<Double_t> &s_13s) {
    // Read text file containing oscillation parameters
    std::ifstream in;
    in.open(runParameterFileName.c_str());
    //std::cout << "opening file: " << runParameterFileName.c_str() << std::endl;

    std::fill(d_21s.begin(), d_21s.end(), 0.);
    std::fill(s_12s.begin(), s_12s.end(), 0.);
    std::fill(s_13s.begin(), s_13s.end(), 0.);

    std::string d_21, s_12, s_13;
    ULong64_t line_no = 0;

    // read until end of file.
    // format of file: 'd_21,s_12,s_13\n'
    while(in.good()){
        std::getline(in,d_21,',');
        std::getline(in,s_12,',');
        std::getline(in,s_13,'\n');

        if (line_no>0){ //skip csv header
            if (strcmp(d_21.c_str(),"")!=0) {
                d_21s.push_back(atof(d_21.c_str()));
                s_12s.push_back(atof(s_12.c_str()));
                s_13s.push_back(atof(s_13.c_str()));

                //std::cout << "v: d_21: " << d_21s[line_no-1] << ", s_12: " << s_12s[line_no-1] << ", s_13: " << s_13s[line_no-1] << std::endl; //debug check ('-1' for header)
            }
        }
        line_no++;
    }
    in.close();

    // // print out read parameters
    // for (size_t i=0; i<(size_t)d_21s.size(); i++)
        // printf("i:%llu, d_21:%.5f, s_12:%.5f, s_13:%.5f\n", i, d_21s[i], s_12s[i], s_13s[i]);
}

BinnedED LHFit_initialise(BinnedED **spectra_ke_pdf, BinnedED **spectra_ev_pdf, const std::string &in_path, const std::string &data_path, ULong64_t flux_data){
    //
    // Load pdf's for reactor types
    // Load data pdf
    //

    printf("Begin init--------------------------------------\n");
    printf("LHFit_initialise...\n");

    char name[1000];
    ULong64_t flux_mc = 100;

    // setup spectra filenames
    const std::string pwr_file = std::string("combinedpwr");//std::string("ohi");//std::string("combinedpwr");
    const std::string phwr_file = std::string("combinedpwr");//std::string("ohi");//std::string("combinedphwr");

    // ke filenames
    printf("Using unoscillated spectrum for:\n");
    sprintf(name, "%sflux%llu/%s_flux%llu_day360_pass1_cleanround1_ke_oxsx.root", in_path.c_str(), flux_mc, pwr_file.c_str(), flux_mc);
    printf("\ttype PWR & BWR: %s\n", name);
    const std::string pwr_unosc_ke_path = std::string(name);
    sprintf(name, "%sflux%llu/%s_flux%llu_day360_pass1_cleanround1_ke_oxsx.root", in_path.c_str(), flux_mc, phwr_file.c_str(), flux_mc);
    printf("\ttype PHWR: %s\n", name);
    const std::string phwr_unosc_ke_path = std::string(name);

    // ev filename
    sprintf(name, "%sflux%llu/%s_flux%llu_day360_pass1_cleanround1_prompt_oxsx.root", in_path.c_str(), flux_mc, pwr_file.c_str(), flux_mc);
    printf("\ttype PWR & BWR: %s\n", name);
    const std::string pwr_unosc_ev_path = std::string(name);
    sprintf(name, "%sflux%llu/%s_flux%llu_day360_pass1_cleanround1_prompt_oxsx.root", in_path.c_str(), flux_mc, phwr_file.c_str(), flux_mc);
    printf("\ttype PHWR: %s\n", name);
    const std::string phwr_unosc_ev_path = std::string(name);

    // setup (oscillated) data filename
    sprintf(name, "%s", data_path.c_str());
    printf("Loading data spectrum: %s\n\n", name);

    // setup ntuple
    ObsSet data_rep(0);

    // set up binning
    AxisCollection axes;
    Double_t e_min = 0.425*2; //2; //*6 for kamland paper #1, *2 for paper #2
    Double_t e_max = 0.425*19; //8;
    Int_t n_bins = 17; //(8-2)*10; //13 for paper #1, 17 for paper #2
    axes.AddAxis(BinAxis("mc_neutrino_energy", e_min, e_max, n_bins));

    // load (oscillated) data ntuple
    BinnedED data_set_pdf("dataSetPdf", axes);
    data_set_pdf.SetObservables(data_rep);
    ROOTNtuple data_ntp(data_path, "nt");
    for(ULong64_t i = 0; i < data_ntp.GetNEntries(); i++)
        data_set_pdf.Fill(data_ntp.GetEntry(i));

    // for kamland reproduction change:
    data_set_pdf.Scale(1./flux_data);
    data_set_pdf.Normalise();
    data_set_pdf.Scale(1.27844983315467830e+003);

    printf("Loading data: %s (osc)integral:%.3f\n", data_path.c_str(), data_set_pdf.Integral());


    // load unoscillated spectra
    // pwr spectrum
    // KE branch
    ROOTNtuple pwr_spectrum_ke_ntp(pwr_unosc_ke_path, "nt");
    sprintf(name, "pwr_spectrum_ke_pdf");
    spectra_ke_pdf[0] = new BinnedED(name, axes);
    spectra_ke_pdf[0]->SetObservables(0);
    for(ULong64_t i = 0; i < pwr_spectrum_ke_ntp.GetNEntries(); i++)
        spectra_ke_pdf[0]->Fill(pwr_spectrum_ke_ntp.GetEntry(i));
    spectra_ke_pdf[0]->Scale(1./flux_mc);
    spectra_ke_pdf[0]->Normalise();

    // EV branch
    ROOTNtuple pwr_spectrum_ev_ntp(pwr_unosc_ev_path, "nt");
    sprintf(name, "pwr_spectrum_ev_pdf");
    spectra_ev_pdf[0] = new BinnedED(name, axes);
    spectra_ev_pdf[0]->SetObservables(0);
    for(ULong64_t i = 0; i < pwr_spectrum_ev_ntp.GetNEntries(); i++)
        spectra_ev_pdf[0]->Fill(pwr_spectrum_ev_ntp.GetEntry(i));
    spectra_ev_pdf[0]->Scale(1./flux_mc);
    spectra_ev_pdf[0]->Normalise();

    // phwr spectrum
    ROOTNtuple phwr_spectrum_ke_ntp(phwr_unosc_ke_path, "nt");
    sprintf(name, "phwr_spectrum_ke_pdf");
    spectra_ke_pdf[1] = new BinnedED(name, axes);
    spectra_ke_pdf[1]->SetObservables(0);
    for(ULong64_t i = 0; i < phwr_spectrum_ke_ntp.GetNEntries(); i++)
        spectra_ke_pdf[1]->Fill(phwr_spectrum_ke_ntp.GetEntry(i));
    spectra_ke_pdf[1]->Scale(1./flux_mc);
    spectra_ke_pdf[1]->Normalise();

    // phwr spectrum
    ROOTNtuple phwr_spectrum_ev_ntp(phwr_unosc_ev_path, "nt");
    sprintf(name, "phwr_spectrum_ev_pdf");
    spectra_ev_pdf[1] = new BinnedED(name, axes);
    spectra_ev_pdf[1]->SetObservables(0);
    for(ULong64_t i = 0; i < phwr_spectrum_ev_ntp.GetNEntries(); i++)
        spectra_ev_pdf[1]->Fill(phwr_spectrum_ev_ntp.GetEntry(i));
    spectra_ev_pdf[1]->Scale(1./flux_mc);
    spectra_ev_pdf[1]->Normalise();

    printf("End init--------------------------------------\n");
    return data_set_pdf;
}
