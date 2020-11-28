#include <TFile.h>
#include <TMath.h>
#include <string>
#include <TTree.h>
#include <TVector3.h>
#include <iostream>
#include <TObject.h>
#include <math.h>
#include <TNtuple.h>
#include <TChain.h>

void process_cuts(const std::string filename_input, const std::string filename_output){

  // load input file
  char *name = new char[1000];
  sprintf(name, "%s/output", filename_input.c_str());
  TChain *tree_input = new TChain("output");
  tree_input->Add(name);

  TFile *file_output = new TFile(filename_output.c_str(), "RECREATE");
  TTree *tree_output = new TTree("nt","Tagged positron + neutron events");

  Int_t mc_entry;
  Double_t ev_pos_x, ev_pos_y, ev_pos_z;
  Double_t ev_energy, ev_energy_p1=0, ev_energy_p2=0;
  Int_t ev_time_seconds, ev_time_days, ev_time_nanoseconds;
  Int_t ev_nhit;
  Bool_t ev_validity;
  Int_t ev_index, ev_index_p1=0, ev_index_p2=0;
  TString *reactor_core_name=0;
  Double_t mc_quench_i, mc_energy_nu, mc_energy_n, mc_energy_ep;
  ULong64_t ev_clock50;

  // set branches input ntuple
  tree_input->SetBranchAddress("mcIndex", &mc_entry);
  tree_input->SetBranchAddress("mcEdepQuenched", &mc_quench_i);
  tree_input->SetBranchAddress("parentKE1", &mc_energy_nu);
  tree_input->SetBranchAddress("mcke1", &mc_energy_ep);
  tree_input->SetBranchAddress("mcke2", &mc_energy_n);
  tree_input->SetBranchAddress("evIndex", &ev_index);
  tree_input->SetBranchAddress("energy", &ev_energy);
  tree_input->SetBranchAddress("fitValid", &ev_validity);
  tree_input->SetBranchAddress("posx", &ev_pos_x);
  tree_input->SetBranchAddress("posy", &ev_pos_y);
  tree_input->SetBranchAddress("posz", &ev_pos_z);
  tree_input->SetBranchAddress("nhits", &ev_nhit);
  tree_input->SetBranchAddress("uTDays", &ev_time_days);
  tree_input->SetBranchAddress("uTSecs", &ev_time_seconds);
  tree_input->SetBranchAddress("uTNSecs", &ev_time_nanoseconds);
  tree_input->SetBranchAddress("parentMeta1", &reactor_core_name);
  tree_input->SetBranchAddress("clockCount50", &ev_clock50);

  // set branches output pruned ntuple
  tree_output->Branch("entry", &mc_entry);
  tree_output->Branch("mc_quench", &mc_quench_i);
  tree_output->Branch("mc_neutrino_energy", &mc_energy_nu);
  tree_output->Branch("mc_positron_energy", &mc_energy_ep);
  tree_output->Branch("mc_neutron_energy", &mc_energy_n);
  tree_output->Branch("ev_index", &ev_index);
  tree_output->Branch("ev_fit_energy", &ev_energy);
  tree_output->Branch("ev_fit_validity", &ev_validity);
  tree_output->Branch("ev_fit_position_x", &ev_pos_x);
  tree_output->Branch("ev_fit_position_y", &ev_pos_y);
  tree_output->Branch("ev_fit_position_z", &ev_pos_z);
  tree_output->Branch("ev_nhit", &ev_nhit);
  tree_output->Branch("ev_time_days", &ev_time_days);
  tree_output->Branch("ev_time_seconds", &ev_time_seconds);
  tree_output->Branch("ev_time_nanoseconds", &ev_time_nanoseconds);
  tree_output->Branch("reactor_core_name", &reactor_core_name);
  tree_output->Branch("ev_fit_energy_p1", &ev_energy_p1);
  tree_output->Branch("ev_fit_energy_p2", &ev_energy_p2);
  tree_output->Branch("ev_index_p1", &ev_index_p1);
  tree_output->Branch("ev_index_p2", &ev_index_p2);

  ULong64_t n_entries = tree_input->GetEntries();
  if (n_entries != 0){

    for(ULong64_t i = 1; i < n_entries; i++){
      tree_input->GetEntry(i);
	  
	  tree_output->Fill(); // use to fill ntuple
	  
	  printf("core name: %s\n", reactor_core_name->Data()); // print
	
	}
  }
  tree_output->AutoSave();
  file_output->Close();
}

Int_t main(Int_t argc, char *argv[]) {

    if (argc != 3) {
        std::cout<<"Error: 2 arguments expected. Got: "<<argc-1<<std::endl;
        return 1; // return>0 indicates error code
    }
    else {
        const std::string &filename_input = argv[1];
        const std::string &filename_output = argv[2];

        process_cuts(filename_input.c_str(), filename_output.c_str());

        return 0; // completed successfully
    }
}
