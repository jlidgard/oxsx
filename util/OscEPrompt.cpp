#include <TFile.h>
#include <TMath.h>
#include <string>
#include <TNtuple.h>
#include <iostream>
#include <TObject.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <sstream>

double NuSurvProb(double nuE, double baseline, double delmsqr21, double sinsqrtheta12, double sinsqrtheta13) {
  double fSSqr2Theta12 = pow(sin(2.0 * TMath::ASin(sqrt(sinsqrtheta12))), 2.0);
  double fS4 = pow(sinsqrtheta13, 2.0);
  double fC4 = pow(1.0-sinsqrtheta13, 2.0);
  double scale = 1.267e3; // for nuE in [MeV] and baseline in [km]
  double sSqrDmBE = pow(sin(scale * delmsqr21 * baseline / nuE), 2.0);
  double fOscProb = (fC4 * (1.0 - fSSqr2Theta12 * sSqrDmBE) + fS4);
  return fOscProb;
}

bool applyosc;

void OscPromptE_Evindex (const std::string infile, const std::string outfile, int nhit1Min, int nhit1Max, int nhit2Min, int nhit2Max, double E1Min, double E1Max, double E2Min, double E2Max, double deltaT,  double delmsqr21, double sinsqrtheta12, double sinsqrtheta13) {
  TFile *fin = TFile::Open(infile.c_str());
  TNtuple *T = (TNtuple *) fin->Get("nt");
 
  std::stringstream sout;
  if (applyosc)
    sout<<outfile<<"ds21_"<<delmsqr21<<"_ss12_"<<sinsqrtheta12<<"_ss13_"<<sinsqrtheta13<<".root";
  else
    sout<<outfile<<".root";
 
  float ReactorDistance;
  float MCentryb4, MCentry,nextMCentry;
  float days,nextdays;
  float secs;
  float nsecs;
  float nhit,nextnhit;
  float evindex, nextevindex;
  double timeD, DeltaR;
  float Fitted, nextFitted;
  float Energy,nextEnergy,parke,part1ke,part2ke;
  float posX,posY,posZ,nextposX,nextposY,nextposZ;
  float EVcount = 0;
  //TVector3 VecR,nextVecR;
  double VecRMag,nextVecRMag;
  double time, nexttime;
  
  int mccount = 0;
  int numsimmed = 1;
  int EV01 = 0;
  int EV02 = 0;
  int EV03 = 0;
  
  T->SetBranchAddress("MCentry", &MCentry);
  T->SetBranchAddress("Evindex", &evindex);
  T->SetBranchAddress("ParKE", &parke);
  T->SetBranchAddress("Part1KE", &part1ke);
  T->SetBranchAddress("Part2KE", &part2ke);
  T->SetBranchAddress("Days", &days);
  T->SetBranchAddress("Seconds", &secs);
  T->SetBranchAddress("Nanoseconds", &nsecs);
  T->SetBranchAddress("nhits", &nhit);
  T->SetBranchAddress("Energy", &Energy);
  T->SetBranchAddress("posx", &posX);
  T->SetBranchAddress("posy", &posY);
  T->SetBranchAddress("posz", &posZ);
  T->SetBranchAddress("posz", &posZ);
  T->SetBranchAddress("Fitted", &Fitted);
  T->SetBranchAddress("ReactorDistance", &ReactorDistance);

  TH1D deltaTimeEVindex("deltaTEVindex","deltaTEVindex (ns)",100,0,1000000);
  TH1D deltaREVindex("deltaREVindex","deltaREVindex (mm)",300,0,3000);
  TH1D E1EVindex("E1EVindex", "E1EVindex (MeV)", 300, 0, 9);
  TH1D E2EVindex("E2EVindex", "E2EVindex (MeV)", 300, 0, 9);
  TH1D Nhit1("Nhit1", "Nhit1", 350, 0, 3500);
  TH1D Nhit2("Nhit2", "Nhit2", 200, 0, 2000);
  TH1D KEPar("KEPar", "KEPar", 100, 0, 8);
  TH1D KEPart1("KEPart1", "KEPart1 (MeV)", 100, 0, 8);
  TH1D KEPart2("KEPart2", "KEPart2 (MeV)", 200, 0, 2);
  TH1D nubarKE_posKE("nubarKE_posKE", "nubarKE_posKE (MeV)", 300, 0, 3);
  TH1D nubarKE_E1("nubarKE_E1", "nubarKE_E1 (MeV)", 300, 0, 3);

  TH1D deltaTime("deltaT","deltaT (ns)",100,0,1000000);
  TH1D deltaR("deltaR","deltaR (mm)",300,0,3000);
  TH1D EPrompt("EPrompt", "EPrompt (MeV)", 300, 0, 9);
  TH1D ELate("E2Late", "E2Late (MeV)", 300, 0, 9);
  TH1D NhitPrompt("NhitPrompt", "NhitPrompt", 350, 0, 3500);
  TH1D NhitLate("NhitLate", "NhitLate", 200, 0, 2000);
  
  TH2D h2("antinuKEvsEPrompt","antinuKEvsEPrompt",300,0,9,300,0,9);

  double survprob;
  TRandom3 *r1 = new TRandom3();
  r1->SetSeed(0);    
  
  ///////////////////////////////////
  ///////   Using EVindex   /////////
  ///////////////////////////////////
  
  //get first event
  T->GetEntry(0);
  MCentryb4 = MCentry;  
  //go through entries and collect triggered evs which correspond to a MC event.
  for(int i = 1; i < T->GetEntries(); i++){
  //for(int i = 1; i < 50; i++){
    T->GetEntry(i);
    //mark when passed group of evs with same mcindex:
    if(abs(MCentry - MCentryb4)>= 1 && evindex == 0){ 
      numsimmed += 1;
      T->GetEntry(i-1);
      EVcount = evindex + 1; //Get Ev count
      
      size_t j = 0;  // j index for potential positron event
      bool pairfound = false;
      //move j index within MC index group, move onto next group when j reaches second to last event
      while (j < EVcount - 1){
	int k = 1; // k index for potential neutron event
	float time_diff_condition = 0;
	//move j index within MC index group, move on to next group if nothing within deltaT or no more events to look at
	while((time_diff_condition < deltaT) && (k < EVcount - j)){
	  T->GetEntry(i-EVcount+j+k); // neutron event
	  nextdays = days;
	  nextMCentry = MCentry;
	  nextFitted = Fitted;
	  nexttime = nsecs + (secs * pow(10, 9)) + (days * 24 * 3600 * pow(10, 9));
	  nextnhit = nhit;
	  nextevindex = evindex;
	  nextEnergy = Energy;
	  nextposX = posX;
	  nextposY = posY;
	  nextposZ = posZ;
	  nextVecRMag = sqrt((pow(posX,2))+(pow(posY,2))+(pow(posZ,2)));
	  //nextVecR = TVector3(posX,posY,posZ);
	  
	  T->GetEntry(i-EVcount+j);
	  bool goodpair = false;
	  timeD = 0;
	  if (Fitted > 0 && nextFitted > 0){
	    //VecR = TVector3(posX,posY,posZ);
	    VecRMag = sqrt((pow(posX,2))+(pow(posY,2))+(pow(posZ,2)));
	    time = nsecs + (secs * pow(10, 9)) + (days * 24 * 3600 * pow(10, 9));
	    timeD = std::fabs(nexttime - time);
	    //DeltaR = (VecR-nextVecR).Mag();
	    DeltaR = sqrt((pow(posX-nextposX,2))+(pow(posY-nextposY,2))+(pow(posZ-nextposZ,2)));
	    if(timeD > 400 && timeD < deltaT){
	      if (VecRMag < 5700){
		if (nextVecRMag < 5700){
		  if (DeltaR < 1500){
		    if (E1Min <= Energy && Energy <= E1Max){
		      if (E2Min <= nextEnergy &&nextEnergy <= E2Max){
			if (nhit1Min <= nhit  && nhit <= nhit1Max){
			  if (nhit2Min <= nextnhit && nextnhit <= nhit2Max){
			    if (applyosc)
			      survprob = NuSurvProb(parke,ReactorDistance,delmsqr21,sinsqrtheta12,sinsqrtheta13);
			    else
			      survprob = 100;
			    if (survprob > r1->Rndm()){
			      mccount += 1;
			      if (evindex == 0 && nextevindex == 1)
				EV01 += 1;
			      if (evindex == 0 && nextevindex == 2)
				EV02 += 1;
			      if (evindex == 0 && nextevindex == 3)
				EV03 += 1;
			      if (evindex == 0)
				h2.Fill(Energy,parke);
			      deltaTimeEVindex.Fill(timeD);
			      deltaREVindex.Fill(DeltaR);
			      E1EVindex.Fill(Energy);
			      E2EVindex.Fill(nextEnergy);
			      Nhit1.Fill(nhit);
			      Nhit2.Fill(nextnhit);
			      KEPar.Fill(parke);
			      KEPart1.Fill(part1ke);
			      KEPart2.Fill(part2ke);
			      nubarKE_posKE.Fill(parke-part1ke);
			      nubarKE_E1.Fill(parke-Energy);
			      
			      goodpair = true;  //pair passed quality cuts, if not continue search in group
			      pairfound = true; // a pair was found
			      k += 100; // cancel sub search
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  } // if pair failed quality cuts, continue sub group search
	  if (!goodpair){
	    k += 1;
	    if (Fitted <= 0){ //if positron/j_index not a valid fit, move onto next sub group
	      k += 100;
	    }else{ //else if good, check that it satisfies time diff, if not move onto next sub group
	      time_diff_condition = timeD;
	    }
	  }
	}
	if (pairfound){
	  j += 100;
	}else{
	  j += 1;
	}
      } 
    }
    MCentryb4 = MCentry;
  }
  
  int lastentry = (T->GetEntries() - 1);
  T->GetEntry(lastentry);
  EVcount = evindex + 1; //Get Ev count
      
  size_t j = 0;  // j index for potential positron event
  bool pairfound = false;
  //move j index within MC index group, move onto next group when j reaches second to last event
  while (j < EVcount - 1){
    int k = 1; // k index for potential neutron event
    float time_diff_condition = 0;
    //move j index within MC index group, move on to next group if nothing within deltaT or no more events to look at
    while((time_diff_condition < deltaT) && (k < EVcount - j)){
      T->GetEntry(lastentry-EVcount+j+k+1); // neutron event
      nextFitted = Fitted;
      nexttime = nsecs + (secs * pow(10, 9)) + (days * 24 * 3600 * pow(10, 9));
      nextnhit = nhit;
      nextevindex = evindex;
      nextEnergy = Energy;
      //nextVecR = TVector3(posX,posY,posZ);
      nextposX = posX;
      nextposY = posY;
      nextposZ = posZ;
      nextVecRMag = sqrt((pow(posX,2))+(pow(posY,2))+(pow(posZ,2)));
      
      T->GetEntry(lastentry+1-EVcount+j);
      bool goodpair = false;
      if (Fitted > 0 && nextFitted > 0){
	//VecR = TVector3(posX,posY,posZ);
	VecRMag = sqrt((pow(posX,2))+(pow(posY,2))+(pow(posZ,2)));
	time = nsecs + (secs * pow(10, 9)) + (days * 24 * 3600 * pow(10, 9));
	timeD = std::fabs(nexttime - time);
	//DeltaR = (VecR-nextVecR).Mag();
	DeltaR = sqrt((pow(posX-nextposX,2))+(pow(posY-nextposY,2))+(pow(posZ-nextposZ,2)));
	if(timeD > 400 && timeD < deltaT){
	  if (VecRMag < 5700){
	    if (nextVecRMag < 5700){
	      if (DeltaR < 1500){
		if (E1Min <= Energy && Energy <= E1Max){
		  if (E2Min <= nextEnergy &&nextEnergy <= E2Max){
		    if (nhit1Min <= nhit  && nhit <= nhit1Max){
		      if (nhit2Min <= nextnhit && nextnhit <= nhit2Max){
			if (applyosc)
			  survprob = NuSurvProb(parke,ReactorDistance,delmsqr21,sinsqrtheta12,sinsqrtheta13);
			else
			  survprob = 100;
			if (survprob > r1->Rndm()){
			  mccount += 1;
			  
			  if (evindex == 0 && nextevindex == 1)
			    EV01 += 1;
			  if (evindex == 0 && nextevindex == 2)
			    EV02 += 1;
			  if (evindex == 0 && nextevindex == 3)
			    EV03 += 1;
			  if (evindex == 0)
			    h2.Fill(Energy,parke);
			  deltaTimeEVindex.Fill(timeD);
			  deltaREVindex.Fill(DeltaR);
			  E1EVindex.Fill(Energy);
			  E2EVindex.Fill(nextEnergy);
			  Nhit1.Fill(nhit);
			  Nhit2.Fill(nextnhit);
			  KEPar.Fill(parke);
			  KEPart1.Fill(part1ke);
			  KEPart2.Fill(part2ke);
			  nubarKE_posKE.Fill(parke-part1ke);
			  nubarKE_E1.Fill(parke-Energy);
		      	
			  goodpair = true;  //pair passed quality cuts, if not continue search in group
			  pairfound = true; // a pair was found
			  k += 100; // cancel sub search
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      } // if pair failed quality cuts, continue sub group search
      if (!goodpair){
	k += 1;
	if (Fitted <= 0){ //if positron/j_index not a valid fit, move onto next sub group
	  k += 100;
	}else{ //else if good, check that it satisfies time diff, if not move onto next sub group
	  time_diff_condition = timeD;
	}
      }
    }
    if (pairfound){
      j += 100;
    }else{
      j += 1;
    }
  } 
  
  std::cout<<" EVindex 0_1: "<<EV01<<std::endl;
  std::cout<<" EVindex 0_2: "<<EV02<<std::endl;
  std::cout<<" EVindex 0_3: "<<EV03<<std::endl;
  
  std::cout<<"\n number of antinus simmed: "<<numsimmed<<std::endl;
  std::cout<<"MC partner events within deltaT: "<<mccount<<std::endl;

  std::stringstream ssout;
  ssout<<outfile<<"ds21_"<<delmsqr21<<"_ss12_"<<sinsqrtheta12<<"_ss13_"<<sinsqrtheta13<<".root";
  TFile fout(ssout.str().c_str(), "RECREATE");

  deltaTimeEVindex.Write();
  deltaREVindex.Write();
  E1EVindex.Write();
  E2EVindex.Write();
  Nhit1.Write();
  Nhit2.Write();
  KEPar.Write();
  KEPart1.Write();
  KEPart2.Write();
  nubarKE_posKE.Write();
  nubarKE_E1.Write();

  h2.Write();
  
  fout.Close();
}

void OscPrunedNtuple (const std::string infile, const std::string outfile,  double delmsqr21, double sinsqrtheta12, double sinsqrtheta13) {
  
  TFile *fin = TFile::Open(infile.c_str());
  TNtuple *T = (TNtuple *) fin->Get("nt");
  
  float parke,ReactorDistance,MCentry,MCentryb4,evindex;
  
  T->SetBranchAddress("MCentry", &MCentry);
  T->SetBranchAddress("Evindex", &evindex);
  T->SetBranchAddress("ParKE", &parke);
  T->SetBranchAddress("ReactorDistance", &ReactorDistance);

  std::stringstream ssout;
  ssout<<outfile<<"ds21_"<<delmsqr21<<"_ss12_"<<sinsqrtheta12<<"_ss13_"<<sinsqrtheta13<<".root";
  TFile *fout = new TFile(ssout.str().c_str(),"RECREATE");
  TTree *outT = T->CloneTree(0);
  
  ///////////////////////////////////
  /////   Osc using EVindex   ///////
  ///////////////////////////////////
  double survprob;
  TRandom3 *r1 = new TRandom3();
  r1->SetSeed(0);    
  
  bool survived = false;
  bool survivedb4 = false;

  //get first event
  T->GetEntry(0);
  MCentryb4 = MCentry;  

  survprob = NuSurvProb(parke, ReactorDistance, delmsqr21, sinsqrtheta12, sinsqrtheta13);  
  MCentryb4 = MCentry;

  if (survprob > r1->Rndm()){
    outT->Fill();
    survivedb4 = true;
  }else
    survivedb4 = false;
  
  //go through entries and collect triggered evs which correspond to a MC event.
  for(int i = 1; i < T->GetEntries(); i++){
    T->GetEntry(i);
    
    //if(abs(MCentry - MCentryb4)>= 1 && evindex == 0){ 
    if (MCentry != MCentryb4){
      survprob = NuSurvProb(parke, ReactorDistance, delmsqr21, sinsqrtheta12, sinsqrtheta13);  
      if (survprob > r1->Rndm()){
	survived = true;
	outT->Fill();
      }else
	survived = false; 	
      
      survivedb4 = survived;
    }
    else if (MCentry == MCentryb4 && survivedb4)
      outT->Fill();
    
    MCentryb4 = MCentry;
  }
  
  outT->Print();
  outT->AutoSave();
  fin->Close();
  fout->Close();
}


void PromptE_Data (const std::string infile, const std::string outfile, int nhit1Min, int nhit1Max, int nhit2Min, int nhit2Max, double E1Min, double E1Max, double E2Min, double E2Max, double deltaT,  double delmsqr21, double sinsqrtheta12, double sinsqrtheta13) {
  
  std::stringstream ssout;

  if (applyosc){
    OscPrunedNtuple(infile,outfile,delmsqr21,sinsqrtheta12, sinsqrtheta13);  
    
    ssout<<outfile<<"ds21_"<<delmsqr21<<"_ss12_"<<sinsqrtheta12<<"_ss13_"<<sinsqrtheta13<<".root"; 
  }
  else{
    ssout<<infile<<".root"; 
  }
  TFile *fin = TFile::Open(ssout.str().c_str());
  TNtuple *T = (TNtuple *) fin->Get("nt");

  float ReactorDistance;
  float days,nextdays;
  float secs;
  float nsecs;
  float nhit,nextnhit;
  double timeD, DeltaR;
  float Fitted, nextFitted;
  float Energy,nextEnergy,parke,part1ke,part2ke;
  float posX,posY,posZ,nextposX,nextposY,nextposZ;
  //TVector3 VecR,nextVecR;
  double VecRMag,nextVecRMag;
  double time, nexttime;
  
  int evcount = 0;
  
  T->SetBranchAddress("ParKE", &parke);
  T->SetBranchAddress("Part1KE", &part1ke);
  T->SetBranchAddress("Part2KE", &part2ke);
  T->SetBranchAddress("Days", &days);
  T->SetBranchAddress("Seconds", &secs);
  T->SetBranchAddress("Nanoseconds", &nsecs);
  T->SetBranchAddress("nhits", &nhit);
  T->SetBranchAddress("Energy", &Energy);
  T->SetBranchAddress("posx", &posX);
  T->SetBranchAddress("posy", &posY);
  T->SetBranchAddress("posz", &posZ);
  T->SetBranchAddress("Fitted", &Fitted);
  T->SetBranchAddress("ReactorDistance", &ReactorDistance);

  TH1D deltaTime("deltaT","deltaT (ns)",100,0,1000000);
  TH1D deltaR("deltaR","deltaR (mm)",300,0,3000);
  TH1D EPrompt("EPrompt", "EPrompt (MeV)", 300, 0, 9);
  TH1D ELate("E2Late", "E2Late (MeV)", 300, 0, 9);
  TH1D NhitPrompt("NhitPrompt", "NhitPrompt", 350, 0, 3500);
  TH1D NhitLate("NhitLate", "NhitLate", 200, 0, 2000);
  
  std::cout<<"\n \n Data: "<<std::endl;
  ///////////////////////////////////
  ///////      For Data     /////////
  ///////////////////////////////////
  int i = 0;
  std::cout<<"numentries: "<<T->GetEntries()<<std::endl;
  while (i < T->GetEntries()-1){
    int j = 1;
    float time_diff_condition = 0;
    bool pairfound = false;
    while(time_diff_condition < deltaT){
      T->GetEntry(i+j);
      nextdays = days;
      nextFitted = Fitted;
      nexttime = nsecs + (secs * pow(10, 9)) + (days * 24 * 3600 * pow(10, 9));
      nextnhit = nhit;
      nextEnergy = Energy;
      //nextVecR = TVector3(posX,posY,posZ);
      nextposX = posX;
      nextposY = posY;
      nextposZ = posZ;
      nextVecRMag = sqrt((pow(posX,2))+(pow(posY,2))+(pow(posZ,2)));
      
      T->GetEntry(i);
      bool pairgood = false;
      timeD = 0;
      if (Fitted > 0 && nextFitted > 0){
	//VecR = TVector3(posX,posY,posZ);
	VecRMag = sqrt((pow(posX,2))+(pow(posY,2))+(pow(posZ,2)));
	time = nsecs + (secs * pow(10, 9)) + (days * 24 * 3600 * pow(10, 9));
	timeD = std::fabs(nexttime - time);
	//DeltaR = (VecR-nextVecR).Mag();
	DeltaR = sqrt((pow(posX-nextposX,2))+(pow(posY-nextposY,2))+(pow(posZ-nextposZ,2)));
	if(timeD > 400 && timeD < deltaT){
	  if (VecRMag < 5700){
	    if (nextVecRMag < 5700){
	      if (DeltaR < 1500){
	      	if (E1Min <= Energy && Energy <= E1Max){
		  if (E2Min <= nextEnergy &&nextEnergy <= E2Max){
		    if (nhit1Min <= nhit  && nhit <= nhit1Max){
		      if (nhit2Min <= nextnhit && nextnhit <= nhit2Max){
			evcount += 1;
			deltaTime.Fill(timeD);
			deltaR.Fill(DeltaR);
			EPrompt.Fill(Energy);
			ELate.Fill(nextEnergy);
			NhitPrompt.Fill(nhit);
			NhitLate.Fill(nextnhit);
			
			pairgood = true;
			pairfound = true;
			time_diff_condition = deltaT + 1e10;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      if (!pairgood){
	j += 1;
	if (Fitted <= 0){
	  pairfound = false;
	  time_diff_condition = deltaT + 1e10;
	}
	else{
	  time_diff_condition = timeD;
	}
	
      }
    }
    
    if (pairfound)
      i += j+1;
    else
      i += 1; 
  }
  
  std::cout<<"EV partner events within deltaT: "<<evcount<<std::endl;
  
  TFile fout(ssout.str().c_str(), "RECREATE");

  deltaTime.Write();
  deltaR.Write();
  EPrompt.Write();
  ELate.Write();
  NhitPrompt.Write();
  NhitLate.Write();
  
  fout.Close();
}

int main(int argc, char *argv[])
{
  if (argc < 16){
    std::cout<<"15 (12) arguments expected: \n 1: input ntuple \n 2: \
output file with hists  (UP TO .root !!!) \n  \n 3: nhit1Min \n 4: nhit1Max\n 5: nhit2Min \n 6: nhit2Max\n 7: E1Min \n 8: E1Max\n 9: E2Min \n 10: E2Max \n 11: deltaT (in ns) \n (12: Data or MC) \n 12: delmsqr21 \n 13: sinsqrtheta12 \n 14:sinsqrtheta13 \n \n 15: Data or MC "<<std::endl;
  }else{
    //std::stringstream argParser;
    const std::string &infile = argv[1];
    const std::string &outfile = argv[2];
    int nhit1Min = atoi(argv[3]);
    int nhit1Max = atoi(argv[4]);
    int nhit2Min = atoi(argv[5]);
    int nhit2Max = atoi(argv[6]);
    double E1Min = atof(argv[7]);
    double E1Max = atof(argv[8]);
    double E2Min = atof(argv[9]);
    double E2Max = atof(argv[10]);
    double deltaT = atof(argv[11]);
    
    double delmsqr21, sinsqrtheta12, sinsqrtheta13;
    std::string DataMC;
    if (argc == 13){
      delmsqr21 = 0;
      sinsqrtheta12 = 0;
      sinsqrtheta13 = 0;

      DataMC = argv[12];
    }
    else{
      delmsqr21 = atof(argv[12]);
      sinsqrtheta12 = atof(argv[13]);
      sinsqrtheta13 = atof(argv[14]);
      
      DataMC= argv[15];
    }
    
    double PromptRmax = 5700;
    double LateRmax = 5700;

    if (DataMC == "MC")
      OscPromptE_Evindex(infile,outfile,nhit1Min,nhit1Max,nhit2Min,nhit2Max,E1Min,E1Max,E2Min,E2Max,deltaT,delmsqr21,sinsqrtheta12,sinsqrtheta13);
    else if (DataMC == "Data")
      PromptE_Data(infile,outfile,nhit1Min,nhit1Max,nhit2Min,nhit2Max,E1Min,E1Max,E2Min,E2Max,deltaT,delmsqr21,sinsqrtheta12,sinsqrtheta13);
    else
      std::cout<<"\n Looking at 'Data' or 'MC'";
  }
}
