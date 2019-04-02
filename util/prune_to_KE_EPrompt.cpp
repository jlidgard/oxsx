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

void OscPromptE_Evindex (const std::string infile, const std::string fileout, int nhit1Min, int nhit1Max, int nhit2Min, int nhit2Max, double E1Min, double E1Max, double E2Min, double E2Max, double deltaT, double PromptRmax, double LateRmax, double delmsqr21, double sinsqrtheta12, double sinsqrtheta13, bool KEflag) {
  TFile *fin = TFile::Open(infile.c_str());
  TNtuple *T = (TNtuple *) fin->Get("nt");

  bool applyosc;
  std::stringstream sout;
  if (delmsqr21 > 0 || sinsqrtheta12 > 0 || sinsqrtheta13 >0){
    if (KEflag)
      sout<<fileout<<"KEds21_"<<delmsqr21<<"_ss12_"<<sinsqrtheta12<<"_ss13_"<<sinsqrtheta13<<"_oxsx.root";
    else
      sout<<fileout<<"E1ds21_"<<delmsqr21<<"_ss12_"<<sinsqrtheta12<<"_ss13_"<<sinsqrtheta13<<"_oxsx.root";

    applyosc = true;
  }
  else{
    if (KEflag)
      sout<<fileout<<"KE_oxsx.root";
    else
      sout<<fileout<<"E1_oxsx.root";

    applyosc = false;
  }
  TFile *fout = new TFile(sout.str().c_str(),"RECREATE");
  TNtuple* ntout;
  if (KEflag)
    ntout = new TNtuple("nt","nt","ParKE");
  else
    ntout = new TNtuple("nt","nt","E1");
  
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
  float posX,posY,posZ;
  float EVcount = 0;
  TVector3 VecR,nextVecR;
  double time, nexttime;
  
  int mccount = 0;
  int numsimmed = 1;
  int EV01 = 0;
  int EV02 = 0;
  int EV03 = 0;
  //int EV02 = 0;
  
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
  T->SetBranchAddress("Fitted", &Fitted);
  T->SetBranchAddress("ReactorDistance", &ReactorDistance);

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
	  nextVecR = TVector3(posX,posY,posZ);
	  
	  T->GetEntry(i-EVcount+j);
	  bool goodpair = false;
	  timeD = 0;
	  if (Fitted > 0 && nextFitted > 0){
	    VecR = TVector3(posX,posY,posZ);
	    time = nsecs + (secs * pow(10, 9)) + (days * 24 * 3600 * pow(10, 9));
	    timeD = std::fabs(nexttime - time);
	    DeltaR = (VecR-nextVecR).Mag();
	    if(timeD > 400 && timeD < deltaT){
	      if (VecR.Mag() < PromptRmax){
		if (nextVecR.Mag() < LateRmax){
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

			      if (KEflag)
				ntout->Fill(parke);
			      else
				ntout->Fill(Energy);
  
			      if (evindex == 0 && nextevindex == 1)
				EV01 += 1;
			      if (evindex == 0 && nextevindex == 2)
				EV02 += 1;
			      if (evindex == 0 && nextevindex == 3)
				EV03 += 1;
			    
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
      nextVecR = TVector3(posX,posY,posZ);
      
      T->GetEntry(lastentry+1-EVcount+j);
      bool goodpair = false;
      if (Fitted > 0 && nextFitted > 0){
	VecR = TVector3(posX,posY,posZ);
	time = nsecs + (secs * pow(10, 9)) + (days * 24 * 3600 * pow(10, 9));
	timeD = std::fabs(nexttime - time);
	DeltaR = (VecR-nextVecR).Mag();
	if(timeD > 400 && timeD < deltaT){
	  if (VecR.Mag() < PromptRmax){
	    if (nextVecR.Mag() < LateRmax){
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

			  if (KEflag)
			    ntout->Fill(parke);
			  else
			    ntout->Fill(Energy);
  
			  if (evindex == 0 && nextevindex == 1)
			    EV01 += 1;
			  if (evindex == 0 && nextevindex == 2)
			    EV02 += 1;
			  if (evindex == 0 && nextevindex == 3)
			    EV03 += 1;
			
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
  
  ntout->Write();
  fin->Close();
  fout->Close();
}



int main(int argc, char *argv[])
{
  if (argc < 13){
    std::cout<<"14 (12) arguments expected: \n 1: input pruned ntuple \n 2: \
output ntuple with single entry EPrompt (UP TO .root!!)\n  \n 3: nhit1Min \n 4: nhit1Max\n 5: nhit2Min \n 6: nhit2Max\n 7: E1Min \n 8: E1Max\n 9: E2Min \n 10: E2Max \n 11: deltaT (in ns) \n (12: KE or E1) \n  \n 12: delmsqr21 \n 13: sinsqrtheta12 \n 14:sinsqrtheta13 \n \n 15: KE or E1"<<std::endl;
  }else{
    
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
    std::string KEorData;

    if (argc == 13){
      delmsqr21 = 0;
      sinsqrtheta12 = 0;
      sinsqrtheta13 = 0;
      
      KEorData = argv[12];
    }
    else{
      delmsqr21 = atof(argv[12]);
      sinsqrtheta12 = atof(argv[13]);
      sinsqrtheta13 = atof(argv[14]);
      
      KEorData = argv[15];
    }
    double PromptRmax = 5700;
    double LateRmax = 5700;
    if (KEorData == "KE")
      OscPromptE_Evindex(infile,outfile,nhit1Min,nhit1Max,nhit2Min,nhit2Max,E1Min,E1Max,E2Min,E2Max,deltaT,PromptRmax,LateRmax,delmsqr21,sinsqrtheta12,sinsqrtheta13,true);
    else if (KEorData == "E1")
      OscPromptE_Evindex(infile,outfile,nhit1Min,nhit1Max,nhit2Min,nhit2Max,E1Min,E1Max,E2Min,E2Max,deltaT,PromptRmax,LateRmax,delmsqr21,sinsqrtheta12,sinsqrtheta13,false);
    else 
      std::cout<<"specificy if 'KE' or 'E1' (EPrompt) fake data is desired"<<std::endl;

  }
}
