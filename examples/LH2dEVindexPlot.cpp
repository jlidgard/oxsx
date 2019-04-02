// A simple fit in energy for signal and a background
#include <string>
#include <vector>
#include <math.h>
#include <Rand.h>
#include <fstream>
#include <iostream>

#include <TH1D.h>
#include <TH2D.h>
#include <TPaveText.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <ROOTNtuple.h>
#include <TFile.h>
#include <TMath.h>
#include <string>
#include <TNtuple.h>
#include <iostream>
#include <TObject.h>
#include <TRandom3.h>
#include <TVector3.h>

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
#include <TH1D.h>
#include <ROOTMultiPlot.h>
#include <TRandom3.h>

// flag which is used for repeating fits that weren't good
bool badfit = true;

int unoscPdfint;
int oscPdfint;

//Surv Probablity function for antinu of KE nuE and distance traveled of baseline.
double probability(double nuE, double baseline, double delmsqr21, double sinsqrtheta12, double sinsqrtheta13) {
  double fSSqr2Theta12 = pow(sin(2.0 * TMath::ASin(sqrt(sinsqrtheta12))), 2.0);
  double fS4 = pow(sinsqrtheta13, 2.0);
  double fC4 = pow(1.0-sinsqrtheta13, 2.0);
  double scale = 1.267e3; // for nuE in [MeV] and baseline in [km]
  double sSqrDmBE = pow(sin(scale * delmsqr21 * baseline / nuE), 2.0);
  double fOscProb = (fC4 * (1.0 - fSSqr2Theta12 * sSqrDmBE) + fS4);
  return fOscProb;
}

//Coincidence tagging - function to apply cuts to isolate positron and neutron events
// At present, looking at the MC event number and EV index number to ensure no pairs are missed
void OscPromptE_EVindex (const std::string infile, const std::string outfile, int nhit1Min, int nhit1Max, int nhit2Min, int nhit2Max, double E1Min, double E1Max, double E2Min, double E2Max, double deltaT,double PromptRmax,double LateRmax,  double delmsqr21, double sinsqrtheta12, double sinsqrtheta13, double Distance) {
  unoscPdfint = 0;
  oscPdfint = 0;
  TFile *fin = TFile::Open(infile.c_str());
  TNtuple *T = (TNtuple *) fin->Get("nt");

  TFile *fout = new TFile(outfile.c_str(), "RECREATE");
  TNtuple* ntout = new TNtuple("nt","nt","E1");

  //entries used for coincidence tagging
  float MCentryb4, MCentry,nextMCentry;
  float days;
  float secs;
  float nsecs;
  float nhit,nextnhit;
  float evindex, nextevindex;
  double timeD, DeltaR;
  float Fitted, nextFitted;
  float Energy,nextEnergy,parke,part1ke,part2ke;
  float posX,posY,posZ,nextposX,nextposY,nextposZ;
  float EVcount = 0;
  //TVector3 VecR, nextVecR;  //TVector3 doesn't seem to work in oxsx, maybe need to add something to the compiler?
  double VecRMag,nextVecRMag;
  double time, nexttime;
  
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
  //T->SetBranchAddress("ReactorDistance", &ReactorDistance); // Average distance taken from info file, not event by event (though it can be)

  double survprob;
  TRandom3 *r1 = new TRandom3();
  r1->SetSeed(0);    

  //get first event MC entry, go through each MC event (MC index), then loop through triggered events from that MC event (EV index) looking for the proton and neutron
  // loop through triggered events until get to the last EV index for a given MC index, or until deltaT between chosen events > deltaT Cut
  T->GetEntry(0);
  MCentryb4 = MCentry;  
  //go through entries and collect triggered evs which correspond to a MC event.
  for(int i = 1; i < T->GetEntries(); i++){
    T->GetEntry(i);
    //mark when passed group of evs with same mcindex:
    if(abs(MCentry - MCentryb4)>= 1 && evindex == 0){ 
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
	  nextMCentry = MCentry;
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
	    //VecR = TVector3(posX,posY,posZ);
	    if(timeD > 400 && timeD < deltaT){
	      if (VecRMag < PromptRmax){
		if (nextVecRMag < LateRmax){
		  if (DeltaR < 1500){
		    if (E1Min <= Energy && Energy <= E1Max){
		      if (E2Min <= nextEnergy &&nextEnergy <= E2Max){
			if (nhit1Min <= nhit  && nhit <= nhit1Max){
			  if (nhit2Min <= nextnhit && nextnhit <= nhit2Max){
			    survprob = probability(parke, Distance, delmsqr21, sinsqrtheta12, sinsqrtheta13);
			    unoscPdfint +=1;
			    if (survprob > r1->Rndm()){
			      oscPdfint +=1;
			      ntout->Fill(Energy);
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
	if (pairfound){ // move onto next mc group once pair found
	  j += 100;
	}else{
	  j += 1;
	}
      } 
    }
    MCentryb4 = MCentry;
  }
  // need to repeat process for the very last MC event, can neglet for high stats
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
      if (Fitted > 0 && nextFitted > 0){  //valid fit
	//VecR = TVector3(posX,posY,posZ);
	VecRMag = sqrt((pow(posX,2))+(pow(posY,2))+(pow(posZ,2)));
	time = nsecs + (secs * pow(10, 9)) + (days * 24 * 3600 * pow(10, 9));
	timeD = std::fabs(nexttime - time);
	//DeltaR = (VecR-nextVecR).Mag();
	DeltaR = sqrt((pow(posX-nextposX,2))+(pow(posY-nextposY,2))+(pow(posZ-nextposZ,2)));
	if(timeD > 400 && timeD < deltaT){
	  if (VecRMag < PromptRmax){
	    if (nextVecRMag < LateRmax){
	      if (DeltaR < 1500){
		if (E1Min <= Energy && Energy <= E1Max){
		  if (E2Min <= nextEnergy &&nextEnergy <= E2Max){
		    if (nhit1Min <= nhit  && nhit <= nhit1Max){
		      if (nhit2Min <= nextnhit && nextnhit <= nhit2Max){
			survprob = probability(parke, Distance, delmsqr21, sinsqrtheta12, sinsqrtheta13);
			unoscPdfint +=1;
			if (survprob > r1->Rndm()){
			  oscPdfint += 1;
			  ntout->Fill(Energy);
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
  
  ntout->Write();
  delete fin;
  delete fout;
}
  

void readInfoFile(const std::string &runInfoFileName, std::vector<std::string> &reactorNames, std::vector<double> &distances, std::vector<std::string> &reactorTypes, std::vector<int> &nCores, std::vector<double> &powers ) {
  // Read couchDB run-info text file
  std::ifstream in;
  in.open(runInfoFileName.c_str());
  std::cout << "opening file: " << runInfoFileName.c_str() << std::endl;

  std::fill(reactorNames.begin(), reactorNames.end(), "");
  std::fill(distances.begin(), distances.end(), 0.);
  std::fill(reactorTypes.begin(), reactorTypes.end(), "");
  std::fill(nCores.begin(), nCores.end(), 0);
  std::fill(powers.begin(), powers.end(), 0.);

  std::string reactorName,distance,reactorType,nCore,power;
  int lineNo = 0;

  // read until end of file.
  while(in.good()){
    std::getline(in,reactorName,',');
    std::getline(in,distance,',');
    std::getline(in,reactorType,',');
    std::getline(in,nCore,',');
    std::getline(in,power,'\n');

    if (lineNo>0){ //skip csv header
      if (strcmp(reactorName.c_str(),"")!=0) {
	reactorNames.push_back(reactorName);
	distances.push_back(atof(distance.c_str()));
	reactorTypes.push_back(reactorType.c_str());
	nCores.push_back(atoi(nCore.c_str()));
	powers.push_back(atof(power.c_str()));

	std::cout << "v: reactorName: " << reactorNames[lineNo-1] << ", distance: " << distances[lineNo-1] << ", reactorType: " << reactorTypes[lineNo-1] << ", nCore: " << nCores[lineNo-1] << ", power: " << powers[lineNo-1] << std::endl; //debug check ('-1' for header)
      }
    }
    lineNo++;
  }
  in.close();
}

double lhmax = -1000000000;
double lhmin = 1000000000;

//double LHFit(const std::string UnOscfile, const std::string dataFile, AxisCollection axes, BinnedED dataSetPdf, const std::string tempFile, int numPdfs, std::vector<std::string> reactorNames, std::vector<double> reactorDistances,  double fluxfrac, int numbins, double Emin, double Emax, double Normmin, double Normmax, double d21fix, double s12fix,  int nhit1Min, int nhit1Max, int nhit2Min, int nhit2Max, double E1Min, double E1Max, double E2Min, double E2Max, double deltaT,double PromptRmax,double LateRmax, int Pdfint){
double LHFit(const std::string UnOscfile, const std::string dataFile, AxisCollection axes, BinnedED dataSetPdf, const std::string tempFile, int numPdfs, std::vector<std::string> reactorNames, std::vector<double> reactorDistances, int numbins, double Emin, double Emax, double Normmin, double Normmax, double d21fix, double s12fix,  int nhit1Min, int nhit1Max, int nhit2Min, int nhit2Max, double E1Min, double E1Max, double E2Min, double E2Max, double deltaT,double PromptRmax,double LateRmax, std::vector<double> datanorms){
  char name[100];
  TRandom3 *r1 = new TRandom3();
  r1->SetSeed(0);

  ////////////////////
  // 1. Set Up PDFs //
  ////////////////////
  
  
  ParameterDict minima;
  ParameterDict maxima;
  ParameterDict initialval;
  ParameterDict initialerr;
 
  BinnedNLLH lhFunction;
  lhFunction.SetBufferAsOverflow(true);
  int Buff = 1;
  lhFunction.SetBuffer(0,Buff,Buff);  
  lhFunction.SetDataDist(dataSetPdf); // initialise withe the data set
  double rand = r1->Rndm();

  for (int i = 0; i< numPdfs; i++){
    OscPromptE_EVindex(UnOscfile, tempFile,nhit1Min,nhit1Max,nhit2Min,nhit2Max,E1Min,E1Max,E2Min,E2Max,deltaT,PromptRmax,LateRmax,d21fix,s12fix,0.0215,reactorDistances[i]);   //oscillate each input reactor KE pdfs -> outputs a pruned ntuple with an oscillated EPrompt positron spectrum called tempFile, to be used below to fill a reactorPdf
    
    BinnedED *reactorPdf = new BinnedED(reactorNames[i],axes);
    reactorPdf->SetObservables(0);
    
    ROOTNtuple reactorNtp(tempFile.c_str(), "nt");
    for(size_t j = 0; j < reactorNtp.GetNEntries(); j++)
      reactorPdf->Fill(reactorNtp.GetEntry(j));
    std::cout<<"pre osc Int: "<<unoscPdfint<<std::endl;   //integral of input unoscillated pdfs
    //int postoscint = (int)(reactorPdf->Integral());
    std::cout<<"post osc Int: "<<oscPdfint<<std::endl;  //integral of pdf after KE->EPrompt + oscillation
    double osc_loss = oscPdfint/(double)unoscPdfint;      // calculate loss in pdf integral
    std::cout<<" osc loss factor: "<<osc_loss<<std::endl;
    reactorPdf->Normalise();

    // Setting optimisation limits
    sprintf(name,"%s_norm",reactorNames[i].c_str());
    minima[name] = 0;//Normmin;
    maxima[name] = 500;//fluxfrac*20000;   // CHANGE THIS DEPENDING ON SPECIFIC EXAMPLE
    initialval[name] = (rand*(maxima[name]-minima[name]))+minima[name];
    initialerr[name] = 0.5*initialval[name];
    
    lhFunction.AddDist(*reactorPdf);
    double constraint = datanorms[i] * osc_loss; //vector of unoscillated norm contraints made in main() Jeff you can probs improve this
    lhFunction.SetConstraint(name, constraint, 0.1*constraint);  // I just randomly picked 10% for contraint uncertainty 
  }
  
  ////////////
  // 4. Fit //
  ////////////
  Minuit min;
  min.SetMethod("Migrad");
  min.SetMaxCalls(10000000);
  min.SetTolerance(0.001);
  min.SetMinima(minima);
  min.SetMaxima(maxima);
  min.SetInitialValues(initialval);
  min.SetInitialErrors(initialerr);
  
  /////////////////////////////////////////////
  ////////        Fit Result        ///////////
  /////////////////////////////////////////////
  FitResult fitResult = min.Optimise(&lhFunction);
  ParameterDict bestFit = fitResult.GetBestFit();

  bool fitValid = fitResult.GetValid();
  if (fitValid){
    fitResult.Print();
    badfit = false;
    std::cout<<"goodfit!"<<std::endl;
  }else{
    std::cout<<"INVALID FIT!! \n Trying Again:"<<std::endl;
    badfit = true;
  }
  
  lhFunction.SetParameters(bestFit);
  double lhval =(-1)*lhFunction.Evaluate();
  if (lhval > lhmax)
    lhmax = lhval;
  if (lhval < lhmin)
    lhmin = lhval;
  printf("\n-------------------------------------------------------------------------------------\n");
  std::cout<<"LH val at best fit: "<<lhval<<std::endl;
  printf("-------------------------------------------------------------------------------------\n");
	
  return lhval;
}

int main(int argc, char *argv[]) {

  if (argc != 17) {
    std::cout<<"16 arguments expected: \n 1: location of UnOsc pruned Ntuple file (needed for tagging/extraction and oscillating of prompt energy event spectrum)  \n 2: location/filename for data spectrum to be fit (pruned ntuple with just one entry, EPrompt) \n 3: location/filename for reactor info file \n 4: filename of temporary TFile needed to oscillate Ntuple (place anywhere, can delete afterwards) \n 5: filename of output plot \n  \n 6: nhit1Min \n 7: nhit1Max\n 8: nhit2Min \n 9: nhit2Max\n 10: E1Min \n 11: E1Max\n 12: E2Min \n 13: E2Max \n 14: deltaT (in ns) \n \n 15: Fake data at d21 \n 16:Fake data at s12"<<std::endl;
  }
  else{
    //desired (flux)/(MC produced flux)
    double fluxfrac = 1;
    
    double s12min = 0.025;
    double s12max = 1.;
    double d21min = 3.4e-5;//6.0e-5;
    double d21max = 6.2e-5;//8.8e-5;
    
    double s12int = 0.025;
    double d21int = 0.05e-5;
    
    /*
    //for quick testing vv
    double s12min = 0.2;
    double s12max = 0.4;
    double d21min = 3.5e-5;//6e-5;
    double d21max = 6.4e-5;//9.0e-5;
    
    double s12int = 0.1;
    double d21int = 0.5e-5;
    */
    //axes range and number of bins, must not include any bins with zero entries!!!
    double Emin = 1.;
    double Emax = 8.;
    int numbins = 21;
    double Normmin = 0;
    double Normmax = 50000; //not used presently

    double Num_s12 = ((s12max-s12min)/s12int)+1;
    double Num_d21 = ((d21max-d21min)/d21int)+1;
    int num_s12 = round(Num_s12);
    int num_d21 = round(Num_d21);

    std::cout<<"num_s12: "<<num_s12<<std::endl;
    std::cout<<"num_d21: "<<num_d21<<std::endl;

    double PromptRmax = 5700;
    double LateRmax = 5700;
    
    TH2D *h2 = new TH2D("h2", "h2", num_s12,s12min,s12max,num_d21,d21min,d21max);
    int k = 0;
    
    std::stringstream argParser;
    const std::string &UnOscfile = argv[1];
    const std::string &dataFile = argv[2];
    const std::string &infoFile = argv[3];
    const std::string &tempFile = argv[4];
    const std::string &outFile = argv[5];
    argParser << argv[5];
    int numexps;
    argParser >> numexps;
    
    int nhit1Min = atoi(argv[6]);
    int nhit1Max = atoi(argv[7]);
    int nhit2Min = atoi(argv[8]);
    int nhit2Max = atoi(argv[9]);
    double E1Min = atof(argv[10]);
    double E1Max = atof(argv[11]);
    double E2Min = atof(argv[12]);
    double E2Max = atof(argv[13]);
    double deltaT = atof(argv[14]);
    
    double FakeDatad21 = atof(argv[15]);
    double FakeDatas12 = atof(argv[16]);

    std::vector<std::string> reactorNames;
    std::vector<double> distances;
    std::vector<std::string> reactorTypes;
    std::vector<int> nCores;
    std::vector<double> powers;
    std::cout<<"starting filling info....... \n"<<std::endl;
    readInfoFile(infoFile, reactorNames, distances, reactorTypes, nCores, powers); // get reactor information
    std::cout<<"\n finished filling info"<<std::endl;
    int numPdfs = reactorNames.size();

    std::cout<<"Number of reactors: "<<numPdfs<<std::endl;
    
    for (int i = 0; i< numPdfs; i++){
      std::cout<<"Reactor name: "<<reactorNames[i]<<" Distance:"<<distances[i]<<std::endl;
    }

    //SET UP DATA PDFS
    // Only interested in first bit of data ntuple
    ObsSet dataRep(0);
    // Set up binning
    AxisCollection axes;
    axes.AddAxis(BinAxis("EPrompt (MeV)", Emin, Emax, numbins));

    //Not using input data at the moment, generating fake data for the 3CAD reactors
    /*
      BinnedED dataSetPdf("dataSetPdf",axes);
      dataSetPdf.SetObservables(dataRep);
      ROOTNtuple dataNtp(dataFile, "nt");
 
      for(size_t i = 0; i < dataNtp.GetNEntries(); i++)
      dataSetPdf.Fill(dataNtp.GetEntry(i));
    */
 
    // Can scale data spectrum to desired length of time running (quite approximately)
    // fluxfrac defined as Desired_Flux/Produced_Data_Flux
    //int dataint = (int)(fluxfrac*dataSetPdf.Integral());
    //dataSetPdf.Normalise();
    //dataSetPdf.Scale(dataint);
  
    // Generating fake data for 3CAD, norms put in by hand below, taken from number of event expected from the 3CAD reactors using flux=1, 1yr data
    BinnedED dataSetPdf("dataSetPdf",axes);
    dataSetPdf.SetObservables(dataRep);
  
    std::vector<double> datanorms;
    datanorms.push_back(93.);
    datanorms.push_back(21.);
    datanorms.push_back(20.);

    for (int i = 0; i< numPdfs; i++){
      OscPromptE_EVindex(UnOscfile, tempFile,nhit1Min,nhit1Max,nhit2Min,nhit2Max,E1Min,E1Max,E2Min,E2Max,deltaT,PromptRmax,LateRmax,FakeDatad21,FakeDatas12,0.0215,distances[i]);   //oscillate each pdf -> outputs a pruned ntuple with an oscillated EPrompt positron spectrum called tempFile, to be used below to fill a reactorPdf
    
      BinnedED *dataPdf = new BinnedED(reactorNames[i],axes);
      dataPdf->SetObservables(0);
    
      ROOTNtuple dataNtp(tempFile.c_str(), "nt");
      for(size_t j = 0; j < dataNtp.GetNEntries(); j++)
	dataPdf->Fill(dataNtp.GetEntry(j));

      std::cout<<"pre osc Int: "<<unoscPdfint<<std::endl;   //integral of input unoscillated pdfs
      //int postoscint = (int)(reactorPdf->Integral());
      std::cout<<"post osc Int: "<<oscPdfint<<std::endl;  //integral of pdf after KE->EPrompt + oscillation
      double osc_loss = oscPdfint/(double)unoscPdfint;      // calculate loss in pdf integral
      std::cout<<" osc loss factor: "<<osc_loss<<std::endl;
      
      dataPdf->Normalise();
      dataPdf->Scale(datanorms[i]*osc_loss);
  
      dataSetPdf.Add(*dataPdf,1);
    } 
    double LHmean = 0;
    for (int i=1; i < num_s12+1; i++){
      for (int j=1; j < num_d21+1; j++){
	double s12 = s12min+((i-1)*s12int);
	double d21 = d21min+((j-1)*d21int); 
	//printf("\n-x------------------------------------------------------------------------------------\n");
	k+=1;
	std::cout<<"point: "<<k<<" of "<<num_s12*num_d21<<std::endl;
	std::cout<<"s12: "<<s12<<"    (from "<<s12min<<" to "<<s12max<<")"<<std::endl;
	std::cout<<"d21: "<<d21<<"    (from "<<d21min<<" to "<<d21max<<")"<<std::endl;
	std::cout<<""<<std::endl;
	std::cout<<"bin x: "<<i<<"    (from "<<1<<" to "<<num_s12<<")"<<std::endl;
	std::cout<<"bin y: "<<j<<"    (from "<<1<<" to "<<num_d21<<")"<<std::endl;
	double LHval;
	int fitattempts = 0;
	badfit = true;
	while (badfit){
	  fitattempts += 1;
	  std::cout<<"fit attempt: "<<fitattempts<<std::endl;
	  //LHval = LHFit(UnOscfile,dataFile,axes,dataSetPdf,tempFile,numPdfs,reactorNames,distances,fluxfrac,numbins,Emin,Emax,Normmin,Normmax,d21,s12,nhit1Min,nhit1Max,nhit2Min,nhit2Max,E1Min,E1Max,E2Min,E2Max,deltaT,PromptRmax,LateRmax,Pdfint);
	  LHval = LHFit(UnOscfile,dataFile,axes,dataSetPdf,tempFile,numPdfs,reactorNames,distances,numbins,Emin,Emax,Normmin,Normmax,d21,s12,nhit1Min,nhit1Max,nhit2Min,nhit2Max,E1Min,E1Max,E2Min,E2Max,deltaT,PromptRmax,LateRmax,datanorms);
	}
	LHmean += LHval;	
	h2->SetBinContent(i,j,LHval);
	//printf("-------------------------------------------------------------------------------------\n");
      }
    }
    LHmean = LHmean/(num_s12*num_d21);
    
    TH1D DataHist;
    DataHist = DistTools::ToTH1D(dataSetPdf);
    DataHist.SetName("DataHist");

    h2->GetXaxis()->SetTitle("(sintheta12)^2");
    h2->GetXaxis()->SetTitleSize(0.05);
    h2->GetXaxis()->SetTitleOffset(0.9);
    h2->GetXaxis()->CenterTitle();
    h2->GetYaxis()->SetTitle("(delm21)^2 (eV^2)");
    h2->GetYaxis()->SetTitleSize(0.05);
    h2->GetYaxis()->SetTitleOffset(1);
    h2->GetYaxis()->CenterTitle();
    h2->SetLabelSize(0.05,"xyz");
    h2->SetTitle("Best Fit LH value vs Osc Parameters");

    std::cout<<"LH MAX VAL:   "<<lhmax<<std::endl;
    std::cout<<"LH MIN VAL:   "<<lhmin<<std::endl;
    std::stringstream ssout;
    ssout<<outFile<<"ds21_"<<FakeDatad21<<"s12_"<<FakeDatas12<<".root";
    TFile *fileOut = new TFile(ssout.str().c_str(),"RECREATE");
    h2->Write();
    DataHist.Write();
    fileOut->Close();
  }
}
