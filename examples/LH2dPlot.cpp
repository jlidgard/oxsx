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

double LHFit(const std::string UnOscfile, const std::string dataFile, int numPdfs, std::vector<double> reactorDistances,  double flux, int numbins, double Emin, double Emax, double Normmin, double Normmax, double d21fix, double s12fix, double s13fix){
  char name[100];

  TRandom3 *r1 = new TRandom3();
  r1->SetSeed(0);

  ////////////////////
  // 1. Set Up PDFs //
  ////////////////////

  // Only interested in first bit of data ntuple
  ObsSet dataRep(0);

  // Set up binning
  AxisCollection axes;
  axes.AddAxis(BinAxis("ParKE", Emin, Emax, numbins));

  BinnedED dataSetPdf("dataSetPdf",axes);
  dataSetPdf.SetObservables(dataRep);
  ROOTNtuple dataNtp(dataFile, "nt");

  for(size_t i = 0; i < dataNtp.GetNEntries(); i++)
    dataSetPdf.Fill(dataNtp.GetEntry(i));

  // flux defined as Desired_Flux/Produced_Data_Flux
  int dataint = (int)(flux*dataSetPdf.Integral());
  dataSetPdf.Normalise();
  dataSetPdf.Scale(dataint);

  ROOTNtuple reactorNtp(UnOscfile, "nt");
  //NuOsc *reactorSystematic;
  ParameterDict minima;
  ParameterDict maxima;
  ParameterDict initialval;
  ParameterDict initialerr;
 
  BinnedNLLH lhFunction;
  lhFunction.SetBufferAsOverflow(true);
  int Buff = 5;
  lhFunction.SetBuffer(0,Buff,Buff);
  lhFunction.SetDataDist(dataSetPdf); // initialise withe the data set
  
  BinnedED * reactorPdf0 = new BinnedED(name,axes);
  reactorPdf0->SetObservables(0);
  for(size_t j = 0; j < reactorNtp.GetNEntries(); j++)
    reactorPdf0->Fill(reactorNtp.GetEntry(j));
  reactorPdf0->Normalise();

  // loop over all reactor pdfs
  for (int i = 0; i< numPdfs; i++){
    sprintf(name,"ReactorPdf%d",i);
    BinnedED * reactorPdf = new BinnedED(name,axes);
    reactorPdf->SetObservables(0);
    
    NuOsc reactorSystematic("reactorSystematic");
    reactorSystematic.SetFunction(new SurvProb(d21fix,s12fix,s13fix,reactorDistances[i]));
    reactorSystematic.SetAxes(axes);
    reactorSystematic.SetTransformationObs(dataRep);
    reactorSystematic.SetDistributionObs(dataRep);
    reactorSystematic.Construct();
    
    reactorPdf->Add(reactorSystematic(*reactorPdf0),1);
    reactorPdf->Normalise();
    
    // Setting optimisation limits
    sprintf(name,"ReactorPdf%d_norm",i);
    minima[name] = 0;
    maxima[name] = 100000;
    initialval[name] = 50000;
    initialerr[name] = 0.1*initialval[name];

    lhFunction.AddDist(*reactorPdf);
  }
  
  ////////////
  // 4. Fit //
  ////////////
  Minuit min;
  min.SetMethod("Migrad");
  min.SetMaxCalls(1000000);
  min.SetMinima(minima);
  min.SetMaxima(maxima);
  min.SetInitialValues(initialval);
  min.SetInitialErrors(initialerr);

  /////////////////////////////////////////////
  ////////        Fit Result        ///////////
  /////////////////////////////////////////////

  FitResult fitResult = min.Optimise(&lhFunction);
  ParameterDict bestFit = fitResult.GetBestFit();
  fitResult.Print();

  lhFunction.SetParameters(bestFit);
  double lhval =(-1)*lhFunction.Evaluate();
  if (lhval > lhmax)
    lhmax = lhval;
  if (lhval < lhmin)
    lhmin = lhval;
  std::cout<<"LH val at best fit: "<<lhval<<std::endl;
  return lhval;
}

int main(int argc, char *argv[]) {

  if (argc != 5) {
    std::cout<<"4 arguments expected: \n 1: location of UnOscfile \n 2: location/filename\
 for data spectrum to fit \n 3: location/filename for reactor info file \n 4: filename of output\
 plot"<<std::endl;
  }
  else{
    int flux = 1;
    
    double s12min = 0.025;
    double s12max = 1.;
    double d21min = 6.0e-5;
    double d21max = 8.8e-5;
    
    double s12int = 0.025;
    double d21int = 0.05e-5;
    
    double Emin = 2;
    double Emax = 8;
    int numbins = 30;
    double Normmin = 0;
    double Normmax = 100000;

    double Num_s12 = ((s12max-s12min)/s12int)+1;
    double Num_d21 = ((d21max-d21min)/d21int)+1;
    int num_s12 = round(Num_s12);
    int num_d21 = round(Num_d21);

    std::cout<<"num_s12: "<<num_s12<<std::endl;
    std::cout<<"num_d21: "<<num_d21<<std::endl;
    
    TH2D *h2 = new TH2D("h2", "h2", num_s12,s12min,s12max,num_d21,d21min,d21max);
    int k = 0;
    
    std::stringstream argParser;
    const std::string &UnOscfile = argv[1];
    const std::string &dataFile = argv[2];
    const std::string &infoFile = argv[3];
    const std::string &outFile = argv[4];
    argParser << argv[5];
    int numexps;
    argParser >> numexps;

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
    for (int i=1; i < num_s12+1; i++){
      for (int j=1; j < num_d21+1; j++){
	double s12 = s12min+((i-1)*s12int);
	double d21 = d21min+((j-1)*d21int); 
	printf("\n-------------------------------------------------------------------------------------\n");
	k+=1;
	std::cout<<"point: "<<k<<" of "<<num_s12*num_d21<<std::endl;
	std::cout<<"s12: "<<s12<<"    (from "<<s12min<<" to "<<s12max<<")"<<std::endl;
	std::cout<<"d21: "<<d21<<"    (from "<<d21min<<" to "<<d21max<<")"<<std::endl;
	std::cout<<""<<std::endl;
	std::cout<<"bin x: "<<i<<"    (from "<<1<<" to "<<num_s12<<")"<<std::endl;
	std::cout<<"bin y: "<<j<<"    (from "<<1<<" to "<<num_d21<<")"<<std::endl;
	h2->SetBinContent(i,j,LHFit(UnOscfile,dataFile,numPdfs,distances,flux,numbins,Emin,Emax,Normmin,Normmax,d21,s12,0.0215));
	printf("-------------------------------------------------------------------------------------\n");
      }
    }
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
    TFile *fileOut = new TFile(outFile.c_str(),"RECREATE");
    h2->Write();
    fileOut->Close();
  }
}
