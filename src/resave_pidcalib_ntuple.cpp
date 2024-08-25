// Author: Manuel Franco Sevilla
// Description: Resave PIDCalib ntuple in 6.24 so that it can be opened in ROOT 5 

// Standard headers
#include <iostream>
#include <string>
#include <vector>
#include <chrono>

// ROOT headers
#include <TFile.h>
#include <TString.h>
#include <TChain.h>
#include <TKey.h>
#include <TList.h>

// Third-party libraries
#include <cxxopts.hpp>

using namespace std;
using Clock = chrono::steady_clock;

//////////
// Main //
//////////

int main(int argc, char **argv) {
  cxxopts::Options argOpts("gen_postfit_weights_ntuple",
                           "Resave PIDCalib ntuple in 6.24 so that it can be opened in ROOT 5.");

  // clang-format off
  argOpts.add_options()
    ("h,help", "print help")
    ("i,inputName", "input ntuple", cxxopts::value<string>())
    ("o,outputFolder", "specify output folder", cxxopts::value<string>()->default_value("gen/"))
    ;
  // clang-format on

  auto parsedArgs = argOpts.parse(argc, argv);
  if (parsedArgs.count("help")) {
    cout << argOpts.help() << endl;
    return 0;
  }

  TString inputName   = parsedArgs["inputName"].as<string>();
  TString outputFolder = parsedArgs["outputFolder"].as<string>();
  TString treeName = "DecayTree";
  
  // Branches needed to run uBDT
  TString prefix = "probe_Brunel_ANNTraining_";
  TString prefix2 = "probe_Brunel_";
  auto varsD = vector<TString>{
    prefix+"TrackChi2PerDof", prefix+"TrackNumDof", prefix+"TrackGhostProb",
    prefix+"TrackFitMatchChi2", prefix+"TrackFitVeloChi2", prefix+"TrackFitVeloNDoF",
    prefix+"TrackFitTChi2", prefix+"TrackFitTNDoF",
    prefix+"RichDLLe", prefix+"RichDLLmu", prefix+"RichDLLk", prefix+"RichDLLp", prefix+"RichDLLbt",
    prefix+"MuonLLBkg", prefix+"MuonLLMu", prefix+"MuonNShared",
    prefix+"InAccEcal", prefix+"EcalPIDe", prefix+"EcalPIDmu",
    prefix+"InAccHcal", prefix+"HcalPIDe", prefix+"HcalPIDmu",
    prefix+"InAccPrs", prefix+"PrsPIDe",
    prefix+"InAccBrem", prefix+"BremPIDe",
    prefix+"VeloCharge",
    prefix+"TrackP", prefix+"TrackPt",
    "runNumber", "eventNumber"
  };
  auto varsI = vector<TString>{
    prefix2+"RICH1GasUsed", prefix2+"RICH2GasUsed",
    prefix2+"RICHThresholdMu", prefix2+"RICHThresholdKa"
  };
  auto varsB = vector<TString>{
    prefix2 + "isMuonTight"
  };

  // Opening file with input trees
  TFile vFile(inputName);
  TList* list = vFile.GetListOfKeys() ;
  if (!list) {
    cout<< "No keys found in "<< inputName <<"  --- EXITING"<<endl;
    return 1;
  }
  TIter next(list);
  TKey* key;
  TObject* obj;

  // Opening output file
  TString outName = inputName;
  outName.Remove(0, outName.Last('/'));
  outName = outputFolder + outName;
  outName.ReplaceAll(".root", "_resaved.root");
  TFile outFile(outName, "recreate");
  
  // Loop over each folder and its DecayTree inside
  auto start_time = Clock::now();
  while ( (key = static_cast<TKey*>(next())) ) {
    obj = key->ReadObj();
    TString objType = obj->IsA()->GetName();
    if(!objType.Contains("TDirectoryFile")) continue; // Check that it is a folder
    TString folderName = obj->GetName();
    if(folderName.Contains("Luminosity")) continue; // Check that it is not the lumi tree

    // Output file
    TChain ntpin(folderName + "/" + treeName);
    ntpin.Add(inputName);
    outFile.mkdir(folderName);
    outFile.cd(folderName);
    TTree ntpout(treeName, treeName);
  
    // Branches for input and output ntuples
    vector<double> brsD(varsD.size(), -999.);
    for(unsigned ind=0; ind<brsD.size(); ind++){
      ntpin.SetBranchAddress(varsD[ind], &brsD[ind]);
      ntpout.Branch(varsD[ind], &brsD[ind], varsD[ind]+"/D");
    }
    vector<int> brsI(varsI.size(), -999.);
    for(unsigned ind=0; ind<brsI.size(); ind++){
      ntpin.SetBranchAddress(varsI[ind], &brsI[ind]);
      ntpout.Branch(varsI[ind], &brsI[ind], varsI[ind]+"/I");
    }
    bool isMuonTight(false);
    ntpin.SetBranchAddress(varsB[0], &isMuonTight);
    ntpout.Branch(varsB[0], &isMuonTight, varsB[0]+"/O");
    

    // Looping over entries
    long nEntries = ntpin.GetEntries();
    cout<<"Running over "<<inputName<<"/"<<folderName<<" with "<<nEntries<<" entries"<<endl;

    bool debug = false;
    if(debug) nEntries = 10;
    for(long entry=0; entry < nEntries; entry++){
      if(debug) cout<<entry<<"/"<<nEntries<<endl;
      else if(entry%1000000==0) cout<<entry<<"/"<<nEntries<<endl;
      ntpin.GetEntry(entry);

      ntpout.Fill();
    } // for over tree events
    
    ntpout.Write();
  } // while over trees

  outFile.Close();
  double seconds = chrono::duration<double>(Clock::now()-start_time).count();
  cout<<"Took "<<seconds<<" seconds"<<endl;
  cout<<endl<<" rut "<<outName<<endl<<endl;
  
  return 0;
}
