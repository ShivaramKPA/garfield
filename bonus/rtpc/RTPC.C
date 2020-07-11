//////////////////////////////////////////
//
// Simulation of the BONuS RTPC
//
// Nate Dzbenski
// ndzbe001@odu.edu
//
// Gabriel Charles
// gcharles@odu.edu
//
// 22 June 2020
//
///////////////////////////////////////////

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <map>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoPcon.h>
#include <TGeoHalfSpace.h>
#include <TGeoMatrix.h>
#include <TGeoCompositeShape.h>

#include "GeometryRoot.hh"
#include "ViewGeometry.hh"
#include "ComponentAnsys123.hh"
#include "ViewField.hh"
#include "MediumMagboltz.hh"
#include "Sensor.hh"
#include "AvalancheMicroscopic.hh"
#include "AvalancheMC.hh"
#include "AvalancheMC.hh"
#include "Random.hh"
#include "Plotting.hh"
#include "MediumMagboltz.hh"
#include "FundamentalConstants.hh"
#include "ComponentAnalyticField.hh"
#include "ViewCell.hh"
#include "SolidBox.hh"
#include "SolidTube.hh"
#include "Sensor.hh"
#include "GeometrySimple.hh"
#include "TrackHeed.hh"
#include "ComponentElmer.hh"
#include "ViewFEMesh.hh"
#include "ComponentVoxel.hh"


using namespace Garfield;
using namespace std;

char tvec_name[15];
char phivec_name[15];


// create the file to store output
TFile *oFile = new TFile("rtpc_out.root","RECREATE","RTPC output");

// create trees to store info
TTree *tTree = new TTree("tTree","Tree of drift times");
TTree *pTree = new TTree("pTree","Tree of drift angles");

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);
    
    if ( !argv[1] ){
        cout << "You need to execute the program with the number of events to simulate, with: " << endl;
        cout << "./garf_rtpc [number of events]" << endl;
        cout << "Where [number of events] is the number of events you want to simulate." << endl;
        exit(EXIT_FAILURE);
    }
    
//______________________________________________________________________________________________
//___________________________________________Variables__________________________________________
//______________________________________________________________________________________________

  //const double Bx = 0.; // Tesla
  //const double By = 0.; // Tesla
  //const double Bz = 4.; // Tesla

  // Set the initial position [cm] and starting time [ns].
  double x0 = 3.1, y0 = 0.0, z0 = -19.0, t0 = 0.;
  // Set the initial energy [eV].
  double e0 = 36.;
  // Set the initial direction (x, y, z).
  // In case of a null vector, the direction is randomized.
  double dx0 = 0., dy0 = 0., dz0 = 0.;
  // Calculate an electron avalanche
  int ne, ni;
  int ne_tot;
  // Electron information after the avalanche
  Double_t x1, y1, z1, t1, e1;
  Int_t status;
    
//______________________________________________________________________________________________
//________________________________________ Canvas and plots ____________________________________
//______________________________________________________________________________________________

    //ViewField *viewfield = new ViewField();

    vector<double> drift_times[9][9];
    vector<double> drift_angles[9][9];
    
    double zeros[9] = {0,0,0,0,0,0,0,0,0};
    double the_rs[9] = {3.1,3.5,4.0,4.5,5.0,5.5,6.0,6.5,6.9};
    double the_zs[9] = {-19.0,-15.0,-10.0,-5.0,0.0,5.0,10.0,15.0,19.0};
    
//_____________________________________________________________________________________________
//___________________________________________ Openings ________________________________________
//_____________________________________________________________________________________________

  // Setup the gas.
  MediumMagboltz* gas = new MediumMagboltz();
//  gas->LoadGasFile("gasFiles/He_100_DME_00.gas");
  gas->SetComposition("He",80.,"CO2",20.);
  gas->SetTemperature(293.);
  gas->SetPressure(760.);
  gas->EnableDrift();                           // Allow for drifting in this medium
  gas->PrintGas();

  // Import an Elmer-created LEM and the weighting field for the readout electrode.
  ComponentElmer * elm = new ComponentElmer("RTPC/mesh.header","RTPC/mesh.elements","RTPC/mesh.nodes",
  "RTPC/dielectrics.dat","RTPC/RTPC.result","cm");
  elm->SetMedium(0,gas);
  // max B-field is 5.22 T
  elm->LoadMagneticField("Fieldmaps/solenoid_map_may2019.dat", 0.778);


//______________________________________________________________________________________________
//_____________________________________________ Code ___________________________________________
//______________________________________________________________________________________________


  // Assemble a Sensor object 
  Sensor* sensor = new Sensor(); 
  sensor->SetArea(-7.0,-7.0,-20.0,7.0,7.0,20.0);
  // Calculate the electric field using the Component object cmp
  sensor->AddComponent(elm);
  //sensor->AddComponent(bfield);

  /*TCanvas *c_e = new TCanvas("Cell","Cell");
  ViewDrift* v_e = new ViewDrift();
  v_e->SetCanvas(c_e);
  v_e->SetArea(-7.0,-7.0,-20.0,7.0,7.0,20.0);*/
    

  // Evaluate the number of electrons in the avalanche
  AvalancheMicroscopic* aval = new AvalancheMicroscopic(); // did not get it to work with AvalancheMC()
  //AvalancheMC* aval = new AvalancheMC();
    aval->SetSensor(sensor);
  //aval->EnablePlotting(v_e);
  aval->EnableMagneticField();
    
    // begin the data creation
    for(int z_i=0; z_i < 9; z_i++){
        z0 = the_zs[z_i];
        
        for(int r_i=0; r_i < 9; r_i++){
            // Set the initial position [cm] and starting time [ns].
            x0 = the_rs[r_i];
            y0 = 0.0;
            t0 = 0.0;
            
            cout << "z = " << z0 << ", r = " << x0 << endl;
            
            sprintf(tvec_name,"bTd_%i_%i",z_i,r_i);
            sprintf(phivec_name,"bPhid_%i_%i",z_i,r_i);
            
            // create a branch for each vector
            tTree->Branch(tvec_name, &drift_times[z_i][r_i]);
            pTree->Branch(phivec_name, &drift_angles[z_i][r_i]);
            
            for(int eve=0;eve<atoi(argv[1]);eve++){
                ne_tot=0;
                cout << "Event number: " << eve << endl;
                aval->AvalancheElectron(x0, y0, z0, t0, e0, dx0, dy0, dz0);
                
                // Get the number of electrons and ions in the avalanche.
                aval->GetAvalancheSize(ne, ni);
                ne_tot+=ne;
                
                if(0<aval->GetNumberOfElectronEndpoints()){
                    for(int nava=0;nava<aval->GetNumberOfElectronEndpoints();nava++){
                        aval->GetElectronEndpoint(nava, x0, y0, z0, t0, e0, x1, y1, z1, t1, e1, status);
                        
                        drift_times[z_i][r_i].push_back(t1);
                        if(x1 != 0)drift_angles[z_i][r_i].push_back(TMath::ATan2(y1,x1));
                    }
                }
            }
        }
    }
    tTree->Fill();
    pTree->Fill();
    
    oFile->Write();
    oFile->Close();
    
    
  /*viewfield->SetComponent(bfield);
  viewfield->SetSensor(sensor);
  viewfield->SetCanvas((TCanvas*)c_field->cd());
    viewfield->PlotProfile(0., 0., 20.0, 0., 0., -20.0);
  //viewfield->SetWeightingFieldRange(0.0, 10000.0);
    
    ViewFEMesh* meshView = new ViewFEMesh();
    meshView->SetCanvas(c_field);
    // Set the component.
    meshView->SetComponent(elm);
    // Set the viewing plane.
    meshView->SetPlane(0, 0, -1, 0, 0, 1);
    meshView->SetFillMesh(true);
    //meshView->SetViewDrift(driftView);
    meshView->SetArea(-8, -8, -20, 8, 8, 20);
    c_field->cd();
    meshView->Plot(); */
    
    
//  Field plot
  //c_field->cd();
  //viewfield->PlotContour();
    
  app.Run(kFALSE);
    return 0;

}
