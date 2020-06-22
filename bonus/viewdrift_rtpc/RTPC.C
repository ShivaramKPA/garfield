//////////////////////////////////////////
//
// Drift View of the BONuS RTPC
//
// Nate Dzbenski
// ndzbe001@odu.edu
//
// Gabriel Charles
// gcharles@odu.edu
//
// 13 Feb 2020
//
///////////////////////////////////////////

#include <iostream>
#include <stdio.h>

#include <map>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TMath.h>
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


int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);


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

    TCanvas* c_dLine = new TCanvas("c_dLine","Drift Lines",800,600);
    
    ViewDrift* driftView = new ViewDrift();
    driftView->SetCanvas(c_dLine);
    
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
  elm->LoadMagneticField("Fieldmaps/solenoid_map_may2019.dat", 1.0);



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
    
    driftView->SetArea(0.0,0.0,-20.0,7.0,7.0,20.0);
  // Evaluate the number of electrons in the avalanche
  AvalancheMicroscopic* aval = new AvalancheMicroscopic(); // did not get it to work with AvalancheMC()
  //AvalancheMC* aval = new AvalancheMC();
    aval->SetSensor(sensor);
  //aval->EnablePlotting(v_e);
  aval->EnableMagneticField();
    
    aval->EnablePlotting(driftView);
    
        for(int r_i=0; r_i < 9; r_i++){
            // Set the initial position [cm] and starting time [ns].
            x0 = 3.0 + 0.5*r_i,
            y0 = 0.0,
            z0 = 0.0,
            t0 = 0.0;
            if(r_i == 0){x0 = 3.1;}
            if(r_i == 8){x0 = 6.9;}
            
            cout << "r = " << x0 << endl;
            
            for(int eve=0;eve<1;eve++){
                cout << "Event number: " << eve << endl;
                aval->AvalancheElectron(x0, y0, z0, t0, e0, dx0, dy0, dz0);
                
                
                if(0<aval->GetNumberOfElectronEndpoints()){
                    for(int nava=0;nava<aval->GetNumberOfElectronEndpoints();nava++){
                        aval->GetElectronEndpoint(nava, x0, y0, z0, t0, e0, x1, y1, z1, t1, e1, status);
                        
                    }
                }
            }
            
            c_dLine->cd();
            driftView->Plot();
            
            
        }
    
    c_dLine->SaveAs("drift_lines.png");
    
  app.Run(kTRUE);

}
