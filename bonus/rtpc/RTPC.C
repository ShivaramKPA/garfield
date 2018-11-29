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
// 08 Nov 2018
//
///////////////////////////////////////////

#include <iostream>
#include <stdio.h>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
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

  //TCanvas *c_field = new TCanvas("c_field", "c_field", 800, 600);
  TCanvas *c_driftT = new TCanvas("c_driftT","Drift time",800,600);
  //TCanvas *c_energy = new TCanvas("c_energy","Energy Loss",800,600);
  //TCanvas *c_driftV = new TCanvas("c_driftV","Drift Velocity",800,600);
  TCanvas *c_phi = new TCanvas("c_phi","Lorentz angle",800,600);

  ViewField *viewfield = new ViewField();


    TH1D *h_driftT = new TH1D("h_driftT","Drit Time [ns] -3500 V (He_80_CO2_20)",50,0,0);
    TH1D *h_energy = new TH1D("h_energy","Energy Loss -3500 V (He_80_CO2_20)",50,0,0);
    //TH1D *h_driftV = new TH1D("h_driftV","Drift Velocity",30,0,0);
    TH1D *h_phi = new TH1D("h_phi","Drift angle [rad] -3500 V (He_80_CO2_20)",50,0,0);


  TF1 *gausfit = new TF1("gausfit","gaus",0,8000);


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
  elm->LoadMagneticField("Fieldmaps/solenoid_map.dat", 1.0);



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
    
    
  for(int eve=0;eve<1;eve++){
    ne_tot=0;
    cout << "Event number: " << eve << endl;
    aval->AvalancheElectron(x0, y0, z0, t0, e0, dx0, dy0, dz0);
      
    // Get the number of electrons and ions in the avalanche.
    aval->GetAvalancheSize(ne, ni);
    ne_tot+=ne;
      
    if(0<aval->GetNumberOfElectronEndpoints()){
      for(int nava=0;nava<aval->GetNumberOfElectronEndpoints();nava++){
        aval->GetElectronEndpoint(nava, x0, y0, z0, t0, e0, x1, y1, z1, t1, e1, status);
          
        h_driftT->Fill(t1);
          h_energy->Fill(TMath::Abs(e0-e1));
        if(x1 != 0) h_phi->Fill(TMath::ATan2(y1,x1));
      }
    }
    //if(0<ne_tot) h_size->Fill(ne_tot);
   //   v_e->Plot(); 
  }
    cout << "----------------------------------------------------------" << endl;
	cout << "Start position x = " << x0 << ", z = " << z0 << endl;
    cout << "----------------------------------------------------------" << endl;

//______________________________________________________________________________________________
//___________________________________________ Display __________________________________________
//______________________________________________________________________________________________


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
    
 /* c_energy->cd();
    h_energy->Draw();
    h_energy->Fit("gausfit","Q");
    h_energy->GetXaxis()->SetTitle("Energy Loss [eV]");
cout << "Energy Loss" << endl;
cout << "Mean = " << gausfit->GetParameter(1) << ", Sigma = " << gausfit->GetParameter(2) << endl;
    c_energy->SaveAs("figs/eLoss.png");*/

  c_driftT->cd();
    h_driftT->Draw();
    h_driftT->Fit("gausfit","Q");
    h_driftT->GetXaxis()->SetTitle("t_{d} [ns]");
cout << "Drift Time" << endl;
cout << "Mean = " << gausfit->GetParameter(1) << ", Sigma = " << gausfit->GetParameter(2) << ", Error = " << gausfit->GetParError(2) << endl;
    c_driftT->SaveAs("figs/drift_time_19z.png");

  c_phi->cd();
    h_phi->Draw();
    h_phi->Fit("gausfit","Q");
    h_phi->GetXaxis()->SetTitle("#phi_{d} [rad]");
cout << "Drift angle" << endl;
cout << "Mean = " << gausfit->GetParameter(1) << ", Sigma = " << gausfit->GetParameter(2) << ", Error = " << gausfit->GetParError(2) << endl;
    c_phi->SaveAs("figs/drift_angle_19z.png");

  //c_pos->cd();
  //  h_pos->Draw();



  app.Run(kTRUE);

}
