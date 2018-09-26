#include <iostream>

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

using namespace Garfield;
using namespace std;

int main(int argc, char * argv[]) {
    
    TApplication app("app", &argc, argv);
    
    plottingEngine.SetDefaultStyle();
    
    // Canvas and Plotting
    TCanvas* c_field = new TCanvas("c_field", "c_field", 800, 600);
    TCanvas* c_driftV = new TCanvas("c_driftV","Drift Velocity",800,600);
    
    ViewField *viewfield = new ViewField();
    
    TH1D *h_driftV = new TH1D("h_driftV","Drift Velocity",30, 0, 0);
    
    TF1 *gausfit = new TF1("gausfit","gaus",0,8000);
    
    // Variables used to describe the geometry
    const double dWire = 0.05;//50.e-4;
    const double rAnode = 1.46;
    const double x_max = 8.;
    const double y_max = 8.;
    const double z_max = 4.;
    
    // Variables used to describe the field
    const double vCath = 10000.;
    const double vAnode = 0.;
    
    // Set the initial position [cm] and starting time [ns].
    double x0 = 0.5;
    double y0 = 0.1; // -(sqrt(rAnode * rAnode - x0 * x0) - 0.002);
    double z0 = 0.;
    double t0 = 0.;
    
    // Magnetic field
    const double Bx = 0.; // Tesla
    const double By = 0.; // Tesla
    const double Bz = 0.; // Tesla
    
    // Electric Field [V/cm]
    const double Ex = 10000.;
    const double Ey = 10000.;
    const double Ez = 0.;
    
    // Drift velocity variables
    Double_t vx, vy, vz;
    
    // Set the initial energy [eV].
    double e0 = 5.e3;
    // Set the initial direction (x, y, z).
    // In case of a null vector, the direction is randomized.
    double dx0 = 0., dy0 = 0., dz0 = 0.;
    // Calculate an electron avalanche
    int ne, ni;
    int ne_tot;
    // Electron information after the avalanche
    Double_t x1, y1, z1, t1, e1;
    Int_t status;
    
    // Setup the gas.
    MediumMagboltz* gas = new MediumMagboltz();
    //  gas->LoadGasFile("gasFiles/He_100_DME_00.gas");
    gas->SetComposition("He",80.,"CO2",20.);
    gas->SetTemperature(293.);
    gas->SetPressure(760.);
    gas->EnableDrift();                           // Allow for drifting in this medium
    gas->PrintGas();
    
    // Build the geometry.
    GeometrySimple* geo = new GeometrySimple();
    SolidBox* box = new SolidBox(0., 0., 0., x_max/2., y_max/2., z_max/2.);
    geo->AddSolid(box, gas);
    
    // Make a component with analytic electric field.
    ComponentAnalyticField* comp = new ComponentAnalyticField();
    comp->SetGeometry(geo);
    comp->AddWire(0., 0., dWire, vCath, "c", 100., 50., 19.3);
    //comp->AddWire(-0.1, 0., dWire, vAnode, "a", 100., 50., 19.3);
    comp->AddTube(rAnode, vAnode, 0, "a");
    comp->AddReadout("c");
    comp->AddReadout("a");
    
    // Make a sensor.
    Sensor* sensor = new Sensor();
    // Calculate the electric field using the Component object comp
    sensor->AddComponent(comp);
    // Request signal calculation for the electrode named "s",
    // using the weighting field provided by the Component object comp
    sensor->AddElectrode(comp, "c");
    sensor->AddElectrode(comp, "a");
    sensor->SetArea(-x_max/2., -y_max/2., -z_max/2., x_max/2., y_max/2., z_max/2.);
    sensor->ClearSignal();
    
    // Evaluate the number of electrons in the avalanche
    AvalancheMicroscopic* aval = new AvalancheMicroscopic(); // did not get it to work with AvalancheMC()
    //AvalancheMC* aval = new AvalancheMC();
    aval->SetSensor(sensor);
    // Switch on signal calculation.
    //aval->EnableSignalCalculation();
    // Do the drift line calculation in time steps of 50 ps.
    //aval->SetTimeSteps(0.05);
    
    // Start simulation
    sensor->ClearSignal();
    sensor->NewSignal();

    for(int eve=0;eve<10;eve++){
        ne_tot=0;
        cout << "Event number: " << eve << endl;
        aval->AvalancheElectron(x0, y0, z0, t0, e0, dx0, dy0, dz0);
        
        // Get the number of electrons and ions in the avalanche.
        //aval->GetAvalancheSize(ne, ni);
        //ne_tot+=ne;
        
        aval->DriftElectron(x0, y0, z0, t0, e0);
        
        
        if(0<aval->GetNumberOfElectronEndpoints()){
            for(int nava=0;nava<aval->GetNumberOfElectronEndpoints();nava++){
                aval->GetElectronEndpoint(nava, x0, y0, z0, t0, e0, x1, y1, z1, t1, e1, status);
                
                gas->ElectronVelocity(Ex, Ey, Ez, Bx, By, Bz, vx, vy, vz);
                
                h_driftV->Fill(vx);
                
                cout << "Here's vx:" << vx << ", vy: " << vy << ", vz: " << vz << endl;
            }
        }
    }
    
    
    viewfield->SetComponent(comp);
    viewfield->SetSensor(sensor);
    viewfield->SetCanvas((TCanvas*)c_field->cd());
    viewfield->SetWeightingFieldRange(0.0, 10000.0);
    
    //  Field plot
    c_field->cd();
    viewfield->PlotContour();
    
    // Plot Drift Velocity
    c_driftV->cd();
    h_driftV->Draw();
    h_driftV->Fit("gausfit","Q");
    h_driftV->GetXaxis()->SetTitle("Drift Velocity [cm/ns]");
    
    app.Run(kTRUE);
    
}

