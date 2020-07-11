//////////////////////////////////////////
//
// Simulation of the BONuS RTPC GEM regions
//
// Nate Dzbenski
// ndzbe001@odu.edu
//
// Gabriel Charles
// gcharles@odu.edu
//
// 08 July 2020
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

char thisto_name[15];
char phihisto_name[15];

char thisto_title[15];
char phihisto_title[15];

char rleg_name[15];
char phileg_name[15];

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);
    
    if ( !argv[1] || !argv[2] ){
        cout << "You need to execute the program with a region specification and number of events to simulate with: " << endl;
        cout << "./garf_rtpc [region number] [number of events]" << endl;
        cout << "Where [region number] is 1, 2, or 3 and [number of events] is the number of events you want to simulate." << endl;
        exit(EXIT_FAILURE);
    }
    
    int region = atoi(argv[1]);
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

    TCanvas *c_driftT = new TCanvas("c_driftT","Drift time",1100,800);
    c_driftT->Divide(3,3);
    TCanvas *c_phi = new TCanvas("c_phi","Lorentz angle",1100,800);
    c_phi->Divide(3,3);
    
    TCanvas *c_phi_vs_z = new TCanvas("c2","phi_d vs r",800,600);
    TCanvas *c_t_vs_z = new TCanvas("c3","t_d vs r",800,600);

    ViewField *viewfield = new ViewField();

    // Array of histos
    TH1D *h_driftT[9];
    TH1D *h_phi[9];
    
    // Array of values
    double phi_vs_z[9];
    double t_vs_z[9];
    double st_vs_z[9];
    double sphi_vs_z[9];
    double zeros[9] = {0,0,0,0,0,0,0,0,0};
    double the_zs[9] = {-19.0,-15.0,-10.0,-5.0,0.0,5.0,10.0,15.0,19.0};
    
    TF1 *gausfit = new TF1("gausfit","gaus",0.0,6000);
    TF1 *gausfit2 = new TF1("gausfit2","gaus",0.0,1.0);
    
    TLegend *leg = new TLegend(0.7,0.4,0.9,0.88);
    TLegend *leg2 = new TLegend(0.7,0.4,0.9,0.88);
    
    // z loop
    for(int i=0; i < 9; i++){
             sprintf(thisto_name,"Td_%i",i);
             sprintf(phihisto_name,"Phid_%i",i);
            
            h_driftT[i] = new TH1D(thisto_name, "Drit Time [ns] -3500 V (He_80_CO2_20)", 50, 0, 0);
            h_phi[i] = new TH1D(phihisto_name, "Drift angle [rad] -3500 V (He_80_CO2_20)", 50, 0, 0);
    }


//_____________________________________________________________________________________________
//___________________________________________ Openings ________________________________________
//_____________________________________________________________________________________________

  // Setup the gas.
  MediumMagboltz* gas = new MediumMagboltz();
  gas->SetComposition("He",85.5,"CO2",14.5);
  gas->SetTemperature(293.);
  gas->SetPressure(760.);
    gas->EnableDrift();                           // Allow for drifting in this medium
  gas->PrintGas();
    
    float gem_x;
    
    // Assemble a Sensor object
    Sensor* sensor = new Sensor();
    sensor->SetArea(-8.0,-8.0,-20.0,8.0,8.0,20.0);
    
    if (region == 1){
        gem_x = 7.02;
        
        ComponentElmer * elm_trans1 = new ComponentElmer("trans1/mesh.header","trans1/mesh.elements","trans1/mesh.nodes",
                                                         "trans1/dielectrics.dat","trans1/trans1.result","cm");
        elm_trans1->SetMedium(0,gas);
        elm_trans1->LoadMagneticField("Fieldmaps/solenoid_map_may2019.dat", 0.7893);
        
        sensor->AddComponent(elm_trans1);
    }
    else if (region == 2){
        gem_x = 7.32;
        
        ComponentElmer * elm_trans2 = new ComponentElmer("trans2/mesh.header","trans2/mesh.elements","trans2/mesh.nodes",
                                                         "trans2/dielectrics.dat","trans2/trans2.result","cm");
        elm_trans2->SetMedium(0,gas);
        elm_trans2->LoadMagneticField("Fieldmaps/solenoid_map_may2019.dat", 0.7893);
        
        sensor->AddComponent(elm_trans2);
    }
    
    else if (region == 3){
        gem_x = 7.62;
        
        ComponentElmer * elm_trans3 = new ComponentElmer("trans3/mesh.header","trans3/mesh.elements","trans3/mesh.nodes",
                                                         "trans3/dielectrics.dat","trans3/trans3.result","cm");
        elm_trans3->SetMedium(0,gas);
        elm_trans3->LoadMagneticField("Fieldmaps/solenoid_map_may2019.dat", 0.7893);
        
        sensor->AddComponent(elm_trans3);
    }
    else{
        cout << "Invalid region" << endl;
        cout << "Enter a value for region 1, 2, or 3." << endl;
        exit(EXIT_FAILURE);
    }
    
  // Evaluate the number of electrons in the avalanche
  AvalancheMicroscopic* aval = new AvalancheMicroscopic(); // did not get it to work with AvalancheMC()
    aval->SetSensor(sensor);
  aval->EnableMagneticField();
    
    // begin the data analysis
    for(int z_i=0; z_i < 9; z_i++){
        
        z0 = the_zs[z_i];
        
        cout << "z = " << z0 << endl;
        
            int r_i = 0;
            x0 = gem_x;
            y0 = 0.0;
            t0 = 0.0;
            
            cout << "r = " << x0 << endl;
            
            for(int eve=0;eve<atoi(argv[2]);eve++){
                ne_tot=0;
                cout << "Event number: " << eve << endl;
                aval->AvalancheElectron(x0, y0, z0, t0, e0, dx0, dy0, dz0);
                
                // Get the number of electrons and ions in the avalanche.
                aval->GetAvalancheSize(ne, ni);
                ne_tot+=ne;
                
                if(0<aval->GetNumberOfElectronEndpoints()){
                    for(int nava=0;nava<aval->GetNumberOfElectronEndpoints();nava++){
                        aval->GetElectronEndpoint(nava, x0, y0, z0, t0, e0, x1, y1, z1, t1, e1, status);
                        
                        if(t1>5) h_driftT[z_i]->Fill(t1);
                        if(x1 > 0.01) h_phi[z_i]->Fill(TMath::ATan2(y1,x1));
                    }
                }
            }
        
            sprintf(thisto_title,"Drift Time [z=%f, r=%f]",z0,x0);
            sprintf(phihisto_title,"Drift Angle [z=%f, r=%f]",z0,x0);
            
            c_driftT->cd(z_i+1);
            h_driftT[z_i]->Draw();
            gausfit->SetParameter(0, h_driftT[z_i]->GetMaximum());
            gausfit->SetParameter(1, h_driftT[z_i]->GetMean());
            gausfit->SetParameter(2, h_driftT[z_i]->GetRMS());
            h_driftT[z_i]->Fit("gausfit","B");
            h_driftT[z_i]->SetTitle(thisto_title);
            h_driftT[z_i]->GetXaxis()->SetTitle("t_{d} [ns]");
            cout << "Drift Time max: " << h_driftT[z_i]->GetMaximum() << endl;
            
            t_vs_z[z_i] = gausfit->GetParameter(1);
            st_vs_z[z_i] = gausfit->GetParameter(2);
            
            cout << "Mean = " << gausfit->GetParameter(1) << ", Sigma = " << gausfit->GetParameter(2) << ", Error = " << gausfit->GetParError(2) << endl;
            c_driftT->Update();
        
            c_phi->cd(z_i+1);
            h_phi[z_i]->Draw();
            gausfit2->SetParameter(0, h_phi[z_i]->GetMaximum());
            gausfit2->SetParameter(1, h_phi[z_i]->GetMean());
            gausfit2->SetParameter(2, h_phi[z_i]->GetRMS());
            h_phi[z_i]->Fit("gausfit2","BR");
            h_phi[z_i]->SetTitle(phihisto_title);
            h_phi[z_i]->GetXaxis()->SetTitle("#phi_{d} [rad]");
            cout << "Drift angle: " << h_phi[z_i]->GetMaximum() << endl;
            
            phi_vs_z[z_i] = gausfit2->GetParameter(1);
            sphi_vs_z[z_i] = gausfit2->GetParameter(2);
            
            cout << "Mean = " << gausfit2->GetParameter(1) << ", Sigma = " << gausfit2->GetParameter(2) << ", Error = " << gausfit2->GetParError(2) << endl;
            c_phi->Update();
        
    }
    
    
    //______________________________________________________________________________________________
    //___________________________________________ Display __________________________________________
    //______________________________________________________________________________________________

    
    // Graphs of t and phi vs z
    TGraphErrors *gr_phi_vs_z = new TGraphErrors(9,the_zs,phi_vs_z,zeros,sphi_vs_z);
    TGraphErrors *gr_t_vs_z = new TGraphErrors(9,the_zs,t_vs_z,zeros,st_vs_z);
    
    c_phi_vs_z->cd();
    gr_phi_vs_z->Draw("AP");
    gr_phi_vs_z->GetXaxis()->SetTitle("r [cm]");
    gr_phi_vs_z->GetYaxis()->SetTitle("#phi_{d} [rad]");
    gr_phi_vs_z->GetYaxis()->SetTitleOffset(1.3);
    c_phi_vs_z->Update();
    
    c_t_vs_z->cd();
    gr_t_vs_z->Draw("AP");
    gr_t_vs_z->GetXaxis()->SetTitle("r [cm]");
    gr_t_vs_z->GetYaxis()->SetTitle("t_{d} [ns]");
    gr_t_vs_z->GetYaxis()->SetTitleOffset(1.3);
    c_t_vs_z->Update();
    
    
    c_driftT->SaveAs("figs/drift_times.png");
    c_phi->SaveAs("figs/drift_angles.png");
    c_phi_vs_z->SaveAs("figs/phi_vs_r.png");
    c_t_vs_z->SaveAs("figs/t_vs_r.png");
 
  app.Run(kTRUE);

}
