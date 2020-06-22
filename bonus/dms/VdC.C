#include <iostream>

#include <iostream>
#include <stdio.h>
#include <string>
#include <random>

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

int num_events = 50;

int main(int argc, char * argv[]) {
    
    TApplication app("app", &argc, argv);
    
    plottingEngine.SetDefaultStyle();
    
    // Canvas and Plotting
    TH1D *h_dTime = new TH1D("h_dTime","Drift Time",200, 0, 4000);
    h_dTime->SetTitle("Drift time for 2160 V");
    TH1D *h_dTime2 = new TH1D("h_dTime2","Drift Time 2",200, 0, 4000);
    TH1D *h_numEl = new TH1D("h_numEl","Num e-",30, 0, 0);
    
    //TCanvas* c_efield = new TCanvas("c_efield", "c_efield", 800, 600);
    //TCanvas* c_vfield = new TCanvas("c_vfield", "c_vfield", 800, 600);
    //TCanvas* c_vprof = new TCanvas("c_vprofile", "c_vprof", 800, 600);
    TCanvas* c_prof = new TCanvas("c_profile", "c_prof", 800, 600);
    TCanvas* c_dLine = new TCanvas("c_dLine","Drift Lines",800,600);
    TCanvas* c_dTime = new TCanvas("c_dTime","Drift Time",800,600);
    TCanvas* c_numEl = new TCanvas("c_numEl","Num e-",800,600);
    
    TF1 *gausfit = new TF1("gausfit","gaus",800,1200);
    TF1 *gausfit2 = new TF1("gausfit2","gaus",2000,2600);
    
    ViewDrift* driftView = new ViewDrift();
    ViewField *viewfield = new ViewField();
    
    driftView->SetCanvas(c_dLine);
    
    char name[10];
    char viewprofile[10];
    const int n_fields=1;
    
    for(int q=0; q<n_fields;q++){
        sprintf(viewprofile,"field%i",q);
        
        ViewField *viewprofile = new ViewField();
        
        // Variables used to describe the geometry
        const int n_rows = 10;
        const int n_layers = 3;
        
        const double s_x = 0.6;
        const double s_y = 0.6;
        const double L_x = 10.0;
        const double L_y = 6.0;
        const double L_z = 4.;
        const double b = 1.6;
        const double s_1 = 0.6;
        const double s_2 = 0.6;
        const double L_cath = 3.2;
        
        const double x_0 = 0.5;
        
        const double r_a = 0.33;
        const double wall_d = 0.1;
        const double wall_space = 0.0001;
        const int n_wires = 29;
        const double c = 0.25;
        
        double L_e = s_1+((n_rows-1)*s_x)+s_2;
        double dist = 0.;
        
        // Variables to optimize
        //const double dWire =  0.2; //0.2 50.e-4+(100.*q)0.1;
        const double dWire =  0.2;
        const double dCath = 0.2;
        const double dAnode = 0.03; //0.1 25.e-4;
        
        // Variables used to describe the field
        const double vCath = -8000.0;
        
        // Magnetic field
        const double Bx = 0.; // Tesla
        const double By = 0.; // Tesla
        const double Bz = 0.; // Tesla
        
        // Electric Field [V/cm]
        const double Ex = (vCath-0.)/L_e; //-vAnode)/L_e;
        const double Ey = 0.;
        const double Ez = 0.;
        
        //const double vAnode = Ex*(-(r_a+dWire));
        const double vAnode = 2160;
        
        // Drift velocity variables
        Double_t vx, vy, vz;
        
        // Set the initial energy [eV].
        //double e0 = 73;
        // Set the initial direction (x, y, z).
        // In case of a null vector, the direction is randomized.
        double dx0 = 0., dy0 = 0., dz0 = 0.;
        // Calculate an electron avalanche
        int ne, ni;
        int ne_tot;
        
        Int_t status;
        
        // Setup the gas.
        MediumMagboltz* gas = new MediumMagboltz();
        gas->SetComposition("He",80.,"CO2",20.);
        //gas->SetComposition("Ar",90.,"CO2",10.);
        gas->SetTemperature(293.);
        gas->SetPressure(760.);
        gas->EnableDrift();                           // Allow for drifting in this medium
        gas->PrintGas();
        
        // Build the geometry.
        GeometrySimple* geo = new GeometrySimple();
        SolidBox* box = new SolidBox(L_x/2., L_y/2., L_z/2., L_x/2., L_y/2., L_z/2.);
        geo->AddSolid(box, gas);
        
        // Make a component with analytic electric field.
        ComponentAnalyticField* comp = new ComponentAnalyticField();
        comp->SetGeometry(geo);
        
        // Make a sensor.
        Sensor* sensor = new Sensor();
        // Calculate the electric field using the Component object comp
        sensor->AddComponent(comp);
        // Request signal calculation for the electrode named "s",
        // using the weighting field provided by the Component object comp
        
        // Create grounded planes at the edges of the box
        comp->AddPlaneX(0.,0.,"x_min");
        comp->AddPlaneX(L_x,0.,"x_max");
        comp->AddPlaneY(0.,0.,"y_min");
        comp->AddPlaneY(L_y,0.,"y_max");

        comp->AddReadout("x_min");
        comp->AddReadout("x_max");
        comp->AddReadout("y_min");
        comp->AddReadout("y_max");
        
        sensor->AddElectrode(comp, "x_min");
        sensor->AddElectrode(comp, "x_max");
        sensor->AddElectrode(comp, "y_min");
        sensor->AddElectrode(comp, "y_max");
        
        // Build grounded wall between anode and field-shaping wires
        for(int k=0; k<2; k++){
            // for the bottom
            if(k==0){
                for(int i=0; i<n_wires; i++){
                    sprintf(name,"wallB%i",i);
                    
                    dist = i*(wall_d+wall_space);
                    cout << " -------------------------------------------------------- " << endl;
                    cout << "Iteration: " << i << ", Distance: " << dist << endl;
                    cout << " -------------------------------------------------------- " << endl;
                    
                    comp->AddWire(x_0+2.*r_a, dist+0.1001, wall_d, 0., name, 100., 50., 19.3);
                    
                    comp->AddReadout(name);
                    sensor->AddElectrode(comp, name);
                }
            }
            // For the top
            else if(k==1){
                for(int i=0; i<n_wires; i++){
                    sprintf(name,"wallT%i",i);
                    
                    dist = L_y/2.+c/2.+(i*(wall_d+wall_space));
                    
                    cout << " -------------------------------------------------------- " << endl;
                    cout << "Iteration: " << i << ", Distance: " << dist << endl;
                    cout << " -------------------------------------------------------- " << endl;
                    
                    comp->AddWire(x_0+2.*r_a, dist, wall_d, 0., name, 100., 50., 19.3);
                    
                    comp->AddReadout(name);
                    sensor->AddElectrode(comp, name);
                }
                
            }
        }
        
        // Build enclosure around anode
        for(int k=0; k<3; k++){
            // bottom of enclosure
            if(k==0){
                for(int i=0; i<int((2.*r_a)/wall_d); i++){
                    sprintf(name,"encB%i",i);
                    
                    dist = x_0+2.*r_a-(i*(wall_d+wall_space));
                    
                    comp->AddWire(dist-wall_d, L_y/2.-r_a-2.*wall_space, wall_d, 0., name, 100., 50., 19.3);
                    
                    comp->AddReadout(name);
                    sensor->AddElectrode(comp, name);
                }
            }
            // side of enclosure
            else if(k==1){
                for(int i=0; i<6; i++){
                    sprintf(name,"encS%i",i);
                    
                    dist = L_y/2.-r_a-wall_space+wall_d+(i*(wall_d+wall_space));
                    
                    comp->AddWire(x_0+2.*r_a-((int((2.*r_a)/wall_d))*(wall_d+wall_space)), dist, wall_d, 0., name, 100., 50., 19.3);
                    
                    comp->AddReadout(name);
                    sensor->AddElectrode(comp, name);
                }
            }
            //top of enclosure
            else if(k==2){
                for(int i=0; i<int((2.*r_a)/wall_d); i++){
                    sprintf(name,"encT%i",i);
                    
                    dist = x_0+2.*r_a-(i*(wall_d+wall_space));
                    
                    comp->AddWire(dist-wall_d, L_y/2.+r_a+3.*wall_space, wall_d, 0., name, 100., 50., 19.3);
                    
                    comp->AddReadout(name);
                    sensor->AddElectrode(comp, name);
                }
            }
        }
        
        // Build the field wires
        for(int k=0; k<2; k++){
            if(k==0){
                for(int i=0; i<n_layers; i++){
                    for(int j=0; j<n_rows; j++){
                        sprintf(name,"cT%i%i",i,j);
                        
                        dist = x_0 + wall_d + 2.*r_a + s_1 + (j*s_x);
                        
                        double v_electrode = Ex*(s_1 + (j*s_x));
                        
                        comp->AddWire(dist, L_y/2.+b/2.+(i*s_y), dWire, v_electrode, name, 100., 50., 19.3);
                        
                        comp->AddReadout(name);
                        sensor->AddElectrode(comp, name);
                        
                    }
                }
            }
            else if(k==1){
                for(int i=0; i<n_layers; i++){
                    for(int j=0; j<n_rows; j++){
                        sprintf(name,"cB%i%i",i,j);
                        
                        dist = x_0 + wall_d + 2.*r_a + s_1 + (j*s_x);
                        
                        double v_electrode = Ex*(s_1 + (j*s_x));
                        comp->AddWire(dist, L_y/2.-b/2.-(i*s_y), dWire, v_electrode, name, 100., 50., 19.3);
                        
                        comp->AddReadout(name);
                        sensor->AddElectrode(comp, name);
                        
                    }
                }
            }
        }
        
        // Build the cathode
        dist = x_0 + 2.*r_a + wall_d + s_1 + (s_x*(n_rows-1)) + s_2;
        comp->AddWire(dist, L_y/2., dCath, vCath, "cath0", 100., 50., 19.3);
        comp->AddReadout("cath0");
        sensor->AddElectrode(comp, "cath0");
        for(int k=0; k<2; k++){
            for(int i=1; i<int(L_cath/(2.*dCath)); i++){
                if(k==0){
                    sprintf(name,"cathT%i",i);
                    comp->AddWire(dist, L_y/2.+((dCath+wall_space)*i), dCath, vCath, name, 100., 50., 19.3);
                    
                    comp->AddReadout(name);
                    sensor->AddElectrode(comp, name);
                }
                else if(k==1){
                    sprintf(name,"cathB%i",i);
                    
                    comp->AddWire(dist, L_y/2.-((dCath+wall_space)*i), dCath, vCath, name, 100., 50., 19.3);
                    
                    comp->AddReadout(name);
                    sensor->AddElectrode(comp, name);
                }
                
            }
        }
        
        // Try a different way to build the cathode
        //SolidBox* cathode = new SolidBox(L_x/2., L_y/2., L_z/2., L_x/2., L_y/2., L_z/2.);

        // Anode
        comp->AddWire(x_0 + r_a + wall_d/2., L_y/2, dAnode, vAnode, "a", 100., 50., 19.3);
        comp->AddReadout("a");
        sensor->AddElectrode(comp, "a");
        
        sensor->SetArea(0., 0., 0., L_x, L_y, L_z);
        driftView->SetArea(0., 0., 0., L_x, L_y, L_z);
        sensor->ClearSignal();
        
        // Evaluate the number of electrons in the avalanche
        AvalancheMicroscopic* aval = new AvalancheMicroscopic(); // did not get it to work with AvalancheMC()
        //AvalancheMC* aval = new AvalancheMC();
        aval->SetSensor(sensor);
        //aval->SetSensor(driftView);
        // Switch on signal calculation.
        //aval->EnableSignalCalculation();
        // Do the drift line calculation in time steps of 50 ps.
        //aval->SetTimeSteps(0.05);
        
        AvalancheMicroscopic* aval2 = new AvalancheMicroscopic();
        aval2->SetSensor(sensor);
        
        aval->EnablePlotting(driftView);
        aval2->EnablePlotting(driftView);
        
        cout << "V_anode = " << vAnode << endl;
        cout << "V_cathode = " << vCath << endl;
        cout << "Ex = " << Ex << endl;
        
        for(int eve=0;eve<num_events;eve++){
            // Start simulation
            sensor->ClearSignal();
            sensor->NewSignal();
            
            // Set the initial position [cm] and starting time [ns].
            double x0 = 2.9 + ((rand()%10)/10.0)*0.6;   // 3.2 middle
            double y0 = 3.0; // -(sqrt(rAnode * rAnode - x0 * x0) - 0.002);
            double z0 = 2.0;
            double t0 = 0.0;
            double e0 = 73.0;
            
            double x0_2 = 6.9 + ((rand()%10)/10.0)*0.6; // 7.2 middle
            double y0_2 = 3.0; // -(sqrt(rAnode * rAnode - x0 * x0) - 0.002);
            double z0_2 = 2.0;
            double t0_2 = 0.0;
            double e0_2 = 73.0;
            
            // Electron information after the avalanche
            double x1 = 0.0;
            double y1 = 0.0;
            double z1 = 0.0;
            double t1 = 0.0;
            double e1 = 0.0;
            double x1_2 = 0.0;
            double y1_2 = 0.0;
            double z1_2 = 0.0;
            double t1_2 = 0.0;
            double e1_2 = 0.0;
            
            ne_tot=0;
            
            cout << "Event number: " << eve << endl;
            cout << "x0 near: " << x0 << ", x0 far: " << x0_2 << endl;
            
            aval->AvalancheElectron(x0, y0, z0, t0, e0, 0.0, 0.0, 0.0);
            aval2->AvalancheElectron(x0_2, y0_2, z0_2, t0_2, e0_2, 0.0, 0.0, 0.0);
            
            // Get the number of electrons and ions in the avalanche.
            aval->GetAvalancheSize(ne, ni);
            ne_tot+=ne;
            
            //aval->DriftElectron(x0, y0, z0, t0, e0);
            
            h_numEl->Fill(ne);
            
            if(0<aval->GetNumberOfElectronEndpoints()){
                for(int nava=0;nava<aval->GetNumberOfElectronEndpoints();nava++){
                    
                    aval->GetElectronEndpoint(nava, x0, y0, z0, t0, e0, x1, y1, z1, t1, e1, status);
                    
                    //gas->ElectronVelocity(Ex, Ey, Ez, Bx, By, Bz, vx, vy, vz);
                    cout << "Time = " << t1 << endl;
                    if(t1 > 20){ h_dTime->Fill(t1); }
                    
                    //cout << "Here's vx:" << vx << ", vy: " << vy << ", vz: " << vz << endl;
                }
            }
            
            if(0<aval2->GetNumberOfElectronEndpoints()){
                for(int nava=0;nava<aval2->GetNumberOfElectronEndpoints();nava++){
                    x1 = y1 = z1 = t1 = e1 = 0;
                    aval2->GetElectronEndpoint(nava, x0_2, y0_2, z0_2, t0_2, e0_2, x1_2, y1_2, z1_2, t1_2, e1_2, status);
                    cout << "Time 2 = " << t1_2 << endl;
                    if(t1_2 > 20){ h_dTime2->Fill(t1_2); }
                }
            }
        }
        
        c_dLine->cd();
        driftView->Plot();
        
        viewprofile->SetComponent(comp);
        viewprofile->SetSensor(sensor);
        viewprofile->SetWeightingFieldRange(0.0, 3000.);
        viewprofile->SetElectricFieldRange(0.,3000.);
        
        //  Profile plot
        /*viewprofile->SetCanvas((TCanvas*)c_vprof->cd());
        c_vprof->cd();
        viewprofile->PlotProfile(0.,L_y/2.,L_z/2.,L_x,L_y/2.,L_z/2.,"v");*/
        
        viewprofile->SetCanvas((TCanvas*)c_prof->cd());
        c_prof->cd();
        viewprofile->PlotProfile(0.,L_y/2.,L_z/2.,L_x,L_y/2.,L_z/2.,"e");
        
        // Plot Drift time
        c_dTime->cd();
        h_dTime->Draw();
        h_dTime2->Draw("same");
        
        gausfit->SetParameters(1, h_dTime->GetMean());
        gausfit->SetParameters(2, h_dTime->GetRMS());
        gausfit2->SetParameters(1, h_dTime2->GetMean());
        gausfit2->SetParameters(2, h_dTime2->GetRMS());
        
        h_dTime->Fit("gausfit","B","R");
        h_dTime2->Fit("gausfit2","B","R");
        h_dTime->GetXaxis()->SetTitle("Drift time [ns]");
        c_dTime->SaveAs("t1_histo.png");
        
        c_numEl->cd();
        h_numEl->Draw();
        h_numEl->Fit("gausfit","Q");
        h_numEl->GetXaxis()->SetTitle("Num e- in avalanche");

        /*viewfield->SetComponent(comp);
        viewfield->SetSensor(sensor);
        viewfield->SetWeightingFieldRange(0.0, 3000.);
        viewfield->SetElectricFieldRange(0.,3000.);*/
        
    }
    
    
    
    
    //  V-Contour plot
    /*viewfield->SetCanvas((TCanvas*)c_vfield->cd());
    c_vfield->cd();
    viewfield->PlotContour("v");
    //viewfield->PlotProfile(0.,L_y/2.,L_z/2.,10.,L_y/2.,L_z/2.,"field");
    c_vfield->SetTitle("Potential Map");
    
    //  E-field plot
    viewfield->SetCanvas((TCanvas*)c_efield->cd());
    viewfield->SetElectricFieldRange(0.,3000.);
    c_efield->cd();
    viewfield->PlotContour("field");*/

    
    app.Run(kTRUE);
    
}

