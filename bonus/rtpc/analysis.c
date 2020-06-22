//////////////////////////////////////////
//
// Analysis of the BONuS RTPC
//
// Nate Dzbenski
// ndzbe001@odu.edu
//
// 16 May 2020
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
#include <TFile.h>
#include <TTree.h>

using namespace std;

char thisto_name[15];
char phihisto_name[15];

char tvec_name[15];
char phivec_name[15];

char thisto_title[15];
char phihisto_title[15];

char rleg_name[15];
char phileg_name[15];

int analysis() {

    double x0, z0 = 0;
    
//______________________________________________________________________________________________
//________________________________________ Canvas and plots ____________________________________
//______________________________________________________________________________________________

    /*TCanvas *c_driftT = new TCanvas("c_driftT","Drift time",1100,800);
    c_driftT->Divide(9,9);
    TCanvas *can_phi = new TCanvas("can_phi","Lorentz angle",1100,800);
    can_phi->Divide(9,9);*/
    
    //TCanvas *c_dt = new TCanvas("c_dt","Drift time",800,600);
    //TCanvas *c_dphi = new TCanvas("c_dphi","Drift Angle",800,600);
    
    //TCanvas *c_t_vs_phi = new TCanvas("c1","t_d vs phi_d",800,600);
    //TCanvas *c_sphi_vs_r = new TCanvas("c4","sphi_d vs r",800,600);
    //TCanvas *c_st_vs_r= new TCanvas("c5","st_d vs r",800,600);
    
    // Array of histos
    TH1D *h_driftT[9][9];
    TH1D *h_phi[9][9];
    
    // Array of Graphs
    TGraphErrors *gr_phi_vs_r[9];
    TGraphErrors *gr_t_vs_r[9];
    
    // Array of values
    double phi_vs_r[9];
    double t_vs_r[9];
    double st_vs_r[9];
    double sphi_vs_r[9];
    double a_t[9], b_t[9], c_t[9], a_phi[9], b_phi[9], c_phi[9], err_a_t[9], err_b_t[9], err_c_t[9], err_a_phi[9], err_b_phi[9], err_c_phi[9];
    double zeros[9] = {0,0,0,0,0,0,0,0,0};
    double the_rs[9] = {3.1,3.5,4.0,4.5,5.0,5.5,6.0,6.5,6.9};
    double the_zs[9] = {-19.0,-15.0,-10.0,-5.0,0.0,5.0,10.0,15.0,19.0};
    
    //TH1D *h_driftT = new TH1D("h_driftT","Drit Time [ns] -3500 V (He_80_CO2_20)",50,0,0);
    //TH1D *h_phi = new TH1D("h_phi","Drift angle [rad] -3500 V (He_80_CO2_20)",50,0,0);
    TH1D *h_energy = new TH1D("h_energy","Energy Loss -3500 V (He_80_CO2_20)",50,0,0);
    //TH1D *h_driftV = new TH1D("h_driftV","Drift Velocity",30,0,0);
    
    TF1 *gausfit = new TF1("gausfit","gaus",0.0,6000);
    TF1 *gausfit2 = new TF1("gausfit2","gaus",0.0,1.0);
    
    // Old equations
    //TF1  *f_t = new TF1("f_t","[0]*(7-x)+[1]*(7-x)*(7-x)",0,7);
    //TF1  *f_phi = new TF1("f_phi","[0]*(7-x)+[1]*(7-x)*(7-x)",0,7);
    
    // New fit equations
    //TF1  *f_t = new TF1("f_t","[0]+([1]/(2.0*[2]))*([2]-40.0+sqrt(4.0*[2]*(49-x*x)+pow((40.0-[2]),2)))",3.0,7.0);
    TF1  *f_t = new TF1("f_t","[0]+[1]*(0.5+(-40.0+sqrt(4.0*(0.1+[2]*[2])*(49-x*x)+pow((40.0-(0.1+[2]*[2])),2)))/(2*(0.1+[2]*[2])))",3.0,7.0);
    TF1  *f_phi = new TF1("f_phi","[0] + [1]*log(7.0/x) + [2]*((1.0/(x*x)) - (1.0/(49.0)))",0,3.14);
    
    // OLD 4th order fit
    //TF1  *fit_ap = new TF1("fit_ap","[0]*(x*x*x*x)+[1]*(x*x*x)+[2]*(x*x)+[3]*x+[4]",-20,20);
    //TF1  *fit_bp = new TF1("fit_bp","[0]*(x*x*x*x)+[1]*(x*x*x)+[2]*(x*x)+[3]*x+[4]",-20,20);
    //TF1  *fit_at = new TF1("fit_at","[0]*(x*x*x*x)+[1]*(x*x*x)+[2]*(x*x)+[3]*x+[4]",-20,20);
    //TF1  *fit_bt = new TF1("fit_bt","[0]*(x*x*x*x)+[1]*(x*x*x)+[2]*(x*x)+[3]*x+[4]",-20,20);
    
    // New parameter fits
    TF1  *fit_ap = new TF1("fit_ap","[0]+[1]*(x*x)+[2]*(x*x*x*x)",-20,20);
    TF1  *fit_bp = new TF1("fit_bp","[0]+[1]*(x*x)+[2]*(x*x*x*x)",-20,20);
    TF1  *fit_cp = new TF1("fit_cp","[0]+[1]*(x*x)+[2]*(x*x*x*x)",-20,20);
    TF1  *fit_at = new TF1("fit_at","[0]+[1]*(x*x)+[2]*(x*x*x*x)",-20,20);
    TF1  *fit_bt = new TF1("fit_bt","[0]+[1]*(x*x)+[2]*(x*x*x*x)",-20,20);
    TF1  *fit_ct = new TF1("fit_ct","[0]+[1]*(x*x)+[2]*(x*x*x*x)",-20,20);
    
    TMultiGraph *mgr_t = new TMultiGraph();
    TMultiGraph *mgr_phi = new TMultiGraph();
    
    TLegend *leg = new TLegend(0.7,0.4,0.9,0.88);
    TLegend *leg2 = new TLegend(0.7,0.4,0.9,0.88);
    
//_____________________________________________________________________________________________
//___________________________________________ Openings ________________________________________
//_____________________________________________________________________________________________
    TFile *iFile = TFile::Open("rtpc_out.root","READ");
    
    if (!iFile) { return 0; }
    
    TTree *tTree = (TTree*)iFile->Get("tTree");
    TTree *pTree = (TTree*)iFile->Get("pTree");
    
    
//______________________________________________________________________________________________
//_____________________________________________ Code ___________________________________________
//______________________________________________________________________________________________
    int r_i, z_i;
    
    for(int num = 0; num < 81; num++){
        if( num < 9 ) { z_i = 0; r_i = num; }
        else{
            r_i = num%9;
            z_i = int (num/9);
        }
        // Set the initial position [cm]
        x0 = the_rs[r_i];
        z0 = the_zs[z_i];
        
        vector<double> *drift_time = 0;
        vector<double> *drift_angle = 0;
        
        cout << "z = " << z0 << ", r = " << x0 << endl;
        
        sprintf(tvec_name,"bTd_%i_%i",z_i,r_i);
        sprintf(phivec_name,"bPhid_%i_%i",z_i,r_i);
        
        TBranch *bt = 0;
        TBranch *bphi = 0;
        pTree->SetBranchAddress(phivec_name,&drift_angle,&bphi);
        tTree->SetBranchAddress(tvec_name,&drift_time,&bt);
        
        bt->GetEntry(0);
        bphi->GetEntry(0);
    
        sprintf(thisto_name,"Td_%i_%i",z_i,r_i);
        sprintf(phihisto_name,"Phid_%i_%i",z_i,r_i);
        
        // create histos for this r,z bin
        h_driftT[z_i][r_i] = new TH1D(thisto_name, "Drit Time [ns] -3500 V (He_80_CO2_20)", 50, 0, 0);
        h_phi[z_i][r_i] = new TH1D(phihisto_name, "Drift angle [rad] -3500 V (He_80_CO2_20)", 50, 0, 0);
        
        sprintf(thisto_title,"Drift Time [z=%f, r=%f]",z0,x0);
        sprintf(phihisto_title,"Drift Angle [z=%f, r=%f]",z0,x0);
        
        for(int c = 0; c < drift_time->size(); c++){
            if(drift_time->at(c) > 10.0) { h_driftT[z_i][r_i]->Fill(drift_time->at(c)); }
        }
        for(int c = 0; c < drift_angle->size(); c++){
            if(drift_angle->at(c) > 0.005) { h_phi[z_i][r_i]->Fill(drift_angle->at(c)); }
        }
    }
    
    int num = 0;
    // begin the data analysis
    for( z_i = 0; z_i < 9; z_i++){
        for(r_i = 0; r_i < 9; r_i++)
        {
            x0 = the_rs[r_i];
            z0 = the_zs[z_i];
            
            sprintf(thisto_title,"Drift Time [z=%f, r=%f]",z0,x0);
            sprintf(phihisto_title,"Drift Angle [z=%f, r=%f]",z0,x0);
            
            //c_driftT->cd(num+1);
            //h_driftT[z_i][r_i]->Draw();
            
            gausfit->SetParameter(0, h_driftT[z_i][r_i]->GetMaximum());
            gausfit->SetParameter(1, h_driftT[z_i][r_i]->GetMean());
            gausfit->SetParameter(2, h_driftT[z_i][r_i]->GetRMS());
            h_driftT[z_i][r_i]->Fit("gausfit","B");
            h_driftT[z_i][r_i]->SetTitle(thisto_title);
            h_driftT[z_i][r_i]->GetXaxis()->SetTitle("t_{d} [ns]");
            
            t_vs_r[r_i] = gausfit->GetParameter(1);
            st_vs_r[r_i] = gausfit->GetParameter(2);
            
            //can_phi->cd(num+1);
            //h_phi[z_i][r_i]->Draw();
            gausfit2->SetParameter(0, h_phi[z_i][r_i]->GetMaximum());
            gausfit2->SetParameter(1, h_phi[z_i][r_i]->GetMean());
            gausfit2->SetParameter(2, h_phi[z_i][r_i]->GetRMS());
            h_phi[z_i][r_i]->Fit("gausfit2","BR");
            h_phi[z_i][r_i]->SetTitle(phihisto_title);
            h_phi[z_i][r_i]->GetXaxis()->SetTitle("#phi_{d} [rad]");
            
            phi_vs_r[r_i] = gausfit2->GetParameter(1);
            sphi_vs_r[r_i] = gausfit2->GetParameter(2);
            
            //can_phi->Update();
            
            /*if( num == 40 ){
                c_dt->cd(0);
                h_driftT[z_i][r_i]->Draw();
                h_driftT[z_i][r_i]->SetTitle(thisto_title);
                
                c_dphi->cd(0);
                h_phi[z_i][r_i]->Draw();
                h_phi[z_i][r_i]->SetTitle(phihisto_title);
                
            }*/
            num++;
        }
        gr_phi_vs_r[z_i] = new TGraphErrors(9,the_rs,phi_vs_r,zeros,sphi_vs_r);
        gr_phi_vs_r[z_i]->SetMarkerColor(1+z_i);
        gr_phi_vs_r[z_i]->SetMarkerStyle(21);
        gr_phi_vs_r[z_i]->SetMarkerSize(1);
        gr_phi_vs_r[z_i]->Fit("f_phi");
        
        gr_t_vs_r[z_i] = new TGraphErrors(9,the_rs,t_vs_r,zeros,st_vs_r);
        gr_t_vs_r[z_i]->SetMarkerColor(1+z_i);
        gr_t_vs_r[z_i]->SetMarkerStyle(21);
        gr_t_vs_r[z_i]->SetMarkerSize(1);
        gr_t_vs_r[z_i]->Fit("f_t");
        
        a_t[z_i] = f_t->GetParameter(0);
        b_t[z_i] = f_t->GetParameter(1);
        c_t[z_i] = f_t->GetParameter(2);
        a_phi[z_i] = f_phi->GetParameter(0);
        b_phi[z_i] = f_phi->GetParameter(1);
        c_phi[z_i] = f_phi->GetParameter(2);
        err_a_t[z_i] = f_t->GetParError(0)/2.0;
        err_b_t[z_i] = f_t->GetParError(1)/2.0;
        err_c_t[z_i] = f_t->GetParError(2)/2.0;
        err_a_phi[z_i] = f_phi->GetParError(0)/2.0;
        err_b_phi[z_i] = f_phi->GetParError(1)/2.0;
        err_c_phi[z_i] = f_phi->GetParError(2)/2.0;
        
        sprintf(rleg_name,"z = %.1f cm", z0);
        leg->AddEntry(gr_t_vs_r[z_i],rleg_name,"lep");
        
        sprintf(phileg_name,"z = %.1f cm", z0);
        leg2->AddEntry(gr_phi_vs_r[z_i],phileg_name,"lep");
        
        mgr_t->Add(gr_t_vs_r[z_i]);
        mgr_phi->Add(gr_phi_vs_r[z_i]);
    }
    
    //______________________________________________________________________________________________
    //___________________________________________ Display __________________________________________
    //______________________________________________________________________________________________

    TGraphErrors *gr_at_vs_z = new TGraphErrors(9,the_zs,a_t,zeros,err_a_t);
    TGraphErrors *gr_bt_vs_z = new TGraphErrors(9,the_zs,b_t,zeros,err_b_t);
    TGraphErrors *gr_ct_vs_z = new TGraphErrors(9,the_zs,c_t,zeros,err_c_t);
    TGraphErrors *gr_aphi_vs_z = new TGraphErrors(9,the_zs,a_phi,zeros,err_a_phi);
    TGraphErrors *gr_bphi_vs_z = new TGraphErrors(9,the_zs,b_phi,zeros,err_b_phi);
    TGraphErrors *gr_cphi_vs_z = new TGraphErrors(9,the_zs,c_phi,zeros,err_c_phi);

    
    //TCanvas *c_sphi_vs_r = new TCanvas("c4","sphi_d vs r",800,600);
    //TCanvas *c_st_vs_r= new TCanvas("c5","st_d vs r",800,600);
    
    TCanvas *c_param = new TCanvas("c_param","c_param",1200,800);
    c_param->Divide(3,2);
    c_param->cd(1);
    gr_at_vs_z->Draw("AP");
    gr_at_vs_z->SetTitle("a_t vs. z");
    gr_at_vs_z->GetXaxis()->SetTitle("z [cm]");
    gr_at_vs_z->GetYaxis()->SetTitle("a_t [ns]");
    gr_at_vs_z->GetXaxis()->SetTitleSize(0.05);
    gr_at_vs_z->GetYaxis()->SetTitleSize(0.05);
    gr_at_vs_z->SetMarkerStyle(21);
    gr_at_vs_z->SetMarkerSize(1);
    //gr_at_vs_z->GetYaxis()->SetTitleOffset(1.3);
    gr_at_vs_z->Fit("fit_at");
    c_param->cd(2);
    gr_bt_vs_z->Draw("AP");
    gr_bt_vs_z->SetTitle("b_t vs. z");
    gr_bt_vs_z->GetXaxis()->SetTitle("z [cm]");
    gr_bt_vs_z->GetYaxis()->SetTitle("b_t [ns]");
    gr_bt_vs_z->GetXaxis()->SetTitleSize(0.05);
    gr_bt_vs_z->GetYaxis()->SetTitleSize(0.05);
    gr_bt_vs_z->SetMarkerStyle(21);
    gr_bt_vs_z->SetMarkerSize(1);
    //gr_bt_vs_z->GetYaxis()->SetTitleOffset(1.3);
    gr_bt_vs_z->Fit("fit_bt");
    c_param->cd(3);
    gr_ct_vs_z->Draw("AP");
    gr_ct_vs_z->SetTitle("c_t vs. z");
    gr_ct_vs_z->GetXaxis()->SetTitle("z [cm]");
    gr_ct_vs_z->GetYaxis()->SetTitle("c_t");
    gr_ct_vs_z->GetXaxis()->SetTitleSize(0.05);
    gr_ct_vs_z->GetYaxis()->SetTitleSize(0.05);
    gr_ct_vs_z->SetMarkerStyle(21);
    gr_ct_vs_z->SetMarkerSize(1);
    //gr_ct_vs_z->GetYaxis()->SetTitleOffset(1.3);
    gr_ct_vs_z->Fit("fit_ct");
    c_param->cd(4);
    gr_aphi_vs_z->Draw("AP");
    gr_aphi_vs_z->SetTitle("a_phi vs. z");
    gr_aphi_vs_z->GetXaxis()->SetTitle("z [cm]");
    gr_aphi_vs_z->GetYaxis()->SetTitle("a_phi [rad]");
    gr_aphi_vs_z->GetXaxis()->SetTitleSize(0.05);
    gr_aphi_vs_z->GetYaxis()->SetTitleSize(0.05);
    gr_aphi_vs_z->SetMarkerStyle(21);
    gr_aphi_vs_z->SetMarkerSize(1);
    //gr_aphi_vs_z->GetYaxis()->SetTitleOffset(1.3);
    gr_aphi_vs_z->Fit("fit_bp");
    c_param->cd(5);
    gr_bphi_vs_z->Draw("AP");
    gr_bphi_vs_z->SetTitle("b_phi vs. z");
    gr_bphi_vs_z->GetXaxis()->SetTitle("z [cm]");
    gr_bphi_vs_z->GetYaxis()->SetTitle("b_phi [rad]");
    gr_bphi_vs_z->GetXaxis()->SetTitleSize(0.05);
    gr_bphi_vs_z->GetYaxis()->SetTitleSize(0.05);
    gr_bphi_vs_z->SetMarkerStyle(21);
    gr_bphi_vs_z->SetMarkerSize(1);
    //gr_bphi_vs_z->GetYaxis()->SetTitleOffset(1.3);
    gr_bphi_vs_z->Fit("fit_bp");
    c_param->cd(6);
    gr_cphi_vs_z->Draw("AP");
    gr_cphi_vs_z->SetTitle("c_phi vs. z");
    gr_cphi_vs_z->GetXaxis()->SetTitle("z [cm]");
    gr_cphi_vs_z->GetYaxis()->SetTitle("c_phi");
    gr_cphi_vs_z->GetXaxis()->SetTitleSize(0.05);
    gr_cphi_vs_z->GetYaxis()->SetTitleSize(0.05);
    gr_cphi_vs_z->SetMarkerStyle(21);
    gr_cphi_vs_z->SetMarkerSize(1);
    //gr_cphi_vs_z->GetYaxis()->SetTitleOffset(1.3);
    gr_cphi_vs_z->Fit("fit_cp");
    
    
    TCanvas *c_phi_vs_r = new TCanvas("c2","phi_d vs r",800,600);
    c_phi_vs_r->cd();
    mgr_phi->Draw("AP");
    mgr_phi->GetXaxis()->SetTitle("r [cm]");
    mgr_phi->GetYaxis()->SetTitle("#phi_{d} [rad]");
    mgr_phi->GetYaxis()->SetTitleOffset(1.3);
    leg2->Draw();
    //c_phi_vs_r->Update();
    
    TCanvas *c_t_vs_r = new TCanvas("c3","t_d vs r",800,600);
    c_t_vs_r->cd();
    mgr_t->Draw("AP");
    mgr_t->GetXaxis()->SetTitle("r [cm]");
    mgr_t->GetYaxis()->SetTitle("t_{d} [ns]");
    mgr_t->GetYaxis()->SetTitleOffset(1.3);
    leg->Draw();
    //c_t_vs_r->Update();
    
    //c_driftT->SaveAs("figs/drift_times.png");
    //can_phi->SaveAs("figs/drift_angles.png");
    c_phi_vs_r->SaveAs("figs/phi_vs_r.png");
    c_t_vs_r->SaveAs("figs/t_vs_r.png");
    c_param->SaveAs("figs/parameters.png");
    //c_dt->SaveAs("figs/td_histo.png");
    //c_dphi->SaveAs("figs/phid_histo.png");
    
    cout << "=================== Here they are, the parameters =====================" << endl;
    cout << "a1_t = " << fit_at->GetParameter(0) << endl;
    cout << "a2_t = " << fit_at->GetParameter(1) << endl;
    cout << "a3_t = " << fit_at->GetParameter(2) << endl;
    cout << "***************************************" << endl;
    cout << "b1_t = " << fit_bt->GetParameter(0) << endl;
    cout << "b2_t = " << fit_bt->GetParameter(1) << endl;
    cout << "b3_t = " << fit_bt->GetParameter(2) << endl;
    cout << "***************************************" << endl;
    cout << "c1_t = " << fit_ct->GetParameter(0) << endl;
    cout << "c2_t = " << fit_ct->GetParameter(1) << endl;
    cout << "c3_t = " << fit_ct->GetParameter(2) << endl;
    cout << "***************************************" << endl;
    cout << "a1_phi = " << fit_ap->GetParameter(0) << endl;
    cout << "a2_phi = " << fit_ap->GetParameter(1) << endl;
    cout << "a3_phi = " << fit_ap->GetParameter(2) << endl;
    cout << "***************************************" << endl;
    cout << "b1_phi = " << fit_bp->GetParameter(0) << endl;
    cout << "b2_phi = " << fit_bp->GetParameter(1) << endl;
    cout << "b3_phi = " << fit_bp->GetParameter(2) << endl;
    cout << "***************************************" << endl;
    cout << "b1_phi = " << fit_cp->GetParameter(0) << endl;
    cout << "b2_phi = " << fit_cp->GetParameter(1) << endl;
    cout << "b3_phi = " << fit_cp->GetParameter(2) << endl;
    cout << "======================================================================" << endl;
    
    //for(int i = 0; i < 9; i++){ cout << the_rs[i] << ", " << t_vs_r[i] << endl; }
    

    return 0;
}
