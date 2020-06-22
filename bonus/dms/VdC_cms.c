#include <iostream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TMath.h>

#include "ViewSignal.hh"
#include "ComponentAnalyticField.hh"
#include "MediumMagboltz.hh"
#include "SolidTube.hh"
#include "SolidBox.hh"
#include "GeometrySimple.hh"
#include "Sensor.hh"
#include "AvalancheMC.hh"
#include "AvalancheMicroscopic.hh"
#include "TrackHeed.hh"
#include "Random.hh"
#include "Plotting.hh"

using namespace Garfield;



int main() {
    
    
    double L_y=7.2
    double cath_y=3.2
    double d_wires=0.01
    double r_anode=0.5
    double s_x=0.6
    double s_y=0.6
    double s_1=0.583
    double L_ydrft=1.6
    double c=0.25
    double d_anode=0.005
    double d_fieldform=0.2
    double L_c=0.5
    double N_rows=10
    double cath_x=0.5
    double sdriftgr=6.3
    double sdriftkl=2.1
    double s_2=0.596
    double L_x={s_1}+({N_rows}-1)*{s_x}+{s_2}+{cath_x}+{L_c}
    double th_wall=0.1
    double dist_2_wires={d_wires}+0.00001
    
    // Electronics
    double V_anode=1800
    double V_cathode=-11000
    double rn=2000000
    double r1=2000000
    double r2=2000000
    double abs_R=r1+r2+({N_rows}-1)*rn
    
    // Wire Numbers
    double anztrennla=entier(({L_y}/2-c/2-th_wall/2-sqrt({((th_wall/2+dist_2_wires/2)^2-((th_wall-d_wires)/2)^2)}))/{dist_2_wires})
    double anzkathode=entier({cath_y/2}/{dist_2_wires})
    double anzanodkam=entier(2*r_anode/{dist_2_wires})
    double anzkathbr=entier({cath_x}/{dist_2_wires})
    
    ////////////////////////
    &CELL
    ////////////////////////
    OPTIONS ISOMETRIC
    OPTIONS LAYOUT
    CELL-IDENTIFIER "s_2:{s_2}cm"
    plane y={L_y/2}, V=0
    plane y=-{L_y/2}, V=0
    plane x={r_anode+th_wall+L_x}, V=0
    plane x=-{r_anode}, V=0
    rows
    // Anode Wire
    a 1 {d_anode} 0 0 {V_anode}
    // Field Shaping Wires
    b 3 {d_fieldform} {r_anode}+{th_wall}+{s_1}+({N_rows}-1)*{s_x} {L_ydrft}/2+i*{s_y} {V_cathode}-{r1}*{V_cathode}/{abs_R}
    u 3 {d_fieldform} {r_anode}+{th_wall}+{s_1}+({N_rows}-1)*{s_x}-{L_ydrft}/2-i*{s_y} {V_cathode}-{r1}*{V_cathode}/{abs_R}
    d 3 {d_fieldform} {r_anode}+{th_wall}+{s_1}+({N_rows}-2)*{s_x} {L_ydrft}/2+i*{s_y} {V_cathode}-{r1}*{V_cathode}/{abs_R}-{r2}*{V_cathode}/{abs_R}
    e 3 {d_fieldform} {r_anode}+{th_wall}+{s_1}+({N_rows}-2)*{s_x}-{L_ydrft}/2-i*{s_y} ...
    {V_cathode}-{r1}*{V_cathode}/{abs_R}-{r2}*{V_cathode}/{abs_R}
    f {N_rows}-2 {d_fieldform} {r_anode}+{th_wall}+{s_1}+i*{s_x}{L_ydrft}/2+2*{s_y} (i+1)*{rn}*{V_cathode}/{abs_R}
    g {N_rows}-2 {d_fieldform} {r_anode}+{th_wall}+{s_1}+i*{s_x}{L_ydrft}/2+{s_y} (i+1)*{rn}*{V_cathode}/{abs_R}
    h {N_rows}-2 {d_fieldform} {r_anode}+{th_wall}+{s_1}+i*{s_x}{L_ydrft}/2 (i+1)*{rn}*{V_cathode}/{abs_R}
    j {N_rows}-2 {d_fieldform} {r_anode}+{th_wall}+{s_1}+i*{s_x}-{L_ydrft}/2 (i+1)*{rn}*{V_cathode}/{abs_R}
    k {N_rows}-2 {d_fieldform} {r_anode}+{th_wall}+{s_1}+i*{s_x}-{L_ydrft}/2-{s_y} (i+1)*{rn}*{V_cathode}/{abs_R}
    l {N_rows}-2 {d_fieldform} {r_anode}+{th_wall}+{s_1}+i*{s_x}-{L_ydrft}/2-2*{s_y} (i+1)*{rn}*{V_cathode}/{abs_R}
    
    //* Thickness of the Wall //*
    m entier({anztrennla/4}) {d_wires} {r_anode}+{d_wires}/2{c}/2+{th_wall/2}+sqrt({((th_wall/2+dist_2_wires/2)^2-((th_wall-d_wires)/2)^2)})+i*{dist_2_wires} 0
    m {anztrennla} {d_wires} {r_anode}+{th_wall}-{d_wires}/2{c}/2+{th_wall/2}+sqrt({((th_wall/2+dist_2_wires/2)^2-((th_wall-d_wires)/2)^2)})+i*{dist_2_wires} 0
    m 1 {th_wall} {r_anode+th_wall/2} {c/2+th_wall/2} 0
    n entier({anztrennla/4}) {d_wires} {r_anode}+{d_wires}/2-({c}/2+{th_wall/2}+sqrt({((th_wall/2+dist_2_wires/2)^2-((th_wall-d_wires)/2)^2)})+i*{dist_2_wires}) 0
    n {anztrennla} {d_wires} {r_anode}+{th_wall}-{d_wires}/2-({c}/2+{th_wall/2}+sqrt({((th_wall/2+dist_2_wires/2)^2-((th_wall-d_wires)/2)^2)})+i*{dist_2_wires}) 0
    n 1 {th_wall} {r_anode+th_wall/2} -{c/2+th_wall/2} 0
    o {anzkathode} {d_wires}{r_anode}+{th_wall}+{s_1}+({N_rows}-1)*{s_x}+{s_2}+{d_wires/2}{d_wires/2}+i*{dist_2_wires} {V_cathode}
    
    //* Cathode //*
    t {anzkathode} {d_wires}{r_anode}+{th_wall}+{s_1}+({N_rows}-1)*{s_x}+{s_2}+{d_wires/2}-{d_wires/2}-i*{dist_2_wires} {V_cathode}
    v {anzkathode} {d_wires}{r_anode}+{th_wall}+{s_1}+({N_rows}-1)*{s_x}+{s_2}+{cath_x}-{d_wires/2}{d_wires/2}+i*{dist_2_wires} {V_cathode}
    w {anzkathode} {d_wires}{r_anode}+{th_wall}+{s_1}+({N_rows}-1)*{s_x}+{s_2}+{cath_x}-{d_wires/2}-{d_wires/2}-i*{dist_2_wires} {V_cathode}
    x 1 {cath_x-2*d_wires-0.00001}{r_anode}+{th_wall}+{s_1}+({N_rows}-1)*{s_x}+{s_2}+{cath_x/2}{cath_y/2} {V_cathode}
    y 1 {cath_x-2*d_wires-0.00001}{r_anode}+{th_wall}+{s_1}+({N_rows}-1)*{s_x}+{s_2}+{cath_x/2}-{cath_y/2} {V_cathode}
    
    //* Anode Chamber //*
    q {anzanodkam} {d_wires} {r_anode}+{d_wires}/2-(i+1)*{dist_2_wires}{r_anode}+{d_wires}/2 0}
    r {anzanodkam} {d_wires} {r_anode}+{d_wires}/2-(i+1)*{dist_2_wires}-{r_anode}-{d_wires}/2 0

////////////////////////**
&FIELD
////////////////////////**
track 1.6 0 6.8 0
plot-field graph ex
&STOP

