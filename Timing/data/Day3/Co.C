#include <iostream>
#include <fstream>
#include "../Coincidences.C"
using namespace std;

Coincidences* Start(){
    Coincidences* C = new Coincidences("CO_night_20_24.root",3);
    C->SetAlphaBeta(0,0.0901,-26);
    C->SetAlphaBeta(1,0.0794,-63);
    C->SetAlphaBeta(2,0.001477,0);
    C->CalibrateHisto();
    
    return C;
}

double FWHM(TH1F* h, double e_min, double e_max, bool draw = false){
    float centroid = 0., width_sx = 0., width_dx = 0.;
    float max_bin_content = 0.;
    int max_i = 0;
    for(unsigned int i = 1; i < h->GetXaxis()->GetNbins(); i++){
        if(h->GetXaxis()->GetBinCenter(i) < e_min)
            continue;
        if(h->GetXaxis()->GetBinCenter(i) > e_max)
            break;
        if(h->GetBinContent(i) > max_bin_content){
            max_bin_content = h->GetBinContent(i);
            centroid = h->GetXaxis()->GetBinCenter(i);
            max_i = i;
        }            
    }
    
    for(unsigned int i = 0; i < max_i; i++){
        if(h->GetXaxis()->GetBinCenter(i) < e_min)
            continue;
        if(h->GetBinContent(i) > max_bin_content*0.5){
            width_sx = centroid - h->GetXaxis()->GetBinCenter(i);
            break;
        }
    }
    
    for(unsigned int i = max_i; i < h->GetXaxis()->GetNbins(); i++){
        if(h->GetXaxis()->GetBinCenter(i) > e_max)
            break;
        if(h->GetBinContent(i) < max_bin_content*0.5){
            width_dx = h->GetXaxis()->GetBinCenter(i) - centroid;
            break;
        }
    }
    
    if(draw){
        h->Draw();
        TLine* l = new TLine(centroid - width_sx,0.5*max_bin_content,centroid + width_dx,0.5*max_bin_content);
        l->Draw("SAME");
    }
    
    return width_sx + width_dx;
}


Coincidences* Co(){
    Coincidences* C = new Coincidences("CO_night_20_24.root",3);
    C->SetAlphaBeta(0,0.0901,-26);
    C->SetAlphaBeta(1,0.0794,-63);
    C->SetAlphaBeta(2,0.001477,0);
    C->CalibrateHisto();

    auto channel = new vector<int>{0,1};
    auto e_min = new vector<double>{0,0};
    auto e_max = new vector<double>{DBL_MAX,DBL_MAX};
    

    TH1F* h;
    double fwhm = 0.;
    double t_min = 10., t_max = 30.; //to be set
    ofstream lthr("lower_thr.txt");
    if(!lthr){
        cout << "Failed to open lower_thr.txt" << endl;
        return 0x0;
    }
    ofstream win("window.txt");
    if(!win){
        cout << "Failed to open window.txt" << endl;
        return 0x0;
    }
    ofstream win_n("window_n.txt");
    ofstream lthr_n("lower_thr_n.txt");
    
    for(int i=100; i <= 1000; i += 100) {
	e_min->at(0) = i;
	e_min->at(1) = i;
	h = C->Filter(channel,e_min)->GetHistogram(2);
    h->Rebin(4);
    fwhm = FWHM(h,t_min,t_max);
	cout<<i<<endl;
    lthr << i << '\t' << fwhm << "\t0\t" << h->GetXaxis()->GetBinWidth(1)/2.45 << endl;
    lthr_n << i << '\t' << h->GetEntries() << endl;
    }

    for(int i=50; i <= 1000; i += 100){ 
        e_min->at(0) = i;
        e_min->at(1) = i;
        e_max->at(0) = i + 100;
        e_max->at(1) = i + 100;
        h = C->Filter(channel,e_min,e_max)->GetHistogram(2);
        h->Rebin(4);
        fwhm = FWHM(h,t_min,t_max);
        cout<<i<<endl;
        win << i << '\t' << fwhm << "\t0\t" << h->GetXaxis()->GetBinWidth(1)/2.45 << endl;
        win_n << i << '\t' << h->GetEntries() << endl;
    }
    
    lthr.close();
    win.close();


    return C;
}



