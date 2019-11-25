#include <iostream>
#include <fstream>
#include "../Coincidences.C"
using namespace std;

Coincidences* Co(){
    Coincidences* C = new Coincidences("CO_night_20_24.root",3);
    C->SetAlphaBeta(0,0.0901,-26);
    C->SetAlphaBeta(1,0.0794,-63);
    C->SetAlphaBeta(2,0.001477,0);
    C->CalibrateHisto();

    auto channel = new vector<int>{0,1};
    auto e_min = new vector<double>{0,0};
    auto e_max = new vector<double>{DBL_MAX,DBL_MAX};
    

    Coincidences* D;
    double fwhm = 0.;
    double t_min = 0., t_max = 0.; //to be set
    
    for ( int i=1; i <= 10; i++ ) {
	e_min->at(0) = 100*i;
	e_min->at(1) = 100*i;
	D = C->Filter(channel,e_min);
    D->GetHistogram(2)->Rebin(2);
    fwhm = FWHM(D->GetHistogram(2),t_min,t_max)
	cout<<100*i<<endl;
    }

    for ( int i=1; i <= 10; i++ ) { 
        e_min->at(0) = 100*i;
        e_min->at(1) = 100*i;
	e_max->at(0) = 100*(i + 1);
	e_max->at(1) = 100*(i + 1);
        D = C->Filter(channel,e_min,e_max);
    }


    return C;
}

double FWHM(TH1F* h, double e_min, double e_max){
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
    
    return width_sx + width_dx;
}
