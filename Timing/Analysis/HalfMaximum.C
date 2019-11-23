#include <iostream>
using namespace std;

void HalfMaximum(TH1F* h, int ch_min, int ch_max){
    float centroid = 0., width = 0.;
    float max_bin_content = 0.;
    int max_i = 0;
    for(unsigned int i = 1; i < h->GetXaxis()->GetNbins(); i++){
        if(h->GetXaxis()->GetBinCenter(i) < ch_min)
            continue;
        if(h->GetXaxis()->GetBinCenter(i) > ch_max)
            break;
        if(h->GetBinContent(i) > max_bin_content){
            max_bin_content = h->GetBinContent(i);
            centroid = h->GetXaxis()->GetBinCenter(i);
            max_i = i;
        }            
    }
    for(unsigned int i = max_i; i < h->GetXaxis()->GetNbins(); i++){
        if(h->GetXaxis()->GetBinCenter(i) > ch_max)
            break;
        if(h->GetBinContent(i) < max_bin_content*0.5){
            width = h->GetXaxis()->GetBinCenter(i) - centroid;
            break;
        }
    }
    
    // Drawing lines on histogram
    TLine* lc = new TLine(centroid,0,centroid,max_bin_content*1.08);
    lc->Draw("SAME");
    TLine* lw = new TLine(centroid + width,0,centroid + width,max_bin_content*0.5*1.4);
    lw->Draw("SAME");
    
    cout << "Centroid, width = " << centroid << ", " << width << endl;
    cout << "Max bin content at centroid: " << max_bin_content << endl;
    cout << "Bin width = " << h->GetBinWidth(max_i) << endl;
    
    return;
}
