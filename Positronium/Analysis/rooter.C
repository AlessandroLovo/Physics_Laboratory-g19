#include <fstream>
#include <iostream>
#include <math.h>

using namespace std;

TH1F* energy_spectrum(const char *filename, int chan, int numbins){
    TH1F* h = new TH1F(filename,"simulated energy spectrum",numbins,0,511);
    
    if(chan < 0 || chan > 2){
        cout << "Invalid channel" << endl;
        return 0;
    }
    
    ifstream in(filename);
    if(!in){
        cout << "Cannot open file " << filename << endl;
        return 0;
    }
    
    int i = 0;
    float e = 0.;
    while(in >> e){
        if(i % 3 == chan)
            h->Fill(e);
        i++;
    }
    
    cout << i << " data" << endl;
    
    return h;
}

TH2F* angular_correlation(const char *filename, int numbins){
    TH2F* h2 = new TH2F(filename,"simulated angular correlation",numbins,0,M_PI,numbins,0,M_PI);
    
    ifstream in(filename);
    if(!in){
        cout << "Cannot open file " << filename << endl;
        return 0;
    }
    
    int i = 0;
    float t1 = 0., t2 = 0., t3 = 0.;
    while(in >> t1){
        in >> t2;
        in >> t3;
        h2->Fill(t1,t2);
        
        i++;
    }
    
    cout << i << " data" << endl;
    
    return h2;
}

TH2F* energy_correlation(const char *filename, int numbins){
    TH2F* h2 = new TH2F(filename,"simulated energy correlation",numbins,0,511,numbins,0,511);
    
    ifstream in(filename);
    if(!in){
        cout << "Cannot open file " << filename << endl;
        return 0;
    }
    
    int i = 0;
    float t1 = 0., t2 = 0., t3 = 0.;
    while(in >> t1){
        in >> t2;
        in >> t3;
        h2->Fill(t1,t2);
        
        i++;
    }
    
    cout << i << " data" << endl;
    
    return h2;
}
