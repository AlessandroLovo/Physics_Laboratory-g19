#include <iostream>
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

    for ( int i=100; i <= 1000; i++ ) {
	e_min->at(0) = i;
	e_min->at(1) = i;
	C->Filter(channel,e_min)->GetHistogram(2)->Draw();
	cout<<i<<endl;
    }

    for ( int i=100; i <= 1000; i++ ) { 
        e_min->at(0) = i;
        e_min->at(1) = i;
	e_max->at(0) = i+100;
	e_max->at(1) = i+100;
        C->Filter(channel,e_min,e_max)->GetHistogram(2)->Draw();
    }


    return C;
}

