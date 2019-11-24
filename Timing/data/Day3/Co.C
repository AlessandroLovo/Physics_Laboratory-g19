#include <iostream>
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
