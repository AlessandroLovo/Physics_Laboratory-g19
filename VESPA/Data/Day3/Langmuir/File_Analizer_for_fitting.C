#include <iostream>
#include <vector>

void an(char* file){
	ifstream in;
	in.open(file);
	in.ignore(10000,'\n');
	TGraph* g=new TGraph();
	double i,j,k,l;
	while(in >> i >> j >> k >> l)
		g->SetPoint(g->GetN(),k,l);
	TF1* LangmuirFourParam = new TF1("LangmuirFourParam","[Is]*(1.+[R]*(x-[Vf]))*(exp((x-[Vf])/[Te])-1.)",-30,30);
	LangmuirFourParam -> SetParameter("Is", 0.028);
	LangmuirFourParam -> SetParameter("R", -0.14);
	LangmuirFourParam -> SetParameter("Vf", -2.);
	LangmuirFourParam -> SetParameter("Te", 0.68);
	g->Draw();
	g->Fit("LangmuirFourParam");
	
}