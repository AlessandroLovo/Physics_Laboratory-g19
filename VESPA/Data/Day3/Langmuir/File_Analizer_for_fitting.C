#include <iostream>
#include <vector>
#include <math.h>

bool fit_in_costructor = true;

class LangFit{
	public:
	LangFit(char* f){
		file = f;
		ifstream in;
		in.open(Form("/Users/andreagrossutti/Documents/GitHub/Physics_Laboratory-g19/VESPA/Data/Day3/Langmuir/DatiAcquisiti/%s.txt",file));
		in.ignore(10000,'\n');
		g=new TGraph();
		double i,j,k,l;
		while(in >> i >> j >> k >> l)
			g->SetPoint(g->GetN(),k,l);
		g->Draw();
		Is_fit();
		Fit();
	}

	void Is_fit(){
		TF1* Linear = new TF1("Linear","[Is]*x+[q]",-19,-2);
		Linear -> SetParameter("Is", 0.070);
		Linear -> SetParameter("q", 0);
		g->Fit("Linear","RS");
	}
	
	void DrawFitFunc(){
		TF1* f = g->GetFunction("LangmuirFourParam");
		f->SetRange(f->GetXmin(),f->GetXmax()+2);
		f->Draw("SAME");
	}

	void DrawGraph(){
		g->GetXaxis()->SetTitle("Potential [V]");
		g->GetYaxis()->SetTitle("Current [mA]");
		g->Draw();
	}

	void Fit(){
		TF1* f = g->GetFunction("Linear");
		Is = f->GetParameter(0);
		Is_err = f->GetParError(0);
		TF1* LangmuirFourParam = new TF1("LangmuirFourParam","[Is]*(1.+[R]*(x-[Vf]))*(exp((x-[Vf])/[Te])-1.)",-19,2);
		LangmuirFourParam -> SetParameter(0, Is);
		LangmuirFourParam -> SetParameter("R", -0.14);
		LangmuirFourParam -> SetParameter("Vf", -2.);
		LangmuirFourParam -> SetParameter("Te", 0.68);
		g->Fit("LangmuirFourParam","RS");
	}

	void DerivatedData(){
		TF1* f = g->GetFunction("LangmuirFourParam");
		Is = f->GetParameter(0);
		Is_err = f->GetParError(0);
		double Is_rel = pow(Is_err/Is,2);
		R = f->GetParameter(1);
		R_err = f->GetParError(1);
		double R_rel = pow(R_err/R,2);
		Te = f->GetParameter(2);
		Te_err = f->GetParError(2);
		double Te_rel = pow(Te_err/Te,2);
		Vf = f->GetParameter(3);
		Vf_err = f->GetParError(3);
		double Vf_rel = pow(Vf_err/Vf,2);
		double Rshunt = 100.;
		double area = 30.; //sq. mm
		double area_rel = pow(1./area,2);
		double e = 1.6022e-19;//C
		double e_rel = pow(.0001e-19/e,2);
		double mp = 1.66054e-27;//KG
		double mp_rel = pow(.00001e-27/mp,2);
		double me = 9.10938356e-31;//KG
		double me_rel = pow(.00000001e-31/me,2);
		double argon_am = 39.948;//uma
		double argon_am_rel=pow(.001/argon_am,2);
		double mi = mp * argon_am;
		double mi_err = mi*sqrt(mp_rel + argon_am_rel);
		double mi_rel = pow(mi_err/mi,2);
		double pi = 3.14159265358979323846;
		double alpha = 0.5*(log(mi/(2*pi*me))+1.);
		double alpha_err = 0.5*sqrt(mi_rel+me_rel)/(mi/(2*pi*me));
		double alpha_rel = pow(alpha_err/alpha,2);
		cs = sqrt(e*Te/mi);
		cs_err = cs/2*sqrt(e_rel+Te_rel+mi_rel); 
		double cs_rel = pow(cs_err/cs,2);
		n = 2.0*Is/(e*cs*area*1e-3);
		n_err = n*sqrt(Is_rel+e_rel+cs_rel+area_rel);
		//n = 2.*Is/(e*cs*area*1e-6)/1000;
		Vp = Vf + alpha*Te;
		Vp_err = sqrt(Vf_rel+pow(alpha*Te,2)*(alpha_rel+Te_rel));
	}

	void DrawResults(){
		DerivatedData();
		ofstream out;
		out.open("/Users/andreagrossutti/Documents/GitHub/Physics_Laboratory-g19/VESPA/Data/Day3/Langmuir/Langmuir_Fit_Results.txt", ofstream::out | ofstream::app);
		out << file << '\t';
		out << Is << '\t' << Is_err << '\t';
		out << R << '\t' << R_err << '\t';
		out << Te << '\t' << Te_err << '\t';
		out << Vf << '\t' << Vf_err << '\t';
		out << cs << '\t'<< cs_err << '\t';
		out << n << '\t'<< n_err << '\t';
		out << Vp <<'\t'<< Vp_err << endl;
		cout << n << '\t'<< Vp << endl;
	}

	TGraph* GetGraph(){
		return g;
	}

	private:
	char* file;
	TGraph* g;
	double Is;
	double Is_err;
	double Te;
	double Te_err;
	double Vf;
	double Vf_err;
	double R;
	double R_err;
	double cs;
	double cs_err;
	double n;
	double n_err;
	double Vp;
	double Vp_err;
};

void MultGraph_sameCurr(){
	fit_in_costructor = false;
	auto l = new LangFit("dx6");
	auto g = l->GetGraph();
	g->SetLineColor(kBlack);
	g->SetName("sx6");
	l->DrawGraph();

	l = new LangFit("dx2");
	g = l->GetGraph();
	g->SetLineColor(kGreen);
	g->SetName("sx2");
	g->Draw("SAME");

	l = new LangFit("dx3");
	g = l->GetGraph();
	g->SetLineColor(kRed);
	g->SetName("sx3");
	g->Draw("SAME");

	l = new LangFit("dx4");
	g = l->GetGraph();
	g->SetLineColor(kBlue);
	g->SetName("sx4");
	g->Draw("SAME");

	l = new LangFit("dx5");
	g = l->GetGraph();
	g->SetLineColor(kViolet);
	g->SetName("sx5");
	g->Draw("SAME");

	auto legend = new TLegend(0.1,0.7,0.48,0.9);
   legend->AddEntry("sx2","20V","l");
   legend->AddEntry("sx3","30V","l");
   legend->AddEntry("sx4","40V","l");
   legend->AddEntry("sx5","50V","l");
   legend->AddEntry("sx6","60V","l");
   legend->Draw();
}

void MultGraph_sameVolt(){
	fit_in_costructor = false;
	auto l = new LangFit("sx7");
	auto g = l->GetGraph();
	g->SetLineColor(kBlack);
	g->SetName("dx7");
	l->DrawGraph();

	l = new LangFit("sx4");
	g = l->GetGraph();
	g->SetLineColor(kBlue);
	g->SetName("dx4");
	g->Draw("SAME");

	l = new LangFit("sx8");
	g = l->GetGraph();
	g->SetLineColor(kRed);
	g->SetName("dx8");
	g->Draw("SAME");


	auto legend = new TLegend(0.1,0.7,0.48,0.9);
   legend->AddEntry("dx8","6.1A","l");
   legend->AddEntry("dx4","6.5A","l");
   legend->AddEntry("dx7","6.7A","l");
   legend->Draw();

}

void MultGraph_changePressure(){
	fit_in_costructor = false;
	auto l = new LangFit("sx10");
	auto g = l->GetGraph();
	g->SetLineColor(kRed);
	g->SetName("dx10");
	l->DrawGraph();

	l = new LangFit("sx4");
	g = l->GetGraph();
	g->SetLineColor(kBlue);
	g->SetName("dx4");
	g->Draw("SAME");

	l = new LangFit("sx9");
	g = l->GetGraph();
	g->SetLineColor(kBlack);
	g->SetName("dx9");
	g->Draw("SAME");


	auto legend = new TLegend(0.1,0.7,0.48,0.9);
   legend->AddEntry("dx9","0.282 muBar","l");
   legend->AddEntry("dx4","2.82 muBar","l");
   legend->AddEntry("dx10","28.2 muBar","l");
   legend->Draw();

}





