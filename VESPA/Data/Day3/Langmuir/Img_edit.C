#include <iostream>
#include <vector>
#include <string>
//side: 0=right(filement), 1=left(pump)
TGraphErrors* ParamGraph(int side, int x_par, int y_par){
	vector <int> N;
	vector <double> x;
	vector <double> y;
	vector <double> ex;
	vector <double> ey;
	vector <double> p;//pressure
	switch (x_par){
		case 0:{
			double vn [] = {2,3,4,5,6};
			double vx [] = {20,30,40,50,60};
			double vex [] = {1,1,1,1,1};
			N.assign(vn,vn+5);
			x.assign(vx,vx+5);
			ex.assign(vex,vex+5);
			break;
		}
		case 1:{
			double vn [] = {2,3,4,5,6};
			if(!side){
				double vx []={323,667,721,756,801};
				x.assign(vx,vx+5);
			}
			else{
				double vx []={323,700,717,762,804};
				x.assign(vx,vx+5);
			}
			double vex []= {1,1,1,1,1};
			N.assign(vn,vn+5);
			ex.assign(vex,vex+5);
			break;
		}
		case 2:{
			double vn []= {8,4,7};
			double vx []= {6.1,6.5,6.7};
			double vex []= {0.1,0.1,0.1};
			N.assign(vn,vn+3);
			x.assign(vx,vx+3);
			ex.assign(vex,vex+3);
			break;
		}
		case 3:{
			double vn []= {9,4,10};
			double vx []= {0.282,2.84,29.5};
			double vex []= {0.001,0.01,0.1};
			N.assign(vn,vn+3);
			x.assign(vx,vx+3);
			ex.assign(vex,vex+3);
			break;
		}
	}
	for(int i = 0; i < N.size(); i++){
		ifstream in;
		in.open("/Users/andreagrossutti/Documents/GitHub/Physics_Laboratory-g19/VESPA/Data/Day3/Langmuir/Langmuir_Fit_Results.txt");
		string s,file;
		int n_line = (!side) ? N[i] + 1 : N[i] + 10;
		for(int j = 0; j < n_line; j++)
			getline(in,s);
		stringstream ss;
		ss.str(s);
		ss >> file;
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
		ss >> Is >> Is_err >> R >> R_err >> Te >> Te_err >> Vf >> Vf_err >> cs >> cs_err >> n >> n_err >> Vp >> Vp_err;
		switch(y_par){
			case 0:
				y.push_back(Te);
				ey.push_back(Te_err);
				break;
			case 1:
				y.push_back(n);
				ey.push_back(n_err);
				break;
			case 2:
				y.push_back(Vp);
				ey.push_back(Vp_err);
				break;
			case 3:
				y.push_back(cs);
				ey.push_back(cs_err);
				break;
			case 4:
				y.push_back(Vf);
				ey.push_back(Vf_err);
				break;
			case 5:
				double p_dx [] = {2.7,2.94,2.84,2.84,2.84,2.84,2.84,2.84,0.282,29.5};
				double p_sx [] = {2.7,2.94,2.84,2.84,2.84,2.84,2.84,2.84,0.282,29.5};
				double err_p [] = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.001,0.1};
				double p = (!side) ? p_dx[N[i]-1] : p_sx[N[i]-1];
				double ep = err_p[N[i]-1];
				double n0=p*1.0e-6*1.0e+5/(1.380649e-23)/(273.15+19);
				double en0=n0*sqrt(pow(ep/p,2)+pow(3/(273.15+19),2));
				double ern=pow(n_err/n,2);
				double f=n/(n+n0);
				double ef=f*sqrt(ern+(en0*en0+n_err*n_err)/pow(n+n0,2));
				y.push_back(f);
				ey.push_back(ef);
				break;
		}
	}
	TGraphErrors* g = new TGraphErrors(x.size());
	for(int i=0 ; i<N.size(); i++){
		g->SetPoint(i,x[i],y[i]);
		g->SetPointError(i,ex[i],ey[i]);
	}
	return g;
}


// TGraph* DependenceOnDischargePolVolt(int side, int par){
// 	double v []={20,30,40,50,60};
// 	double T_dx []={0.619492,0.727488,0.738858,0.69059,0.671489};
// 	double T_sx []={0.745883,0.535847,0.505854,0.476674,0.530029};
// 	double n_dx []={1.12E+16,2.70E+16,3.72E+16,4.39E+16,5.63E+16};
// 	double n_sx []={5.22E+15,8.55E+15,8.66E+15,1.05E+16,1.18E+16};
// 	double Vp_dx []={1.3416,1.79968,2.01706,1.88958,2.0838};
// 	double Vp_sx []={1.80092,1.13899,1.03999,1.02233,1.12923};
// 	TGraph* g;
// 	switch (par){
// 		case 0:{
// 			g = !((bool) side) ? new TGraph (5, v, T_dx) : new TGraph (5, v, T_sx);
// 			g->GetYaxis()->SetRangeUser(0.4, 0.9);
// 			break;
// 		};
// 		case 1:{
// 			g = !((bool) side) ? new TGraph (5, v, n_dx) : new TGraph (5, v, n_sx);
// 			g->GetYaxis()->SetRangeUser(0, 5.0E+16);
// 			break;
// 		};
// 		case 2:{
// 			g = !((bool) side) ? new TGraph (5, v, Vp_dx) : new TGraph (5, v, Vp_sx);
// 			g->GetYaxis()->SetRangeUser(0.75, 2.2);
// 			break;
// 		};
// 	}
// 	return g;
// }

void ParamGraph(int x_par, int y_par){
	auto g = ParamGraph(0, x_par, y_par);
	g->SetLineColor(kRed);
	g->SetName("right");

	auto g1 = ParamGraph(1, x_par, y_par);
 	g1->SetLineColor(kBlue);
 	g1->SetName("left");

	switch(x_par){
		case 0:{
			g->GetXaxis()->SetTitle("Discharge Polarization Voltage [V]");
			break;
		}
		case 1:{
			g->GetXaxis()->SetTitle("Discharge Current [mA]");
			break;
		}
		case 2:{
			g->GetXaxis()->SetTitle("Filament Current [A]");
			break;
		}
		case 3:{
			g->GetXaxis()->SetTitle("Pressure [muBar]");
			break;
		}
	}

	switch(y_par){
 		case 0:{
	 		g->GetYaxis()->SetTitle("Te [eV]");
	 		break;
 		}
 		case 1:{
 			g->GetYaxis()->SetTitle("n [m^{-3}]");
 			break;
 		}
 		case 2:{
 			g->GetYaxis()->SetTitle("Vp [V]");
 			break;
 		}
 		case 3:{
 			g->GetYaxis()->SetTitle("c_s [m s^{-1}]");
 			break;
 		}
 		case 4:{
 			g->GetYaxis()->SetTitle("Vf [V]");
 			break;
 		}
 		case 5:{
 			g->GetYaxis()->SetTitle("f []");
 			break;
 		}
 	}

 	double YMax = max(TMath::MaxElement(g->GetN(),g->GetY()),TMath::MaxElement(g1->GetN(),g1->GetY()));
 	double YMin = min(TMath::MinElement(g->GetN(),g->GetY()),TMath::MinElement(g1->GetN(),g1->GetY()));

 	g->GetYaxis()->SetRangeUser(YMin*0.8, YMax*1.2);
 	if(y_par==5){
 		//g->GetYaxis()->SetRangeUser(0,1.0e-3);
 	}
 	g->SetTitle("");

 	
 	

 	g->Draw();
 	g1->Draw("LP");

 	g->SetMarkerStyle(20);
 	g->SetMarkerSize(0.8);
 	g1->SetMarkerStyle(20);
 	g1->SetMarkerSize(0.8);

	if(y_par==5){
		TF1* OneLine = new TF1("OneLine","1.",-100,1000);
		OneLine->Draw("same");
	}

 	auto legend = new TLegend(0.1,0.7,0.3,0.85);
   	legend->AddEntry("right","Filament Side","lep");
   	legend->AddEntry("left","Pump Side","lep");
   	legend->Draw();

}

// void DischargePolVolt(int par){
// 	auto g = DependenceOnDischargePolVolt(0,par);
// 	g->SetLineColor(kViolet);
//  	g->SetName("right");

//  	g->GetXaxis()->SetTitle("Discharge Polarization Voltage [V]");
//  	switch(par){
//  		case 0:{
// 	 		g->GetYaxis()->SetTitle("Te [eV]");
// 	 		break;
//  		}
//  		case 1:{
//  			g->GetYaxis()->SetTitle("n [m^(-3)]");
//  			break;
//  		}
//  		case 2:{
//  			g->GetYaxis()->SetTitle("Vp [V]");
//  			break;
//  		}
//  	}

//  	g->Draw();

//  	g = DependenceOnDischargePolVolt(1,par);
// 	g->SetLineColor(kRed);
//  	g->SetName("left");
//  	g->Draw("SAME");

//  	auto legend = new TLegend(0.1,0.7,0.48,0.9);
//  	legend->AddEntry("right","Filament Side","l");
//  	legend->AddEntry("left","Pump Side","l");
//  	legend->Draw();
// }

// TGraph* DependenceOnDischargePolCurr(int side, int par){
// 	double i_dx []={323,700,721,762,801};
// 	double i_sx []={323,667,717,756,804};
// 	double T_dx []={0.5101,0.6179,0.5901,0.5995,0.5481};
// 	double T_sx []={0.8566,0.5152,0.4471,0.4125,0.5254};
// 	double n_dx []={9.43E+15,2.40E+16,3.23E+16,4.18E+16,4.87E+16};
// 	double n_sx []={5.45E+15,7.80E+15,7.77E+15,9.54E+15,1.19E+16};
// 	double Vp_dx []={0.89,1.48,1.65,1.71,1.76};
// 	double Vp_sx []={2.12,0.93,0.80,0.82,1.12};
// 	// double T_dx []={0.619492,0.727488,0.738858,0.69059,0.671489};
// 	// double T_sx []={0.745883,0.535847,0.505854,0.476674,0.530029};
// 	// double n_dx []={1.12E+16,2.70E+16,3.72E+16,4.39E+16,5.63E+16};
// 	// double n_sx []={5.22E+15,8.55E+15,8.66E+15,1.05E+16,1.18E+16};
// 	// double Vp_dx []={1.3416,1.79968,2.01706,1.88958,2.0838};
// 	// double Vp_sx []={1.80092,1.13899,1.03999,1.02233,1.12923};
// 	TGraph* g;
// 	switch (par){
// 		case 0:{
// 			g = !((bool) side) ? new TGraph (5, i_dx, T_dx) : new TGraph (5, i_sx, T_sx);
// 			g->GetYaxis()->SetRangeUser(0.4, 0.9);
// 			break;
// 		};
// 		case 1:{
// 			g = !((bool) side) ? new TGraph (5, i_dx, n_dx) : new TGraph (5, i_sx, n_sx);
// 			g->GetYaxis()->SetRangeUser(0, 5.0E+16);
// 			break;
// 		};
// 		case 2:{
// 			g = !((bool) side) ? new TGraph (5, i_dx, Vp_dx) : new TGraph (5, i_sx, Vp_sx);
// 			g->GetYaxis()->SetRangeUser(0.75, 2.2);
// 			break;
// 		};
// 	}
// 	return g;
// }

// void DischargePolCurr(int par){
// 	auto g = DependenceOnDischargePolCurr(0,par);
// 	g->SetLineColor(kViolet);
//  	g->SetName("right");

//  	g->GetXaxis()->SetTitle("Discharge Current [mA]");
//  	switch(par){
//  		case 0:{
// 	 		g->GetYaxis()->SetTitle("Te [eV]");
// 	 		break;
//  		}
//  		case 1:{
//  			g->GetYaxis()->SetTitle("n [m^(-3)]");
//  			break;
//  		}
//  		case 2:{
//  			g->GetYaxis()->SetTitle("Vp [V]");
//  			break;
//  		}
//  	}

//  	g->Draw();

//  	g = DependenceOnDischargePolCurr(1,par);
// 	g->SetLineColor(kRed);
//  	g->SetName("left");
//  	g->Draw("SAME");

//  	auto legend = new TLegend(0.1,0.7,0.48,0.9);
//  	legend->AddEntry("right","Filament Side","l");
//  	legend->AddEntry("left","Pump Side","l");
//  	legend->Draw();
// }



// TGraph* DependenceOnFilamentCurr(int side, int par){
// 	double v []={6.1,6.5,6.7};
// 	double T_dx []={0.2834,0.5901,0.8145};
// 	double T_sx []={0.21,0.4471,0.6184};
// 	double n_dx []={1.04E+16,3.23E+16,5.92E+16};
// 	double n_sx []={3.92E+15,7.77E+15,1.40E+16};
// 	double Vp_dx []={0.87,1.65,1.91};
// 	double Vp_sx []={0.48,0.8,1.09};
// 	TGraph* g;
// 	switch (par){
// 		case 0:{
// 			g = !((bool) side) ? new TGraph (5, v, T_dx) : new TGraph (5, v, T_sx);
// 			g->GetYaxis()->SetRangeUser(0.2, 0.9);
// 			break;
// 		};
// 		case 1:{
// 			g = !((bool) side) ? new TGraph (5, v, n_dx) : new TGraph (5, v, n_sx);
// 			g->GetYaxis()->SetRangeUser(0, 8.0E+16);
// 			break;
// 		};
// 		case 2:{
// 			g = !((bool) side) ? new TGraph (5, v, Vp_dx) : new TGraph (5, v, Vp_sx);
// 			g->GetYaxis()->SetRangeUser(0.3, 2.2);
// 			break;
// 		};
// 	}
// 	return g;
// }

// void FilamentCurr(int par){
// 	auto g = DependenceOnDischargePolVolt(0,par);
// 	g->SetLineColor(kViolet);
//  	g->SetName("right");

//  	g->GetXaxis()->SetTitle("Discharge Polarization Voltage [V]");
//  	switch(par){
//  		case 0:{
// 	 		g->GetYaxis()->SetTitle("Te [eV]");
// 	 		break;
//  		}
//  		case 1:{
//  			g->GetYaxis()->SetTitle("n [m^(-3)]");
//  			break;
//  		}
//  		case 2:{
//  			g->GetYaxis()->SetTitle("Vp [V]");
//  			break;
//  		}
//  	}

//  	g->Draw();

//  	g = DependenceOnDischargePolVolt(1,par);
// 	g->SetLineColor(kRed);
//  	g->SetName("left");
//  	g->Draw("SAME");

//  	auto legend = new TLegend(0.1,0.7,0.48,0.9);
//  	legend->AddEntry("right","Filament Side","l");
//  	legend->AddEntry("left","Pump Side","l");
//  	legend->Draw();
// }

// TGraph* DependenceOnPressure(int side, int par){
// 	double v []={0.282,2.84,29.5};
// 	double T_dx []={2.24551,0.5901,0.5189};
// 	double T_sx []={0.907499,0.4471,6.2641};
// 	double n_dx []={4.23E+15,3.23E+16,3.46E+16};
// 	double n_sx []={2.65E+15,7.77E+15,3.41E+14};
// 	double Vp_dx []={4.65265,1.65,1.1};
// 	double Vp_sx []={1.72722,0.8,17.5403};
// 	TGraph* g;
// 	switch (par){
// 		case 0:{
// 			g = !((bool) side) ? new TGraph (3, v, T_dx) : new TGraph (3, v, T_sx);
// 			g->GetYaxis()->SetRangeUser(0.2, 8);
// 			break;
// 		};
// 		case 1:{
// 			g = !((bool) side) ? new TGraph (3, v, n_dx) : new TGraph (3, v, n_sx);
// 			g->GetYaxis()->SetRangeUser(0, 4.0E+16);
// 			break;
// 		};
// 		case 2:{
// 			g = !((bool) side) ? new TGraph (3, v, Vp_dx) : new TGraph (3, v, Vp_sx);
// 			g->GetYaxis()->SetRangeUser(0.3, 20);
// 			break;
// 		};
// 	}
// 	return g;
// }

// void Pressure(int par){
// 	auto g = DependenceOnPressure(0,par);
// 	g->SetLineColor(kViolet);
//  	g->SetName("right");

//  	g->GetXaxis()->SetTitle("Pressure in the chamber [muBar]");
//  	switch(par){
//  		case 0:{
// 	 		g->GetYaxis()->SetTitle("Te [eV]");
// 	 		break;
//  		}
//  		case 1:{
//  			g->GetYaxis()->SetTitle("n [m^(-3)]");
//  			break;
//  		}
//  		case 2:{
//  			g->GetYaxis()->SetTitle("Vp [V]");
//  			break;
//  		}
//  	}

//  	g->Draw();

//  	g = DependenceOnPressure(1,par);
// 	g->SetLineColor(kRed);
//  	g->SetName("left");
//  	g->Draw("SAME");

//  	auto legend = new TLegend(0.1,0.7,0.48,0.9);
//  	legend->AddEntry("right","Filament Side","l");
//  	legend->AddEntry("left","Pump Side","l");
//  	legend->Draw();
// }
