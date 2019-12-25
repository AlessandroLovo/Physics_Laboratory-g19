#include <iostream>

TGraph* DependenceOnDischargePolVolt(int side, int par){
	double v []={20,30,40,50,60};
	// double T_dx []={0.5101,0.6179,0.5901,0.5995,0.5481};
	// double T_sx []={0.8566,0.5152,0.4471,0.4125,0.5254};
	// double n_dx []={9.43E+15,2.40E+16,3.23E+16,4.18E+16,4.87E+16};
	// double n_sx []={5.45E+15,7.80E+15,7.77E+15,9.54E+15,1.19E+16};
	// double Vp_dx []={0.89,1.48,1.65,1.71,1.76};
	// double Vp_sx []={2.12,0.93,0.80,0.82,1.12};
	double T_dx []={0.619492,0.727488,0.738858,0.69059,0.671489};
	double T_sx []={0.745883,0.535847,0.505854,0.476674,0.530029};
	double n_dx []={1.12E+16,2.70E+16,3.72E+16,4.39E+16,5.63E+16};
	double n_sx []={5.22E+15,8.55E+15,8.66E+15,1.05E+16,1.18E+16};
	double Vp_dx []={1.3416,1.79968,2.01706,1.88958,2.0838};
	double Vp_sx []={1.80092,1.13899,1.03999,1.02233,1.12923};
	TGraph* g;
	switch (par){
		case 0:{
			g = !((bool) side) ? new TGraph (5, v, T_dx) : new TGraph (5, v, T_sx);
			g->GetYaxis()->SetRangeUser(0.4, 0.9);
			break;
		};
		case 1:{
			g = !((bool) side) ? new TGraph (5, v, n_dx) : new TGraph (5, v, n_sx);
			g->GetYaxis()->SetRangeUser(0, 5.0E+16);
			break;
		};
		case 2:{
			g = !((bool) side) ? new TGraph (5, v, Vp_dx) : new TGraph (5, v, Vp_sx);
			g->GetYaxis()->SetRangeUser(0.75, 2.2);
			break;
		};
	}
	return g;
}

void DischargePolVolt(int par){
	auto g = DependenceOnDischargePolVolt(0,par);
	g->SetLineColor(kViolet);
 	g->SetName("right");

 	g->GetXaxis()->SetTitle("Discharge Polarization Voltage [V]");
 	switch(par){
 		case 0:{
	 		g->GetYaxis()->SetTitle("Te [eV]");
	 		break;
 		}
 		case 1:{
 			g->GetYaxis()->SetTitle("n [m^(-3)]");
 			break;
 		}
 		case 2:{
 			g->GetYaxis()->SetTitle("Vp [V]");
 			break;
 		}
 	}

 	g->Draw();

 	g = DependenceOnDischargePolVolt(1,par);
	g->SetLineColor(kRed);
 	g->SetName("left");
 	g->Draw("SAME");

 	auto legend = new TLegend(0.1,0.7,0.48,0.9);
 	legend->AddEntry("right","Filament Side","l");
 	legend->AddEntry("left","Pump Side","l");
 	legend->Draw();
}

TGraph* DependenceOnDischargePolCurr(int side, int par){
	double i_dx []={323,700,721,762,801};
	double i_sx []={323,667,717,756,804};
	double T_dx []={0.5101,0.6179,0.5901,0.5995,0.5481};
	double T_sx []={0.8566,0.5152,0.4471,0.4125,0.5254};
	double n_dx []={9.43E+15,2.40E+16,3.23E+16,4.18E+16,4.87E+16};
	double n_sx []={5.45E+15,7.80E+15,7.77E+15,9.54E+15,1.19E+16};
	double Vp_dx []={0.89,1.48,1.65,1.71,1.76};
	double Vp_sx []={2.12,0.93,0.80,0.82,1.12};
	// double T_dx []={0.619492,0.727488,0.738858,0.69059,0.671489};
	// double T_sx []={0.745883,0.535847,0.505854,0.476674,0.530029};
	// double n_dx []={1.12E+16,2.70E+16,3.72E+16,4.39E+16,5.63E+16};
	// double n_sx []={5.22E+15,8.55E+15,8.66E+15,1.05E+16,1.18E+16};
	// double Vp_dx []={1.3416,1.79968,2.01706,1.88958,2.0838};
	// double Vp_sx []={1.80092,1.13899,1.03999,1.02233,1.12923};
	TGraph* g;
	switch (par){
		case 0:{
			g = !((bool) side) ? new TGraph (5, i_dx, T_dx) : new TGraph (5, i_sx, T_sx);
			g->GetYaxis()->SetRangeUser(0.4, 0.9);
			break;
		};
		case 1:{
			g = !((bool) side) ? new TGraph (5, i_dx, n_dx) : new TGraph (5, i_sx, n_sx);
			g->GetYaxis()->SetRangeUser(0, 5.0E+16);
			break;
		};
		case 2:{
			g = !((bool) side) ? new TGraph (5, i_dx, Vp_dx) : new TGraph (5, i_sx, Vp_sx);
			g->GetYaxis()->SetRangeUser(0.75, 2.2);
			break;
		};
	}
	return g;
}

void DischargePolCurr(int par){
	auto g = DependenceOnDischargePolCurr(0,par);
	g->SetLineColor(kViolet);
 	g->SetName("right");

 	g->GetXaxis()->SetTitle("Discharge Current [mA]");
 	switch(par){
 		case 0:{
	 		g->GetYaxis()->SetTitle("Te [eV]");
	 		break;
 		}
 		case 1:{
 			g->GetYaxis()->SetTitle("n [m^(-3)]");
 			break;
 		}
 		case 2:{
 			g->GetYaxis()->SetTitle("Vp [V]");
 			break;
 		}
 	}

 	g->Draw();

 	g = DependenceOnDischargePolCurr(1,par);
	g->SetLineColor(kRed);
 	g->SetName("left");
 	g->Draw("SAME");

 	auto legend = new TLegend(0.1,0.7,0.48,0.9);
 	legend->AddEntry("right","Filament Side","l");
 	legend->AddEntry("left","Pump Side","l");
 	legend->Draw();
}



TGraph* DependenceOnFilamentCurr(int side, int par){
	double v []={6.1,6.5,6.7};
	double T_dx []={0.2834,0.5901,0.8145};
	double T_sx []={0.21,0.4471,0.6184};
	double n_dx []={1.04E+16,3.23E+16,5.92E+16};
	double n_sx []={3.92E+15,7.77E+15,1.40E+16};
	double Vp_dx []={0.87,1.65,1.91};
	double Vp_sx []={0.48,0.8,1.09};
	TGraph* g;
	switch (par){
		case 0:{
			g = !((bool) side) ? new TGraph (5, v, T_dx) : new TGraph (5, v, T_sx);
			g->GetYaxis()->SetRangeUser(0.2, 0.9);
			break;
		};
		case 1:{
			g = !((bool) side) ? new TGraph (5, v, n_dx) : new TGraph (5, v, n_sx);
			g->GetYaxis()->SetRangeUser(0, 8.0E+16);
			break;
		};
		case 2:{
			g = !((bool) side) ? new TGraph (5, v, Vp_dx) : new TGraph (5, v, Vp_sx);
			g->GetYaxis()->SetRangeUser(0.3, 2.2);
			break;
		};
	}
	return g;
}

void FilamentCurr(int par){
	auto g = DependenceOnDischargePolVolt(0,par);
	g->SetLineColor(kViolet);
 	g->SetName("right");

 	g->GetXaxis()->SetTitle("Discharge Polarization Voltage [V]");
 	switch(par){
 		case 0:{
	 		g->GetYaxis()->SetTitle("Te [eV]");
	 		break;
 		}
 		case 1:{
 			g->GetYaxis()->SetTitle("n [m^(-3)]");
 			break;
 		}
 		case 2:{
 			g->GetYaxis()->SetTitle("Vp [V]");
 			break;
 		}
 	}

 	g->Draw();

 	g = DependenceOnDischargePolVolt(1,par);
	g->SetLineColor(kRed);
 	g->SetName("left");
 	g->Draw("SAME");

 	auto legend = new TLegend(0.1,0.7,0.48,0.9);
 	legend->AddEntry("right","Filament Side","l");
 	legend->AddEntry("left","Pump Side","l");
 	legend->Draw();
}

TGraph* DependenceOnPressure(int side, int par){
	double v []={0.282,2.84,29.5};
	double T_dx []={2.24551,0.5901,0.5189};
	double T_sx []={0.907499,0.4471,6.2641};
	double n_dx []={4.23E+15,3.23E+16,3.46E+16};
	double n_sx []={2.65E+15,7.77E+15,3.41E+14};
	double Vp_dx []={4.65265,1.65,1.1};
	double Vp_sx []={1.72722,0.8,17.5403};
	TGraph* g;
	switch (par){
		case 0:{
			g = !((bool) side) ? new TGraph (3, v, T_dx) : new TGraph (3, v, T_sx);
			g->GetYaxis()->SetRangeUser(0.2, 8);
			break;
		};
		case 1:{
			g = !((bool) side) ? new TGraph (3, v, n_dx) : new TGraph (3, v, n_sx);
			g->GetYaxis()->SetRangeUser(0, 4.0E+16);
			break;
		};
		case 2:{
			g = !((bool) side) ? new TGraph (3, v, Vp_dx) : new TGraph (3, v, Vp_sx);
			g->GetYaxis()->SetRangeUser(0.3, 20);
			break;
		};
	}
	return g;
}

void Pressure(int par){
	auto g = DependenceOnPressure(0,par);
	g->SetLineColor(kViolet);
 	g->SetName("right");

 	g->GetXaxis()->SetTitle("Pressure in the chamber [muBar]");
 	switch(par){
 		case 0:{
	 		g->GetYaxis()->SetTitle("Te [eV]");
	 		break;
 		}
 		case 1:{
 			g->GetYaxis()->SetTitle("n [m^(-3)]");
 			break;
 		}
 		case 2:{
 			g->GetYaxis()->SetTitle("Vp [V]");
 			break;
 		}
 	}

 	g->Draw();

 	g = DependenceOnPressure(1,par);
	g->SetLineColor(kRed);
 	g->SetName("left");
 	g->Draw("SAME");

 	auto legend = new TLegend(0.1,0.7,0.48,0.9);
 	legend->AddEntry("right","Filament Side","l");
 	legend->AddEntry("left","Pump Side","l");
 	legend->Draw();
}
