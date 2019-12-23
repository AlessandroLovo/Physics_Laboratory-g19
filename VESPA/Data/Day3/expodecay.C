
void read(TGraphErrors* tge, string title) {
	ifstream in(title);
	double a,b;
	double moltiplic = 1000;
	in.ignore(1000000,'\n');
	while( in>>a>>b ) {
		tge->SetPoint( tge->GetN(), (32.5-a)/100 * moltiplic, b );
		tge->SetPointError( tge->GetN()-1, 0.001 * moltiplic, 0.05*b );
	}
}

void simpleplot() {
	TGraphErrors* tge3 = new TGraphErrors();
	read(tge3,"expodecay_3.txt");
	tge3->Draw();


	tge3->SetTitle("Wave amplitude");
	tge3->GetXaxis()->SetTitle("Position [m]");
	tge3->GetYaxis()->SetTitle("Amplitude [V]");
	
}

void plot() {
//	TGraph* tg1 = new TGraph("expodecay_1.txt");
//	TGraph* tg2 = new TGraph("expodecay_2.txt");
//	TGraph* tg3 = new TGraph("expodecay_3.txt");
	TGraphErrors* tge1 = new TGraphErrors();
	TGraphErrors* tge3 = new TGraphErrors();
	
	read(tge1,"expodecay_1.txt");
	read(tge3,"expodecay_3.txt");
/*
	TLegend* tl = new TLegend(10,5);

	tg1->SetMarkerColor(1);
	tg3->SetMarkerColor(3);
	tl->AddEntry(tg1,"1.0\\cdot10^{-3}mbar","*");
	tl->AddEntry(tg3,"3.3\\cdot10^{-3}mbar","*");
*/
	tge3->Draw();
//	tg2->Draw("SAME");
	tge1->Draw("SAME");
//	tl->Draw();
}

void FitAndPlot() {
	TGraph2DErrors* tge = new TGraph2DErrors();
	TGraphErrors* tge1 = new TGraphErrors();
	TGraphErrors* tge3 = new TGraphErrors();
	double a,b;

	ifstream in1("expodecay_1.txt");
	in1.ignore(1000000,'\n');
	while( in1>>a>>b ) {
		tge->SetPoint( tge->GetN(), a, 0, b );
		tge->SetPointError( tge->GetN()-1, 0.01, 0, 0.05*b );

		tge1->SetPoint( tge1->GetN(), a, b );
		tge1->SetPointError( tge1->GetN()-1, 0.01, 0.05*b );
	}

	ifstream in2("expodecay_3.txt");
	in2.ignore(1000000,'\n');
	while( in2>>a>>b ) {
		tge->SetPoint( tge->GetN(), a, 1, b );
		tge->SetPointError( tge->GetN()-1, 0.01, 0, 0.05*b );

		tge3->SetPoint( tge3->GetN(), a, b );
		tge3->SetPointError( tge3->GetN()-1, 0.01, 0.05*b );
	}

	TF2* func = new TF2("Func","( (1-y)*[c1] + y*[c3] ) * e ^ ( [a] * x )", 25,35, 0,1);

	tge->Draw();
	func->SetParameters(3.2e6,3.2e6,0.5);
	tge->Fit(func,"N0");

/*
	TF1* func1 = new TF1("Func1","[c1] * e^([a]*x)",30,35);
	TF1* func3 = new TF1("Func3","[c3] * e^([a]*x)",30,35);


	func1->SetParameters( func->GetParameter(0), func->GetParameter(2));
	func3->SetParameters( func->GetParameter(1), func->GetParameter(2));

	tge3->Draw();
	tge1->Draw("SAME");
	func1->Draw();
	func3->Draw();*/
}
