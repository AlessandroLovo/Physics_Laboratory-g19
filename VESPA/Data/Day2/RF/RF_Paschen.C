void plot() {
	TGraphErrors* tg = new TGraphErrors();
	ifstream in("RF_Paschen.txt");
	in.ignore(10000,'\n');
	double a,b,i;
	while ( in>>a>>b>>i ) {
		tg->SetPoint( tg->GetN(), a,b );
		tg->SetPointError( tg->GetN()-1, 0.05*a, 0.01*b );
	}
	tg->Draw();
	tg->SetTitle("Breakdown voltage - RF discharge");
	tg->GetXaxis()->SetTitle("Pressure [#mubar]");
	tg->GetYaxis()->SetTitle("Peak-to-peak [V]");
	gPad->SetLogx();
}
