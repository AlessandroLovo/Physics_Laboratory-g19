void plot() {
	TGraphErrors* tg = new TGraphErrors();
	ifstream in("RF_Resonance.txt");
	in.ignore(10000,'\n');
	double a,b;
	while ( in>>a>>b ) {
		tg->SetPoint( tg->GetN(), a,b );
		tg->SetPointError( tg->GetN()-1, 0.005*a, 0.01*b );
	}
	tg->Draw();
}
