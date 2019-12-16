void plot_I137V( string filename ) {
	TGraphErrors* tg = new TGraphErrors();
	ifstream in(filename);
	double a,b;
	while( in>>a>>b ) {
		a = pow(a, 13.0/7.0);
		tg->SetPoint( tg->GetN(), a, b);
		tg->SetPointError( tg->GetN()-1, 0.3, 0.1);
	}
	tg->SetTitle(filename.c_str());
	tg->GetXaxis()->SetTitle("I^{13/7} [A^{13/7}]");
	tg->GetYaxis()->SetTitle("V [V]");
	tg->Draw();
}

