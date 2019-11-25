void plot_CFTD_simulations () {
	ifstream in("CFTD_simulations.txt");
	TGraphErrors* tge25 = new TGraphErrors();
	TGraphErrors* tge50 = new TGraphErrors();
	TGraphErrors* tge75 = new TGraphErrors();
	tge25->SetMarkerColor(1);
	tge50->SetMarkerColor(2);
	tge75->SetMarkerColor(3);
	double a,b,c,d;
	in.ignore(10000,'\n');
	while (in>>a>>b>>c>>d) {
		TGraphErrors* tg;
		if ( a < 0.3 ) tg = tge25;
		else if ( a < 0.6 ) tg = tge50;
		else tg = tge75;
		tg->SetPoint( tg->GetN(), b, c);
		tg->SetPointError ( tg->GetN() - 1, 0, d);
	}
	tge25->Draw();
	if ( tge50->GetN() > 0 ) tge50->Draw("SAME");
	if ( tge75->GetN() > 0 ) tge75->Draw("SAME");
}
