void plot_CFTD_simulations (int id = 0, bool same = false) {
	char* filename[] = {"CFTD_simulations.txt","CFTD_simulations_2.txt"};
	ifstream in(filename[id]);
	TGraphErrors* tge25 = new TGraphErrors();
	TGraphErrors* tge50 = new TGraphErrors();
	TGraphErrors* tge75 = new TGraphErrors();
	tge25->SetLineColor(1 + same*3);
	tge50->SetLineColor(2 + same*3);
	tge75->SetLineColor(3 + same*3);
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
	tge25->Draw( ( same ? "SAME P" : "AP") );
	if ( tge50->GetN() > 0 ) tge50->Draw("SAME P");
	if ( tge75->GetN() > 0 ) tge75->Draw("SAME P");

	tge25->GetXaxis()->SetTitle("Delay [ns]");
	tge25->GetYaxis()->SetTitle("FWHM [ps]");
	tge25->SetTitle(Form("CFTD simulation, dataset %d",id+1));
	TLegend* tl = new TLegend(2,5);
	tl->AddEntry(tge25,"Fraction: 25%","l");
	tl->AddEntry(tge50,"Fraction: 50%","l");
	tl->AddEntry(tge75,"Fraction: 75%","l");
	tl->Draw();
}
