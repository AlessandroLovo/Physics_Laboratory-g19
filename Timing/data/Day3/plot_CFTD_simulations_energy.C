void plot_CFTD_simulations_energy (int what, int id = 0) {

	// what: 0->Mean, 1->FWHM, 2->Kurtosis

	//if ( what < 0 || what > 2 ) return;
	char* filename[] = {"CFTD_energythresh_1_2D.txt","CFTD_energythresh_2_2D.txt"};
	ifstream in(filename[id]);
    TGraphErrors* tge_notop = new TGraphErrors();
    TGraphErrors* tge_sitop = new TGraphErrors();

	double a,b,c,d,e,f,g,h,point,error, min = 100000, max = -1000000;
	in.ignore(10000,'\n');
	while (in>>a>>b>>c>>d>>e>>f>>g>>h) {
		if(what == 0) { point = c; error = d; }
		else if(what == 1) { point =  e; error =  f; }
		else if(what == 2) { point = g; error = h; }

		TGraphErrors* tge = ( b > 550 ? tge_notop : tge_sitop );
		
		if ( point > max ) max = point;
		if ( point < min ) min = point;
		tge->SetPoint( tge->GetN(), a, point);
		tge->SetPointError ( tge->GetN() - 1, 0, error);
	}

	if ( what == 0) tge_notop->SetTitle("Mean");
	else if ( what == 1) {
		tge_notop->SetTitle("#sigma [ns]");
	}
	else if ( what == 2) tge_notop->SetTitle("Kurtosis");

	tge_notop->GetXaxis()->SetTitle("LET [keV]");
	tge_notop->GetYaxis()->SetTitle("#sigma [ns]");
	tge_notop->GetYaxis()->SetRangeUser( min - 0.1 * (max-min) , max + 0.1 * (max-min) );

	tge_notop->SetLineColor(2);

	TLegend* tl = new TLegend(2.5,20);
	tl->AddEntry(tge_sitop,"Upper Energy Threshold (UET) = LET + 100 keV","l");
	tl->AddEntry(tge_notop,"Only Lower Energy Threshold (LET)","l");

	tge_notop->Draw();
	tge_sitop->Draw("SAME");
	tl->Draw();
}
