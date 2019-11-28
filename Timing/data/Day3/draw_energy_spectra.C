void draw_energy_spectra() {
	auto wf = new Waveform("Digital_CFTD.root", 0);
	auto hist = wf->GetEnergyHisto();

	hist->GetXaxis()->SetTitle("Energy [keV]");
	hist->GetYaxis()->SetTitle("Counts");
	hist->SetTitle("A-Posteriori energy spectra");

	hist->Draw();

	// calibration

	double m = 0.093;
	double q = -282;

	double xmax = hist->GetXaxis()->GetXmax();
	double xmin = hist->GetXaxis()->GetXmin();
	hist->GetXaxis()->SetLimits( xmin * m + q, xmax * m + q );

	hist->GetYaxis()->SetRangeUser(0, 11000);

	TLine* tl1 = new TLine( 283.6, 8000, 283.6, 11000 );
	tl1->SetLineColor(2);
	tl1->Draw();


	TLine* tl2 = new TLine( 358.6, 6250, 358.6, 3250 );
	tl2->SetLineColor(2);
	tl2->Draw();
}
