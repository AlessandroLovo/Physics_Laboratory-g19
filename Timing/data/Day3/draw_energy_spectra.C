void draw_energy_spectra(int id = 0) {
	auto wf = new Waveform("Digital_CFTD.root", id);
	auto hist = wf->GetEnergyHisto();

	hist->GetXaxis()->SetTitle("Energy [a.u.]");
	hist->GetYaxis()->SetTitle("Counts");
	hist->SetTitle("A-Posteriori energy spectra");

	hist->Draw();

	// calibration

	double m[] = { 0.08199, 0.06579 };
	double q[] = { -208.5, -106.8 };

	double xmax = hist->GetXaxis()->GetXmax();
	double xmin = hist->GetXaxis()->GetXmin();
	hist->GetXaxis()->SetLimits( xmin * m[id] + q[id], xmax * m[id] + q[id] );

	hist->GetYaxis()->SetRangeUser(0, 11000);

	TLine* tl1 = new TLine( 283.6, 8000, 283.6, 11000 );
	tl1->SetLineColor(2);
	tl1->Draw();


	TLine* tl2 = new TLine( 358.6, 6250, 358.6, 3250 );
	tl2->SetLineColor(2);
	tl2->Draw();

}
