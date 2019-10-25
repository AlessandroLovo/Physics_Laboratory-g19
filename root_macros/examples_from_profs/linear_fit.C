void fit(float peak0, float peak1) {

	float m=1., q=0.; // default values for calibration (m=1, q=0 means "no calibration")

	float na_bin[2];
	na_bin[0] = peak0; // Sodium photopeak energy in bins
	na_bin[1] = peak1; // Sodium photopeak energy in bins

	float na_kev[2];
	na_kev[0] = 511.; // Sodium photopeak energy in keV
	na_kev[1] = 1275.; // Sodium photopeak energy in keV

	// Define a TGraphErrors which contains the experimental points in x
	// and the theoretical points in y; I semplify by not taking into account errors.
	TGraphErrors *graphErr = new TGraphErrors(2,na_bin,na_kev);

	// Define a function which fits the points
	TF1 *fitfun = new TF1("calfitfun1","pol1",na_bin[0],na_bin[1]);

	// Graphical options to draw the fit
	graphErr->SetMarkerStyle(21);
	graphErr->SetLineColor(2);
	graphErr->Draw("ALP");

	// Fit
	graphErr->Fit(fitfun);

	// Get the parameters
	m = fitfun->GetParameter(1);
	q = fitfun->GetParameter(0);

	cout << "Slope = " << m << ", Intercept = " << q << endl;

}

