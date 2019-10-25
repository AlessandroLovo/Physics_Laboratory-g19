#include "gethisto.C"

void gaussian_fit(char *file_na) {

	// Canvas
	TCanvas *c[3];
	for (int i=0; i<3; i++)
	   c[i] = new TCanvas(Form("c%i",i));

	// Get a histo with gaussian peaks
	TH1F *h_na = getHistoFromTree(file_na,1026,0,16384);
	c[0]->cd();
	h_na->Draw();
	
	// Get the two peaks of Na
	TSpectrum *s = new TSpectrum(30);
	int nPeaks;
	float *xPeaks;
	nPeaks = s->Search(h_na,1,"goff",0.12);
	xPeaks = s->GetPositionX();

	// I fit the two peaks with a Gaussian to get better values of centroids and sigmas

	int npoints;
	TH1F *h_fitted[2];
	float left_mk, right_mk, x_fit_min, x_fit_max, chisqr, mean, sigma;
	apertura_sx = 250; // arbitrary values. Per un calcolo raffinato si può fare un secondo
	apertura_dx = 450; // fit dove questi valori corrispondono alle sigma del primo fit...
	for (int p=0; p<2; p++) {
		c[p+1]->cd(); // change canvas
		h_fitted[p] = (TH1F*)h_na->Clone(); // Clone on a new TH1F to show separately the fits
		x_fit_min = xPeaks[p] - left_mk; // fit left margin
		x_fit_max = xPeaks[p] + right_mk; // fit right margin
		// fit function
		TF1 *fit = new TF1("fit","gaus",x_fit_min,x_fit_max);
		fit->SetParameters(1,xPeaks[p]); // Fix the centroid from the already known x
		fit->SetNpx(1000); // Number of points to define the function
		// fit
		h_fitted[p]->Fit("fit","","",x_fit_min,x_fit_max);
		npoints = (fit->GetNumberFitPoints())-3;
		if (npoints>0) 
		    chisqr = fit->GetChisquare()/npoints;
		else
		    chisqr = 9999;
	//	const = fit->GetParameter(0);
		mean = fit->GetParameter(1);
		sigma = fit->GetParameter(2);
		// Results
		cout << "Fit the peak n." << p << ": with TSpectrum x = " << xPeaks[p] << endl;
		cout << "                    with the fit    x = " << mean << " and sigma = " << sigma << endl;
		fit->Delete();
	}

	// http://root.cern.ch/viewvc/trunk/tutorials/fit/multifit.C?view=markup

}

