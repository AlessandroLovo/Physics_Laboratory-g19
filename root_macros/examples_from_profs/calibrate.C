#include "gethisto.C"

void CalibrateHisto(TH1F *h_uncal, float m, float q) { //Re-scaling of axis, as in the slides

	int max_bin = h_uncal->GetNbinsX(); // This method returns the number of bins in x of the histogram
	float max_kev = h_uncal->GetBinCenter(max_bin)*m + q;
	h_uncal->GetXaxis()->SetLimits(q,max_kev);
	if (m!=1 && q!=0) //This means that I actually changed the calibration!
	    h_uncal->SetXTitle("keV");

}

void calibrate(char *file_na, char *file_bg) {

	// canvas
	TCanvas *c0 = new TCanvas("c0");

	// read the histos
	TH1F *h_na = getHistoFromTree(file_na,514,0,16384);
	TH1F *h_bg = getHistoFromTree(file_bg,514,0,16384);
	h_na->Add(h_bg,-1);

	// search the peaks
	TSpectrum *s = new TSpectrum(30);
	int nPeaks;
	float *na_bin;
	nPeaks = s->Search(h_sodio,1,"",0.1);
 	na_bin = s->GetPositionX();

	// fit the peaks with the Na energies
	float na_kev[2]; 
	na_kev[0] = 511.; 
	na_kev[1] = 1275.;
	TGraphErrors *graphErr = new TGraphErrors(nPeaks,na_bin,na_kev);
	TF1 *fitfun = new TF1("calfitfun1","pol1",na_bin[0],na_bin[1]);
	graphErr->Fit(fitfun);

	// get the fit parameters
	float m = fitfun->GetParameter(1);
	float q = fitfun->GetParameter(0);

	// call the function to calibrate the histogram
	CalibrateHisto(h_na,m,q);

	h_na->Draw();

}

