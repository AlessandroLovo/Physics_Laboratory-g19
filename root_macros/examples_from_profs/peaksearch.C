#include "gethisto.C"

void peaksearch(char *file_na, char *file_bg) {

	// creare il canvas
	TCanvas *c0 = new TCanvas("c0");
	c0->Divide(2,2); // to draw more histos on a Canvas (see slides)

	TH1F *h_na = getHistoFromTree(file_na,1026,0,16384);
	TH1F *h_bg = getHistoFromTree(file_bg,1026,0,16384);
	TH1F *h_subtr = h_na->Clone();
	h_subtr->Add(h_bg,-1);
	TH1F *h_peaks = h_subtr->Clone();

	h_na->SetTitle("Measured spectrum");
	h_bg->SetTitle("Measure of background");
	h_subtr->SetTitle("Subtracted background");
	h_peaks->SetTitle("Peak search");

	c0->cd(1);
	h_na->Draw();
	c0->cd(2);
	h_bg->Draw();
	c0->cd(3);
	h_subtr->Draw();
	c0->cd(4);
	h_peaks->Draw();

	TSpectrum *s = new TSpectrum(30);
	int nPeaks;
	float *xPeaks;
	int sigma = 1;
	double minratio = 0.12;	// minimum ratio between a peak and the main peak

	nPeaks = s->Search(h_peaks,sigma,"",minratio);

	xPeaks = s->GetPositionX();
	for (int p = 0; p<nPeaks; p++) {
		cout << "Peak #" << p << " @ channel " << xPeaks[p] << endl;	
	}
	
	//Print the result as png image
	c0->Print(Form("peaks_sigma%i.png",sigma));

}

