//These methods are useful to open canvas with the right setting of the fit pannel

//(Maybe) TO DO: add constraint in parameters, for example p[0]>0 for exponential amplitude


void DetectorCalibration(TH1F* h){

	TF1* gaussCostantNoise = new TF1("gaussCostantNoise","gaus(0)+[3]*x+[4]",0, 32768);
	h->Fit("gaussCostantNoise", "R");

}

void TacCalibration(TH1F* h){

	TF1* gausNoNoise = new TF1("gausNoNoise","gaus(0)",0, 32768);
	h->Fit("gausNoNoise","R");

}

void DecayRate(TH1F* h){

	TF1* expDecay = new TF1("expDecay","[0]*exp([1]*(x-[2]))+[3]",0, 32768);
	h->Fit("expDecay", "R");

}

void PeaksToCoefficents(int ch1275, int ch511){

	double alpha=(1275-511)/(ch1275-ch511);
	double beta=1275-alpha*ch1275;
	cout << "Channel corresponding to 1275 keV: " << ch1275 << endl;
	cout << "Channel corresponding to 511 keV: " << ch511 << endl;
	cout << "Energy(keV) = alpha * channel + beta" << endl;
	cout << "alpha = " << alpha << endl;
	cout << "beta = " << beta << endl;
}

int EnergyToChannel