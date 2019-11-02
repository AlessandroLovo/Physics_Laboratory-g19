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