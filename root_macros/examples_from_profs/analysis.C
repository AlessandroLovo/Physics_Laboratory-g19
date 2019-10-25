struct slimport_data_t {
	UInt_t	timestamp;
	UInt_t	baseline;
	UShort_t	qshort;
	UShort_t	qlong;
	UShort_t	pur;
	UShort_t	scope[144];
};


void analysis(char *name_file) {

	// variables
	slimport_data_t indata;
	TFile *infile = new TFile(name_file);
	TTree *intree = (TTree*)infile->Get("datatree");
	TBranch *inbranch = intree->GetBranch("ACQ_ch0");
	inbranch->SetAddress(&indata.timestamp);
	histo2D = new TH2F("hs_2D","PSD plot",256,0,12288,64,0,0.4);
	float ratio;

	// histogram filling
	for (int i=0; i<inbranch->GetEntries(); i++) {
		inbranch->GetEntry(i);
		ratio = (float)(indata.qlong-indata.qshort)/(indata.qlong);
		if (ratio>0.15 && indata.qlong>2000) {
			histo2D->Fill( indata.qlong, ratio );
		}
	}

	TCanvas *c0 = new TCanvas("c0");
	histo2D->Draw("colz");

}
