struct slimport_data_t {
	ULong64_t	timetag; //time stamp
	UInt_t		baseline;
	UShort_t	qshort; //integration with shorter time
	UShort_t	qlong; //integration with longer time
	UShort_t	pur;
	UShort_t	samples[4096];
};


class Histograms{

	public: 

	//Costructor
	Histograms(const char *name_file){

		original_data_file = name_file;
		TFile* infile = new TFile(name_file);
		TTree* intree = (TTree*)infile->Get("acq_tree_0");

		for(int chan=0; chan<4; chan++){

			slimport_data_t indata;
			TBranch *inbranch = (TBranch*)intree->GetBranch(Form("acq_ch%d",chan));
			inbranch->SetAddress(&indata.timetag);
			hist [chan] = new TH1F(Form("hist_ch%d",chan),"",4096,0,32768);

			for (int i=0; i<inbranch->GetEntries(); i++) {
				inbranch->GetEntry(i);
				hist [chan]->Fill(indata.qlong);
			}
		}
	}

	//Get Histogram of channel with index "chan"
	TH1F* GetHistogram(int chan){
		if (chan < 0 || chan > 3)
			return 0x0;
		return hist[chan];
	}

	//Get path to the file with datas of the histograms
	const char* GetDataFile(){
		return original_data_file;
	}

	//Calibrate histogram through a relation E = alpha * channel + beta
	void CalibrateHisto(int chan, double alpha, double beta) {
		if (chan < 0 || chan > 3)
			return;
		TAxis *axis = hist[chan]->GetXaxis();
		axis->SetLimits(axis->GetXmin()*alpha+beta, axis->GetXmax()*alpha+beta);
	}

	private:

	const char* original_data_file;
	TH1F* hist [4];
};
