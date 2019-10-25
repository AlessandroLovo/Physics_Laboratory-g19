struct slimport_data_t {
	ULong64_t	timetag; //time stamp
	UInt_t		baseline;
	UShort_t	qshort; //integration with shorter time
	UShort_t	qlong; //integration with longer time
	UShort_t	pur;
	UShort_t	samples[4096];
};

class MyBranch : public TBranch {
	public:
		Int_t baskets_number(){
			return fWriteBasket;
		} 
};



class Histograms{

	public: 

	Histograms(const char *name_file){

		original_file = name_file;
		TFile* infile = new TFile(name_file);
		TTree* intree = (TTree*)infile->Get("acq_tree_0");

		for(int chan=0; chan<4; chan++){

			slimport_data_t indata;
			MyBranch *inbranch = (MyBranch*)intree->GetBranch(Form("acq_ch%d",chan));
			inbranch->SetAddress(&indata.timetag);
			int numBins = inbranch->baskets_number();

			hist [chan] = new TH1F(Form("hist_ch%d",chan),"",numBins,0,16384);


			for (int i=0; i<inbranch->GetEntries(); i++) {
				inbranch->GetEntry(i);
				hist [chan]->Fill(indata.qlong);
			}

		}
	}

	TH1F* GetHistogram(int chan){
		if (chan < 0 || chan > 3)
			return 0x0;
		return hist[chan];
	}

	private:

	const char* original_file;

	TH1F* hist [4];





};

