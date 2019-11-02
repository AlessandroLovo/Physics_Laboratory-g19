#include <iostream>
using namespace std;

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

			alpha[chan]=beta[chan]=-1.0;

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

	void CalibrationPeaksChannels(int chan, double ch1275, double ch511){
		if (chan < 0 || chan > 2)
			return;
		alpha[chan]=(1275.0-511.0)/(ch1275-ch511);
		beta[chan]=1275.0-alpha[chan]*ch1275;
		cout << "Channel corresponding to 1275 keV: " << ch1275 << endl;
		cout << "Channel corresponding to 511 keV: " << ch511 << endl;
		cout << "Energy(keV) = alpha * channel + beta" << endl;
		cout << "alpha = " << alpha[chan] << endl;
		cout << "beta = " << beta[chan] << endl;
	}

	void SetAlphaBeta(int chan, double a, double b){
		if (chan < 0 || chan > 2){
			cout<<"Invalid histogram index"<<endl;
			return;
		}
		alpha[chan] = a;
		beta[chan] = b;
	}

	void ChannelToEnergy(int chan, int c){
		if (chan < 0 || chan > 2)
			return;
		if(alpha[chan]==-1)
			cout << "Detector not calibrated" << endl;
		else
			cout << "Channel " << c << " corresponds to " << alpha[chan]*c+beta[chan] << " keV" << endl;
	}

	void EnergyToChannel(int chan, double e){
		if (chan < 0 || chan > 2)
			return;
		if(alpha[chan]==-1.0)
			cout << "Detector not calibrated" << endl;
		else{
			cout << e << "keV corresponds to channel " << (int)(e-beta[chan])/alpha[chan] << endl;
		}		
	}


	private:

	const char* original_data_file;
	TH1F* hist [4];
	double alpha[4], beta[4];
};
