


// NEW EDITION WITH NEW CONFIG OF CHANNELS
// (not all the metods are updated)




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

		for(int chan=0; chan< 3; chan++){

			alpha[chan]=beta[chan]=-1.0;

			slimport_data_t indata;
			TBranch *inbranch = (TBranch*)intree->GetBranch(Form("acq_ch%d",chan));
			inbranch->SetAddress(&indata.timetag);
			hist [chan] = new TH1F(Form("hist_ch%d",chan),"",8192,0,65536);

			for (int i=0; i<inbranch->GetEntries(); i++) {
				inbranch->GetEntry(i);
				hist [chan]->Fill(indata.qlong);
			}
		}
	}

	void AddData(const char *name_file){
		TFile* infile = new TFile(name_file);
		TTree* intree = (TTree*)infile->Get("acq_tree_0");

		for(int chan=0; chan< 3; chan++){

			slimport_data_t indata;
			TBranch *inbranch = (TBranch*)intree->GetBranch(Form("acq_ch%d",chan));
			inbranch->SetAddress(&indata.timetag);

			for (int i=0; i<inbranch->GetEntries(); i++) {
				inbranch->GetEntry(i);
				hist [chan]->Fill(indata.qlong);
			}
		}
	}

	TH1F* GetSumHistogram(int nchan) {
		
		if(nchan < 1 || nchan > 3){
			cout << "Invalid num_detect parameter" << endl;
			return NULL;
		}
		for(int i = 0; i < nchan; i++){
			if(alpha[i] == -1){
				cout << "First calibrate histograms" << endl;
				return NULL;
			}
		}

		TFile* infile = new TFile(original_data_file);
		TTree* intree = (TTree*)infile->Get("acq_tree_0");
		slimport_data_t indata[4];
		TBranch *inbranch[4];
		int curr_index[4];

		for(int chan=0; chan<4; chan++){
			inbranch[chan] = (TBranch*)intree->GetBranch(Form("acq_ch%d",chan));
			inbranch[chan]->SetAddress(&indata[chan].timetag);
			curr_index[chan]=0;		
		}

		TH1F* histsum  = new TH1F("HistSum","HistSum",10000,0,5000);
		
		int maxTimeChan=0;
		//double prog=0;
		//double step = inbranch[0]->GetEntries()/100;

		while (curr_index[0]<inbranch[0]->GetEntries() && curr_index[1]<inbranch[1]->GetEntries() && curr_index[2]<inbranch[2]->GetEntries() && curr_index[3]<inbranch[3]->GetEntries()){

			//if(curr_index[0]>prog+step){
			//	prog+=step;
			//	cout << (int)prog/step << " % \t->\t " << (int)prog << "/" << (int)step*100 << "\n";
			//}
			
			for(int c=0;c<4;c++)
				inbranch[c]->GetEntry(curr_index[c]);

			bool found_coinc=true;

			for(int c=1;c<4;c++)
				found_coinc &= (indata[0].timetag == indata[c].timetag);

			for(int c=1;c<4;c++){
				if(indata[0].timetag < indata[c].timetag){
					curr_index[0]++;
					break;
				}
				if(indata[0].timetag > indata[c].timetag)
					curr_index[c]++;
			}
			
			if(!found_coinc)
				continue;

			if(indata[3].qlong>maxTimeChan)
				maxTimeChan = indata[3].qlong;

			double e_sum = 0.0;

			for(int c = 0; c < nchan; c++){
				//if(GetChannelToEnergyValue(c, indata[c].qlong)<100 || GetChannelToEnergyValue(c, indata[c].qlong)>430)
				//	e_sum=-10000000;
				e_sum += GetChannelToEnergyValue(c, indata[c].qlong);
			}

			histsum->Fill(e_sum);

			for(int c=0; c<4; c++)
				curr_index[c]++;

		}

		cout << "Max entry for channel 3 (TAC): " << maxTimeChan << endl;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		return histsum;
	}

	//Get Histogram of channel with index "chan"
	TH1F* GetHistogram(int chan){
		if (chan < 0 || chan > 2)
			return 0x0;
		return hist[chan];
	}

	//Get path to the file with datas of the histograms
	const char* GetDataFile(){
		return original_data_file;
	}

	//Calibrate histogram through a relation E (or t) = alpha * channel + beta
	void CalibrateHisto(int chan = -1, double a = 0, double b = 0) {
		if (chan == -1) {
			CalibrateHisto(0);
			CalibrateHisto(1);
			CalibrateHisto(2);
		}
		if (chan < 0 || chan > 3)
			return;
		if (a == 0 ) { a = alpha[chan]; b = beta[chan]; }
		alpha[chan]=a;
		beta[chan]=b;
		TAxis *axis = hist[chan]->GetXaxis();
		axis->SetLimits(axis->GetXmin()*a+b, axis->GetXmax()*a+b);
		//hist[chan]->Scale(1/a); //to ensure keeping the proportions
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
	
	int e2ch(int chan, double e){
        if (chan < 0 || chan > 2)
			return;
		if(alpha[chan]==-1.0){
			cout << "Detector not calibrated" << endl;
            return -1;
        }
        else
            return (int)(e-beta[chan])/alpha[chan];
    }

	//Filter the spectra taking entries with sum of energy of first "num_detect" detectors equal to "energysum"
	//Accept a window cenetered in energy_sum equal to "window" 
	void SpectraFiltering(double energy_sum, double window, int num_detect, double lowbound, double upbound){
		
		if(num_detect < 1 || num_detect > 3){
			cout << "Invalid num_detect parameter" << endl;
			return;
		}
		if(energy_sum < 0 || window < 0){
			cout << "energy_sum and window must be positive" << endl;
            return;
		}
		for(int i = 0; i < num_detect; i++){
			if(alpha[i] == -1){
				cout << "First calibrate histograms" << endl;
				return;
			}
		}

		TFile* infile = new TFile(original_data_file);
		TTree* intree = (TTree*)infile->Get("acq_tree_0");
		slimport_data_t indata[4];
		TBranch *inbranch[4];
		int curr_index[4];

		for(int chan=0; chan<4; chan++){
			inbranch[chan] = (TBranch*)intree->GetBranch(Form("acq_ch%d",chan));
			inbranch[chan]->SetAddress(&indata[chan].timetag);
			hist [chan] = new TH1F(Form("hist_ch%d",chan),"",8192,0,65536);
			curr_index[chan]=0;		
		}

		
		int maxTimeChan=0;
		//double prog=0;
		//double step = inbranch[0]->GetEntries()/100;

		while (curr_index[0]<inbranch[0]->GetEntries() && curr_index[1]<inbranch[1]->GetEntries() && curr_index[2]<inbranch[2]->GetEntries() && curr_index[3]<inbranch[3]->GetEntries()){

			//if(curr_index[0]>prog+step){
			//	prog+=step;
			//	cout << (int)prog/step << " % \t->\t " << (int)prog << "/" << (int)step*100 << "\n";
			//}
			
			for(int c=0;c<4;c++)
				inbranch[c]->GetEntry(curr_index[c]);

			bool found_coinc=true;

			for(int c=1;c<4;c++)
				found_coinc &= (indata[0].timetag == indata[c].timetag);

			for(int c=1;c<4;c++){
				if(indata[0].timetag < indata[c].timetag){
					curr_index[0]++;
					break;
				}
				if(indata[0].timetag > indata[c].timetag)
					curr_index[c]++;
			}
			
			if(!found_coinc)
				continue;

			if(indata[3].qlong>maxTimeChan)
				maxTimeChan = indata[3].qlong;

			double e_sum = 0.0;

			for(int c = 0; c < num_detect; c++){
				if(GetChannelToEnergyValue(c, indata[c].qlong)<lowbound || GetChannelToEnergyValue(c, indata[c].qlong)>upbound)
					e_sum=-10000000;
				e_sum += GetChannelToEnergyValue(c, indata[c].qlong);
			}

			//If sum of energies are near "energy_sum" i can insert entries in histogram
			if(energy_sum - window/2 < e_sum && energy_sum + window/2 > e_sum){
				for(int c=0; c < 4; c++)
					hist [c]->Fill(indata[c].qlong);
			}

			for(int c=0; c<4; c++)
				curr_index[c]++;

		}

		cout << "Max entry for channel 3 (TAC): " << maxTimeChan << endl;
	}

	void ChannelFiltering(int minChan, int maxChan, int num_detect){
		
		if(num_detect < 0 || num_detect > 1){
			cout << "Invalid num_detect parameter" << endl;
			return;
		}
		if(minChan < 0 || maxChan < 0){
			cout << "minChan and maxChan must be positive" << endl;
		}

		TFile* infile = new TFile(original_data_file);
		TTree* intree = (TTree*)infile->Get("acq_tree_0");
		slimport_data_t indata[3];
		TBranch *inbranch[3];
		int curr_index[3];

		for(int chan=0; chan<3; chan++){
			inbranch[chan] = (TBranch*)intree->GetBranch(Form("acq_ch%d",chan));
			inbranch[chan]->SetAddress(&indata[chan].timetag);
			hist [chan] = new TH1F(Form("hist_ch%d",chan),"",8192,0,65536);
			curr_index[chan]=0;		
		}

		while (curr_index[0]<inbranch[0]->GetEntries() && curr_index[1]<inbranch[1]->GetEntries() && curr_index[2]<inbranch[2]->GetEntries()){
	
			for(int c=0;c<3;c++)
				inbranch[c]->GetEntry(curr_index[c]);

			bool found_coinc=true;

			for(int c=1;c<3;c++)
				found_coinc &= (indata[0].timetag == indata[c].timetag);

			for(int c=1;c<3;c++){
				if(indata[0].timetag < indata[c].timetag){
					curr_index[0]++;
					break;
				}
				if(indata[0].timetag > indata[c].timetag)
					curr_index[c]++;
			}
			
			if(!found_coinc)
				continue;

			if(minChan < indata[num_detect].qlong && maxChan > indata[num_detect].qlong){
				for(int c=0; c < 3; c++)
					hist [c]->Fill(indata[c].qlong);
			}

			for(int c=0; c<3; c++)
				curr_index[c]++;

		}
	}
	private:

	//Faster version of ChannelToEnergy function. Before use it verify that alpha, beta, and chan are currectly initialized.
	double GetChannelToEnergyValue(int chan, int c){
		return alpha[chan]*c+beta[chan];
	}

	const char* original_data_file;
	TH1F* hist [3];
	double alpha[3], beta[3];
};
