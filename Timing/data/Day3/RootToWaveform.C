#include <iostream>
#include <vector>
using namespace std;

struct slimport_data_t {
	ULong64_t	timetag; //time stamp
	UInt_t		baseline;
	UShort_t	qshort; //integration with shorter time
	UShort_t	qlong; //integration with longer time
	UShort_t	pur;
	UShort_t	samples[360];
};

struct waveform_event{
	UShort_t	value[360];
	TH1F* hist;
	TGraph* gr;
};

struct waveform_channel {
	int 	num_events;
	vector <waveform_event> events;
};

class Waveform{

	public: 

	//Costructor
	Waveform(const char *name_file, int ch){

		if (ch < 0 || ch > 1 )
			return 0x0;

		original_data_file = name_file;
		TFile* infile = new TFile(name_file);
		TTree* intree = (TTree*)infile->Get("acq_tree_0");
		TBranch *inbranch = (TBranch*)intree->GetBranch(Form("acq_ch%d",ch));
		slimport_data_t indata;
		inbranch->SetAddress(&indata.timetag);
		channel = ch;
		wc.num_events = inbranch->GetEntries();
		wc.events.resize(inbranch->GetEntries());		
		int step = inbranch->GetEntries() / 100;
		int prog = 0;
		sum_events = new TH2F(Form("Sum_events"),"", 360,0,360,1000,0,1000);
		
		for (int i=0; i<100000; i++) {
		//for (int i=0; i<inbranch->GetEntries(); i++) {
			if(i == prog){
				cout<<((float)i)/step<<"%"<<endl;
				prog += step;
			}

			inbranch->GetEntry(i);
			TH1F* h = new TH1F(Form("Sample_%d", i),"", 360, 0, 360);
			Int_t x[360], y[360];
	
			for(int j=0; j<360; j++){
				wc.events[i].value[j] = indata.samples[j];
				x[j]=j;
				y[j]=indata.samples[j];
				h->Fill(j,indata.samples[j]);
				sum_events->Fill(j,indata.samples[j]);
			}
			wc.events[i].gr = new TGraph(360,x,y);
			wc.events[i].hist = h;
		}
	}

	void DrawSumEvents(){ sum_events->Draw("colz"); }

	TH2F* GetSumHisto(){ return sum_events;}

	TGraph* GetGraph(int event){ return wc.events[event].gr; }

	waveform_event GetWaveform(int event){ return wc.events[event]; }

	TH1F* GetWaveformHisto(int event){ return wc.events[event].hist; }

	const char* GetDataFile(){ return original_data_file; }

	int GetChannel(){ return channel; }

	private:
	int 				channel;
	const char* 		original_data_file;
	TH2F*				sum_events;
	waveform_channel 	wc;
};
