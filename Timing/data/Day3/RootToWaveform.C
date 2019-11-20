#include <iostream>
#include <vector>
using namespace std;

const int entries_analized = 100000;
const bool analize_all_entries = true;

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
		num_events = entries_analized;
		if(analize_all_entries)
			num_events = inbranch->GetEntries();
		events.resize(num_events);		
		int step = num_events / 100;
		int prog = 0;
		sum_events = new TH2F(Form("Sum_events"),"", 360,0,360,256,0,1024);

		for (int i=0; i<num_events; i++) {
			inbranch->GetEntry(i);
			for(int j=0; j<360; j++){
				events[i].value[j] = indata.samples[j];
				sum_events->Fill(j,indata.samples[j]);
			}
		}
	}

	void DrawSumEvents(){ sum_events->Draw("colz"); }

	TH2F* GetSumHisto(){ return sum_events;}

	TGraph* GetGraph(int event){
		Int_t x[360], y[360];
		for(int i=0;i<360;i++){
			x[i]=i;
			y[i]=events[event].value[i];
		}
		return new TGraph(360,x,y);	}

	waveform_event GetEvent(int event){ return events[event]; }

	const char* GetDataFile(){ return original_data_file; }

	int GetChannel(){ return channel; }

	//Floor: 923.1 +- 1.8 -> Mean of noise before peak in sum of all events

	private:
	int 				channel;
	int 				num_events;
	const char* 		original_data_file;
	TH2F*				sum_events;
	vector <waveform_event> events;
};
