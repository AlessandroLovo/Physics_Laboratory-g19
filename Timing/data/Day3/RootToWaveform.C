#include <iostream>
#include <vector>
using namespace std;

//This is useful to improve speed in loading heavy files. 
//Anyway in this version it tooks less than 30 sec to load 500Mb of datas (in my Mac u.u ) 
const int entries_analized = 50000;
const bool analize_all_entries = false;

//CFTD
const float attenuation_fraction = 0.25f;
const float delay_in_ns = 5.0f;

//WAVEFORM CLASS
const bool remove_cutted_energies = true;
const bool make_sum_histo_in_costructor = true;
const bool remove_zero_events_in_costructor = true;

//TWO_CHANNEL_WAVEFORM CLASS
const bool check_coincidences_in_costructor = true;

const int baseline[2] = {923, 925};

struct slimport_data_t {
	ULong64_t	timetag; //time stamp
	UInt_t		baseline;
	UShort_t	qshort; //integration with shorter time
	UShort_t	qlong; //integration with longer time
	UShort_t	pur;
	UShort_t	samples[360];
};

struct waveform_event{
	ULong64_t	timetag;
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
		int step = num_events / 100;
		int prog = 0;

		if(remove_zero_events_in_costructor){
			int new_num_events = 0;
			for (int i=0; i<num_events; i++) {
				bool no_zero_events = true;
				inbranch->GetEntry(i);
				waveform_event we;
				for(int j=0; j<360 && no_zero_events; j++){
					we.value[j] = indata.samples[j];
					no_zero_events &= (indata.samples[j] != 0);
				}
				if(no_zero_events){
					events.push_back(we);
					new_num_events++;
				}
			}
			num_events = new_num_events;
		}
		else{
			events.resize(num_events);
			for (int i=0; i<num_events; i++) {
				inbranch->GetEntry(i);
				events[i].timetag = indata.timetag;
				for(int j=0; j<360; j++)
					events[i].value[j] = indata.samples[j];
			}
		}
		

		if(make_sum_histo_in_costructor) MakeSumHisto();
	}

	Waveform* GetOtherChannel(){ return new Waveform(original_data_file, (int)(channel == 0) );	}

	void RemoveZeroEnergyEvents(){
		int new_num_events = 0;
		vector <waveform_event> filtered_events;
		for (int i=0; i<num_events; i++) {
			bool no_zero_events = true;
			for(int j=0; j<360 && no_zero_events; j++)
				no_zero_events &= (events[i].value[j] != 0);
			if(no_zero_events){
				filtered_events.push_back(events[i]);
				new_num_events++;
			}
		}
		events = filtered_events;
		num_events = new_num_events;
		MakeSumHisto();	}

	TH2F* MakeSumHisto(){
		sum_events = new TH2F(Form("Sum_events"),"", 360,0,360,256,0,1024);
		for (int i=0; i<num_events; i++)
			for(int j=0; j<360; j++)
				sum_events->Fill(j, events[i].value[j]);
		return sum_events; }

	TH2F* GetSumHisto(){ return sum_events;}

	void DrawSumHisto(){ sum_events->Draw("colz"); }

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

	void CheckCoincidences(Waveform* wf1){
		vector <waveform_event> :: iterator it0 = events.begin(), it1 = wf1->events.begin();
		while(it0 != events.end() && it1 < wf1->events.end()){
			if ((*it0).timetag < (*it1).timetag)
				events.erase(it0);
			else if ((*it0).timetag > (*it1).timetag)
				wf1->events.erase(it1);
			else{
				it0++;
				it1++;
			}
		}
		MakeSumHisto();
		wf1->MakeSumHisto();}

	TGraph* CFTD(int e){
		UShort_t* ev = events[e].value;
		TGraph* gr = GetGraph(e);
		float x[720], y[720];
		for(int i=0; i < 360; i++)
			x[i] = (float)i;
		for(int i=0; i<=(int)delay_in_ns; i++)
			y[i] = 0.0;
		for(int i=(int)delay_in_ns+1; i<360; i++)
			y[i] = ((float)ev[i] - baseline[channel]) * attenuation_fraction + ((float)baseline[channel] - gr->Eval(i-delay_in_ns));
		for(int i=0; i<360; i++)
			x[360+i]=(float)i+delay_in_ns;
		for(int i=0; i <360-((int)delay_in_ns+1); i++)
			y[360+i]=(float)baseline[channel]-ev[i] + (gr->Eval(i+delay_in_ns) - (float)baseline[channel]) * attenuation_fraction;
		for(int i=360-((int)delay_in_ns+1); i<360; i++)
			y[360+i]=0;
		return new TGraph(720, x, y);
	}

	private:
	int 				channel;
	int 				num_events;
	const char* 		original_data_file;
	TH2F*				sum_events;
	vector <waveform_event> events;
};

class TwoChannelWaveform{
	public:

	TwoChannelWaveform(const char *name_file){
		original_data_file = name_file;
		wf[0] = new Waveform(original_data_file, 0);
		wf[1] = new Waveform(original_data_file, 1);
		if (check_coincidences_in_costructor)
			CheckCoincidences();
	}

	Waveform* at(int ch){
		if (ch < 0 || ch > 1 )
			return 0x0;
		return wf[ch];}

	void CheckCoincidences(){
		wf[0]->CheckCoincidences(wf[1]);
	}

	private:
	Waveform* 		wf [2];
	const char* 	original_data_file;
};
