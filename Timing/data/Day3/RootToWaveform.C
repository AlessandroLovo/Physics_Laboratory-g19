#include <iostream>
#include <vector>
using namespace std;

const int entries_analized = 500; //000;
const bool analize_all_entries = false;

//CFTD
float attenuation_fraction = 0.2f;
float delay_in_ns_ch0 = 5.0f;
float delay_in_ns_ch1 = 5.0f;
const int steps_in_zero_crossing_binary_search = 13; //-> precision < 1ps
const int steps_in_delay_optimization = 15;

//WAVEFORM CLASS
const bool remove_cutted_energies = true;
const bool make_sum_histo_in_costructor = false;
const bool remove_zero_events_in_costructor = true;
const int  aligment_channel = 100;

//TWO_CHANNEL_WAVEFORM CLASS
const bool energy_filter_in_costructor = true;
const bool check_coincidences_in_costructor = true;
const bool calculate_time_distribution_in_each_time_distr_hist_request = false;
const double delay_introduced_in_ns = 150;

//ENERGY FILTERING WINDOW
const int energy_window_center_channel[2] = {6000, 6000};
const int energy_window_half_size[2] = {100, 100};
const double e_cal_m[2] = {0.08199, 0.06579};
const double e_cal_q[2] = {-208.5, -106.8};

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
	UShort_t	energy;
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
		

		if(make_sum_histo_in_costructor) MakeSumHisto();}

	~Waveform() {
		if ( sum_events != nullptr ) delete sum_events;
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
		num_events = new_num_events;}

	TH2F* MakeSumHisto(){
		sum_events = new TH2F(Form("Sum_events_%d",channel),"", 360,0,360,256,0,1024);
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

	TGraph* GetGraphWithoutBaseline(int event) {
		TGraph* tg = GetGraph(event);
		for( int i=0 ; i < tg->GetN(); i++) {
			double x,y;
			tg->GetPoint(i, x,y);
			tg->SetPoint(i, x,y-(double)baseline[channel]);
		}
		return tg;
	}

	const char* GetDataFile(){ return original_data_file; }

	int GetChannel(){ return channel; }

	void CheckCoincidences(Waveform* wf1){
		vector <waveform_event> :: iterator it0 = events.begin(), it1 = wf1->events.begin();
		while(it0 != events.end() && it1 < wf1->events.end()){
			if ((*it0).timetag < (*it1).timetag)
				events.erase(it0);
			else if ((*it0).timetag > (*it1).timetag)
				wf1->events.erase(it1);
			else{ it0++; it1++; }

		}
		//MakeSumHisto();
		//wf1->MakeSumHisto();
	}

	void EnergyFiltering(double elow = -1, double emax = -1){
		//GetEnergyHisto();	
		fast_CalculateEventEnergies();
		int new_num_events = 0;
		vector <waveform_event> filtered_events;
		for (int i=0; i<num_events; i++) {
			if (elow == -1) {
				if(TMath::Abs(events[i].energy - energy_window_center_channel[channel]) < energy_window_half_size[channel]){
					filtered_events.push_back(events[i]);
					new_num_events++;
				}
			}
			else {
				double e_en = events[i].energy * e_cal_m[channel] + e_cal_q[channel];
				if ( e_en >= elow && e_en <= emax ) {
					filtered_events.push_back(events[i]);
					new_num_events++;
				}
			}
		}
		events = filtered_events;
		num_events = new_num_events;
		//GetEnergyHisto();
	}

	TGraph* CFTD(int e){
		UShort_t* ev = events[e].value;
		TGraph* gr = GetGraph(e);
		double delay_in_ns = (channel==0 ? delay_in_ns_ch0 : delay_in_ns_ch1);
		double x[720], y[720];
		for(int i=0; i < 360; i++)
			x[i] = (double)i;
		for(int i=0; i<=(int)delay_in_ns; i++)
			y[i] = 0.0;
		for(int i=(int)delay_in_ns+1; i<360; i++)
			y[i] = ((double)ev[i] - baseline[channel]) * attenuation_fraction + ((double)baseline[channel] - gr->Eval(i-delay_in_ns));
		for(int i=0; i<360; i++)
			x[360+i]=(double)i+delay_in_ns;
		for(int i=0; i <360-((int)delay_in_ns+1); i++)
			y[360+i]=(double)baseline[channel]-ev[i] + (gr->Eval(i+delay_in_ns) - (double)baseline[channel]) * attenuation_fraction;
		for(int i=360-((int)delay_in_ns+1); i<360; i++)
			y[360+i]=0.0;
		delete gr;
		return new TGraph(720, x, y);
	}

	TGraph* CFTD(int e, double CFTD_delay){
		UShort_t* ev = events[e].value;
		TGraph* gr = GetGraph(e);
		double x[720], y[720];
		for(int i=0; i < 360; i++)
			x[i] = (double)i;
		for(int i=0; i<=(int)CFTD_delay; i++)
			y[i] = 0.0;
		for(int i=(int)CFTD_delay+1; i<360; i++)
			y[i] = ((double)ev[i] - baseline[channel]) * attenuation_fraction + ((double)baseline[channel] - gr->Eval(i-CFTD_delay));
		for(int i=0; i<360; i++)
			x[360+i]=(double)i+CFTD_delay;
		for(int i=0; i <360-((int)CFTD_delay+1); i++)
			y[360+i]=(double)baseline[channel]-ev[i] + (gr->Eval(i+CFTD_delay) - (double)baseline[channel]) * attenuation_fraction;
		for(int i=360-((int)CFTD_delay+1); i<360; i++)
			y[360+i]=0.0;
		return new TGraph(720, x, y); }

	double ZeroCrossing(int e, double CFTD_delay){
		TGraph* gr = CFTD(e, CFTD_delay);
		double* x = gr->GetX();
		double* y = gr->GetY();
		double xmin, xmax, ymin=0, ymax=0;
		for(int i=0; i<720; i++){
			if(y[i]<ymin){xmin=x[i];ymin=y[i];}
			if(y[i]>ymax){xmax=x[i];ymax=y[i];}
		}
		double r = xmax, l = xmin;
		for(int i=0;i<steps_in_zero_crossing_binary_search;i++){
			if(gr->Eval((l+r)/2)<0) l=(l+r)/2;
			else r=(l+r)/2;
			//cout << r-l <<endl;
		}
		return (l+r)/2;}

	double ZeroCrossing(int e){
		TGraph* gr = CFTD(e);
		double* x = gr->GetX();
		double* y = gr->GetY();
		double xmin, xmax, ymin=0, ymax=0;
		for(int i=0; i<720; i++){
			if(y[i]<ymin){xmin=x[i];ymin=y[i];}
			if(y[i]>ymax){xmax=x[i];ymax=y[i];}
		}
		double r = xmax, l = xmin;
		for(int i=0;i<steps_in_zero_crossing_binary_search;i++){
			if(gr->Eval((l+r)/2)<0) l=(l+r)/2;
			else r=(l+r)/2;
			//cout << r-l <<endl;
		}
		delete gr;
		return (l+r)/2;
	}

	TH1F* ZeroCrossingDistribution(){
		TH1F* h = new TH1F("Zero Crossing Distr.","", 360*1000/4, 0,360*1000);
		for(int i=0;i<events.size();i++) h->Fill(ZeroCrossing(i)*1000);
		h->Draw();
		return h;}

	TH1F* ZeroCrossingDistribution(double delay_in_ns){
		TH1F* h = new TH1F("Zero Crossing Distr.","", 360*1000/4, 0,360*1000);
		for(int i=0;i<events.size();i++) h->Fill(ZeroCrossing(i, delay_in_ns)*1000);
		//h->Draw();
		return h;}

	//Try to find optimal delay aligning (minimal) peaks of raw samples to the value "aligment_channel"
	//and then minimizing variation of crossing points in CFTD from "aligment_channel".
	TGraph* OptimalCFTDDelay(double start, double stop){
		EnergyFiltering();
		AlignEventPeaks();
		double x[steps_in_delay_optimization], y[steps_in_delay_optimization];
		for(int i=0; i<steps_in_delay_optimization; i++){
			x[i] = start + (stop-start)/(steps_in_delay_optimization-1)*i;
			y[i]=0;
			for(int j=0; j < events.size(); j++) y[i] += TMath::Power(ZeroCrossing(j, x[i])-aligment_channel,2);
			cout << "Delay = " << x[i] << "ns, sum of square variations from peaks alignment = " << y[i] << endl;
		}
		TGraph* gr = new TGraph(steps_in_delay_optimization,x,y);
		gr->Draw();
		return gr;
	}

	TH1F* GetEnergyHisto(){
		TH1F* h = new TH1F("Distribution_of_energies","",300, 0, 20000);
		for(int i=0; i<num_events;i++){
			h->Fill(CalculateEventEnergy(i));
		}
		//h->Draw();
		return h;}

	int GetNumEvents(){return num_events;}

	void AlignEventPeaks(){
		for(vector<waveform_event>::iterator it = events.begin(); it != events.end(); it++){
			int xmin=0, ymin=1000;
			for(int i=0;i<360;i++) if((*it).value[i]<ymin) {xmin=i; ymin=(*it).value[i];}
			int delta = aligment_channel - xmin;
			if(delta > 0) for(int i = 360-1; i > delta; i--) (*it).value[i]=(*it).value[i-delta];
			else if (delta < 0) for(int i=0; i < 360+delta;i++) (*it).value[i]=(*it).value[i-delta];
		}
	}

	private:

	int CalculateEventEnergy(int e){
		int x[360], y[360];
		for(int i=0; i < 360; i++){
			x[i]=i; y[i]=baseline[channel]-events[e].value[i]+5;
		}
		y[0]=y[360-1]=0;
		//(new TGraph(360, x, y))->Draw();
		events[e].energy = ((new TGraph(360, x, y))->Integral());
		return events[e].energy;}

	void fast_CalculateEventEnergies() {
		for (int e=0; e< events.size(); e++) {
			double sum = 0;
			for(int i=0; i < 360; i++) sum += baseline[channel]-events[e].value[i]+5;
			events[e].energy = sum;
		}
	}

	int 				channel;
	int 				num_events;
	const char* 		original_data_file;
	TH2F*				sum_events = nullptr;
	vector <waveform_event> events;
};





class TwoChannelWaveform{
	public:

	TwoChannelWaveform(const char *name_file, double att_frac = -1, double delay = -1, double elow = -1, double etop = -1){
		if( att_frac != -1 ) attenuation_fraction = att_frac;
		if( delay != -1 ) { delay_in_ns_ch0 = delay; delay_in_ns_ch1 = delay; }
		original_data_file = name_file;
		wf[0] = new Waveform(original_data_file, 0);
		wf[1] = new Waveform(original_data_file, 1);
		if (energy_filter_in_costructor || elow != -1)
			EnergyFiltering(elow, etop);
		if (check_coincidences_in_costructor)
			CheckCoincidences();
		}

	~TwoChannelWaveform() {
		delete wf[0];
		delete wf[1];
	}

	Waveform* at(int ch){
		if (ch < 0 || ch > 1 )
			return 0x0;
		return wf[ch];}

	void CheckCoincidences(){ wf[0]->CheckCoincidences(wf[1]); }

	void CalculateTimeDistribution(){
		if ( !check_coincidences_in_costructor ) CheckCoincidences();
		int num_events = wf[0]->GetNumEvents();
		//delta_time_distribution.resize(num_events);

		// To have something like a progress bar
		int step = num_events / 20;
		cout << setw(3) << 0 <<"%" <<flush;


		for(int i=0; i < num_events; i++){
			delta_time_distribution . push_back ( wf[0]->ZeroCrossing(i) - wf[1]->ZeroCrossing(i) + delay_introduced_in_ns);
			//delta_time_distribution[i]=wf[0]->ZeroCrossing(i) - wf[1]->ZeroCrossing(i) + delay_introduced_in_ns;

			if( i % step == 0 ) cout << "\b\b\b\b" << setw(3) << i/step*5 <<"%" <<flush;

		}
	}

	TH1F* GetTimeDistrHisto(ofstream* out = NULL, double elow = -1, double etop = -1){
		if( calculate_time_distribution_in_each_time_distr_hist_request || delta_time_distribution.size() < 1 )
			CalculateTimeDistribution();
		TH1F* h = new TH1F(Form("Time_distribution_in_ps_%f_%f", attenuation_fraction, delay_in_ns_ch0 ), "", 5000,0, 360*1000); // 5000 (experimentally found binning that permit noise neglecting keeping sufficient resolution)
		for(int i=0;i<delta_time_distribution.size();i++)
			h->Fill(int(delta_time_distribution[i]*1000));
		//h->Draw();
				// Print FWHM with error
		// Find Maximum
		double max_half = h->GetMaximum() / 2;
		double max_bin = h->GetMaximumBin();
		int n_bin = h->GetNbinsX();
		double last_below_before = -1, first_above_before, last_above_after = -1, first_below_after;

		for (int i=0; i<n_bin; i++) {
			double c = h->GetBinContent(i);
			if ( c > max_half && last_below_before == -1 ) last_below_before = i-1;
			if ( c < max_half && i < max_bin ) first_above_before = i+1;
			if ( c < max_half && i > max_bin && last_above_after == -1 ) last_above_after = i-1;
			if ( c > max_half ) first_below_after = i+1;
		}

		last_below_before = h->GetBinCenter( last_below_before );
		first_above_before = h->GetBinCenter( first_above_before );
		last_above_after = h->GetBinCenter( last_above_after );
		first_below_after = h->GetBinCenter( first_below_after );
		
		double FWHM = ( first_below_after + last_above_after ) / 2 - ( first_above_before + last_below_before ) / 2 ;
		double FWHM_err = sqrt( pow( first_above_before - last_below_before ,2) + pow( first_below_after - last_above_after ,2) ) / sqrt(12);
		double mean = h->GetMean();
		double mean_sigma = h->GetMeanError();
		double kurt = h->GetKurtosis();
		double kurt_sigma = h->GetKurtosis(11);
		
		/*
		// COMPUTE WITHOUT FUNCTIONS
		mean = 0;
		for (int i=0; i<n_bin; i++) {
			mean += ( h->GetBinCenter(i) * h->GetBinContent(i) );
		}
		mean /= h->GetEntries();
		double sum2 = 0, sum4 = 0;
		for (int i=0; i<n_bin; i++) {
			sum2 += ( pow(h->GetBinCenter(i) - mean, 2) * h->GetBinContent(i) );
			sum4 += ( pow(h->GetBinCenter(i) - mean, 4) * h->GetBinContent(i) );
		}
		sum2 /= h->GetEntries();
		sum4 /= h->GetEntries();
		kurt = sum4 / sum2 / sum2 - 3;

						   cout << endl << attenuation_fraction << '\t' << delay_in_ns_ch0 << '\t' << mean << '\t' << mean_sigma << '\t' << FWHM << '\t' << FWHM_err << '\t' << kurt << '\t' << kurt_sigma <<endl;
		*/
		/*
			COMPUTE DIRECTLY FROM ARRAY
		mean = 0;
		for (int i=0; i< delta_time_distribution.size(); i++) {
			mean += delta_time_distribution[i]*1000;
		}
		mean /= delta_time_distribution.size();
		sum2 = 0, sum4 = 0;
		for (int i=0; i<delta_time_distribution.size(); i++) {
			sum2 += ( pow( delta_time_distribution[i] *1000 - mean, 2) );
			sum4 += ( pow( delta_time_distribution[i] *1000 - mean, 4) );
		}
		sum2 /= delta_time_distribution.size();
		sum4 /= delta_time_distribution.size();
		kurt = sum4 / sum2 / sum2 - 3;
						   cout << endl << attenuation_fraction << '\t' << delay_in_ns_ch0 << '\t' << mean << '\t' << mean_sigma << '\t' << FWHM << '\t' << FWHM_err << '\t' << kurt << '\t' << kurt_sigma <<endl;
		cout << delta_time_distribution.size();	
		*/

		if(elow != -1) {
		 		  	   cout << endl << elow << '\t' << etop << '\t' << mean << '\t' << mean_sigma << '\t' << FWHM << '\t' << FWHM_err << '\t' << kurt << '\t' << kurt_sigma <<endl;
			if ( out != NULL ) *out << elow << '\t' << etop << '\t' << mean << '\t' << mean_sigma << '\t' << FWHM << '\t' << FWHM_err << '\t' << kurt << '\t' << kurt_sigma <<endl;
		
		} else {
					   cout << endl << attenuation_fraction << '\t' << delay_in_ns_ch0 << '\t' << mean << '\t' << mean_sigma << '\t' << FWHM << '\t' << FWHM_err << '\t' << kurt << '\t' << kurt_sigma <<endl;
			if ( out != NULL ) *out << attenuation_fraction << '\t' << delay_in_ns_ch0 << '\t' << mean << '\t' << mean_sigma << '\t' << FWHM << '\t' << FWHM_err << '\t' << kurt << '\t' << kurt_sigma <<endl;
		}
		return h;
	}

	void GetTimeDistrHisto_withstd(ofstream* out = NULL, double elow = -1, double etop = -1) {
		if( calculate_time_distribution_in_each_time_distr_hist_request || delta_time_distribution.size() < 1 )
			CalculateTimeDistribution();

		double sum = 0, sum2 = 0, t;
		int n = delta_time_distribution.size();
		for(int i = 0; i < n; i++) {
			t = delta_time_distribution[i];
			sum += t;
			sum2+= (t*t);
		}

		double mean = sum / n;
		double mean_sigma = ( sum2 / n - (mean*mean) ) / sqrt(n);

		sum = 0;
		sum2 = 0;
		for(int i = 0; i < n; i++) {
			t = ( delta_time_distribution[i] - mean ) * ( delta_time_distribution[i] - mean );
			sum += t;
			sum2+= (t*t);
		}
		double mu2 = sqrt ( sum / n );
		double mu2_sigma = ( sum2 / n - (mu2*mu2) ) / sqrt(n);
		double std = sqrt ( mu2 );
		double std_sigma = mu2_sigma / 2 / std;


		sum = 0;
		sum2 = 0;
		for(int i = 0; i < n; i++) {
			t = ( delta_time_distribution[i] - mean ) * ( delta_time_distribution[i] - mean ) * ( delta_time_distribution[i] - mean ) * ( delta_time_distribution[i] - mean );
			sum += t;
			sum2+= (t*t);
		}
		double mu4 = sum / n;
		double mu4_sigma = ( sum2 / n - (mu4*mu4) ) / sqrt(n);	

		double kurt = 3 - ( mu4 / std / std );
		double kurt_sigma = sqrt( mu4_sigma / std / std * mu4_sigma / std / std + 2 * mu4 * std_sigma / std / std / std * 2 * mu4 * std_sigma / std / std / std );

		if(elow != -1) {
		 		  	   cout << endl << elow << '\t' << etop << '\t' << mean << '\t' << mean_sigma << '\t' << std << '\t' << std_sigma << '\t' << kurt << '\t' << kurt_sigma <<endl;
			if ( out != NULL ) *out << elow << '\t' << etop << '\t' << mean << '\t' << mean_sigma << '\t' << std << '\t' << std_sigma << '\t' << kurt << '\t' << kurt_sigma <<endl;
		
		} else {
					   cout << endl << attenuation_fraction << '\t' << delay_in_ns_ch0 << '\t' << mean << '\t' << mean_sigma << '\t' << std << '\t' << std_sigma << '\t' << kurt << '\t' << kurt_sigma <<endl;
			if ( out != NULL ) *out << attenuation_fraction << '\t' << delay_in_ns_ch0 << '\t' << mean << '\t' << mean_sigma << '\t' << std << '\t' << std_sigma << '\t' << kurt << '\t' << kurt_sigma <<endl;
		}
	}

	void EnergyFiltering(double elow, double etop){
		wf[0]->EnergyFiltering(elow, etop);
		wf[1]->EnergyFiltering(elow, etop);}

	private:
	vector<double>	delta_time_distribution;
	Waveform* 		wf [2];
	const char* 	original_data_file;
};

void simulateCFTD(int id = -1) {
	if (id == -1) { simulateCFTD(0); simulateCFTD(1); return; }
	const char* outfilename[] = {"CFTD_simulations_1_2D.txt","CFTD_simulations_2_2D.txt"};
	const char* sourcename[] = {"Digital_CFTD.root", "Digital_CFTD_2.root"};
	vector<double> fracs{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
	vector<double> delays{1 , 1.25 , 1.5 , 1.75 , 2 , 2.25 , 2.5 , 2.75 , 3 , 3.25 , 3.5 , 3.75 , 4 , 4.25 , 4.5 , 4.75 , 5 , 5.25 , 5.5 , 5.75 , 6 , 6.25 , 6.5 , 6.75 , 7 , 7.25 , 7.5 , 7.75 , 8 , 8.25 , 8.5 , 8.75 , 9 , 9.25 , 9.5 , 9.75 , 10};
	int i = 0;
	ofstream out(outfilename[id]);
	out << "Frac\tDelay\tMean\tMean_sigma\tFWHM\tFWHM_sigma\tKurtosis\tKurtosis_sigma"<<endl;
	for ( double f : fracs)
		for ( double d : delays) {
			i++;
			cout<<setw(5)<<f<<" - "<<setw(5)<<d<<" - "<<setw(2)<<i<<" of "<<setw(2)<<fracs.size()*delays.size()<<" - Loading data; "<<flush;
			auto tcw = new TwoChannelWaveform(sourcename[id], f, d);
			cout<<"Analizing data: "<<flush;
			auto tdh = tcw->GetTimeDistrHisto(&out);
			delete tdh;
			delete tcw;
		}
}

void simulateCFTD_energythresh(int id = -1, bool std = false) {
	if (id == -1) { simulateCFTD_energythresh(0, std); simulateCFTD_energythresh(1, std); return; }
	const char* outfilename[] = {"CFTD_energythresh_1_2D.txt","CFTD_energythresh_2_2D.txt"};
	const char* sourcename[] = {"Digital_CFTD.root", "Digital_CFTD_2.root"};
	vector<double> energy_low{ 50, 100, 150, 200, 250, 300, 350,  50, 100, 150, 200, 250 };
	vector<double> energy_top{600, 600, 600, 600, 600, 600, 600, 150, 250, 350, 450, 550 };
	ofstream out(outfilename[id]);
	if ( std ) out << "ELow\tETop\tMean\tMean_sigma\tVar\tVar_sigma\tKurtosis\tKurtosis_sigma"<<endl;
	else out << "ELow\tETop\tMean\tMean_sigma\tFWHM\tFWHM_sigma\tKurtosis\tKurtosis_sigma"<<endl;
	for( int i=0; i<energy_low.size(); i++) {
		cout<<setw(5)<<energy_low[i]<<" - "<<energy_top[i]<<" - "<<setw(2)<<i<<" of "<<setw(2)<<energy_low.size()<<" - Loading data; "<<flush;
		auto tcw = new TwoChannelWaveform(sourcename[id], -1, -1, energy_low[i], energy_top[i]);
		cout<<"Analizing data: "<<flush;
		if (std) tcw->GetTimeDistrHisto_withstd(&out,energy_low[i],energy_top[i]);
		else { auto tdh = tcw->GetTimeDistrHisto(&out,energy_low[i],energy_top[i]); delete tdh; }
		delete tcw;
	}
}