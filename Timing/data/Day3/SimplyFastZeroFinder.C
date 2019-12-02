const int max_events = 50000;

#ifndef slimport
#define slimport
struct slimport_data_t {
	ULong64_t	timetag; //time stamp
	UInt_t		baseline;
	UShort_t	qshort; //integration with shorter time
	UShort_t	qlong; //integration with longer time
	UShort_t	pur;
	UShort_t	samples[360];
};
#endif

struct statistic
{
    double value;
    double error;
};


class Event {
public:

    Event (ULong64_t timetag, UShort_t samples[360]) {
        this->timetag = timetag;
        baseline = 0;
        double sum = 0;
        zero_sample=false;
        for( int i = 0; i < 360; i++) {
            this->samples[i] = samples[i];
            if ( i < 15 ) baseline += samples[i];
            zero_sample |= (samples[i]==0);
            sum += samples[i];
        }
        baseline /= 15.0;
        energy = baseline*360.0 - sum;
    }

    double GetZeroCrossing(double frac,int delay) {
        double CFTD[360];
        for(int i=0; i < 360; i++) CFTD[i] = - frac * (baseline - samples[i]);
        for(int i=delay; i < 360; i++) CFTD[i] += (baseline - samples[i - delay]);

        // Find minimum and maximum crossing
        int min = 0, max = 0;
        for(int i=1; i<360; i++) {
            if ( CFTD[i] < CFTD[min]) min = i;
            if ( CFTD[i] > CFTD[max]) max = i;
        }
        // Find zero crossing
        while ( CFTD[min+1] < 0 ) min++;
        while ( CFTD[max-1] > 0 ) max--;
        return (double)max - CFTD[max] * ( (double)max - (double)min ) / ( CFTD[max] - CFTD[min] );
    }

    long timetag;
    short samples[360];
    double baseline;
    double energy;
    bool zero_sample;
};

class SimplyFastZeroFinder {
public:
    SimplyFastZeroFinder(const char *name_file) {
        TFile* infile = new TFile(name_file);
		TTree* intree = (TTree*)infile->Get("acq_tree_0");
		slimport_data_t indata;
        int num_events;

        // CHANNEL 0
		TBranch *inbranch = (TBranch*)intree->GetBranch("acq_ch0");
		inbranch->SetAddress(&indata.timetag);
        num_events = min( (int)inbranch->GetEntries(), max_events);
        for (int i=0; i<num_events; i++) {
				inbranch->GetEntry(i);
                Event* e = new Event(indata.timetag, indata.samples);
                if (!(e->zero_sample))
                    ch0.push_back( e );
		}

        // CHANNEL 1
		inbranch = (TBranch*)intree->GetBranch("acq_ch1");
		inbranch->SetAddress(&indata.timetag);
        num_events = min( (int)inbranch->GetEntries(), max_events);
        for (int i=0; i<num_events; i++) {
				inbranch->GetEntry(i);
                Event* e = new Event(indata.timetag, indata.samples);
                if (!(e->zero_sample))
                    ch1.push_back( e );
		}

        // CHECK COINCIDENCES
        cout<<"Event before:\t"<<ch0.size()<<'\t'<<ch1.size()<<endl;
        int i0=0, i1=0;
        while ( i0 < ch0.size() && i1 < ch1.size() ) {
            if ( ch0[i0]->timetag < ch1[i1]->timetag ) {
                ch0.erase(ch0.begin() + i0);
            }
            else if ( ch0[i0]->timetag > ch1[i1]->timetag ) {
                ch1.erase(ch1.begin() + i1);
            }
            else { i0++; i1++; }
        }
        ch0.resize(min(ch0.size(),ch1.size()));
        ch1.resize(min(ch0.size(),ch1.size()));
        cout<<"Event after:\t"<<ch0.size()<<'\t'<<ch1.size()<<endl;

        enabled = new vector<bool>(ch0.size(), true);
    }

    vector<double>* getWidths(double frac, int delay) {
        vector<double>* width = new vector<double>();
        for( int i=0; i<ch0.size(); i++ ) {
            if( enabled->at(i) )
                width->push_back ( ch0[i]->GetZeroCrossing(frac,delay) - ch1[i]->GetZeroCrossing(frac,delay) );
        }
        return width;
    }

    void getProprieties(vector<double>* vec, statistic &mean, statistic &FWHM, statistic &kurt, bool delete_vec = false) {
        
        double sum = 0;
		int n = vec->size();
		for(int i = 0; i < n; i++) {
			sum += vec->at(i);
		}

		mean.value = sum / n;

		double sum2 = 0.0, sum4 = 0.0;
		for(int i = 0; i < n; i++) {
             
            
			sum2 += ( vec->at(i) - mean.value ) * ( vec->at(i) - mean.value );
			sum4 += ( vec->at(i) - mean.value ) * ( vec->at(i) - mean.value ) * ( vec->at(i) - mean.value ) * ( vec->at(i) - mean.value );
		}

		mean.error = sqrt ( sum2 ) / n;

        if( mean.value < 20 && mean.value > 5 ) {
            TH1F* histo = getHisto( vec );
            FWHM = this->FWHM( histo );
            delete histo;
        } else { FWHM.value = 100; FWHM.error = 100; }
        /*
		statistic mu2, mu4;
        mu2.value = sum2 / n;
		mu4.value = sum4 / n;

		double sum2e = 0.0, sum4e = 0.0, t;
		for(int i = 0; i < n; i++) {
			t = ( vec->at(i) - mean.value ) * ( vec->at(i) - mean.value );
			sum2e = ( t - mu2.value )*( t - mu2.value );
			
			t = ( vec->at(i) - mean.value ) * ( vec->at(i) - mean.value ) * ( vec->at(i) - mean.value ) * ( vec->at(i) - mean.value );
			sum4e = ( t - mu4.value )*( t - mu4.value );
		}
		mu2.error = sqrt ( sum2e ) / sqrt(n);
		mu4.error = sqrt ( sum4e ) / sqrt(n);
		
		std.value = sqrt ( mu2.value );
		std.error = mu2.error / 2.0 / sqrt ( mu2.value );
		kurt.value = mu4.value / mu2.value / mu2.value - 3.0;
		kurt.error = sqrt( mu4.error / mu2.value / mu2.value * mu4.error / mu2.value / mu2.value + 2.0 *  mu4.value * mu2.error / mu2.value / mu2.value / mu2.value * 2 * mu4.value * mu2.error / mu2.value / mu2.value / mu2.value );
        */
        if(delete_vec) delete vec;
    }

    TH1F* getHisto(vector<double>* vec) {
        TH1F* hist = new TH1F ( "hist", "hist", 2000, 0, 20);
        for (int i=0; i < vec->size(); i++) hist->Fill( vec->at(i) );
        return hist;
    }

    TH1F* getEnergyHisto(int ch) {
        TH1F* hist = new TH1F ( Form("Energy_%d",ch), Form("Energy_%d",ch), 1000, 0, 100000);
        vector<Event*>* vec = ( ch == 0 ? &ch0 : &ch1 );
        for (int i=0; i < vec->size(); i++) hist->Fill( vec->at(i)->energy );
        return hist;
    }

    void calibrateHisto(TH1F* histo, double m, double q) {
        auto xax = histo->GetXaxis();
        xax->SetLimits( (double)xax->GetXmin() * m + q, (double)xax->GetXmax() * m + q );
    }

    void energyFilter(double elow, double etop) {
        int n = ch0.size();
        for (int i=0; i < n; i++) {
            double e0 = ch0[i]->energy * m[0] + q[0];
            double e1 = ch1[i]->energy * m[1] + q[1];
            if ( e0 < etop && e0 > elow && e1 < etop && e0 > elow ) enabled->at(i) = true;
            else enabled->at(i) = false;
        }
    }

    statistic FWHM(TH1F* h, double e_min = 5, double e_max = 20, bool draw = false){
        float centroid = 0., width_sx = 0., width_dx = 0.;
        float max_bin_content = 0.;
        int max_i = 0;
        for(unsigned int i = 1; i < h->GetXaxis()->GetNbins(); i++){
            if(h->GetXaxis()->GetBinCenter(i) < e_min)
                continue;
            if(h->GetXaxis()->GetBinCenter(i) > e_max)
                break;
            if(h->GetBinContent(i) > max_bin_content){
                max_bin_content = h->GetBinContent(i);
                centroid = h->GetXaxis()->GetBinCenter(i);
                max_i = i;
            }            
        }
        
        for(unsigned int i = 0; i < max_i; i++){
            if(h->GetXaxis()->GetBinCenter(i) < e_min)
                continue;
            if(h->GetBinContent(i) > max_bin_content*0.5){
                width_sx = centroid - h->GetXaxis()->GetBinCenter(i);
                break;
            }
        }
        
        for(unsigned int i = max_i; i < h->GetXaxis()->GetNbins(); i++){
            if(h->GetXaxis()->GetBinCenter(i) > e_max)
                break;
            if(h->GetBinContent(i) < max_bin_content*0.5){
                width_dx = h->GetXaxis()->GetBinCenter(i) - centroid;
                break;
            }
        }
        
        if(draw){
            h->Draw();
            TLine* l = new TLine(centroid - width_sx,0.5*max_bin_content,centroid + width_dx,0.5*max_bin_content);
            l->Draw("SAME");
        }
        statistic FWHM;
        FWHM.value = width_sx + width_dx;
        FWHM.error = h->GetXaxis()->GetBinWidth(5) / 2.45;
        return FWHM;
    }

    vector<Event*> ch0, ch1;
    vector<bool>* enabled;

    double m[2] = { 0.0825, 0.0668 };
    double q[2] = { -21.07, -3.86  };
};


void analysis(int id = -1) {
    if (id == -1) { analysis(0); analysis(1); return; }
	const char* outfilename[] = {"CFTD_simulations_1_2D.txt","CFTD_simulations_2_2D.txt"};
	const char* sourcename[] = {"Digital_CFTD.root", "Digital_CFTD_2.root"};
	vector<double> fracs{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
	vector<int> delays{1 , 2 , 3 , 4 , 5 ,  6 , 7 , 8 , 9 , 10};
	int i = 0;
	ofstream out(outfilename[id]);
	out << "Frac\tDelay\tMean\tMean_sigma\tFWHM\tFWHM_sigma\tKurtosis\tKurtosis_sigma"<<endl;

    statistic mean, FWHM, kurtosis;
    SimplyFastZeroFinder* s = new SimplyFastZeroFinder(sourcename[id]);
	for ( double f : fracs)
		for ( int d : delays) {
			i++;
			s->getProprieties ( s->getWidths(f,d), mean, FWHM, kurtosis, true );
            cout << i << " of " << fracs.size()*delays.size() << '\t' << f << '\t' << d << '\t' << mean.value  << '\t' << mean.error << '\t' << FWHM.value << '\t' << FWHM.error << '\t' << kurtosis.value << '\t' << kurtosis.error << endl;
            out << f << '\t' << d << '\t' << mean.value  << '\t' << mean.error << '\t' << FWHM.value << '\t' << FWHM.error << '\t' << kurtosis.value << '\t' << kurtosis.error << endl;
		}

    delete s;
}

void analysis_energythresh(int id = -1) {
    double fraction = 0.2;
    int delay = 5.0;

	if (id == -1) { analysis_energythresh(0); analysis_energythresh(1); return; }
	const char* outfilename[] = {"CFTD_energythresh_1_2D.txt","CFTD_energythresh_2_2D.txt"};
	const char* sourcename[] = {"Digital_CFTD.root", "Digital_CFTD_2.root"};
	vector<double> energy_low{ 50, 100, 150, 200, 250, 300, 350,  50, 100, 150, 200, 250 };
	vector<double> energy_top{600, 600, 600, 600, 600, 600, 600, 150, 250, 350, 450, 550 };
	ofstream out(outfilename[id]);
	out << "ELow\tETop\tMean\tMean_sigma\tFWHM\tFWHM_sigma\tKurtosis\tKurtosis_sigma"<<endl;

    statistic mean, FWHM, kurtosis;
    SimplyFastZeroFinder* s = new SimplyFastZeroFinder(sourcename[id]);

	for( int i=0; i<energy_low.size(); i++) {

        s->energyFilter( energy_low[i], energy_top[i] );
        s->getProprieties ( s->getWidths(fraction,delay), mean, FWHM, kurtosis, true );

        cout << setw(2) << i << " of " << energy_low.size() << '\t' << energy_low[i] << '\t' << energy_top[i] << '\t' << mean.value  << '\t' << mean.error << '\t' << FWHM.value << '\t' << FWHM.error << '\t' << kurtosis.value << '\t' << kurtosis.error << endl;
        out << energy_low[i] << '\t' << energy_top[i] << '\t' << mean.value  << '\t' << mean.error << '\t' << FWHM.value << '\t' << FWHM.error << '\t' << kurtosis.value << '\t' << kurtosis.error << endl;

	}
}