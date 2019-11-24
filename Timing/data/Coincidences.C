// This macro rearranges data from the VeRDI acquisition into a better format: no more need to look for coincidences

#ifndef Coincidences_h
#define Coincidences_h

#include <iostream>
#include <vector>

using namespace std;

struct slimport_data_t {
	ULong64_t	timetag; //time stamp
	UInt_t		baseline;
	UShort_t	qshort; //integration with shorter time
	UShort_t	qlong; //integration with longer time
	UShort_t	pur;
	UShort_t	samples[4096];
};


struct enhanced_data_t {
    ULong64_t               timetag; //time stamp
    vector<unsigned int>*   qlongs; //vector of qlongs
};

// vector<unsigned short>* a = new vector<unsigned short>(4);

class Coincidences{
    public:
        
        //Constructor
        Coincidences(int channnel_number){
            n_chan = channnel_number;
            original_data_file = "";
            hist = vector<TH1F*>(n_chan);
            alpha = vector<double>(n_chan,-1.0);
            beta = vector<double>(n_chan,-1.0);
            data = new vector<enhanced_data_t*>;
        }
        
        Coincidences(const char* filename, unsigned int channnel_number){
            n_chan = channnel_number;
            original_data_file = filename;
            hist = vector<TH1F*>(n_chan);
            alpha = vector<double>(n_chan,-1.0);
            beta = vector<double>(n_chan,-1.0);
            
            TFile* infile = new TFile(original_data_file);
            TTree* intree = (TTree*)infile->Get("acq_tree_0");
            slimport_data_t indata[n_chan];
            TBranch *inbranch[n_chan];
            int curr_index[n_chan];
            
            data = new vector<enhanced_data_t*>;
            
            for(int chan=0; chan<n_chan; chan++){
                inbranch[chan] = (TBranch*)intree->GetBranch(Form("acq_ch%d",chan));
                inbranch[chan]->SetAddress(&indata[chan].timetag);
    			hist [chan] = new TH1F(Form("hist_ch%d",chan),"",8192,0,65536);
                curr_index[chan]=0;		
            }
            
            data->reserve(inbranch[0]->GetEntries());
            bool exit = false;
            int discarded_events = 0;
            
//             curr_index[0]<inbranch[0]->GetEntries() && curr_index[1]<inbranch[1]->GetEntries() && curr_index[2]<inbranch[2]->GetEntries()
            while(true){
                
                for(int c=0;c<n_chan;c++){
                    if(!(curr_index[c]<inbranch[c]->GetEntries())){
                        exit = true;
                        break;
                    }
                    inbranch[c]->GetEntry(curr_index[c]);
                }
                if(exit)
                    break;

                bool found_coinc=true;

                for(int c=1;c<n_chan;c++)
                    found_coinc &= (indata[0].timetag == indata[c].timetag);

                for(int c=1;c<n_chan;c++){
                    if(indata[0].timetag < indata[c].timetag){
                        curr_index[0]++;
                        break;
                    }
                    if(indata[0].timetag > indata[c].timetag)
                        curr_index[c]++;
                }
                
                if(!found_coinc){
                    discarded_events++;
                    continue;
                }
                
                enhanced_data_t* d = new enhanced_data_t;
                d->timetag = indata[0].timetag;
                d->qlongs = new vector<unsigned int>(n_chan,0);
                for(int c=0;c<n_chan;c++){
                    d->qlongs->at(c) = indata[c].qlong;
                    hist[c]->Fill(indata[c].qlong);
                }
                data->push_back(d);
                
                //update condition
                for(int c=0; c<n_chan; c++)
                    curr_index[c]++;
            }
            
            cout << data->size() << " coincidences, " << discarded_events << " discarded events" << endl;
        }
        
        Coincidences* Copy(){
            Coincidences* C = new Coincidences(n_chan);
            C->original_data_file = original_data_file;
            C->data = new vector<enhanced_data_t*>(*data);
            C->hist = vector<TH1F*>(hist);
            C->alpha = vector<double>(alpha);
            C->beta = vector<double>(beta);
            
            return C;
        }
        
        //Destructor
        ~Coincidences(){
            delete data;
        }
        
        
        
        Coincidences* Filter(vector<int>* chan, vector<double>* e_min, vector<double>* e_max){
            Coincidences* C = Copy();
            
            //reset C
            C->data = new vector<enhanced_data_t*>;
            C->data->reserve(data->size());
            for(int c=0; c<n_chan; c++)
                C->hist[c] = new TH1F(Form("hist_ch%d",c),"",8192,0,65536);
            
            for(enhanced_data_t* d : *data){
                bool valid = true;
                for(int v = 0; v<chan->size(); v++){
                    int c = chan->at(v);
                    if((d->qlongs->at(c) < e2ch(c,e_min->at(v))) || (d->qlongs->at(c) > e2ch(c,e_max->at(v)))){
                        valid = false;
                        break;
                    }
                }
                if(valid){
                    C->data->push_back(d);
                    for(int c=0; c<n_chan; c++)
                        C->hist[c]->Fill(d->qlongs->at(c));
                }
            }
            C->CalibrateHisto();
            
            return C;
        }
        
        
        
        TH1F* GetHistogram(int chan){
            if (chan < 0 || chan >= n_chan)
                return 0x0;
            return hist[chan];
        }
        
        
        double ch2e(int chan, double ch){
            if (chan < 0 || chan >= n_chan)
                return -1;
            if(alpha[chan]==-1.0){
                cout << "Detector not calibrated" << endl;
                return -1;
            }
            else
                return alpha[chan]*ch+beta[chan];
        }
        
        int e2ch(int chan, double e){
            if (chan < 0 || chan >= n_chan)
                return -1;
            if(alpha[chan]==-1.0){
                cout << "Detector not calibrated" << endl;
                return -1;
            }
            else
                return (int)(e-beta[chan])/alpha[chan];
        }
        
        void SetAlphaBeta(int chan, double a, double b){
            if (chan < 0 || chan >= n_chan){
                cout<<"Invalid histogram index"<<endl;
                return;
            }
            alpha[chan] = a;
            beta[chan] = b;
        }
        
        void CalibrateHisto(int chan = -1, double a = 0, double b = 0) {
            if (chan == -1)
                for(int c=1;c<n_chan;c++)
                    CalibrateHisto(c);
            
            if (chan < 0 || chan > n_chan)
                return;
            if (a == 0 ) { a = alpha[chan]; b = beta[chan]; }
            alpha[chan]=a;
            beta[chan]=b;
            TAxis *axis = hist[chan]->GetXaxis();
            axis->SetLimits(axis->GetXmin()*a+b, axis->GetXmax()*a+b);
            //hist[chan]->Scale(1/a); //to ensure keeping the proportions
        }
    
    
    private:
        unsigned int n_chan;
        const char* original_data_file;
        vector<enhanced_data_t*>* data;
        vector<TH1F*> hist;
        vector<double> alpha;
        vector<double> beta;
};

#endif
