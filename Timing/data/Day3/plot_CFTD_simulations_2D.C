TGraph2DErrors* plot_CFTD_simulations_2D (int what, int id = 0, bool same = false) {

	// what: 0->Mean, 1->FWHM, 2->Kurtosis

	//if ( what < 0 || what > 2 ) return;
	char* filename[] = {"CFTD_simulations_1_2D.txt","CFTD_simulations_2_2D.txt"};
	ifstream in(filename[id]);
    TGraph2DErrors* tge = new TGraph2DErrors();

	double a,b,c,d,e,f,g,h,point,error;
	in.ignore(10000,'\n');
	while (in>>a>>b>>c>>d>>e>>f>>g>>h) {
	
		if(what == 0) { point = c; error = d; }
		else if(what == 1) { point = e; error = f; }
		else if(what == 2) { point = g; error = h; }

		tge->SetPoint( tge->GetN(), a, b, point);
		tge->SetPointError ( tge->GetN() - 1, 0, 0, error);
	}
    return tge;
}
