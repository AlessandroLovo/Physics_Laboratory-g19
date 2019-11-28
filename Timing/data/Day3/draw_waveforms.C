void draw_waveforms() {
	int ev = 18;
	
	auto wf = new Waveform("Digital_CFTD.root", 0);

	auto wfcf = wf->CFTD(ev);
	for( int i=0; i < 360; i++)
		wfcf->SetPoint(i,0,0);
	
	auto wfwb = wf->GetGraphWithoutBaseline(ev);
	wfwb->GetYaxis()->SetRangeUser(-300,300);
	wfwb->GetXaxis()->SetRange(90,140);
	
	wfwb->GetXaxis()->SetTitle("Time [ns]");
	wfwb->GetYaxis()->SetTitle("Voltage [a.u.]");
	wfwb->SetTitle("Simulated CFTD");

	double x0 = wf->ZeroCrossing(18);
	TLine* tl = new TLine(x0,-25,x0,25);
	tl->SetLineColor(2);
	
	wfwb->SetLineColor(5);
	wfcf->SetLineColor(7);
	wfwb->SetFillColor(1);

	wfwb->Draw("AL");
	wfcf->Draw("SAME l");
	tl->Draw();
}
