void plot() {
	TGraph* tg1 = new TGraph("output_1.txt");
	TGraph* tg2 = new TGraph("output_2.txt");
	tg1->Draw();
	tg2->Draw("SAME");
}
