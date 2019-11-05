Histograms* ThreePhotonsHistCal(bool filter){
	Histograms* hs = new Histograms(" /Users/andreagrossutti/Documents/GitHub/Physics_Laboratory-g19/Positronium/data/lavoltabuona/night_g19.root ");
	hs->CalibrationPeaksChannels(0, 6.58890e+03, 2.70877e+03);
	hs->CalibrationPeaksChannels(1, 1.17060e+04, 4.82225e+03);
	hs->CalibrationPeaksChannels(2, 1.25830e+04, 5.19850e+03);
	if(filter)
		hs->SpectraFiltering(1022, 150, 3);
	return hs;
}

Histograms* ThreePhotonsHistCal1(bool filter){
	Histograms* hs = new Histograms(" /Users/andreagrossutti/Documents/GitHub/Physics_Laboratory-g19/Positronium/data/lavoltabuona/DET123_calibr.root ");
	hs->CalibrationPeaksChannels(0, 6.58890e+03, 2.70877e+03);
	hs->CalibrationPeaksChannels(1, 1.17060e+04, 4.82225e+03);
	hs->CalibrationPeaksChannels(2, 1.25830e+04, 5.19850e+03);
	if(filter)
		hs->SpectraFiltering(1022, 150, 3);
	return hs;
}

Histograms* TwoPhotonsHistCal(bool filter){
	Histograms* hs = new Histograms(" /Users/andreagrossutti/Documents/GitHub/Physics_Laboratory-g19/Positronium/data/day3/2Photon_acq_4.root ");
	hs->CalibrationPeaksChannels(0, 4.32049e+03, 1.78084e+03);
	hs->CalibrationPeaksChannels(1, 1.27812e+04, 5.26642e+03);
	hs->CalibrationPeaksChannels(2, 4.97904e+03, 2.05268e+03);
	if(filter)
		hs->SpectraFiltering(1022, 300, 2);
	return hs;
}