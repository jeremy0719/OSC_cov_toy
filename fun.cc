// need to slightly vary the parameters to account for PROSPECT detector.

double Energy_Resolution_Function (double E, double a = 6.453e-08, double b = 4.683e-02, double c = 7.848e-03 ) 
{
	double result;
	result = E * TMath::Sqrt( a*a + b*b/E + c*c/TMath::Power(E,2) );
	// =====================================================================
	// cout << result << " << energy resolution " << endl;
	// =====================================================================
	return result;
}

double Energy_Scale_Function (double E , double alpha = -0.149 , double beta = 0.9544 , double tau = 0.6079 )
{	
	double result;
	result = E * beta * ( 1+ alpha * TMath::Exp( - (1.0/tau) * E ) );
	return result;
}

int Gaussian_Distribution_Generator ( double E , double sigma , TH1D* input_hist , double number_of_event )
{
	// by default \pm 5 \sigma 
	// Need to discuss when 5 sigma is too large or too small to prevent binning issues


	int bin_initial = 0 ;
	int bin_final = 0 ;
	int number_of_sigma = 5 ;
	// =====================================================================
	// cout << E << "<< E " << endl;
	// cout << sigma << " << sigma " << endl;
	// =====================================================================

	if (E - number_of_sigma * sigma  >= input_hist->GetBinLowEdge(1) )
	{
		bin_initial = input_hist->FindBin( E -  number_of_sigma * sigma );
	}
	else {
		bin_initial = 1 ;
	}
	if ( E +  number_of_sigma * sigma <= input_hist->GetBinLowEdge( input_hist->GetNbinsX() + 1 ) )
	{
		bin_final = input_hist->FindBin( E +  number_of_sigma * sigma );
	}
	else {
		bin_final = input_hist->GetNbinsX() ;
	}


	double scaling_factor = 0 ;
	for(int i = bin_initial; i <= bin_final ; i++ )
	{
		scaling_factor += number_of_event * TMath::Gaus( input_hist->GetBinCenter(i) , E , sigma ) ; 
		// cout << number_of_event * TMath::Gaus( input_hist->GetBinCenter(i) , E , sigma ) << " << number of gaussian  " << endl;
	}
	// =====================================================================
	// cout << bin_initial << " to " << bin_final << endl ;
	// cout << scaling_factor << " << scaling factor " << endl;
	// =====================================================================
	// double temp = 0 ;
	for(int i = bin_initial; i <= bin_final ; i++ )
	{
		input_hist->Fill( input_hist->GetBinCenter( i ) ,  number_of_event * TMath::Gaus( input_hist->GetBinCenter(i) , E , sigma ) * number_of_event / scaling_factor );
		// temp += number_of_event * TMath::Gaus( input_hist->GetBinCenter(i) , E , sigma ) * number_of_event / scaling_factor ;
	}
	// cout << number_of_event << " << number of events " << endl;
	// cout << temp << " << number of event after gaussian." << endl;

	// if ( bin_final - bin_initial <= 2 ) cout << "bin spread less than 2 " << endl;

	return 0;
}
int fun()
{
	TFile* file = new TFile("detector.root","OPEN");
	// Uncertainty in the histogram; 
	TH1D* ETrue = (TH1D*) file->Get("ETrue101");
	TH1D* ENull = (TH1D*) file->Get("ENull101");

	
	std::cout << ENull->Integral() << " << number of ENull events " << std::endl;
	std::cout << ETrue->Integral() << " << number of ETrue events " << std::endl;
	
	
	// Define the lower and upper limit of the spectrum range
	double ENull_xlow = 0 ;
	double ENull_xup  =	10 ;

	double ETrue_xlow = ETrue->GetBinLowEdge(1);
	double ETrue_xup  =	ETrue->GetBinLowEdge( ETrue->GetNbinsX() + 1 );

	double BWidth = 0.100000;

	int ENull_NBin = ( ENull_xup - ENull_xlow ) / BWidth ;
	TFile* file_test = new TFile("test_result.root","RECREATE");



	// phase sequence details
	// phase 1: energy shift to go from neutrino energy to position energy depostion
	// phase 2: apply energy scale  
	// phase 3: apply energy resolution
	// phase 4: apply gamma escape energy loss 
	// ==========================================================================================================================
	TH1D* ENull_phase1 = new TH1D("phase_1" , "ETrue after energy shift ",  ENull_NBin  , ENull_xlow , ENull_xup  );
	TH1D* ENull_phase2 = new TH1D("phase_2" , "ETrue after energy scale ",  ENull_NBin  , ENull_xlow , ENull_xup  );
	TH1D* ENull_phase3 = new TH1D("phase_3" , "ETrue after energy resolution ",  ENull_NBin  , ENull_xlow , ENull_xup  );
	TH1D* ENull_phase4 = new TH1D("phase_4" , "ETrue after gamma escape ",  ENull_NBin  , ENull_xlow , ENull_xup  );
	// ==========================================================================================================================

	const double neutrino_conversion_energy = 0.800 ; 
	for( int i = 1; i <= ETrue->GetNbinsX(); i++ ) 
	{	
		// phase 1 execution 
		double shifted_position_energy_deposition = ETrue->GetBinCenter(i) - neutrino_conversion_energy ;
		ENull_phase1->Fill( shifted_position_energy_deposition , ETrue->GetBinContent(i) );
		// phase 2 execution with default enregy scale function 
		double scaled_detected_energy =  Energy_Scale_Function ( shifted_position_energy_deposition ); 
		// std::cout << scaled_detected_energy << "    " << ETrue->GetBinContent(i) << std::endl;
		ENull_phase2->Fill( scaled_detected_energy , ETrue->GetBinContent(i) ); 
		// cout << ETrue->GetBinContent(i) << " << ETrue number of event " << endl ;
		// phase 3 execution with default energy resolution function 
		Gaussian_Distribution_Generator( scaled_detected_energy, Energy_Resolution_Function( scaled_detected_energy ), ENull_phase3 , ETrue->GetBinContent(i) );
	}
	// std::cout << ENull_phase3->Integral() << " << number of ENull_phase3 events " << std::endl;
	// phase 4 execution with average energy loss ~ 0.4 Mev 
	for( int i = 1; i <= ENull_phase3->GetNbinsX(); i++ )
	{
		ENull_phase4->Fill( ENull_phase3->GetBinCenter(i) - 0.4 , ENull_phase3->GetBinContent(i) );
	}
	// test section
	int  NBin = ENull->GetNbinsX() ;
	double xlow = ENull->GetBinLowEdge(1) ;
	double xup = ENull->GetBinLowEdge( ENull->GetNbinsX() + 1 );
	auto ENull_new = new TH1D( "ENullNew" , " Rebinning ETrue ",  NBin  , xlow ,  xup  );
	
	for (int i  = 1 ; i <= ENull_phase4->GetNbinsX(); i++ )
	{	
		ENull_new->Fill( ENull_phase4->GetBinCenter(i) , ENull_phase4->GetBinContent(i) );
	}
	
	
	// ENull->Draw("HIST");




	std::cout << ENull_new->Integral() << " << number of ENull_new events " << std::endl;

	file_test->cd();
	ENull_new->Write();
	ENull->Write();

	// to decouple it from the open file directory !!
	ENull->SetDirectory(0);
	ENull_new->SetDirectory(0);

	file_test->Write();
	// file_test->Flush();
	file_test->Close();

	TCanvas* c1 = new TCanvas("c1","title");
	ENull->SetLineColor(1);
	ENull->Draw("HIST");
	ENull_new->SetLineColor(2);
	ENull_new->Scale(1.0/ENull_new->Integral()*ENull->Integral());
	ENull_new->Draw("HIST SAME");

	c1->Print("ENullvsENull_new.pdf");
	return 0;
}

