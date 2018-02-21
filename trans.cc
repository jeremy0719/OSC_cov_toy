{
	TFile* f1 = new TFile("test.root","OPEN");
	TFile* f2 = new TFile("detector.root","RECREATE");
	TH1D* a = (TH1D*) f1->Get("ETrue101");
	TH1D* b = (TH1D*) f1->Get("ENull101");
	f2->cd();
	a->Write();
	b->Write();
	f2->Write();
	f2->Close();
}


