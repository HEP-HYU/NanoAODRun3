void make_DATApileup()
{
    // 2016 pre
    auto outfile16pre = new TFile("PileupDATA_UL16pre.root","RECREATE");
    auto pufile_16pre = new TFile("Pileup_Histograms/PileupHistogram-goldenJSON-13tev-2016-preVFP-69200ub-99bins.root","READ");
    auto h16pre_norm = pufile_16pre->Get("pileup");
    outfile16pre->cd();
    h16pre_norm->Write();
    pufile_16pre->Close();

    auto pufile_16pre_plus = new TFile("Pileup_Histograms/PileupHistogram-goldenJSON-13tev-2016-preVFP-72400ub-99bins.root","READ");
    auto h16pre_plus = pufile_16pre_plus->Get("pileup");
    auto h16pre_plus_clone = h16pre_plus->Clone("pileup_plus");
    outfile16pre->cd();
    h16pre_plus_clone->Write();
    pufile_16pre_plus->Close();
    
    auto pufile_16pre_minus = new TFile("Pileup_Histograms/PileupHistogram-goldenJSON-13tev-2016-preVFP-66000ub-99bins.root","READ");
    auto h16pre_minus = pufile_16pre_minus->Get("pileup");
    auto h16pre_minus_clone = h16pre_minus->Clone("pileup_minus");
    outfile16pre->cd();
    h16pre_minus_clone->Write();
    pufile_16pre_minus->Close();
    outfile16pre->Close();
    

    // 2016 post
    auto outfile16post = new TFile("PileupDATA_UL16post.root","RECREATE");
    auto pufile_16post = new TFile("Pileup_Histograms/PileupHistogram-goldenJSON-13tev-2016-postVFP-69200ub-99bins.root","READ");
    auto h16post_norm = pufile_16post->Get("pileup");
    outfile16post->cd();
    h16post_norm->Write();
    pufile_16post->Close();

    auto pufile_16post_plus = new TFile("Pileup_Histograms/PileupHistogram-goldenJSON-13tev-2016-postVFP-72400ub-99bins.root","READ");
    auto h16post_plus = pufile_16post_plus->Get("pileup");
    auto h16post_plus_clone = h16post_plus->Clone("pileup_plus");
    outfile16post->cd();
    h16post_plus_clone->Write();
    pufile_16post_plus->Close();
    
    auto pufile_16post_minus = new TFile("Pileup_Histograms/PileupHistogram-goldenJSON-13tev-2016-postVFP-66000ub-99bins.root","READ");
    auto h16post_minus = pufile_16post_minus->Get("pileup");
    auto h16post_minus_clone = h16post_minus->Clone("pileup_minus");
    outfile16post->cd();
    h16post_minus_clone->Write();
    pufile_16post_minus->Close();
    outfile16post->Close();
    

    // 2017
    auto outfile17 = new TFile("PileupDATA_UL17.root","RECREATE");
    auto pufile_17 = new TFile("Pileup_Histograms/PileupHistogram-goldenJSON-13tev-2017-69200ub-99bins.root","READ");
    auto h17_norm = pufile_17->Get("pileup");
    outfile17->cd();
    h17_norm->Write();
    pufile_17->Close();

    auto pufile_17_plus = new TFile("Pileup_Histograms/PileupHistogram-goldenJSON-13tev-2017-72400ub-99bins.root","READ");
    auto h17_plus = pufile_17_plus->Get("pileup");
    auto h17_plus_clone = h17_plus->Clone("pileup_plus");
    outfile17->cd();
    h17_plus_clone->Write();
    pufile_17_plus->Close();
    
    auto pufile_17_minus = new TFile("Pileup_Histograms/PileupHistogram-goldenJSON-13tev-2017-66000ub-99bins.root","READ");
    auto h17_minus = pufile_17_minus->Get("pileup");
    auto h17_minus_clone = h17_minus->Clone("pileup_minus");
    outfile17->cd();
    h17_minus_clone->Write();
    pufile_17_minus->Close();
    outfile17->Close();
    

    // 2018
    auto outfile18 = new TFile("PileupDATA_UL18.root","RECREATE");
    auto pufile_18 = new TFile("Pileup_Histograms/PileupHistogram-goldenJSON-13tev-2018-69200ub-99bins.root","READ");
    auto h18_norm = pufile_18->Get("pileup");
    outfile18->cd();
    h18_norm->Write();
    pufile_18->Close();

    auto pufile_18_plus = new TFile("Pileup_Histograms/PileupHistogram-goldenJSON-13tev-2018-72400ub-99bins.root","READ");
    auto h18_plus = pufile_18_plus->Get("pileup");
    auto h18_plus_clone = h18_plus->Clone("pileup_plus");
    outfile18->cd();
    h18_plus_clone->Write();
    pufile_18_plus->Close(); 
    
    auto pufile_18_minus = new TFile("Pileup_Histograms/PileupHistogram-goldenJSON-13tev-2018-66000ub-99bins.root","READ");
    auto h18_minus = pufile_18_minus->Get("pileup");
    auto h18_minus_clone = h18_minus->Clone("pileup_minus");
    outfile18->cd();
    h18_minus_clone->Write();
    pufile_18_minus->Close();
    outfile18->Close();
}
