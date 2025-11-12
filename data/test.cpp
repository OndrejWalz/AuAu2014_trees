void test() {
    TFile* infile = TFile::Open("test_embedding.jetTree.root");
    TTree* jetTree = (TTree*)infile->Get("R0.2/CENT_0_10/JetTree");
    TH1D * NUM = new TH1D("NUM","NUM",100,-40,60);
    TH1D * DEN = new TH1D("DEN","DEN",100,-40,60);

    jetTree->Draw("reco_pt_corr>>DEN","","goff");
    jetTree->Draw("reco_pt_corr>>NUM","reco_trigger_match","goff");

    NUM->Divide(NUM,DEN,1,1,"B");
    TCanvas* c1 = new TCanvas("c1","c1",800,600);
    NUM->SetTitle("Trigger Efficiency vs Reco Jet p_{T}^{corr};p_{T}^{corr} [GeV];Efficiency");
    NUM->Draw("E");
    c1->SaveAs("trigger_efficiency.pdf");
}   