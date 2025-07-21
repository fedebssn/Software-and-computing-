 double mygausexp(double* x, double* p)
    {
    double out;
    if (x[0] <= p[1] - p[3]) {
        // Lato sinistro con decadimento esponenziale
        out = p[0] * TMath::Gaus(p[1] - p[3], p[1], p[2]) * TMath::Exp((x[0] - p[1] + p[3]) * p[3] / (p[2] * p[2]));
    } else if (x[0] <= p[1] + p[3]) {
        // Parte centrale: gaussiana pura
        out = p[0] * TMath::Gaus(x[0], p[1], p[2]);
    } else {
        // Lato destro con decadimento esponenziale
        out = p[0] * TMath::Gaus(p[1] + p[3], p[1], p[2]) * TMath::Exp(-p[3] * (x[0] - p[1] - p[3]) / (p[2] * p[2]));
    }
    return out;
    }
 void readtreeMC(){
    
    // Open the ROOT file containing the TTree
    TFile *file = TFile::Open("mytreeMC1.root");
    if (!file || file->IsZombie()) {
        // If file could not be opened, print error and return
        std::cerr << "Error opening file: " << "mytreeMC1.root" << std::endl;
        return;
    }

    // Load the TTree named "tracks" from the file
    TTree *tree = (TTree*)file->Get("tracks");
    if (!tree) {
        // If the TTree doesn't exist, print error and return
        std::cerr << "Tree 'tracks' not found in file." << std::endl;
        file->Close();
        return;
    }

    // Load MC-generated histograms for different particle species and charges
    TH1D* hGenPiPlus = (TH1D*)file->Get("hPiPlusGenVsPt");     
    TH1D* hGenPrPlus = (TH1D*)file->Get("hPrPlusGenVsPt");
    TH1D* hGenKaPlus = (TH1D*)file->Get("hKaPlusGenVsPt");

    TH1D* hGenPiMinus = (TH1D*)file->Get("hPiMinusGenVsPt;1");     
    TH1D* hGenPrMinus = (TH1D*)file->Get("hPrMinusGenVsPt;1");
    TH1D* hGenKaMinus = (TH1D*)file->Get("hKaMinusGenVsPt;1");

    // Variables to store data from each branch in the TTree
    int pdgcode = 0, isphysicalprimary = -1;
    float pt = 0., eta = 0., beta = 0., p = 0., px = 0., py = 0., pz = 0., sign = 0., tpc = 0., tof = 0.,
          tevent = 0., tpi = 0., tpr = 0., tka = 0., hasTOF = 0.,
          tpcnsigmapi = 0., tpcnsigmapr = 0., tpcnsigmaka = 0.,
          tofnsigmapi = 0., tofnsigmapr = 0., tofnsigmaka = 0.,
          tofnsigmapim = 0., tofnsigmaprm = 0., tofnsigmakam = 0.;

    // Connect variables to branches in the TTree
    tree->SetBranchAddress("pt", &pt);
    tree->SetBranchAddress("eta", &eta);
    tree->SetBranchAddress("beta", &beta);
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("sign", &sign);
    tree->SetBranchAddress("tpcSignal", &tpc);
    tree->SetBranchAddress("tofSignal", &tof);
    tree->SetBranchAddress("tofEvTime", &tevent);
    tree->SetBranchAddress("tofExpTimePi", &tpi);
    tree->SetBranchAddress("tofExpTimePr", &tpr);
    tree->SetBranchAddress("tofExpTimeKa", &tka);
    tree->SetBranchAddress("hasTOF", &hasTOF);
    tree->SetBranchAddress("tpcNSigmaPi", &tpcnsigmapi);
    tree->SetBranchAddress("tpcNSigmaPr", &tpcnsigmapr);
    tree->SetBranchAddress("tpcNSigmaKa", &tpcnsigmaka);
    tree->SetBranchAddress("tofNSigmaPi", &tofnsigmapi);
    tree->SetBranchAddress("tofNSigmaPr", &tofnsigmapr);
    tree->SetBranchAddress("tofNSigmaKa", &tofnsigmaka);
    tree->SetBranchAddress("pdgcode", &pdgcode);
    tree->SetBranchAddress("isphysicalprimary", &isphysicalprimary);

    // Get the number of entries (rows) in the TTree
    Long64_t nEntries = tree->GetEntries();
    std::cout << "Reading " << nEntries << " entries from tree: tracks" << std::endl;


        // Define 1D histograms to store basic track quantities
    TH1D* hPt = new TH1D("hPt", "titolo; assex; assey", 250, 0, 25);  // Transverse momentum
    TH1D* hEta = new TH1D("hEta", "eta distribution; #eta; entries", 200, -1, +1);  // Pseudorapidity
    TH1D* hPx = new TH1D("hPx", "px; px; entries", 250, 0, 25);  // Momentum component px
    TH1D* hPy = new TH1D("hPy", "py; py; entries", 250, 0, 25);  // Momentum component py
    TH1D* hPz = new TH1D("hPz", "pz; pz; entries", 250, 0, 25);  // Momentum component pz
    TH1D* hSign = new TH1D("hSign", "sign; sign; sign", 10, -1, 1);  // Charge sign
    TH1D* hTpc = new TH1D("hTpc", "tpc; tpc; tpc", 200, 0, 150);  // TPC signal (dE/dx)
    TH1D* hP = new TH1D("hP", "p; p; entries", 250, 0, 25);  // Total momentum
    TH1D* hBeta = new TH1D("hBeta", "#beta; beta; entries", 250, -60, 60);  // Beta (v/c)
    TH1D* hTof = new TH1D("hTof", "tof; tof; tof", 200, 0, 20);  // TOF signal
    TH1D* hTevent = new TH1D("hTevent", "tevent; tevent; tevent", 500, 0, 2000);  // TOF event time
    TH1D* hTpi = new TH1D("hTpi", "tpi; tpi; tpi", 500, 0, 20000);  // Expected TOF for pions
    TH1D* hTpr = new TH1D("hTpr", "tpr; tpr; tpr", 500, 0, 60000);  // Expected TOF for protons
    TH1D* hTka = new TH1D("hTka", "tka; tka; tka", 500, 0, 30000);  // Expected TOF for kaons

    // Histograms for TPC PID response (number of sigmas)
    TH1D* hTpcsigpi = new TH1D("htpcnsigmapi", "tpcnsigmapi; tpcnsigmapi; tpcnsigmapi", 500, 0, 10);
    TH1D* hTpcsigpr = new TH1D("htpcnsigmapr", "tpcnsigmapr; tpcnsigmapr; tpcnsigmapr", 500, 0, 10);
    TH1D* hTpcsigka = new TH1D("htpcnsigmaka", "tpcnsigmaka; tpcnsigmaka; tpcnsigmaka", 500, 0, 10);

    
        // Define custom pT bins for different particle types
    float ptbins[] = {0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};
    const int nbins = sizeof(ptbins) / sizeof(float) - 1;

    float ptbinska[] = {0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4};
    const int nbinska = sizeof(ptbinska) / sizeof(float) - 1;

    
    // Define histograms to store efficiencies and denominator counts for each particle species and charge

    // Positive pions
    TH1D* hEffPiPlus = new TH1D("hEffPiPlus", "; p_{T}(GeV/c); #varepsilon", nbins, ptbins);
    TH1D* hDenPiPlus = new TH1D("hDenPiPlus", "; p_{T}(GeV/c); #varepsilon", nbins, ptbins);

    // Positive protons
    TH1D* hEffPrPlus = new TH1D("hEffPrPlus", "; p_{T}(GeV/c); #varepsilon", nbins, ptbins);
    TH1D* hDenPrPlus = new TH1D("hDenPrPlus", "; p_{T}(GeV/c); #varepsilon", nbins, ptbins);

    // Positive kaons
    TH1D* hEffKaPlus = new TH1D("hEffKaPlus", "; p_{T}(GeV/c); #varepsilon", nbinska, ptbinska);
    TH1D* hDenKaPlus = new TH1D("hDenKaPlus", "; p_{T}(GeV/c); #varepsilon", nbinska, ptbinska);

    // Negative particles (same as above but with minus charge)
    TH1D* hEffPiMinus = new TH1D("hEffPiMinus", "; p_{T}(GeV/c); #varepsilon", nbins, ptbins);
    TH1D* hDenPiMinus = new TH1D("hDenPiMinus", "; p_{T}(GeV/c); #varepsilon", nbins, ptbins);
    TH1D* hEffPrMinus = new TH1D("hEffPrMinus", "; p_{T}(GeV/c); #varepsilon", nbins, ptbins);
    TH1D* hDenPrMinus = new TH1D("hDenPrMinus", "; p_{T}(GeV/c); #varepsilon", nbins, ptbins);
    TH1D* hEffKaMinus = new TH1D("hEffKaMinus", "; p_{T}(GeV/c); #varepsilon", nbinska, ptbinska);
    TH1D* hDenKaMinus = new TH1D("hDenKaMinus", "; p_{T}(GeV/c); #varepsilon", nbinska, ptbinska);

    
    // Histograms for raw yields per pT bin for each particle species and charge
    TH1D* hrawpi = new TH1D("hrawpi", "; p_{T}(GeV/c); #frac{1}{N_{ev}} #frac{d^{2}N}{dP_{t}dy}(GeV/c)^{-1}", nbins, ptbins);
    TH1D* hrawpr = new TH1D("hrawpr", "; p_{T}(GeV/c); #frac{1}{N_{ev}} #frac{d^{2}N}{dP_{t}dy}(GeV/c)^{-1}", nbins, ptbins);
    TH1D* hrawka = new TH1D("hrawka", "; p_{T}(GeV/c); #frac{1}{N_{ev}} #frac{d^{2}N}{dP_{t}dy}(GeV/c)^{-1}", nbinska, ptbinska);

    TH1D* hrawpim = new TH1D("hrawpim", "; p_{T}(GeV/c); #frac{1}{N_{ev}} #frac{d^{2}N}{dP_{t}dy}(GeV/c)^{-1}", nbins, ptbins);
    TH1D* hrawprm = new TH1D("hrawprm", "; p_{T}(GeV/c); #frac{1}{N_{ev}} #frac{d^{2}N}{dP_{t}dy}(GeV/c)^{-1}", nbins, ptbins);
    TH1D* hrawkam = new TH1D("hrawkam", "; p_{T}(GeV/c); #frac{1}{N_{ev}} #frac{d^{2}N}{dP_{t}dy}(GeV/c)^{-1}", nbinska, ptbinska);


    // 2D histogram of energy loss vs momentum/charge (dE/dx vs p/Z)
    TH2D* h2energyloss = new TH2D("h2energyloss", " ; p/Z (GeV/c); dE/dx (a.u.)", 400, 0, 10, 300, 0, 250);

    // 2D histogram of beta (v/c) vs total momentum
    TH2D* h2beta = new TH2D("h2beta", " ; p (GeV/c); Beta", 400, 0, 10, 160, 0.3, 1.1);

    // 2D histograms of TOF residuals (measured TOF - expected TOF)
    TH2D* h2pi = new TH2D("h2pi", "#DeltaT ; p (Gev/c); Tpi (ns)", 400, 0, 10, 500, -20000, 20000);  // pions
    TH2D* h2pr = new TH2D("h2pr", "#DeltaT ; p (Gev/c); Tpr (ns)", 400, 0, 10, 500, -20000, 20000);  // protons
    TH2D* h2ka = new TH2D("h2ka", "#DeltaT ; p (Gev/c); Tka (ns)", 400, 0, 10, 500, -20000, 20000);  // kaons

    // 2D histograms for TPC PID in number of sigmas (per particle type)
    TH2D* h2nsigpi = new TH2D("h2nsigpi", "n #sigma ; p_{t} (GeV/c); n #sigma_{pi}", 200, 0, 10, 160, -8, 8);
    TH2D* h2nsigpr = new TH2D("h2nsigpr", "n #sigma ; p_{t} (GeV/c); n #sigma_{pr}", 200, 0, 10, 160, -8, 8);
    TH2D* h2nsigka = new TH2D("h2nsigka", "n #sigma ; p_{t} (GeV/c); n #sigma_{ka}", 200, 0, 10, 160, -8, 8);

    // TOF PID response histograms for positively charged particles
    TH2D* h2nsigtofpi = new TH2D("h2nsigtofpi", "n #sigma tof ; p_{t} (GeV/c); n #sigma_{pi}", 200, 0, 10, 160, -8, 8);
    TH2D* h2nsigtofpr = new TH2D("h2nsigtofpr", "n #sigma tof ; p_{t} (GeV/c); n #sigma_{pr}", 200, 0, 10, 160, -8, 8);
    TH2D* h2nsigtofka = new TH2D("h2nsigtofka", "n #sigma tof ; p_{t} (GeV/c); n #sigma_{ka}", 200, 0, 10, 160, -8, 8);

    // TOF PID response histograms for negatively charged particles
    TH2D* h2nsigtofpim = new TH2D("h2nsigtofpim", "n #sigma tof ; p_{t} (GeV/c); n #sigma_{pi-}", 200, 0, 10, 160, -8, 8);
    TH2D* h2nsigtofprm = new TH2D("h2nsigtofprm", "n #sigma tof ; p_{t} (GeV/c); n #sigma_{apr}", 200, 0, 10, 160, -8, 8);
    TH2D* h2nsigtofkam = new TH2D("h2nsigtofkam", "n #sigma tof ; p_{t} (GeV/c); n #sigma_{ka-}", 200, 0, 10, 160, -8, 8);


        for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);  // Load data for event i from tree

        // Fill histograms with basic track info
        hPt->Fill(pt);
        hEta->Fill(eta);
        hPx->Fill(px);
        hPy->Fill(py);
        hPz->Fill(pz);
        hSign->Fill(sign);
        hTpc->Fill(tpc);

        // Compute total momentum from px, py, pz
        p = TMath::Sqrt(px*px + py*py + pz*pz);
        hP->Fill(p);             // Fill momentum histogram
        hTpc->Fill(tpc);         // Fill again TPC signal
        h2energyloss->Fill(p/sign, tpc);  // Plot dE/dx vs p/Z
        h2beta->Fill(p, beta);   // Plot beta vs momentum

        // If TOF information is available, fill corresponding histograms
        if (hasTOF == 1) {
            hTof->Fill(tof);          // Measured TOF signal
            hTevent->Fill(tevent);    // Event time
            hTpi->Fill(tpi);          // Expected TOF for pion
            hTpr->Fill(tpr);          // Expected TOF for proton
            hTka->Fill(tka);          // Expected TOF for kaon

            // TOF residuals (difference between measured and expected)
            h2pi->Fill(p, tof - tevent - tpi);
            h2pr->Fill(p, tof - tevent - tpr);
            h2ka->Fill(p, tof - tevent - tka);
        }

        // TPC PID responses vs pT
        h2nsigpi->Fill(pt, tpcnsigmapi);
        h2nsigpr->Fill(pt, tpcnsigmapr);
        h2nsigka->Fill(pt, tpcnsigmaka);

                // Skip secondary particles (keep only primary particles)
        if (isphysicalprimary != 1) {
            continue;
        }

        // Particle identification by PDG code
        if (pdgcode == 211) {  // Positive pion
            if (std::abs(tpcnsigmapi) < 3 && sign > 0) {
                h2nsigtofpi->Fill(pt, tofnsigmapi);  // TOF PID response for pi+
            }
        }
        if (pdgcode == 2212) {  // Positive proton
            if (std::abs(tpcnsigmapr) < 3 && sign > 0) {
                h2nsigtofpr->Fill(pt, tofnsigmapr);  // TOF PID response for p+
            }
        }
        if (pdgcode == 321) {  // Positive kaon
            if (std::abs(tpcnsigmaka) < 3 && sign > 0) {
                h2nsigtofka->Fill(pt, tofnsigmaka);  // TOF PID response for K+
            }
        }



                // Check again for physical primary 
        if (isphysicalprimary != 1) {
            continue;
        }

        if (pdgcode == -211) {  // Negative pion
            if (std::abs(tpcnsigmapi) < 3 && sign < 0) {
                h2nsigtofpim->Fill(pt, tofnsigmapi);  // TOF PID response for pi-
            }
        }
        if (pdgcode == -2212) {  // Negative proton
            if (std::abs(tpcnsigmapr) < 3 && sign < 0) {
                h2nsigtofprm->Fill(pt, tofnsigmapr);  // TOF PID response for p-
            }
        }
        if (pdgcode == -321) {  // Negative kaon
            if (std::abs(tpcnsigmaka) < 3 && sign < 0) {
                h2nsigtofkam->Fill(pt, tofnsigmaka);  // TOF PID response for K-
            }
        }


    }
  // Arrays to store histograms and functions for each pt bin
TH1D* hprojnsigtofpi[nbins];
TF1* fgauspi[nbins];
TF1* fgausclonepi[nbins];
TF1* fsignalpi[nbins];

// Loop over pt bins to project the 2D histogram onto Y-axis (tofNSigmaPi)
for (int i = 0; i < nbins; i++)
{
    // Project Y axis from h2nsigtofpi for the i-th pt bin
    hprojnsigtofpi[i] = h2nsigtofpi->ProjectionY(Form("hprojnsigtofpi%i", i), h2nsigtofpi->GetXaxis()->FindBin(ptbins[i]+1E-6), h2nsigtofpi->GetXaxis()->FindBin(ptbins[i+1]-1E-6));
    hprojnsigtofpi[i]->SetTitle(Form("%.1f<pt<%.1f", ptbins[i], ptbins[i+1]));  // Title of histogram with pt range

    // Define the custom Gaussian + exponential function for fitting
    fgauspi[i] = new TF1(Form("fb%i", i), mygausexp, -5, 5, 4);
    fgauspi[i]->SetParameter(0, hprojnsigtofpi[i]->GetMaximum());  // Initial guess for amplitude
    fgauspi[i]->SetParameter(1, 0);                                // Initial guess for mean
    fgauspi[i]->SetParameter(2, 1);                                // Initial guess for sigma
}

// Setup LaTeX for displaying fit parameters on canvas
TLatex *tex = new TLatex();
tex->SetNDC(); tex->SetTextFont(42); tex->SetTextSize(0.06);

// Create canvas for pi+ signal plots and divide into pads
TCanvas* csignalpi = new TCanvas ("Pi Signal","pi signal", 1500, 900);
csignalpi->Divide(4, 4);

for (int i = 0; i < nbins; i++)
{
    csignalpi->cd(i+1);                       // Move to next canvas pad
    hprojnsigtofpi[i]->Draw();               // Draw histogram
    hprojnsigtofpi[i]->Fit(fgauspi[i],"RW"); // Fit histogram with custom function

    // Draw fit parameters on the plot
    tex->DrawLatex(0.15, 0.8, Form("#mu=%.2f", fgauspi[i]->GetParameter(1)));
    tex->DrawLatex(0.15, 0.7, Form("#sigma=%.2f", fgauspi[i]->GetParameter(2)));
    tex->DrawLatex(0.15, 0.6, Form("#tau=%.2f", fgauspi[i]->GetParameter(3)));

    // Create pure Gaussian function with the same parameters
    fsignalpi[i] = new TF1(Form("fsignal%i", i), "gaus", -5, 5);
    fsignalpi[i]->SetParameter(0, fgauspi[i]->GetParameter(0));
    fsignalpi[i]->SetParameter(1, fgauspi[i]->GetParameter(1));
    fsignalpi[i]->SetParameter(2, fgauspi[i]->GetParameter(2));
    fsignalpi[i]->SetLineColor(kGreen+1);
    fsignalpi[i]->Draw("same");  // Overlay the Gaussian

    csignalpi->SaveAs("signalpi.png");  // Save the canvas

    // Fill the raw yield and its error for pi+
    hrawpi->SetBinContent(i+1, hprojnsigtofpi[i]->Integral(hprojnsigtofpi[i]->GetXaxis()->FindBin(-6), hprojnsigtofpi[i]->GetXaxis()->FindBin(6)));
    hrawpi->SetBinError(i+1, TMath::Sqrt(hrawpi->GetBinContent(i+1)));

    // Copy the same content into the efficiency numerator
    hEffPiPlus->SetBinContent(i+1, hrawpi->GetBinContent(i+1));
    hEffPiPlus->SetBinError(i+1, hrawpi->GetBinError(i+1));
}

// Map the generated true particles from the MC into the denominator histogram
double temp;
for(int ibin = 1; ibin < hGenPiPlus->GetNbinsX(); ibin++) {
    temp = hGenPiPlus->GetXaxis()->GetBinCenter(ibin);  // Get bin center
    double content = hGenPiPlus->GetBinContent(ibin);   // Get bin content
    for (int ifill = 0; ifill < content; ifill++) {      // Fill many times according to the content
        hDenPiPlus->Fill(temp);
    }
}

// Calculate efficiency: numerator / denominator
hEffPiPlus->Divide(hDenPiPlus);

// Plot pi+ efficiency
TCanvas* cEffPiPlus = new TCanvas("Eff Pi", "Eff Pi", 1000, 900);
hEffPiPlus->Draw();
cEffPiPlus->SaveAs("efficiencyPi.png");

// Plot numerator vs denominator for pi+
TCanvas* cNumDenpi = new TCanvas("num dem pi", "num dem pi", 1000, 900);
hDenPiPlus->Draw();
hrawpi->Draw("same hist");
hrawpi->SetLineColor(2);

// Add legend
TLegend* legpi = new TLegend(0.4,0.7,0.7,0.9);
legpi->AddEntry(hDenPiPlus, "Denominator", "l");
legpi->AddEntry(hrawpi, "Numerator", "l");
legpi->Draw();
cNumDenpi->SaveAs("numerator pi.png");

// Define 2D histograms for bin-by-bin (pt) distributions of nSigma_TPC vs nSigma_TOF
TH2D* hBatPiPlus[nbins];
for (int i = 0; i < nbins; i++)
{
    hBatPiPlus[i] = new TH2D(Form("htpcnsigmapi%i", i), "; tpcnsigmapi; tofnsigmapi", 200, -10, 10, 200, -10, 10);
}

// Fill these histograms with 1/10 of the total entries (faster analysis)
for (Long64_t i = 0; i < nEntries / 10; ++i)
{
    tree->GetEntry(i);  // Load entry

    for (int j = 0; j < nbins; j++)
    {
        if (pt > ptbins[j] && pt < ptbins[j+1])
            hBatPiPlus[j]->Fill(tpcnsigmapi, tofnsigmapi);  // Fill 2D histo
    }
}

// Draw and save the 2D bin-by-bin histograms
for (int j = 0; j < nbins; j++)
{
    TCanvas* cBinbinPiPlus = new TCanvas(Form("cBinbinPiPlus_%d", j), Form("pt bin %d", j), 800, 700);
    hBatPiPlus[j]->Draw("colz");
    cBinbinPiPlus->SaveAs(Form("binbinpiplus%d.png", j));
}

    /*This part of the code performs signal extraction and efficiency computation for π⁺ (pi plus) particles.

    It uses Gaussian+exponential fits on the nSigma_TOF distributions for each transverse momentum bin.

    Then, it maps Monte Carlo truth-level generated particles (hGenPiPlus) to fill the denominator.

    Efficiency is computed as reconstructed / generated.

    The same logic is repeated later for π⁻ (pi minus). 

    This exact same analysis (signal extraction, fitting, yield calculation, efficiency evaluation, and 2D diagnostics) is repeated analogously for:

    Protons (pdgcode == 2212)

    Anti-protons (pdgcode == -2212)

    Kaons (pdgcode == 321)

    Anti-kaons (pdgcode == -321)

    They use similar histograms: hprojnsigtofpr, hprojnsigtofpim, hprojnsigtofka, hprojnsigtofkam, etc.
    
    Each with their own efficiency and raw yield histograms like hEffPrPlus, hEffKaMinus, etc.*/

    

    TH1D* hprojnsigtofpim[nbins];
    TF1* fgauspim[nbins];
    TF1* fgausclonepim[nbins];
    TF1* fsignalpim[nbins];
   
    for (int i = 0; i < nbins; i++)
    {
        hprojnsigtofpim[i]= h2nsigtofpim->ProjectionY(Form("hprojnsigtofpim%i",i), h2nsigtofpim->GetXaxis()->FindBin(ptbins[i]+1E-6), h2nsigtofpim->GetXaxis()->FindBin(ptbins[i+1]-1E-6));
        hprojnsigtofpim[i]->SetTitle(Form("%.1f<pt<%.1f", ptbins[i], ptbins[i+1]));      //%.1f keep one digit after comma
        fgauspim[i]=new TF1(Form("fb%i", i), mygausexp, -5,5,4);   
        fgauspim[i]->SetParameter(0,hprojnsigtofpim[i]->GetMaximum());
        fgauspim[i]->SetParameter(1,0);
        fgauspim[i]->SetParameter(2,1);
    }

    TLatex *texm = new TLatex();
    texm->SetNDC();
    texm->SetTextFont(42);
    texm->SetTextSize(0.06);
    
    TCanvas* csignalpim = new TCanvas ("Pim Signal","pim signal", 1500, 900);
    csignalpim->Divide(4,4);
    for (int i = 0; i < nbins; i++)
    {
        csignalpim->cd(i+1);
        hprojnsigtofpim[i]->Draw();
        hprojnsigtofpim[i]->Fit(fgauspim[i],"RW");
        texm->DrawLatex(0.15, 0.8, Form("#mu=%.2f", fgauspim[i]->GetParameter(1)));
        texm->DrawLatex(0.15, 0.7, Form("#sigma=%.2f", fgauspim[i]->GetParameter(2)));
        texm->DrawLatex(0.15, 0.6, Form("#tau=%.2f", fgauspim[i]->GetParameter(3)));
        fsignalpim[i]= new TF1(Form("fsignal%i", i), "gaus", -5,5);
        fsignalpim[i]->SetParameter(0,fgauspim[i]->GetParameter(0));
        fsignalpim[i]->SetParameter(1,fgauspim[i]->GetParameter(1));
        fsignalpim[i]->SetParameter(2,fgauspim[i]->GetParameter(2));
        fsignalpim[i]->SetLineColor(kGreen+1);
        fsignalpim[i]->Draw("same");
        csignalpim->SaveAs("signalpim.png");

        hrawpim->SetBinContent(i+1, hprojnsigtofpim[i]->Integral(hprojnsigtofpim[i]->GetXaxis()->FindBin(-6),hprojnsigtofpim[i]->GetXaxis()->FindBin(6)));
        hrawpim->SetBinError(i+1, TMath::Sqrt(hprojnsigtofpim[i]->Integral(hprojnsigtofpim[i]->GetXaxis()->FindBin(-6),hprojnsigtofpim[i]->GetXaxis()->FindBin(6)))); 

        hEffPiMinus->SetBinContent(i+1, hprojnsigtofpim[i]->Integral(hprojnsigtofpim[i]->GetXaxis()->FindBin(-6),hprojnsigtofpim[i]->GetXaxis()->FindBin(6)));
        hEffPiMinus->SetBinError(i+1, TMath::Sqrt(hprojnsigtofpim[i]->Integral(hprojnsigtofpim[i]->GetXaxis()->FindBin(-6),hprojnsigtofpim[i]->GetXaxis()->FindBin(6)))); 
    }

    double tempm;             //ricorda di inizializzare a 0 prima del prossimo ciclo!!

    for(int ibin = 1; ibin<hGenPiMinus->GetNbinsX(); ibin++){
        tempm = hGenPiMinus->GetXaxis()->GetBinCenter(ibin);    //temp ora è il centro dei bin 
        double contentm = hGenPiMinus->GetBinContent(ibin);     //salvato il valore del bin

        for (int ifill = 0; ifill < contentm ; ifill++){      //loop sui bin dell'istogramma da mappare
            hDenPiMinus->Fill(tempm);
        }
    }

    hEffPiMinus->Divide(hDenPiMinus);

    TCanvas* cEffPi= new TCanvas ("Eff Pim", "Eff Pim", 1000,900);
    hEffPiPlus->Draw("hist");
    hEffPiMinus->Draw("same hist");
    hEffPiMinus->SetLineColor(kGreen+1);
    cEffPi->SaveAs("efficiencyPim.png");

    TCanvas* cNumDenpim= new TCanvas ("num dem pim","num dem pim", 1000, 900); 
    hDenPiMinus->Draw();
    hrawpim->Draw("same hist");
    hrawpim->SetLineColor(2);
    
    TLegend* legpim = new TLegend(0.4,0.7,0.7,0.9);
    //legpi->SetHeader();
    legpim->AddEntry(hDenPiMinus,"Denominator","l");
    legpim->AddEntry(hrawpim,"Numerator","l");
    legpim->Draw();
    cNumDenpim->SaveAs("numerator pim.png");

    ////
    TH1D* hprojnsigtofpr[nbins];
    TF1* fgauspr[nbins];
    TF1* fgausclonepr[nbins];
    TF1* fsignalpr[nbins];

    for (int i = 0; i < nbins; i++)
    {
        hprojnsigtofpr[i]= h2nsigtofpr->ProjectionY(Form("hprojnsigtofpr%i",i), h2nsigtofpr->GetXaxis()->FindBin(ptbins[i]+1E-6), h2nsigtofpr->GetXaxis()->FindBin(ptbins[i+1]-1E-6));
        hprojnsigtofpr[i]->SetTitle(Form("%.1f<pt<%.1f", ptbins[i], ptbins[i+1]));      //%.1f keep one digit after comma
        fgauspr[i]=new TF1(Form("fb%i", i), mygausexp, -5,6,4);
        fgauspr[i]->SetParameter(0,hprojnsigtofpr[i]->GetMaximum());
        fgauspr[i]->SetParameter(1,0);
        fgauspr[i]->SetParameter(2,1);
        fgauspr[i]->SetParameter(3,1);

        
    }
    
    TCanvas* csignalpr = new TCanvas ("Pr Signal","pr signal", 1500, 900);
    csignalpr->Divide(4,4);
    for (int i = 0; i < nbins; i++)
    {
        csignalpr->cd(i+1);
        hprojnsigtofpr[i]->Draw();
        hprojnsigtofpr[i]->Fit(fgauspr[i],"RW");
        tex->DrawLatex(0.15, 0.8, Form("#mu=%.2f", fgauspr[i]->GetParameter(1)));
        tex->DrawLatex(0.15, 0.7, Form("#sigma=%.2f", fgauspr[i]->GetParameter(2)));
        tex->DrawLatex(0.15, 0.6, Form("#tau=%.2f", fgauspr[i]->GetParameter(3)));
        fsignalpr[i]= new TF1(Form("fsignal%i", i), "gaus", -5,5);
        fsignalpr[i]->SetParameter(0,fgauspr[i]->GetParameter(0));
        fsignalpr[i]->SetParameter(1,fgauspr[i]->GetParameter(1));
        fsignalpr[i]->SetParameter(2,fgauspr[i]->GetParameter(2));
        fsignalpr[i]->SetLineColor(kGreen+1);
        fsignalpr[i]->Draw("same");
        csignalpr->SaveAs("signalpr.png");

        hrawpr->SetBinContent(i+1, hprojnsigtofpr[i]->Integral(hprojnsigtofpr[i]->GetXaxis()->FindBin(-6),hprojnsigtofpr[i]->GetXaxis()->FindBin(6)));
        hrawpr->SetBinError(i+1, TMath::Sqrt(hprojnsigtofpr[i]->Integral(hprojnsigtofpr[i]->GetXaxis()->FindBin(-6),hprojnsigtofpr[i]->GetXaxis()->FindBin(6)))); 

        hEffPrPlus->SetBinContent(i+1, hprojnsigtofpr[i]->Integral(hprojnsigtofpr[i]->GetXaxis()->FindBin(-6),hprojnsigtofpr[i]->GetXaxis()->FindBin(6)));
        hEffPrPlus->SetBinError(i+1, TMath::Sqrt(hprojnsigtofpr[i]->Integral(hprojnsigtofpr[i]->GetXaxis()->FindBin(-6),hprojnsigtofpr[i]->GetXaxis()->FindBin(6)))); 


    }
    temp= 0;
    for(int ibin = 1; ibin<hGenPrPlus->GetNbinsX(); ibin++){
        temp = hGenPrPlus->GetXaxis()->GetBinCenter(ibin);     
        double content = hGenPrPlus->GetBinContent(ibin);     

        for (int ifill = 0; ifill < content ; ifill++){      
            hDenPrPlus->Fill(temp);
        }
    }

    hEffPrPlus->Divide(hDenPrPlus);

    TCanvas* cEffPrPlus= new TCanvas ("Eff Pr", "Eff Pr", 1000,900);
    hEffPrPlus->Draw();
    cEffPrPlus->SaveAs("efficiencyPr.png");
    TCanvas* cNumDenpr= new TCanvas ("num dem pr","num dem pr", 1000, 900); 
    hDenPrPlus->Draw();
    hrawpr->Draw("same hist");
    hrawpr->SetLineColor(2);

    TLegend* legpr = new TLegend(0.4,0.7,0.7,0.9);
    //legpr->SetHeader();
    legpr->AddEntry(hDenPrPlus,"Denominator","l");
    legpr->AddEntry(hrawka,"Numerator","l");
    legpr->Draw();
    cNumDenpr->SaveAs("numerator pr.png");

    TH1D* hprojnsigtofprm[nbins];
    TF1* fgausprm[nbins];
    TF1* fgauscloneprm[nbins];
    TF1* fsignalprm[nbins];
   
    for (int i = 0; i < nbins; i++)
    {
        hprojnsigtofprm[i]= h2nsigtofprm->ProjectionY(Form("hprojnsigtofprm%i",i), h2nsigtofprm->GetXaxis()->FindBin(ptbins[i]+1E-6), h2nsigtofprm->GetXaxis()->FindBin(ptbins[i+1]-1E-6));
        hprojnsigtofprm[i]->SetTitle(Form("%.1f<pt<%.1f", ptbins[i], ptbins[i+1]));      //%.1f keep one digit after comma
        fgausprm[i]=new TF1(Form("fb%i", i), mygausexp, -5,5,4);   
        fgausprm[i]->SetParameter(0,hprojnsigtofprm[i]->GetMaximum());
        fgausprm[i]->SetParameter(1,0);
        fgausprm[i]->SetParameter(2,1);
    }

    TLatex *texmpr = new TLatex();
    texmpr->SetNDC();
    texmpr->SetTextFont(42);
    texmpr->SetTextSize(0.06);
    
    TCanvas* csignalprm = new TCanvas ("Prm Signal","prm signal", 1500, 900);
    csignalprm->Divide(4,4);
    for (int i = 0; i < nbins; i++)
    {
        csignalprm->cd(i+1);
        hprojnsigtofprm[i]->Draw();
        hprojnsigtofprm[i]->Fit(fgausprm[i],"RW");
        texmpr->DrawLatex(0.15, 0.8, Form("#mu=%.2f", fgausprm[i]->GetParameter(1)));
        texmpr->DrawLatex(0.15, 0.7, Form("#sigma=%.2f", fgausprm[i]->GetParameter(2)));
        texmpr->DrawLatex(0.15, 0.6, Form("#tau=%.2f", fgausprm[i]->GetParameter(3)));
        fsignalprm[i]= new TF1(Form("fsignal%i", i), "gaus", -5,5);
        fsignalprm[i]->SetParameter(0,fgausprm[i]->GetParameter(0));
        fsignalprm[i]->SetParameter(1,fgausprm[i]->GetParameter(1));
        fsignalprm[i]->SetParameter(2,fgausprm[i]->GetParameter(2));
        fsignalprm[i]->SetLineColor(kGreen+1);
        fsignalprm[i]->Draw("same");
        csignalprm->SaveAs("signalprm.png");

        hrawprm->SetBinContent(i+1, hprojnsigtofprm[i]->Integral(hprojnsigtofprm[i]->GetXaxis()->FindBin(-6),hprojnsigtofprm[i]->GetXaxis()->FindBin(6)));
        hrawprm->SetBinError(i+1, TMath::Sqrt(hprojnsigtofprm[i]->Integral(hprojnsigtofprm[i]->GetXaxis()->FindBin(-6),hprojnsigtofprm[i]->GetXaxis()->FindBin(6)))); 

        hEffPrMinus->SetBinContent(i+1, hprojnsigtofprm[i]->Integral(hprojnsigtofprm[i]->GetXaxis()->FindBin(-6),hprojnsigtofprm[i]->GetXaxis()->FindBin(6)));
        hEffPrMinus->SetBinError(i+1, TMath::Sqrt(hprojnsigtofprm[i]->Integral(hprojnsigtofprm[i]->GetXaxis()->FindBin(-6),hprojnsigtofprm[i]->GetXaxis()->FindBin(6)))); 
    }

    tempm=0;             //ricorda di inizializzare a 0 prima del prossimo ciclo!!

    for(int ibin = 1; ibin<hGenPrMinus->GetNbinsX(); ibin++){
        tempm = hGenPrMinus->GetXaxis()->GetBinCenter(ibin);    //temp ora è il centro dei bin 
        double contentm = hGenPrMinus->GetBinContent(ibin);     //salvato il valore del bin

        for (int ifill = 0; ifill < contentm ; ifill++){      //loop sui bin dell'istogramma da mappare
            hDenPrMinus->Fill(tempm);
        }
    }

    hEffPrMinus->Divide(hDenPrMinus);

    TCanvas* cEffPr= new TCanvas ("Eff Prm", "Eff Prm", 1000,900);
    hEffPrPlus->Draw("hist");
    hEffPrMinus->Draw("same hist");
    hEffPrMinus->SetLineColor(kGreen+1);
    cEffPr->SaveAs("efficiencyPrm.png");

    TCanvas* cNumDenprm= new TCanvas ("num dem prm","num dem prm", 1000, 900); 
    hDenPrMinus->Draw();
    hrawprm->Draw("same hist");
    hrawprm->SetLineColor(2);
    
    TLegend* legprm = new TLegend(0.4,0.7,0.7,0.9);
    legprm->AddEntry(hDenPiMinus,"Denominator","l");
    legprm->AddEntry(hrawpim,"Numerator","l");
    legprm->Draw();
    cNumDenprm->SaveAs("numerator prm.png");



    TH1D* hprojnsigtofka[nbins];
    TF1* fgauska[nbins];
    TF1* fgauscloneka[nbins];
    TF1* fsignalka[nbins];
    for (int i = 0; i < nbins; i++)
    {
        hprojnsigtofka[i]= h2nsigtofka->ProjectionY(Form("hprojnsigtofka%i",i), h2nsigtofka->GetXaxis()->FindBin(ptbinska[i]+1E-6), h2nsigtofka->GetXaxis()->FindBin(ptbinska[i+1]-1E-6));
        hprojnsigtofka[i]->SetTitle(Form("%.1f<pt<%.1f", ptbinska[i], ptbinska[i+1]));      //%.1f keep one digit after comma
        fgauska[i]=new TF1(Form("fgaus%i", i),mygausexp,-5,6,4);
        fgauska[i]->SetParameter(0,hprojnsigtofka[i]->GetMaximum());
        fgauska[i]->SetParameter(1,0);
        fgauska[i]->SetParameter(2,1);
       
    }
    
    TCanvas* csignalka = new TCanvas ("Ka Signal","ka signal", 1500, 900);
    csignalka->Divide(5,2);
    for (int i = 0; i < nbins; i++)
    {
        csignalka->cd(i+1);
        hprojnsigtofka[i]->Draw();
        hprojnsigtofka[i]->Fit(fgauska[i],"RW");
        tex->DrawLatex(0.15, 0.8, Form("#mu=%.2f", fgauska[i]->GetParameter(1)));
        tex->DrawLatex(0.15, 0.7, Form("#sigma=%.2f", fgauska[i]->GetParameter(2)));
        tex->DrawLatex(0.15, 0.6, Form("#tau=%.2f", fgauska[i]->GetParameter(3)));
        fsignalka[i]= new TF1(Form("fsignal%i", i), "gaus", -5,5);
        fsignalka[i]->SetParameter(0,fgauska[i]->GetParameter(0));
        fsignalka[i]->SetParameter(1,fgauska[i]->GetParameter(1));
        fsignalka[i]->SetParameter(2,fgauska[i]->GetParameter(2));
        fsignalka[i]->SetLineColor(kGreen+1);
        fsignalka[i]->Draw("same");
        csignalka->SaveAs("signalka.png");

        hrawka->SetBinContent(i+1, hprojnsigtofka[i]->Integral(hprojnsigtofka[i]->GetXaxis()->FindBin(-6),hprojnsigtofka[i]->GetXaxis()->FindBin(6)));
        hrawka->SetBinError(i+1, TMath::Sqrt(hprojnsigtofka[i]->Integral(hprojnsigtofka[i]->GetXaxis()->FindBin(-6),hprojnsigtofka[i]->GetXaxis()->FindBin(6)))); 

        hEffKaPlus->SetBinContent(i+1, hprojnsigtofka[i]->Integral(hprojnsigtofka[i]->GetXaxis()->FindBin(-6),hprojnsigtofka[i]->GetXaxis()->FindBin(6)));
        hEffKaPlus->SetBinError(i+1, TMath::Sqrt(hprojnsigtofka[i]->Integral(hprojnsigtofka[i]->GetXaxis()->FindBin(-6),hprojnsigtofka[i]->GetXaxis()->FindBin(6)))); 

        
    }
    temp= 0;
    for(int ibin = 1; ibin<hGenKaPlus->GetNbinsX(); ibin++){
        temp = hGenKaPlus->GetXaxis()->GetBinCenter(ibin);     
        double content = hGenKaPlus->GetBinContent(ibin);     

        for (int ifill = 0; ifill < content ; ifill++){      
            hDenKaPlus->Fill(temp);
        }
    }

    hEffKaPlus->Divide(hDenKaPlus);

    TCanvas* cEffKaPlus= new TCanvas ("Eff Ka", "Eff Ka", 1000,900);
    hEffKaPlus->Draw();
    cEffKaPlus->SaveAs("efficiencyKa.png");
    

    TCanvas* cNumDenka= new TCanvas ("num dem ka","num dem ka", 1000, 900); 
    hDenKaPlus->Draw();
    hrawka->Draw("same hist");
    hrawka->SetLineColor(2);

    TLegend* legka = new TLegend(0.4,0.7,0.7,0.9);
    legka->AddEntry(hDenKaPlus,"Denominator","l");
    legka->AddEntry(hrawka,"Numerator","l");
    legka->Draw();
    cNumDenka->SaveAs("numerator ka.png");

    TH1D* hprojnsigtofkam[nbins];
    TF1* fgauskam[nbins];
    TF1* fgausclonekam[nbins];
    TF1* fsignalkam[nbins];
   
    for (int i = 0; i < nbins; i++)
    {
        hprojnsigtofkam[i]= h2nsigtofkam->ProjectionY(Form("hprojnsigtofprm%i",i), h2nsigtofkam->GetXaxis()->FindBin(ptbins[i]+1E-6), h2nsigtofkam->GetXaxis()->FindBin(ptbins[i+1]-1E-6));
        hprojnsigtofkam[i]->SetTitle(Form("%.1f<pt<%.1f", ptbins[i], ptbins[i+1]));      //%.1f keep one digit after comma
        fgauskam[i]=new TF1(Form("fb%i", i), mygausexp, -5,5,4);   
        fgauskam[i]->SetParameter(0,hprojnsigtofkam[i]->GetMaximum());
        fgauskam[i]->SetParameter(1,0);
        fgauskam[i]->SetParameter(2,1);
    }

    TLatex *texmka = new TLatex();
    texmka->SetNDC();
    texmka->SetTextFont(42);
    texmka->SetTextSize(0.06);
    
    TCanvas* csignalkam = new TCanvas ("kam Signal","kam signal", 1500, 900);
    csignalkam->Divide(4,4);
    for (int i = 0; i < nbins; i++)
    {
        csignalkam->cd(i+1);
        hprojnsigtofkam[i]->Draw();
        hprojnsigtofkam[i]->Fit(fgauskam[i],"RW");
        texmka->DrawLatex(0.15, 0.8, Form("#mu=%.2f", fgauskam[i]->GetParameter(1)));
        texmka->DrawLatex(0.15, 0.7, Form("#sigma=%.2f", fgauskam[i]->GetParameter(2)));
        texmka->DrawLatex(0.15, 0.6, Form("#tau=%.2f", fgauskam[i]->GetParameter(3)));
        fsignalkam[i]= new TF1(Form("fsignal%i", i), "gaus", -5,5);
        fsignalkam[i]->SetParameter(0,fgauskam[i]->GetParameter(0));
        fsignalkam[i]->SetParameter(1,fgauskam[i]->GetParameter(1));
        fsignalkam[i]->SetParameter(2,fgauskam[i]->GetParameter(2));
        fsignalkam[i]->SetLineColor(kGreen+1);
        fsignalkam[i]->Draw("same");
        csignalkam->SaveAs("signalkam.png");

        hrawkam->SetBinContent(i+1, hprojnsigtofkam[i]->Integral(hprojnsigtofkam[i]->GetXaxis()->FindBin(-6),hprojnsigtofkam[i]->GetXaxis()->FindBin(6)));
        hrawkam->SetBinError(i+1, TMath::Sqrt(hprojnsigtofkam[i]->Integral(hprojnsigtofkam[i]->GetXaxis()->FindBin(-6),hprojnsigtofkam[i]->GetXaxis()->FindBin(6)))); 

        hEffKaMinus->SetBinContent(i+1, hprojnsigtofkam[i]->Integral(hprojnsigtofkam[i]->GetXaxis()->FindBin(-6),hprojnsigtofkam[i]->GetXaxis()->FindBin(6)));
        hEffKaMinus->SetBinError(i+1, TMath::Sqrt(hprojnsigtofkam[i]->Integral(hprojnsigtofkam[i]->GetXaxis()->FindBin(-6),hprojnsigtofkam[i]->GetXaxis()->FindBin(6)))); 
    }

    tempm=0;             //ricorda di inizializzare a 0 prima del prossimo ciclo!!

    for(int ibin = 1; ibin<hGenKaMinus->GetNbinsX(); ibin++){
        tempm = hGenKaMinus->GetXaxis()->GetBinCenter(ibin);    //temp ora è il centro dei bin 
        double contentm = hGenKaMinus->GetBinContent(ibin);     //salvato il valore del bin

        for (int ifill = 0; ifill < contentm ; ifill++){      //loop sui bin dell'istogramma da mappare
            hDenKaMinus->Fill(tempm);
        }
    }

    hEffKaMinus->Divide(hDenKaMinus);

    TCanvas* cEffKa= new TCanvas ("Eff Kam", "Eff Kam", 1000,900);
    hEffKaPlus->Draw("hist");
    hEffKaMinus->Draw("same hist");
    hEffKaMinus->SetLineColor(kGreen+1);
    cEffKa->SaveAs("efficiencyKam.png");

    TCanvas* cNumDenkam= new TCanvas ("num dem kam","num dem kam", 1000, 900); 
    hDenKaMinus->Draw();
    hrawkam->Draw("same hist");
    hrawkam->SetLineColor(2);
    
    TLegend* legkam = new TLegend(0.4,0.7,0.7,0.9);
    legkam->AddEntry(hDenKaMinus,"Denominator","l");
    legkam->AddEntry(hrawkam,"Numerator","l");
    legkam->Draw();
    cNumDenkam->SaveAs("numerator kam.png");




    // Create canvas to draw pT distribution and display it
TCanvas* cpt = new TCanvas("cpt","cpt",1000,900);
hPt->Draw();  // Plot transverse momentum histogram

// Create canvas for eta (pseudorapidity) distribution
TCanvas* ceta = new TCanvas("ceta","ceta",1000,900);
hEta->Draw();  // Plot pseudorapidity histogram

// ---------------------------- TPC π signal (π± from TPC nsigma) ----------------------------

TCanvas* c2sigpi = new TCanvas ("tpcsigpi","tpcsigpi", 1000, 900);
c2sigpi->SetLogz();  // Logarithmic color scale for better visualization
h2nsigpi->SetStats(0);  // Remove statistics box

// Define the pT bin range to perform Y-slice fits (nsigma TPC vs pT)
const int min_binpi = h2nsigpi->GetXaxis()->FindBin(0.2);
const int max_binpi = h2nsigpi->GetXaxis()->FindBin(0.8);

// Define a simple Gaussian function for fitting each Y slice
TF1* pigaus = new TF1("fgaus","gaus",-1.5,+1.5);

// Fit each Y slice with Gaussian and store results (mean, sigma) in slicespi
TObjArray* slicespi = new TObjArray();
h2nsigpi->FitSlicesY(pigaus, min_binpi, max_binpi, 0, "WW,QNR", slicespi);

// Get sigma and mean histograms from the fit slices
TH1D* hSigmapi = (TH1D*)slicespi->At(2);
hSigmapi->SetName("hSigmapi");
TH1D* hMeanpi = (TH1D*)slicespi->At(1);

// Style the curves for visibility
hSigmapi->SetLineColor(2); hSigmapi->SetLineWidth(3);
hMeanpi->SetLineColor(4); hMeanpi->SetLineWidth(3);

// Limit the pT range shown on X axis
h2nsigpi->GetXaxis()->SetRangeUser(0,2);

// Draw original 2D histogram and overlay mean and sigma curves
h2nsigpi->Draw();
hSigmapi->Draw("same");
hMeanpi->Draw("same");

// Save the canvas to image
c2sigpi->SaveAs("nsigpi.png");

// ---------------------------- Raw π+ yield ----------------------------

TCanvas* crawpi = new TCanvas("crawpi","crawpi", 1000, 900);
hrawpi->Draw();  // Draw raw yield for π+
crawpi->SaveAs("raw pi.png");

// ---------------------------- TPC Proton signal ----------------------------

TCanvas* c2sigpr = new TCanvas ("tpcsigpr","tpcsigpr", 1000, 900);
c2sigpr->SetLogz();
h2nsigpr->SetStats(0);

const int min_binpr = h2nsigpr->GetXaxis()->FindBin(0.2);
const int max_binpr = h2nsigpr->GetXaxis()->FindBin(0.8);

TF1* prgaus = new TF1("fgaus","gaus",-1.5,+1.5);

TObjArray* slicespr = new TObjArray();
h2nsigpr->FitSlicesY(prgaus, min_binpr, max_binpr, 0, "WW,QNR", slicespr);

TH1D* hSigmapr = (TH1D*)slicespr->At(2);
hSigmapr->SetName("hSigmapr");
TH1D* hMeanpr = (TH1D*)slicespr->At(1);
hSigmapr->SetLineColor(2); hSigmapr->SetLineWidth(3);
hMeanpr->SetLineColor(4); hMeanpr->SetLineWidth(3);
h2nsigpr->GetXaxis()->SetRangeUser(0,2);
h2nsigpr->Draw();
hSigmapr->Draw("same");
hMeanpr->Draw("same");
c2sigpr->SaveAs("nsigpr.png");

// ---------------------------- Raw Proton yield ----------------------------

TCanvas* crawpr = new TCanvas("crawpr","crawpr", 1000, 900);
hrawpr->Draw();
crawpr->SaveAs("raw pr.png");

// ---------------------------- TPC Kaon signal ----------------------------

TCanvas* c2sigka = new TCanvas ("tpcsigka","tpcsigka", 1000, 900);
c2sigka->SetLogz();
h2nsigka->SetStats(0);

const int min_binka = h2nsigka->GetXaxis()->FindBin(0.2);
const int max_binka = h2nsigka->GetXaxis()->FindBin(0.6);

TF1* kagaus = new TF1("fgaus","gaus",-2,2);

TObjArray* sliceska = new TObjArray();
h2nsigka->FitSlicesY(kagaus, min_binka, max_binka, 0, "WW,QNR", sliceska);

TH1D* hSigmaka = (TH1D*)sliceska->At(2);
TH1D* hMeanka = (TH1D*)sliceska->At(1);
hSigmaka->SetLineColor(2); hSigmaka->SetLineWidth(3);
hMeanka->SetLineColor(3); hMeanka->SetLineWidth(3);
h2nsigka->GetXaxis()->SetRangeUser(0,2);
h2nsigka->Draw();
hSigmaka->Draw("same");
hMeanka->Draw("same");
c2sigka->SaveAs("nsigka.png");

// ---------------------------- Raw Kaon yield ----------------------------

TCanvas* crawka = new TCanvas("crawka","crawka", 1000, 900);
hrawka->Draw();
crawka->SaveAs("raw ka.png");

// ---------------------------- TOF π signal (π± from TOF nsigma) ----------------------------

TCanvas* c2nsigtofpi = new TCanvas ("tofsigpi","tofsigpi", 1000, 900);
c2nsigtofpi->SetLogz();
h2nsigtofpi->SetStats(0);
h2nsigtofpi->Draw("colz");

const int min_binpitof = h2nsigtofpi->GetXaxis()->FindBin(0.4);
const int max_binpitof = h2nsigtofpi->GetXaxis()->FindBin(2);

TF1* pigaustof = new TF1("fgaus","gaus",-1.5,+1.5);

TObjArray* slicespitof = new TObjArray();
h2nsigtofpi->FitSlicesY(pigaus, min_binpitof, max_binpitof, 0, "WW,QNR", slicespitof);

TH1D* hSigmapitof = (TH1D*)slicespitof->At(2);
hSigmapitof->SetName("hSigmapitof");
TH1D* hMeanpitof = (TH1D*)slicespitof->At(1);
hSigmapitof->SetLineColor(2); hSigmapitof->SetLineWidth(3);
hMeanpitof->SetLineColor(4); hMeanpitof->SetLineWidth(3);
h2nsigtofpi->GetXaxis()->SetRangeUser(0,2);
h2nsigtofpi->Draw();
hSigmapitof->Draw("same");
hMeanpitof->Draw("same");
hMeanpitof->SetName("hMeanpi");
c2nsigtofpi->SaveAs("nsigtofpi.png");

// ---------------------------- TOF Proton signal ----------------------------

TCanvas* c2nsigtofpr = new TCanvas ("tofsigpr","tofsigpr", 1000, 900);
h2nsigtofpr->SetStats(0);
h2nsigtofpr->Draw("colz");

const int min_binprtof = h2nsigtofpr->GetXaxis()->FindBin(0.4);
const int max_binprtof = h2nsigtofpr->GetXaxis()->FindBin(2);

TF1* prgaustof = new TF1("fgaus","gaus",-1.5,+1.5);

TObjArray* slicesprtof = new TObjArray();
h2nsigtofpr->FitSlicesY(prgaus, min_binprtof, max_binprtof, 0, "WW,QNR", slicesprtof);

TH1D* hSigmaprtof = (TH1D*)slicesprtof->At(2);
hSigmaprtof->SetName("hSigmaprtof");
TH1D* hMeanprtof = (TH1D*)slicesprtof->At(1);
hSigmaprtof->SetLineColor(2); hSigmaprtof->SetLineWidth(3);
hMeanprtof->SetLineColor(4); hMeanprtof->SetLineWidth(3);
h2nsigtofpr->GetXaxis()->SetRangeUser(0,2);
h2nsigtofpr->Draw();
hSigmaprtof->Draw("same");
hMeanprtof->Draw("same");
hMeanprtof->SetName("hMeanpr");
c2nsigtofpr->SaveAs("nsigtofpr.png");

// ---------------------------- TOF Kaon signal ----------------------------

TCanvas* c2nsigtofka = new TCanvas ("tofsigka","tofsigka", 1000, 900);
h2nsigtofka->SetStats(0);
h2nsigtofka->Draw("colz");

const int min_binkatof = h2nsigtofka->GetXaxis()->FindBin(0.4);
const int max_binkatof = h2nsigtofka->GetXaxis()->FindBin(2);

TF1* kagaustof = new TF1("fgaus","gaus",-2,2);

TObjArray* sliceskatof = new TObjArray();
h2nsigtofka->FitSlicesY(kagaustof, min_binkatof, max_binkatof, 0, "WW,QNR", sliceskatof);

TH1D* hSigmakatof = (TH1D*)sliceskatof->At(2);
TH1D* hMeankatof = (TH1D*)sliceskatof->At(1);
hSigmakatof->SetLineColor(2); hSigmakatof->SetLineWidth(3);
hMeankatof->SetLineColor(3); hMeankatof->SetLineWidth(3);
h2nsigtofka->GetXaxis()->SetRangeUser(0,2);
h2nsigtofka->Draw();
hSigmakatof->Draw("same");
hMeankatof->Draw("same");
hMeankatof->SetName("hMeanpr");
c2nsigtofka->SaveAs("nsigtofka.png");

// ---------------------------- Save Efficiencies to ROOT file ----------------------------

// Save all the efficiency histograms into a ROOT file for later analysis
TFile* write = new TFile ("efficienze.root", "RECREATE");
hEffPiPlus->Write();
hEffPiMinus->Write();
hEffPrPlus->Write();
hEffPrMinus->Write();
hEffKaPlus->Write();
hEffKaMinus->Write();
}
/*Draws and saves histograms for transverse momentum, pseudorapidity, and raw yields.

Performs Gaussian fits on slices of 2D nsigma (from TPC or TOF) vs pT histograms to extract mean and sigma values of particle signals.

Saves plots of these fits for π, K, and p (both TPC and TOF).

Stores efficiency histograms into a ROOT file called "efficienze.root" for π⁺/π⁻, K⁺/K⁻, and p⁺/p⁻.*/