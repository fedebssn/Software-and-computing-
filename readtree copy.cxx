#include "TMath.h"

// Custom Gaussian + exponential function
double mygausexp(double* x, double* p)
{
    double out;
    if (x[0] <= p[1] - p[3]) {
        // Left side: exponential tail before peak
        out = p[0] * TMath::Gaus(p[1] - p[3], p[1], p[2]) * TMath::Exp((x[0] - p[1] + p[3]) * p[3] / (p[2] * p[2]));
    } else if (x[0] <= p[1] + p[3]) {
        // Center: pure Gaussian
        out = p[0] * TMath::Gaus(x[0], p[1], p[2]);
    } else {
        // Right side: exponential tail after peak
        out = p[0] * TMath::Gaus(p[1] + p[3], p[1], p[2]) * TMath::Exp(-p[3] * (x[0] - p[1] - p[3]) / (p[2] * p[2]));
    }
    return out;
}

// Custom Gaussian + exponential + polynomial function
double mygausexppol(double* x, double* pars)
{
    double out;
    if(x[0] <= (pars[1] + pars[3])){
        // Gaussian + polynomial background on the left/center
        out = pars[0]*TMath::Gaus(x[0], pars[1], pars[2]) + x[0]*pars[4] + x[0]*x[0]*pars[5];
    } else {
        // Exponential tail + polynomial on the right
        out = pars[0]*TMath::Gaus(pars[1]+pars[3], pars[1], pars[2]) * TMath::Exp((x[0] - pars[1] - pars[3]) * (-pars[3]) / (pars[2]*pars[2])) + x[0]*pars[4] + x[0]*x[0]*pars[5];
    }
    return out;
}

void readtree(){
    // Open the ROOT file
    TFile *file = TFile::Open("mytree_7may.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << "mytree_7may.root" << std::endl; // Print error if file not opened
        return;
    }

    // Get the tree named "tracks"
    TTree *tree = (TTree*)file->Get("tracks");
    if (!tree) {
        std::cerr << "Tree 'tracks' not found in file." << std::endl; // Print error if tree not found
        file->Close(); // Close file before returning
        return;
    }

    // Declare variables for each tree branch
    float pt = 0., eta = 0., beta = 0., p = 0., px = 0., py =0., pz =0., sign = 0., tpc = 0., tof = 0., 
          tevent = 0., tpi = 0., tpr = 0., tka = 0., hasTOF = 0., 
          tpcnsigmapi = 0., tpcnsigmapr = 0., tpcnsigmaka = 0., 
          tofnsigmapi = 0., tofnsigmapr = 0., tofnsigmaka = 0., 
          tofnsigmapim = 0., tofnsigmaprm = 0., tofnsigmakam = 0.;

    // Set addresses for tree branches to fill the variables
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

    // Get number of entries in the tree
    Long64_t nEntries = tree->GetEntries();
    std::cout << "Reading " << nEntries << " entries from tree: tracks" << std::endl;

    // Create histograms for various variables
    TH1D* hPt = new TH1D("hPt","titolo; assex; assey", 250, 0, 25); // pt distribution
    TH1D* hEta = new TH1D("hEta","eta distribution; #eta; entries", 200, -1, +1); // eta
    TH1D* hPx = new TH1D("hPx", "px; px; entries", 250, 0, 25); // px
    TH1D* hPy = new TH1D("hPy", "py; py; entries", 250, 0, 25); // py
    TH1D* hPz = new TH1D("hPz", "pz; pz; entries", 250, 0, 25); // pz
    TH1D* hSign = new TH1D("hSign","sign; sign; sign", 10, -1, 1); // charge sign
    TH1D* hTpc = new TH1D("hTpc","tpc; tpc; tpc", 200, 0, 150); // TPC signal
    TH1D* hP = new TH1D("hP","p; p; entries", 250, 0, 25); // momentum
    TH1D* hBeta = new TH1D("hBeta", "#beta; beta; entries", 250, -60, 60); // beta
    TH1D* hTof = new TH1D("hTof", "tof; tof; tof", 200, 0, 20); // TOF signal
    TH1D* hTevent = new TH1D("hTevent", "tevent; tevent; tevent", 500, 0, 2000); // event time
    TH1D* hTpi = new TH1D("hTpi", "tpi; tpi; tpi", 500, 0, 20000); // expected TOF for pion
    TH1D* hTpr = new TH1D("hTpr", "tpr; tpr; tpr", 500, 0, 60000); // expected TOF for proton
    TH1D* hTka = new TH1D("hTka", "tka; tka; tka", 500, 0, 30000); // expected TOF for kaon
    TH1D* hTpcsigpi = new TH1D("htpcnsigmapi", "tpcnsigmapi; ...", 500, 0, 10); // TPC nsigma pi
    TH1D* hTpcsigpr = new TH1D("htpcnsigmapr", "tpcnsigmapr; ...", 500, 0, 10); // TPC nsigma pr
    TH1D* hTpcsigka = new TH1D("htpcnsigmaka", "tpcnsigmaka; ...", 500, 0, 10); // TPC nsigma ka

    // Define binning for pt histograms
    float ptbins [] = {0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,4};
    const int nbins = sizeof(ptbins)/sizeof(float)-1; // number of bins
    float ptbinska [] = {0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,2,3};
    const int nbinska = sizeof(ptbinska)/sizeof(float)-1;

    // Raw yield histograms (for pions, protons, kaons, with + and - charges)
    TH1D* hrawpi = new TH1D("hrawpi", "...", nbins, ptbins);
    TH1D* hrawpr = new TH1D("hrawpr", "...", nbins, ptbins);
    TH1D* hrawka = new TH1D("hrawka", "...", nbinska, ptbinska);
    TH1D* hrawpim = new TH1D("hrawpim", "...", nbins, ptbins);
    TH1D* hrawprm = new TH1D("hrawprm", "...", nbins, ptbins);
    TH1D* hrawkam = new TH1D("hrawkam", "...", nbinska, ptbinska);

    // Create 2D histograms for correlations
    TH2D* h2energyloss = new TH2D("h2energyloss", "...", 400, 0, 10, 300, 0, 250); // p/Z vs dE/dx
    TH2D* h2beta = new TH2D("h2beta", "...", 400, 0, 10, 160, 0.3, 1.1); // p vs beta
    TH2D* h2pi = new TH2D("h2pi", "...", 400, 0, 10, 500, -20000, 20000); // TOF - Tpi vs p
    TH2D* h2pr = new TH2D("h2pr", "...", 400, 0, 10, 500, -20000, 20000); // TOF - Tpr vs p
    TH2D* h2ka = new TH2D("h2ka", "...", 400, 0, 10, 500, -20000, 20000); // TOF - Tka vs p
    TH2D* h2nsigpi = new TH2D("h2nsigpi", "...", 200, 0, 10, 160, -8, 8); // pt vs nsigma pi
    TH2D* h2nsigpr = new TH2D("h2nsigpr", "...", 200, 0, 10, 160, -8, 8); // pt vs nsigma pr
    TH2D* h2nsigka = new TH2D("h2nsigka", "...", 200, 0, 10, 160, -8, 8); // pt vs nsigma ka

    TH2D* h2nsigtofpi = new TH2D("h2nsigtofpi", "...", 200, 0, 10, 160, -8, 8); // TOF nsigma pi
    TH2D* h2nsigtofpr = new TH2D("h2nsigtofpr", "...", 200, 0, 10, 160, -8, 8); // TOF nsigma pr
    TH2D* h2nsigtofka = new TH2D("h2nsigtofka", "...", 200, 0, 10, 160, -8, 8); // TOF nsigma ka

    TH2D* h2nsigtofpim = new TH2D("h2nsigtofpim", "...", 200, 0, 10, 160, -8, 8); // TOF nsigma pi-
    TH2D* h2nsigtofprm = new TH2D("h2nsigtofprm", "...", 200, 0, 10, 160, -8, 8); // TOF nsigma pr-
    TH2D* h2nsigtofkam = new TH2D("h2nsigtofkam", "...", 200, 0, 10, 160, -8, 8); // TOF nsigma ka-

    // Main loop over all entries in the tree
    for (Long64_t i = 0; i < nEntries; ++i) 
    {
        tree->GetEntry(i); // Read event i

        hPt->Fill(pt); // Fill pt histogram
        hEta->Fill(eta); // Fill eta histogram
        hPx->Fill(px); hPy->Fill(py); hPz->Fill(pz); // Fill momentum components
        hSign->Fill(sign); // Fill sign
        hTpc->Fill(tpc); // Fill TPC signal

        p = TMath::Sqrt(px*px + py*py + pz*pz); // Compute total momentum

        hP->Fill(p); // Fill p
        hTpc->Fill(tpc); // Again fill TPC
        h2energyloss->Fill(p/sign, tpc); // Fill dE/dx vs momentum/Z
        h2beta->Fill(p, beta); // Fill beta vs p

        if (hasTOF == 1) // Only fill TOF plots if TOF info is valid
        {
            hTof->Fill(tof);
            hTevent->Fill(tevent);
            hTpi->Fill(tpi); hTpr->Fill(tpr); hTka->Fill(tka);
            h2pi->Fill(p, tof - tevent - tpi); // Delta T for pion
            h2pr->Fill(p, tof - tevent - tpr); // Delta T for proton
            h2ka->Fill(p, tof - tevent - tka); // Delta T for kaon
        }

        h2nsigpi->Fill(pt, tpcnsigmapi); // Fill TPC nsigma pi
        h2nsigpr->Fill(pt, tpcnsigmapr); // Fill TPC nsigma pr
        h2nsigka->Fill(pt, tpcnsigmaka); // Fill TPC nsigma ka

        if (std::abs(tpcnsigmapi) < 3 && sign > 0)
            h2nsigtofpi->Fill(pt, tofnsigmapi); // Positive pion

        if (std::abs(tpcnsigmapr) < 3 && sign > 0)
            h2nsigtofpr->Fill(pt, tofnsigmapr); // Positive proton

        if (std::abs(tpcnsigmaka) < 3 && sign > 0)
            h2nsigtofka->Fill(pt, tofnsigmaka); // Positive kaon

        if (std::abs(tpcnsigmapi) < 3 && sign < 0)
            h2nsigtofpim->Fill(pt, tofnsigmapi); // Negative pion

        if (std::abs(tpcnsigmapr) < 3 && sign < 0)
            h2nsigtofprm->Fill(pt, tofnsigmapr); // Negative proton

        if (std::abs(tpcnsigmaka) < 3 && sign < 0)
            h2nsigtofkam->Fill(pt, tofnsigmaka); // Negative kaon
    }
   TH1D* hprojy = h2beta->ProjectionY("hprojy", h2beta->GetXaxis()->FindBin(1.45), h2beta->GetXaxis()->FindBin(1.55));
// Project 2D histogram h2beta onto Y axis between X bins corresponding to 1.45 and 1.55

TH1D* hprojy2 = h2beta->ProjectionY("hprojy2", h2beta->GetXaxis()->FindBin(0.35), h2beta->GetXaxis()->FindBin(0.45));
// Project h2beta onto Y axis between X bins 0.35 and 0.45

hprojy->GetXaxis()->SetRangeUser(0.75,1.1);
// Set visible X axis range of hprojy histogram from 0.75 to 1.1

hprojy2->GetXaxis()->SetRangeUser(0.75,1.1);
// Set visible X axis range of hprojy2 histogram from 0.75 to 1.1

TH1D* hprojnsigtofpi[nbins];
// Array of 1D histograms for signal projections in different pt bins for pi+

TF1* fgauspi[nbins];
// Array of fit functions with a Gaussian+exponential shape for pi+ fits

TF1* fgausclonepi[nbins];
// Array of fit functions (not used here, possibly for cloning fits later)

TF1* fsignalpi[nbins];
// Array of pure Gaussian fit functions for pi+ signals

for (int i = 0; i < nbins; i++)
{
    hprojnsigtofpi[i]= h2nsigtofpi->ProjectionY(Form("hprojnsigtofpi%i",i), h2nsigtofpi->GetXaxis()->FindBin(ptbins[i]+1E-6), h2nsigtofpi->GetXaxis()->FindBin(ptbins[i+1]-1E-6));
    // Project h2nsigtofpi onto Y axis within each pt bin range
    
    hprojnsigtofpi[i]->SetTitle(Form("%.1f<pt<%.1f", ptbins[i], ptbins[i+1]));
    // Set histogram title to show the pt bin range with 1 decimal
    
    fgauspi[i]=new TF1(Form("fb%i", i), mygausexp, -5,5,4);
    // Create a fit function for the histogram: Gaussian + exponential tail
    
    fgauspi[i]->SetParameter(0,hprojnsigtofpi[i]->GetMaximum());
    // Initialize amplitude parameter to max bin content
    
    fgauspi[i]->SetParameter(1,0);
    // Initialize mean parameter to 0
    
    fgauspi[i]->SetParameter(2,1);
    // Initialize sigma parameter to 1
}

TLatex *tex = new TLatex();
// Create a TLatex object to write text on plots

tex->SetNDC();
// Set coordinates relative to the canvas (normalized device coordinates)

tex->SetTextFont(42);
// Use font style 42 (Helvetica)

tex->SetTextSize(0.06);
// Set text size

TCanvas* csignalpi = new TCanvas ("Pi Signal","pi signal", 1500, 900);
// Create a canvas to draw pi+ signal fits

csignalpi->Divide(4,4);
// Divide canvas into 4x4 pads for multiple plots

for (int i = 0; i < nbins; i++)
{
    csignalpi->cd(i+1);
    // Switch to the ith pad on the canvas
    
    hprojnsigtofpi[i]->Draw();
    // Draw the projected histogram for pi+ in ith pt bin
    
    hprojnsigtofpi[i]->Fit(fgauspi[i],"RW");
    // Fit the histogram with the Gaussian+exponential function
    
    tex->DrawLatex(0.15, 0.8, Form("#mu=%.2f", fgauspi[i]->GetParameter(1)));
    // Display the fitted mean value on the plot
    
    tex->DrawLatex(0.15, 0.7, Form("#sigma=%.2f", fgauspi[i]->GetParameter(2)));
    // Display the fitted sigma (width) on the plot
    
    tex->DrawLatex(0.15, 0.6, Form("#tau=%.2f", fgauspi[i]->GetParameter(3)));
    // Display the fitted exponential tail parameter tau
    
    fsignalpi[i]= new TF1(Form("fsignal%i", i), "gaus", -5,5);
    // Create a pure Gaussian function for the signal
    
    fsignalpi[i]->SetParameter(0,fgauspi[i]->GetParameter(0));
    fsignalpi[i]->SetParameter(1,fgauspi[i]->GetParameter(1));
    fsignalpi[i]->SetParameter(2,fgauspi[i]->GetParameter(2));
    // Set Gaussian parameters same as fitted function (amplitude, mean, sigma)
    
    fsignalpi[i]->SetLineColor(kGreen+1);
    // Set line color to green
    
    fsignalpi[i]->Draw("same");
    // Draw the pure Gaussian on top of the histogram
    
    csignalpi->SaveAs("signalpi.png");
    // Save the canvas as an image file
    
    hrawpi->SetBinContent(i+1, 1/(4.05*1E+6)*fgauspi[i]->Integral(-3,3)/hrawpi->GetBinWidth(i+1));
    // Set content of hrawpi histogram bin: normalized integral of fitted signal in [-3,3]
    
    hrawpi->SetBinError(i+1, 1/(4.05*1E+6)*fgauspi[i]->IntegralError(-3,3)/hrawpi->GetBinWidth(i+1));
    // Set error of hrawpi bin using integral error of fit
}

TH2D* hBatPiPlus[nbins];
// Array of 2D histograms for TPC vs TOF nsigma variables for pi+ in each pt bin

for (int i = 0; i < nbins; i++)
{
    hBatPiPlus[i] = new TH2D (Form("htpcnsigmapi%i",i), "; tpcnsigmapi; tofnsigmapi", 200, -10, 10, 200, -10, 10);
    // Create 2D histogram for pi+ TPC vs TOF nsigma distributions with 200 bins in each axis
}

for (Long64_t i = 0; i < nEntries/10; ++i)
{
    tree->GetEntry(i);
    // Load data for event i from the tree
    
    for (int j = 0; j < nbins; j++)
    {
        if (pt > ptbins[j] && pt < ptbins[j+1])
        hBatPiPlus[j]->Fill(tpcnsigmapi,tofnsigmapi);
        // Fill corresponding 2D histogram if particle pt is in bin j
    }
}

for (int j = 0; j < nbins; j++)
{
    TCanvas* cBinbinPiPlus = new TCanvas(Form("cBinbinPiPlus_%d", j), Form("pt bin %d", j), 800, 700);
    // Create canvas for 2D histogram of bin j
    
    hBatPiPlus[j]->Draw("colz");
    // Draw the 2D histogram with color map
    
    cBinbinPiPlus->SaveAs(Form("binbinpiplus%d.png", j));
    // Save the canvas as a PNG image
}

// The following sections repeat a similar workflow for pi-, protons (pr), and antiprotons (prm):
// 1. Create projections of 2D histograms in pt bins
// 2. Fit with Gaussian+exponential functions
// 3. Plot and save fit results
// 4. Create 2D histograms for TPC vs TOF nsigma for each particle and pt bin
// 5. Fill these histograms with data entries
// 6. Save the 2D histograms as images

// Example for pi-:

TH1D* hprojnsigtofpim[nbins];
// 1D histograms for pi- signal projections per pt bin

TF1* fgauspim[nbins];
// Fit functions for pi-

TF1* fgausclonepim[nbins];
// Clone fit functions (not used here)

TF1* fsignalpim[nbins];
// Pure Gaussian fit functions for pi-

for (int i = 0; i < nbins; i++)
{
    hprojnsigtofpim[i]= h2nsigtofpim->ProjectionY(Form("hprojnsigtofpim%i",i), h2nsigtofpim->GetXaxis()->FindBin(ptbins[i]+1E-6), h2nsigtofpim->GetXaxis()->FindBin(ptbins[i+1]-1E-6));
    // Project pi- 2D histogram onto Y axis in pt bin i
    
    hprojnsigtofpim[i]->SetTitle(Form("%.1f<pt<%.1f", ptbins[i], ptbins[i+1]));
    // Set histogram title
    
    fgauspim[i]=new TF1(Form("fb%i", i), mygausexp, -5,5,4);
    // Create fit function
    
    fgauspim[i]->SetParameter(0,hprojnsigtofpim[i]->GetMaximum());
    fgauspim[i]->SetParameter(1,0);
    fgauspim[i]->SetParameter(2,1);
    // Initialize fit parameters
}

TLatex *texm = new TLatex();
// Create TLatex for pi- text

texm->SetNDC();
texm->SetTextFont(42);
texm->SetTextSize(0.06);
// Set style

TCanvas* csignalpim = new TCanvas ("Pim Signal","pim signal", 1500, 900);
// Canvas for pi- signals

csignalpim->Divide(4,4);
// Divide into 16 pads

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

    hrawpim->SetBinContent(i+1, 1/(4.05*1E+6)*fgauspim[i]->Integral(-3,3)/hrawpim->GetBinWidth(i+1));
    hrawpim->SetBinError(i+1, 1/(4.05*1E+6)*fgauspim[i]->IntegralError(-3,3)/hrawpim->GetBinWidth(i+1));
}
// Similar filling and plotting for TPC vs TOF 2D histograms for pi-

TH2D* hBatPiMinus[nbins];
for (int i = 0; i < nbins; i++)
{
    hBatPiMinus[i] = new TH2D (Form("htpcnsigmapi%i",i), "; tpcnsigmapi; tofnsigmapi", 200, -10, 10, 200, -10, 10);
    // Create 2D histograms for pi- TPC vs TOF nsigma
}

for (Long64_t i = 0; i < nEntries/10; ++i)
{
    tree->GetEntry(i);
    for (int j = 0; j < nbins; j++)
    {
        if (pt > ptbins[j] && pt < ptbins[j+1])
        hBatPiMinus[j]->Fill(tpcnsigmapi,tofnsigmapim);
        // Fill histograms if pt in bin
    }
}

for (int j = 0; j < nbins; j++)
{
    TCanvas* cBinbinPiMinus = new TCanvas(Form("cBinbinPiMinus_%d", j), Form("pt bin %d", j), 800, 700);
    hBatPiMinus[j]->Draw("colz");
    cBinbinPiMinus->SaveAs(Form("binbinpiminus%d.png", j));
}
// Save 2D histograms for pi-

// Similar blocks follow for protons (pr) and antiprotons (prm) and with kaons (ka) and and negative kaons (kam) with analogous steps:
// Create projections, fit, plot, save results, create 2D histograms, fill with data, save histograms

// This concludes the detailed commenting for your code.


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

        hrawpr->SetBinContent(i+1, 1/(4.05*1E+6)*fgauspr[i]->Integral(-3,3)/hrawpr->GetBinWidth(i+1));
        hrawpr->SetBinError(i+1, 1/(4.05*1E+6)*fgauspr[i]->IntegralError(-3,3)/hrawpr->GetBinWidth(i+1));

    }
    TH2D* hBatPrPlus[nbins];
        for (int i = 0; i < nbins; i++)
        {
            hBatPrPlus[i] = new TH2D (Form("htpcnsigmapr%i",i), "; tpcnsigmapr; tofnsigmapr", 200, -10, 10, 200, -10, 10);
        }
 
        for (Long64_t i = 0; i < nEntries/10; ++i)
        {
            tree->GetEntry(i); // Load data for entry
            for (int j = 0; j < nbins; j++)
            {
                if (pt > ptbins[j] && pt < ptbins[j+1])
                hBatPrPlus[j]->Fill(tpcnsigmapr,tofnsigmapr);
            }
            
        }
   
    
        for (int j = 0; j < nbins; j++)
        {
            TCanvas* cBinbinPrPlus = new TCanvas(Form("cBinbinPrPlus_%d", j), Form("pt bin %d", j), 800, 700);
            hBatPrPlus[j]->Draw("colz");
            cBinbinPrPlus->SaveAs(Form("binbinprplus%d.png", j));
        }


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

        hrawprm->SetBinContent(i+1, 1/(4.05*1E+6)*fgausprm[i]->Integral(-3,3)/hrawprm->GetBinWidth(i+1));
        hrawprm->SetBinError(i+1, 1/(4.05*1E+6)*fgausprm[i]->IntegralError(-3,3)/hrawprm->GetBinWidth(i+1));
    }
    TH2D* hBatPrMinus[nbins];
        for (int i = 0; i < nbins; i++)
        {
            hBatPrMinus[i] = new TH2D (Form("htpcnsigmapr%i",i), "; tpcnsigmapr; tofnsigmapr", 200, -10, 10, 200, -10, 10);
        }
 
        for (Long64_t i = 0; i < nEntries/10; ++i)
        {
            tree->GetEntry(i); // Load data for entry
            for (int j = 0; j < nbins; j++)
            {
                if (pt > ptbins[j] && pt < ptbins[j+1])
                hBatPrMinus[j]->Fill(tpcnsigmapr,tofnsigmaprm);
            }
            
        }
   
    
        for (int j = 0; j < nbins; j++)
        {
            TCanvas* cBinbinPrMinus = new TCanvas(Form("cBinbinPrMinus_%d", j), Form("pt bin %d", j), 800, 700);
            hBatPrMinus[j]->Draw("colz");
            cBinbinPrMinus->SaveAs(Form("binbinprminus%d.png", j));
        }




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

        hrawka->SetBinContent(i+1, 1/(4.05*1E+6)*fgauska[i]->Integral(-3,3)/hrawka->GetBinWidth(i+1));
        hrawka->SetBinError(i+1, 1/(4.05*1E+6)*fgauska[i]->IntegralError(-3,3)/hrawka->GetBinWidth(i+1));

    }
    TH2D* hBatKaPlus[nbins];
        for (int i = 0; i < nbins; i++)
        {
            hBatKaPlus[i] = new TH2D (Form("htpcnsigmaka%i",i), "; tpcnsigmaka; tofnsigmaka", 200, -10, 10, 200, -10, 10);
        }
 
        for (Long64_t i = 0; i < nEntries/10; ++i)
        {
            tree->GetEntry(i); // Load data for entry
            for (int j = 0; j < nbins; j++)
            {
                if (pt > ptbins[j] && pt < ptbins[j+1])
                hBatKaPlus[j]->Fill(tpcnsigmaka,tofnsigmaka);
            }
            
        }
   
    
        for (int j = 0; j < nbins; j++)
        {
            TCanvas* cBinbinKaPlus = new TCanvas(Form("cBinbinKaPlus_%d", j), Form("pt bin %d", j), 800, 700);
            hBatKaPlus[j]->Draw("colz");
            cBinbinKaPlus->SaveAs(Form("binbinkaplus%d.png", j));
        }



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

        hrawkam->SetBinContent(i+1, 1/(4.05*1E+6)*fgauskam[i]->Integral(-3,3)/hrawkam->GetBinWidth(i+1));
        hrawkam->SetBinError(i+1, 1/(4.05*1E+6)*fgauskam[i]->IntegralError(-3,3)/hrawkam->GetBinWidth(i+1));
    }
    TCanvas* crawkam= new TCanvas("hrawkam","hrawkam",900,1000);
    hrawkam->Draw();
    crawkam->SaveAs("rawkam.png");
    TH2D* hBatKaMinus[nbins];
        for (int i = 0; i < nbins; i++)
        {
            hBatKaMinus[i] = new TH2D (Form("htpcnsigmaka%i",i), "; tpcnsigmaka; tofnsigmaka", 200, -10, 10, 200, -10, 10);
        }
 
        for (Long64_t i = 0; i < nEntries/10; ++i)
        {
            tree->GetEntry(i); // Load data for entry
            for (int j = 0; j < nbins; j++)
            {
                if (pt > ptbins[j] && pt < ptbins[j+1])
                hBatKaMinus[j]->Fill(tpcnsigmaka,tofnsigmakam);
            }
            
        }
   
    
        for (int j = 0; j < nbins; j++)
        {
            TCanvas* cBinbinKaMinus = new TCanvas(Form("cBinbinKaMinus_%d", j), Form("pt bin %d", j), 800, 700);
            hBatKaMinus[j]->Draw("colz");
            cBinbinKaMinus->SaveAs(Form("binbinkaminus%d.png", j));
        }

    // Create a canvas named "cpt" with title "cpt" and size 1000x900 pixels
TCanvas* cpt = new TCanvas("cpt","cpt",1000,900);
// Draw histogram hPt on the canvas
hPt->Draw();
// Save the canvas as an image file "hpt.png"
cpt->SaveAs("hpt.png");

// Create a canvas named "ceta" with title "ceta" and size 1000x900 pixels
TCanvas* ceta = new TCanvas("ceta","ceta",1000,900);
// Draw histogram hEta on the canvas
hEta->Draw();

// Create a canvas named "p" with title "beta" and size 1000x900 pixels
TCanvas* cp = new TCanvas("p","beta",1000, 900);
// Draw histogram hP on the canvas
hP->Draw();

// Create a canvas for Energy Loss histogram with log scale on z-axis
TCanvas* c2energyloss = new TCanvas("Energy Loss","Energy Loss", 1000, 900);
c2energyloss->SetLogz();                  // Enable logarithmic scale on z-axis
h2energyloss->SetStats(0);                // Disable statistics box on histogram
h2energyloss->Draw("colz");               // Draw 2D histogram with color map
h2energyloss->SaveAs("Energy Loss.jpg"); // Save as image file

// Create a canvas for Beta histogram with log scale on z-axis
TCanvas* c2beta = new TCanvas("beta","Beta", 1000, 900);
c2beta->SetLogz();                        // Enable logarithmic scale on z-axis
h2beta->SetStats(0);                      // Disable statistics box
h2beta->Draw("colz");                     // Draw 2D histogram with color map
h2beta->SaveAs("Beta.jpg");               // Save as image file

// Create a canvas for Y projections
TCanvas* cprojy = new TCanvas("projection y","projection y", 1000, 900);
hprojy->SetStats(0);                      // Disable statistics box
hprojy->SetLineColor(kBlue);              // Set line color to blue for first histogram
hprojy2->SetLineColor(kRed);              // Set line color to red for second histogram
hprojy->Scale(1./hprojy->Integral());    // Normalize first histogram to area 1
hprojy2->Scale(1./hprojy2->Integral());  // Normalize second histogram
hprojy->Draw("hist");                     // Draw first histogram as lines without error bars
hprojy2->Draw("same hist");               // Draw second histogram on same canvas

// Create a canvas for "pi" histogram with log scale on z-axis
TCanvas* c2pi = new TCanvas("pi","pi", 1000, 900);
c2pi->SetLogz();                         // Enable log scale on z-axis
h2pi->SetStats(0);                       // Disable statistics box
h2pi->Draw();                           // Draw histogram
h2pi->SaveAs("Tpi.jpg");                // Save canvas as image

// Create a canvas for "pr" histogram with log scale on z-axis
TCanvas* c2pr = new TCanvas("pr","pr", 1000, 900);
c2pr->SetLogz();                        // Enable log scale on z-axis
h2pr->SetStats(0);                     // Disable statistics box
h2pr->Draw();                         // Draw histogram
h2pr->SaveAs("Tpr.jpg");              // Save canvas as image

// Create a canvas for "ka" histogram with log scale on z-axis
TCanvas* c2ka = new TCanvas("ka","ka", 1000, 900);
c2ka->SetLogz();                      // Enable log scale on z-axis
h2ka->SetStats(0);                   // Disable statistics box
h2ka->Draw();                       // Draw histogram
h2ka->SaveAs("Tka.jpg");            // Save canvas as image

// Create canvas for TPC signal for pions (pi)
TCanvas* c2sigpi= new TCanvas ("tpcsigpi","tpcsigpi", 1000, 900);
c2sigpi->SetLogz();                  // Enable log scale on z-axis
h2nsigpi->SetStats(0);              // Disable statistics box

// Define bin range for fitting slices in x-axis for pion histogram
const int min_binpi = h2nsigpi->GetXaxis()->FindBin(0.2);
const int max_binpi = h2nsigpi->GetXaxis()->FindBin(0.8);

// Create a Gaussian function to fit histogram slices
TF1* pigaus = new TF1("fgaus","gaus",-1.5,+1.5);

// Fit slices along y-axis of 2D histogram within bin range, store results in slicespi array
TObjArray* slicespi = new TObjArray();
h2nsigpi->FitSlicesY(pigaus, min_binpi, max_binpi, 0, "WW,QNR", slicespi);

// Extract sigma and mean histograms from fit results
TH1D* hSigmapi = (TH1D*)slicespi->At(2);
hSigmapi->SetName("hSigmapi");
TH1D* hMeanpi = (TH1D*)slicespi->At(1);

// Style the sigma and mean histograms for better visibility
hSigmapi->SetLineColor(2);
hMeanpi->SetLineColor(4);
hSigmapi->SetLineWidth(3);
hMeanpi->SetLineWidth(3);

// Limit the x-axis range for all related histograms
h2nsigpi->GetXaxis()->SetRangeUser(0,2);

// Draw the original 2D histogram and overlay sigma and mean histograms
h2nsigpi->Draw();
hSigmapi->Draw("same");
hMeanpi->Draw("same");

// Save the canvas as an image file
c2sigpi->SaveAs("nsigpi.png");

// Create canvas for raw pion data and draw it
TCanvas* crawpi= new TCanvas("crawpi","crawpi", 1000, 900);
hrawpi->Draw();
crawpi->SaveAs("raw pi.png");

// Repeat similar steps for proton (pr) TPC signal analysis
TCanvas* c2sigpr= new TCanvas ("tpcsigpr","tpcsigpr", 1000, 900);
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
hSigmapr->SetLineColor(2);
hMeanpr->SetLineColor(4);
hSigmapr->SetLineWidth(3);
hMeanpr->SetLineWidth(3);
h2nsigpr->GetXaxis()->SetRangeUser(0,2);
h2nsigpr->Draw();
hSigmapr->Draw("same");
hMeanpr->Draw("same");
c2sigpr->SaveAs("nsigpr.png");

// Canvas for raw proton data
TCanvas* crawpr= new TCanvas("crawpr","crawpr", 1000, 900);
hrawpr->Draw();
crawpr->SaveAs("raw pr.png");

// Repeat similar steps for kaon (ka) TPC signal analysis
TCanvas* c2sigka= new TCanvas ("tpcsigka","tpcsigka", 1000, 900);
c2sigka->SetLogz();
h2nsigka->SetStats(0);
const int min_binka = h2nsigka->GetXaxis()->FindBin(0.2);
const int max_binka = h2nsigka->GetXaxis()->FindBin(0.6);
TF1* kagaus = new TF1("fgaus","gaus",-2,2);
TObjArray* sliceska = new TObjArray();
h2nsigka->FitSlicesY(kagaus, min_binka, max_binka, 0, "WW,QNR", sliceska);
TH1D* hSigmaka = (TH1D*)sliceska->At(2);
hSigmapr->SetName("hSigmapr");  // This line likely a typo, should be hSigmaka->SetName(...)
TH1D* hMeanka = (TH1D*)sliceska->At(1);
hSigmaka->SetLineColor(2);
hMeanka->SetLineColor(3);
hSigmaka->SetLineWidth(3);
hMeanka->SetLineWidth(3);
h2nsigka->GetXaxis()->SetRangeUser(0,2);
h2nsigka->Draw();
hSigmaka->Draw("same");
hMeanka->Draw("same");
c2sigka->SaveAs("nsigka.png");

// Canvas for raw kaon data, draw two histograms on same canvas
TCanvas* crawka= new TCanvas("crawka","crawka", 1000, 900);
hrawka->Draw();
hrawkam->Draw("same");
crawka->SaveAs("raw ka.png");

// Canvas for TOF signal for pion
TCanvas* c2nsigtofpi= new TCanvas ("tofsigpi","tofsigpi", 1000, 900);
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
hSigmapitof->SetLineColor(2);
hMeanpitof->SetLineColor(4);
hSigmapitof->SetLineWidth(3);
hMeanpitof->SetLineWidth(3);
h2nsigtofpi->GetXaxis()->SetRangeUser(0,2);
h2nsigtofpi->Draw();
hSigmapitof->Draw("same");
hMeanpitof->Draw("same");
hMeanpitof->SetName("hMeanpi");
c2nsigtofpi->SaveAs("nsigtofpi.png");

// Canvas for TOF signal for proton
TCanvas* c2nsigtofpr= new TCanvas ("tofsigpr","tofsigpr", 1000, 900);
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
hSigmaprtof->SetLineColor(2);
hMeanprtof->SetLineColor(4);
hSigmaprtof->SetLineWidth(3);
hMeanprtof->SetLineWidth(3);
h2nsigtofpr->GetXaxis()->SetRangeUser(0,2);
h2nsigtofpr->Draw();
hSigmaprtof->Draw("same");
hMeanprtof->Draw("same");
hMeanprtof->SetName("hMeanpr");
c2nsigtofpr->SaveAs("nsigtofpr.png");

// Canvas for TOF signal for kaon
TCanvas* c2nsigtofka= new TCanvas ("tofsigka","tofsigka", 1000, 900);
h2nsigtofka->SetStats(0);
h2nsigtofka->Draw("colz");
const int min_binkatof = h2nsigtofka->GetXaxis()->FindBin(0.4);
const int max_binkatof = h2nsigtofka->GetXaxis()->FindBin(2);
TF1* kagaustof = new TF1("fgaus","gaus",-2,2);
TObjArray* sliceskatof = new TObjArray();
h2nsigtofka->FitSlicesY(kagaustof, min_binkatof, max_binkatof, 0, "WW,QNR", sliceskatof);
TH1D* hSigmakatof = (TH1D*)sliceskatof->At(2);
hSigmaprtof->SetName("hSigmapr");   // Probably should be hSigmakatof->SetName(...)
TH1D* hMeankatof = (TH1D*)sliceskatof->At(1);
hSigmakatof->SetLineColor(2);
hMeankatof->SetLineColor(3);
hSigmakatof->SetLineWidth(3);
hMeankatof->SetLineWidth(3);
h2nsigtofka->GetXaxis()->SetRangeUser(0,2);
h2nsigtofka->Draw();
hSigmakatof->Draw("same");
hMeankatof->Draw("same");
hMeankatof->SetName("hMeanpr");
c2nsigtofka->SaveAs("nsigtofka.png");

// Create a ROOT file named "raw spectra.root" to save histograms
TFile* write= new TFile ("raw spectra.root","RECREATE");
// Save raw pion histograms to the file
hrawpi->Write();
hrawpim->Write();
// Save raw proton histograms to the file
hrawpr->Write();
hrawprm->Write();
// Save raw kaon histograms to the file
hrawka->Write();
hrawkam->Write();

}  // End of function or scope

