#include "myAnaConsts.h"
#include "massfitting.C"
#include <iostream>
#include <memory>

TH1D* hDcaMCAll(TH3*, TH1&, const double&, const std::string&);
void hDcaData(TH3*, std::map<std::string, TH3*>, std::map<std::string, TH3*>, TH1&, const double&, const std::string&);
void dcaHistsOptimal(int mode = 1)
{
   //TGaxis::SetMaxDigits(3);
   gStyle->SetOptStat(0);
   //TFile* f1 = new TFile(Form("%s_hists.root", ana::whichtree[mode].c_str()));
   std::unique_ptr<TFile> f1 = std::unique_ptr<TFile>(new TFile(Form("%s_hists.root", ana::whichtree[mode].c_str())));
   std::map<std::string, TH3*> hDcaVsMassAndMvaPD0;
   std::map<std::string, TH3*> hDcaVsMassAndMvaNPD0;
   hDcaVsMassAndMvaPD0["h_match_unswap"] = (TH3D*) f1->Get("hDcaVsMassAndMvaPD0");
   hDcaVsMassAndMvaPD0["h_match_all"] = (TH3D*) f1->Get("hDcaVsMassAndMvaPD0_All");
   hDcaVsMassAndMvaNPD0["h_match_unswap"] = (TH3D*) f1->Get("hDcaVsMassAndMvaNPD0");
   hDcaVsMassAndMvaNPD0["h_match_all"] = (TH3D*) f1->Get("hDcaVsMassAndMvaNPD0_All");
   TH3D* hDcaVsMassAndMvaDataD0 = (TH3D*) f1->Get("hDcaVsMassAndMvaDataD0");

   TFile f2(Form("%s_dca_hists.root", ana::whichtree[mode].c_str()), "recreate");

   // create Pic to store pictures
   const std::string dirPic = "if [ ! -d \"Pic\" ]; then\n"
                              "    mkdir Pic \n"
                              "fi";
   gSystem->Exec(dirPic.c_str());

   //for(int iMva=0; iMva<15; iMva++){
   for(int iMva=0; iMva<9; iMva++){

      //create Pic/mva to store pictures
      std::string dirMvaPic = Form( "if [ ! -d \"Pic/mva%d\" ]; then\n"
                                 "    mkdir Pic/mva%d \n"
                                 "fi", iMva, iMva);
      gSystem->Exec(dirMvaPic.c_str());

      // DCA distribution of Prompt D0 
      double mvaCut = (double)iMva*ana::mvaStep + ana::mvaMin;
  
      TCanvas cPrompt("cPrompt", "", 550, 450);
      cPrompt.SetLogy();
      TH1D hDcaMCPD0("hDcaMCPD0", "hDCaMCPD0", ana::nuofDca, ana::dcaBin);
      TH1D* hDcaMCPD0_Proj = hDcaMCAll(hDcaVsMassAndMvaPD0["h_match_all"], hDcaMCPD0, mvaCut, "hprompt");
      hDcaMCPD0.Scale(1./hDcaMCPD0.Integral("width"));
      hDcaMCPD0.Draw();
      hDcaMCPD0_Proj->Scale(hDcaMCPD0.GetMaximum()/hDcaMCPD0_Proj->GetMaximum());
      hDcaMCPD0_Proj->SetLineColor(kRed);
      hDcaMCPD0_Proj->Draw("same");
      cPrompt.Print(Form("Pic/mva%d/cPrompt.png", iMva));
      delete hDcaMCPD0_Proj;

      // DCA distribution of Non-Prompt D0

      TCanvas cNonPrompt("cNonPrompt", "", 550, 450);
      cNonPrompt.SetLogy();
      TH1D hDcaMCNPD0("hDcaMCNPD0", "hDCaMCNPD0", ana::nuofDca, ana::dcaBin);
      TH1D* hDcaMCNPD0_Proj = hDcaMCAll(hDcaVsMassAndMvaNPD0["h_match_all"], hDcaMCNPD0, mvaCut, "hnonprompt");
      hDcaMCNPD0.Scale(1./hDcaMCNPD0.Integral("width"));
      hDcaMCNPD0.Draw();
      hDcaMCNPD0_Proj->Scale(hDcaMCNPD0.GetMaximum()/hDcaMCNPD0_Proj->GetMaximum());
      hDcaMCNPD0_Proj->SetLineColor(kRed);
      hDcaMCNPD0_Proj->Draw("same");
      cNonPrompt.Print(Form("Pic/mva%d/cNonPrompt.png", iMva));
      delete hDcaMCNPD0_Proj;

   // DCA distribution of signal + swap, by fitting mass;
      TCanvas cData("cData", "", 550, 450);
      TH1D hDcaDataD0("hDcaDataD0", "hDcaDataD0", ana::nuofDca, ana::dcaBin);
      std::string tmpName(Form("Pic/mva%d/massPerDcaBin", iMva));
      hDcaData(hDcaVsMassAndMvaDataD0, hDcaVsMassAndMvaNPD0, hDcaVsMassAndMvaPD0, hDcaDataD0, mvaCut, tmpName);
      cData.cd();
      cData.SetLogy();
      hDcaDataD0.Scale(1./hDcaDataD0.Integral("width"));
      hDcaDataD0.Draw();
      cData.Print(Form("Pic/mva%d/cData.png", iMva));
      
      // extract invariant mass
      int mvaBinMin = hDcaVsMassAndMvaNPD0["h_match_all"]->GetYaxis()->FindBin(mvaCut+0.1*ana::mvaStep); // 0.1 offset to make sure the correct bin is returned
      int mvaBinMax = hDcaVsMassAndMvaNPD0["h_match_all"]->GetYaxis()->GetNbins()+1;
      int label = mvaCut*100;
      TH1D* hMassNPD0 = hDcaVsMassAndMvaNPD0["h_match_unswap"]->ProjectionX(Form("hMassMCNPD0mva%d", label), mvaBinMin, mvaBinMax);
      TH1D* hMassNPD0All = hDcaVsMassAndMvaNPD0["h_match_all"]->ProjectionX(Form("hMassMCNPD0Allmva%d", label), mvaBinMin, mvaBinMax);
      TH1D* hMassPD0 = hDcaVsMassAndMvaPD0["h_match_unswap"]->ProjectionX(Form("hMassMCPD0mva%d", label), mvaBinMin, mvaBinMax);
      TH1D* hMassPD0All = hDcaVsMassAndMvaPD0["h_match_all"]->ProjectionX(Form("hMassMCPD0Allmva%d", label), mvaBinMin, mvaBinMax);
      TH1D* hMassData = hDcaVsMassAndMvaDataD0->ProjectionX(Form("hMassDataD0mva%d", label), mvaBinMin, mvaBinMax);

      f2.cd();

      hDcaDataD0.Write(Form("hDcaDataD0mva%d", label));
      hDcaMCNPD0.Write(Form("hDcaMCNPD0mva%d", label));
      hDcaMCPD0.Write(Form("hDcaMCPD0mva%d", label));

      hMassNPD0->Write();
      hMassNPD0->Delete();
      hMassNPD0All->Write();
      hMassNPD0All->Delete();
      hMassPD0->Write();
      hMassPD0->Delete();
      hMassPD0All->Write();
      hMassPD0All->Delete();
      hMassData->Write();
      hMassData->Delete();
   }
}

TH1D* hDcaMCAll(TH3* hDcaVsMassAndMva, TH1& hDcaMC, const double& mvaCut, const  std::string& proj_name){
   int mvaBinMin = hDcaVsMassAndMva->GetYaxis()->FindBin(mvaCut+0.1*ana::mvaStep); // 0.1 offset to make sure the correct bin is returned
   int mvaBinMax = hDcaVsMassAndMva->GetYaxis()->GetNbins()+1;
   TH1D* hDcaMCProj = hDcaVsMassAndMva->ProjectionZ(proj_name.c_str(), 1, 60, mvaBinMin, mvaBinMax); // x: mass, y: dca 
   for(int iDca=0; iDca<ana::nuofDca; iDca++){
      int binlw = hDcaMCProj->FindBin(ana::dcaBin[iDca] + hDcaMCProj->GetBinWidth(1)/2.); // binwidth/2 offset to make sure the correct bin is returned
      int binup = hDcaMCProj->FindBin(ana::dcaBin[iDca+1] + hDcaMCProj->GetBinWidth(1)/2.) - 1;  // 1 bin decreased due to the integral is taken within [binlw, binup]
      //std::cout << binlw << std::endl;
      //std::cout << ana::dcaBin[iDca] << std::endl;
      //std::cout << binup << std::endl;
      //std::cout << ana::dcaBin[iDca+1] << std::endl;
      double binwidth = ana::dcaBin[iDca+1] - ana::dcaBin[iDca];
      //std::cout << binwidth << std::endl;
      double err = 0.;
      double yield = hDcaMCProj->IntegralAndError(binlw, binup, err);
      hDcaMC.SetBinContent(iDca+1, yield/binwidth);
      hDcaMC.SetBinError(iDca+1, err/binwidth);
   }
   
   return hDcaMCProj;
}

void hDcaData(TH3* hData, std::map<std::string, TH3*> hNP, std::map<std::string, TH3*> hP, TH1& hDataD0, const double& mvaCut, const std::string& tmpName)
{
   int mvaBinMin = hP["h_match_unswap"]->GetYaxis()->FindBin(mvaCut+0.1*ana::mvaStep); // 0.1 offset to make sure the correct bin is returned
   int mvaBinMax = hP["h_match_unswap"]->GetYaxis()->GetNbins()+1; // overflow of y axis, mva
   int mcDcaBinMin = 0; // underflow of the z axis, the DCA
   int mcDcaBinMax = hP["h_match_unswap"]->GetZaxis()->GetNbins()+1; // overflow of z axis, the DCA

   // x: mass, y: mva, z: dca
   TH1D* hMCPromptMass = hP["h_match_unswap"]->ProjectionX("hPromptMass", mvaBinMin, mvaBinMax, mcDcaBinMin, mcDcaBinMax); // projection, y-mva:[cut, overflow], z-DCA:[underflow, overflow]
   TH1D* hMCPromptMassAll = hP["h_match_all"]->ProjectionX("hPromptMassAll", mvaBinMin, mvaBinMax, mcDcaBinMin, mcDcaBinMax);
   TH1D* hMCNonPromptMass = hNP["h_match_unswap"]->ProjectionX("hNonPromptMass", mvaBinMin, mvaBinMax, mcDcaBinMin, mcDcaBinMax);
   TH1D* hMCNonPromptMassAll = hNP["h_match_all"]->ProjectionX("hNonPromptMassAll", mvaBinMin, mvaBinMax, mcDcaBinMin, mcDcaBinMax);

   TH1D* hTemp;
   TCanvas cTemp("cTemp", "", 550, 450);
   for(int iDca=0; iDca<ana::nuofDca; iDca++){
      int binlw = hData->GetZaxis()->FindBin(ana::dcaBin[iDca] + hData->GetZaxis()->GetBinWidth(1) * 0.2); // 0.2 offset to make sure the correct bin is returned
      int binup = hData->GetZaxis()->FindBin(ana::dcaBin[iDca+1] + hData->GetZaxis()->GetBinWidth(1) * 0.2) -1; // -1 due to the integral is taken within [binlw, binup]
      hTemp = hData->ProjectionX("tmp", mvaBinMin, mvaBinMax, binlw, binup); // projection, y-mva:[cut, overflow], z-DCA:[dca[iDca], dca[iDca+1]]
      hTemp->GetYaxis()->SetRangeUser(0, hTemp->GetMaximum() * 1.3);
      TF1 f = massfitting(hTemp, hMCNonPromptMass, hMCNonPromptMassAll, "iDca");
      TF1 signal("signal", "[0]* [5] * (" "[4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))"
            "+ (1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))" ")", 1.7, 2.0);
      signal.SetLineColor(kOrange-3);
      signal.SetLineWidth(1);
      signal.SetLineStyle(2);
      signal.SetFillColorAlpha(kOrange-3,0.3);
      signal.SetFillStyle(1001);
      for(int ipar=0; ipar<6+1; ipar++){
         signal.FixParameter(ipar, f.GetParameter(ipar));
      }

      TF1 swap("swap", "[0]*((1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])))"
            "+0 *[1]*[2]*[3]*[4]", 1.7, 2.0);
      swap.SetLineColor(kGreen+4);
      swap.SetLineWidth(1);
      swap.SetLineStyle(1);
      swap.SetFillColorAlpha(kGreen+4,0.3);
      swap.SetFillStyle(1001);
      for(int ipar=0; ipar<8+1; ipar++){
         swap.FixParameter(ipar, f.GetParameter(ipar));
      }

      TF1 bkg("bkg", "[9] + [10]*x + [11]*x*x + [12]*x*x*x"
            "+ 0 *[0]*[1]*[2]*[3]*[4]*[5]*[6]*[7]*[8]", 1.7, 2.0);
      bkg.SetLineColor(4);
      bkg.SetLineWidth(1);
      bkg.SetLineStyle(2);
      for(int ipar=0; ipar<12+1; ipar++){
         bkg.FixParameter(ipar, f.GetParameter(ipar));
      }

      signal.Draw("same FC");
      swap.Draw("same FC");
      bkg.Draw("same L");
      cTemp.Print(Form("%s%d.png", tmpName.c_str(), iDca));

      double binWidth = hTemp->GetBinWidth(1) * (ana::dcaBin[iDca+1] - ana::dcaBin[iDca]);
      hDataD0.SetBinContent(iDca+1, f.GetParameter(0)/binWidth);
      hDataD0.SetBinError(iDca+1, f.GetParError(0)/binWidth);

      // clear pointers
      hTemp->Delete();

   }
   hMCPromptMass->Delete(); 
   hMCPromptMass = nullptr;
   hMCPromptMassAll->Delete();
   hMCPromptMassAll = nullptr;
   hMCNonPromptMass->Delete();
   hMCNonPromptMass = nullptr;
   hMCNonPromptMassAll->Delete();
   hMCNonPromptMassAll = nullptr;
}