#include "myAnaConsts.h"
#include "massfitting.C"
#include <iostream>
#include <memory>

void massHistsOptimal(int mode = 0)
{
   TH1::SetDefaultSumw2();

   //TGaxis::SetMaxDigits(3);
   gStyle->SetOptStat(0);
   //TFile* f1 = new TFile(Form("%s_hists.root", ana::whichtree[mode].c_str()));
   std::unique_ptr<TFile> f1 = std::unique_ptr<TFile>(new TFile(Form("hists/%s_hists_%s/%s_hists_pT%.1f-%.1f_y%.1f-%.1f.root",ana::whichtree[mode].c_str(),ana::whichbk[ana::bk_mode].c_str(),ana::whichtree[mode].c_str(),ana::pTMin,ana::pTMax,ana::yMin,ana::yMax)));

   std::cout << Form("ready to write ouput file: hists/%s_hists_%s/%s_hists_pT%.1f-%.1f_y%.1f-%.1f.root",ana::whichtree[mode].c_str(),ana::whichbk[ana::bk_mode].c_str(), ana::whichtree[mode].c_str(),ana::pTMin,ana::pTMax,ana::yMin,ana::yMax) << std::endl;
   std::map<std::string, TH3*> hDcaVsMassAndMvaPD0;
   std::map<std::string, TH3*> hDcaVsMassAndMvaNPD0;

   //char hist_promptDCA_unswap[128];
   //char hist_promptDCA_all[128];
   //char hist_nonPromptDCA_unswap[128];
   //char hist_nonPromptDCA_all[128];
   //char hist_data[128];

   std::map<std::string,int> ptbin_index; 
   ptbin_index.insert(std::pair<std::string,int>("1.5-2.2",1)); 
   ptbin_index.insert(std::pair<std::string,int>("2.2-3.0",2)); 
   ptbin_index.insert(std::pair<std::string,int>("3.0-4.0",3)); 
   ptbin_index.insert(std::pair<std::string,int>("4.0-5.0",4)); 
   ptbin_index.insert(std::pair<std::string,int>("5.0-6.0",5)); 
   ptbin_index.insert(std::pair<std::string,int>("6.0-7.0",6)); 
   ptbin_index.insert(std::pair<std::string,int>("7.0-8.0",7)); 
   ptbin_index.insert(std::pair<std::string,int>("8.0-10.0",8)); 
   ptbin_index.insert(std::pair<std::string,int>("10.0-20.0",9)); 

   std::map<std::string, int> ybin_index; 
   ybin_index.insert(pair<std::string,int>("0.0-0.5",1)); 
   ybin_index.insert(pair<std::string,int>("0.0-1.5",1)); 
   ybin_index.insert(pair<std::string,int>("0.0-2.0",1)); 


   //char pTBinRange[128];   
   //char yBinRange[128];   

   std::string pTBinRange(TString::Format("%.1f-%.1f", ana::pTMin, ana::pTMax)); 
   std::string yBinRange(TString::Format("%.1f-%.1f", ana::yMin, ana::yMax)); 

   cout << "pTBinRange" << ":" << pTBinRange.c_str()<< endl;    
   cout << "yBinRange" << ":" << yBinRange.c_str() << endl;    

   std::string hist_promptDCA_unswap(TString::Format("hDcaVsMassAndMvaPD0_%d_%d",ptbin_index[pTBinRange.c_str()],ybin_index[yBinRange.c_str()]));
   std::string hist_promptDCA_all(TString::Format("hDcaVsMassAndMvaPD0_All_%d_%d",ptbin_index[pTBinRange.c_str()],ybin_index[yBinRange.c_str()]));
   std::string hist_nonPromptDCA_unswap(TString::Format("hDcaVsMassAndMvaNPD0_%d_%d",ptbin_index[pTBinRange.c_str()],ybin_index[yBinRange.c_str()]));
   std::string hist_nonPromptDCA_all(TString::Format("hDcaVsMassAndMvaNPD0_All_%d_%d",ptbin_index[pTBinRange.c_str()],ybin_index[yBinRange.c_str()]));
   std::string hist_data(TString::Format("hDcaVsMassAndMvaDataD0_%d_%d",ptbin_index[pTBinRange.c_str()],ybin_index[yBinRange.c_str()]));

   cout << "hist_promptDCA_unswap" << ":" << hist_promptDCA_unswap.c_str() << endl;
   cout << "hist_promptDCA_all" << ":" << hist_promptDCA_all.c_str() << endl;
   cout << "hist_nonPromptDCA_unswap" << ":" << hist_nonPromptDCA_unswap.c_str() << endl;
   cout << "hist_nonPromptDCA_all" << ":" << hist_nonPromptDCA_all.c_str() << endl;
   cout << "hist_data" << ":" << hist_data.c_str() << endl;
   
   hDcaVsMassAndMvaPD0["h_match_unswap"] = (TH3D*) f1->Get(hist_promptDCA_unswap.c_str());
   hDcaVsMassAndMvaPD0["h_match_all"] = (TH3D*) f1->Get(hist_promptDCA_all.c_str());
   hDcaVsMassAndMvaNPD0["h_match_unswap"] = (TH3D*) f1->Get(hist_nonPromptDCA_unswap.c_str());
   hDcaVsMassAndMvaNPD0["h_match_all"] = (TH3D*) f1->Get(hist_nonPromptDCA_all.c_str());
   TH3D* hDcaVsMassAndMvaDataD0 = (TH3D*) f1->Get(hist_data.c_str());
//++++++++++
 
 
   TFile f2(Form("%s_%s/%s_dca_hists_pT%.1f-%.1f_y%.1f-%.1f.root",ana::whichtree[mode].c_str(),ana::whichbk[ana::bk_mode].c_str(),ana::whichtree[mode].c_str(),ana::pTMin,ana::pTMax,ana::yMin,ana::yMax), "recreate");

   // create Pic to store pictures
   const std::string dirPic = Form("if [ ! -d \"%sPic\" ]; then\n"
                              "    mkdir %sPic \n"
                              "fi", ana::whichtree[mode].c_str(), ana::whichtree[mode].c_str());
   gSystem->Exec(dirPic.c_str());

   const std::string dirPicPtY = Form("if [ ! -d \"%sPic/%s_%s_Pic/pT%.1f-%.1f_y%.1f-%.1f\" ]; then\n"
                              "    mkdir %sPic/%s_%s_Pic/pT%.1f-%.1f_y%.1f-%.1f \n"
                              "fi", ana::whichtree[mode].c_str(),ana::whichtree[mode].c_str(),ana::whichbk[ana::bk_mode].c_str(),ana::pTMin,ana::pTMax,ana::yMin,ana::yMax,
                                    ana::whichtree[mode].c_str(),ana::whichtree[mode].c_str(),ana::whichbk[ana::bk_mode].c_str(),ana::pTMin,ana::pTMax,ana::yMin,ana::yMax);
   gSystem->Exec(dirPicPtY.c_str());

   double significance[100];
   double yields_signal[100];
   double yields_bkg[100];
   double mvaCuts[100];
   
   //To select the appropriate the number of iMva based on different ptbin, ybin. 
   std::map<std::string,int> nuofMva; 
   nuofMva.insert(std::pair<std::string,int>("1.5-2.2_0.0-0.5",50)); 
   nuofMva.insert(std::pair<std::string,int>("1.5-2.2_0.0-1.5",50)); 
   nuofMva.insert(std::pair<std::string,int>("1.5-2.2_0.0-2.0",50)); 
   nuofMva.insert(std::pair<std::string,int>("2.2-3.0_0.0-0.5",50)); 
   nuofMva.insert(std::pair<std::string,int>("2.2-3.0_0.0-1.5",48)); 
   nuofMva.insert(std::pair<std::string,int>("2.2-3.0_0.0-2.0",48)); 
   nuofMva.insert(std::pair<std::string,int>("3.0-4.0_0.0-0.5",50)); 
   nuofMva.insert(std::pair<std::string,int>("3.0-4.0_0.0-1.5",50)); 
   nuofMva.insert(std::pair<std::string,int>("3.0-4.0_0.0-2.0",48)); 
   nuofMva.insert(std::pair<std::string,int>("4.0-5.0_0.0-0.5",48)); 
   nuofMva.insert(std::pair<std::string,int>("4.0-5.0_0.0-1.5",48)); 
   nuofMva.insert(std::pair<std::string,int>("4.0-5.0_0.0-2.0",48)); 
   nuofMva.insert(std::pair<std::string,int>("5.0-6.0_0.0-0.5",48)); 
   nuofMva.insert(std::pair<std::string,int>("5.0-6.0_0.0-1.5",48)); 
   nuofMva.insert(std::pair<std::string,int>("5.0-6.0_0.0-2.0",48)); 
   nuofMva.insert(std::pair<std::string,int>("6.0-7.0_0.0-0.5",45)); 
   nuofMva.insert(std::pair<std::string,int>("6.0-7.0_0.0-1.5",48)); 
   nuofMva.insert(std::pair<std::string,int>("6.0-7.0_0.0-2.0",48)); 
   nuofMva.insert(std::pair<std::string,int>("7.0-8.0_0.0-0.5",48)); 
   nuofMva.insert(std::pair<std::string,int>("7.0-8.0_0.0-1.5",47)); 
   nuofMva.insert(std::pair<std::string,int>("7.0-8.0_0.0-2.0",48)); 
   nuofMva.insert(std::pair<std::string,int>("8.0-10.0_0.0-0.5",48)); 
   nuofMva.insert(std::pair<std::string,int>("8.0-10.0_0.0-1.5",48)); 
   nuofMva.insert(std::pair<std::string,int>("8.0-10.0_0.0-2.0",48)); 
   nuofMva.insert(std::pair<std::string,int>("10.0-20.0_0.0-0.5",48)); 
   nuofMva.insert(std::pair<std::string,int>("10.0-20.0_0.0-1.5",48)); 
   nuofMva.insert(std::pair<std::string,int>("10.0-20.0_0.0-2.0",48)); 
   
   //char pTBinRange_yBinRange[128]; 
   std::string pTBinRange_yBinRange(TString::Format("%.1f-%.1f_%.1f-%.1f",ana::pTMin,ana::pTMax,ana::yMin,ana::yMax));
   cout << "pTBinRange_yBinRange" << " " << ":" << " " << pTBinRange_yBinRange.c_str() << endl; 
   cout << "nuofMva[pTBinRange_yBinRange]" << " " <<":" << nuofMva[pTBinRange_yBinRange] << endl; 
   for(int iMva=0; iMva<nuofMva[pTBinRange_yBinRange]; iMva++){
	   //for(int iMva=2; iMva<4; iMva++){

      //create Pic/mva to store pictures
      std::string dirMvaPic(TString::Format( "if [ ! -d \"%sPic/%s_%s_Pic/pT%.1f-%.1f_y%.1f-%.1f/mva%d\" ]; then\n"
                                 "    mkdir %sPic/%s_%s_Pic/pT%.1f-%.1f_y%.1f-%.1f/mva%d \n"
                                 "fi", ana::whichtree[mode].c_str(),ana::whichtree[mode].c_str(),ana::whichbk[ana::bk_mode].c_str(),ana::pTMin,ana::pTMax,ana::yMin,ana::yMax, iMva, 
                                       ana::whichtree[mode].c_str(),ana::whichtree[mode].c_str(),ana::whichbk[ana::bk_mode].c_str(),ana::pTMin,ana::pTMax,ana::yMin,ana::yMax, iMva));
      gSystem->Exec(dirMvaPic.c_str());

      double mvaCut = (double)iMva*ana::mvaStep + ana::mvaMin;

      // extract invariant mass
      int mvaBinMin = hDcaVsMassAndMvaNPD0["h_match_all"]->GetYaxis()->FindBin(mvaCut+0.1*ana::mvaStep); // 0.1 offset to make sure the correct bin is returned
      cout << "mvaBinMin" << " = " << mvaBinMin << endl; 

      int mvaBinMax = hDcaVsMassAndMvaNPD0["h_match_all"]->GetYaxis()->GetNbins()+1;
      int label = mvaCut*100;
      TH1D* hMassNPD0 = hDcaVsMassAndMvaNPD0["h_match_unswap"]->ProjectionX(Form("hMassMCNPD0mva%d", label), mvaBinMin, mvaBinMax);
      TH1D* hMassNPD0All = hDcaVsMassAndMvaNPD0["h_match_all"]->ProjectionX(Form("hMassMCNPD0Allmva%d", label), mvaBinMin, mvaBinMax);
      TH1D* hMassPD0 = hDcaVsMassAndMvaPD0["h_match_unswap"]->ProjectionX(Form("hMassMCPD0mva%d", label), mvaBinMin, mvaBinMax);
      TH1D* hMassPD0All = hDcaVsMassAndMvaPD0["h_match_all"]->ProjectionX(Form("hMassMCPD0Allmva%d", label), mvaBinMin, mvaBinMax);
      TH1D* hMassData = hDcaVsMassAndMvaDataD0->ProjectionX(Form("hMassDataD0mva%d", label), mvaBinMin, mvaBinMax);

      //hMassNPD0->Scale(1.0/hMassNPD0->GetBinWidth(1));
      //hMassNPD0All->Scale(1.0/hMassNPD0All->GetBinWidth(1));
      //hMassPD0->Scale(1.0/hMassPD0->GetBinWidth(1));
      //hMassPD0All->Scale(1.0/hMassPD0All->GetBinWidth(1));
      //hMassData->Scale(1.0/hMassData->GetBinWidth(1));

      // TFitResultPtr constructed by std::shared_ptr<TFitResult>, be free of it even if you do not delete it
      // fitting the mass
      // rescale the mass
      
      TFitResultPtr fitResultPtr;
      TF1 fMass;
      if(mode == 0) fMass = massfitting(hMassData, hMassPD0, hMassPD0All, TString::Format("fMassMva%d", iMva), fitResultPtr);
      if(mode == 1) fMass = massfitting(hMassData, hMassNPD0, hMassNPD0All, TString::Format("fMassMva%d", iMva), fitResultPtr);
      if(mode == 2) fMass = massfitting(hMassData, hMassNPD0, hMassNPD0All, TString::Format("fMassMva%d", iMva), fitResultPtr);

      // draw the fitting
      drawMassFitting(hMassData, fMass, TString::Format("%sPic/%s_%s_Pic/pT%.1f-%.1f_y%.1f-%.1f/mva%d/MassFittingUnNormalized.png", ana::whichtree[mode].c_str(),ana::whichtree[mode].c_str(),ana::whichbk[ana::bk_mode].c_str(),ana::pTMin,ana::pTMax,ana::yMin,ana::yMax, iMva), 
            TString::Format("%.1f<p_{T}<%.1fGeV, %.1f<y<%.1f",ana::pTMin,ana::pTMax,ana::yMin,ana::yMax), TString::Format("MVA > %.2f", -0.3+0.02*iMva));

      // https://root.cern/doc/v616/TF1Helper_8cxx_source.html and https://root.cern.ch/doc/master/classTF1.html
      // help you understand how the integral error of TF1 is calculated
      // 1-dimensional function have five parameters, x_low, x_min, array of pars, array of covariance matrix, precision
      // array of pars would replace the original set of pars by last fitting unless no set of pars is passed
      // it is very important to pass in the correct set of pars and covariance matrix
      // explicitly pass the parameters to make sure anything won't go wrong
      // be aware that the fixed parameter would bring effect on the evaluation
      
      // get the covariance matrix
      auto covMat = fitResultPtr->GetCovarianceMatrix();

      // define the background function and then fix or set the parameters
      TF1 bkg("bkg", "[9] + [10]*x + [11]*x*x + [12]*x*x*x"
            "+ 0 *[0]*[1]*[2]*[3]*[4]*[5]*[6]*[7]*[8]", 1.7, 2.0);
      for(int ipar=0; ipar<8+1; ipar++){
         bkg.FixParameter(ipar, fMass.GetParameter(ipar));
      }
      for(int ipar=9; ipar<12+1; ipar++){
         bkg.SetParameter(ipar, fMass.GetParameter(ipar));
      }

      // define the signal function and then fix or set the parameters
      TF1 signal("signal", "[0]*[5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6])) + (1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+0 *[7]*[8]*[9]*[10]*[11]*[12]");
      for(int ipar=7; ipar<12+1; ipar++){
         signal.FixParameter(ipar, fMass.GetParameter(ipar));
      }
      for(int ipar=0; ipar<6+1; ipar++){
         signal.SetParameter(ipar, fMass.GetParameter(ipar));
      }

      // calculate the yield significance 
      std::map<std::string, double> yield_signal;
      std::map<std::string, double> yieldError_signal;
      std::map<std::string, double> yield_bkg;
      std::map<std::string, double> yieldError_bkg;

//      yield_signal["peak"] = signal.Integral(ana::peak_min, ana::peak_max);
//      yieldError_signal["peak"] = signal.IntegralError(ana::peak_min, ana::peak_max, 0, covMat.GetMatrixArray());
//      yield_bkg["peak"] = bkg.Integral(ana::peak_min, ana::peak_max);
//      yieldError_bkg["peak"] = bkg.IntegralError(ana::peak_min, ana::peak_max, 0, covMat.GetMatrixArray());
//
//      yields_signal[iMva] = yield_signal["peak"]; 
//      yields_bkg[iMva] = yield_bkg["peak"];
//      significance[iMva] = yield_signal["peak"] / sqrt(yield_signal["peak"] + yield_bkg["peak"]);
//

      significance[iMva] = fMass.GetParameter(0)/fMass.GetParError(0);
      yields_signal[iMva] = fMass.GetParameter(0);
      cout << "iMva" << "=" << " " << iMva << " " << "significance[iMva]" << "=" << significance[iMva] << endl;   
      mvaCuts[iMva] = mvaCut;

      f2.cd();
      fMass.Write();

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

   TCanvas* cc = new TCanvas("cc","cc",600,600);
   TGraph* gr = new TGraph(nuofMva[pTBinRange_yBinRange], mvaCuts, significance);
   gr->SetMarkerStyle(20);
   gr->Draw("AP");

   TCanvas* cc1 = new TCanvas("cc1","cc1",600,600);
   TGraph* gr1 = new TGraph(nuofMva[pTBinRange_yBinRange], mvaCuts, yields_signal);
   gr1->SetMarkerStyle(20);
   gr1->Draw("AP");

   gr1->Write("signalvsmva");
   gr->Write("significancevsmva");

}
