#include <iostream>
#include <map>
#include <vector>
#include <string>

#include "TSystem.h"
#include "TROOT.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TH3D.h"
#include "TMath.h"
#include "TTree.h"
#include "TFileCollection.h"
#include "TCollection.h"
#include "THashList.h"

#include "include/d0data.h"
#include "include/d0mc.h"
#include "include/myAnaConsts.h"

void loopTree(d0tree* , std::map<std::string, TH3*> d0hists[][20], bool isMC);
bool passKinematicalCut(d0tree*);

void readMVAtreeFull
(int mode = 0, std::string mcPD0List = "./list/MCPD0.list", std::string mcNPD0List = "./list/MCNPD0.list",
 std::string dataList = "./list/Data.list")
{
   TH1::SetDefaultSumw2();
   //std::map<std::string, TH1*> prompt;
   //prompt["dca3D"] = new TH1D("hDca3D_pd0", "hDca3D_pd0", 20, 0, 0.04);
   //prompt["pseudo_decayL"] = new TH1D("hDL_pd0", "hpseudoDL_pd0", 40, -0.12, 0.32);
   //std::map<std::string, TH1*> nonprompt;
   //nonprompt["dca3D"] = new TH1D("hDca3D_npd0", "hDca3D_npd0", 20, 0, 0.04);
   //nonprompt["pseudo_decayL"] = new TH1D("hDL_npd0", "hpseudoDL_npd0", 40, -0.12, 0.32);

   /*
   std::map<std::string, TH2*> promptDCA;
   promptDCA["h_match_unswap"] = new TH2D("hDcaVsMassPD0", "hDcaVsMassPD0", 60, 1.7, 2.0,  60, 0, 0.06);
   promptDCA["h_match_all"] = new TH2D("hDcaVsMassPD0_All", "hDcaVsMassPD0_All", 60, 1.7, 2.0,  60, 0, 0.06);

   std::map<std::string, TH2*> nonPromptDCA;
   nonPromptDCA["h_match_unswap"] = new TH2D("hDcaVsMassNPD0", "hDcaVsMassNPD0", 60, 1.7, 2.0,  60, 0, 0.06);
   nonPromptDCA["h_match_all"] = new TH2D("hDcaVsMassNPD0_All", "hDcaVsMassNPD0_All", 60, 1.7, 2.0,  60, 0, 0.06);

   std::map<std::string, TH2*> dataDCA;
   dataDCA["hdata"] = new TH2D("hDcaVsMassDataD0", "hDcaVsMassDataD0", 60, 1.7, 2.0, 60, 0., 0.06);
   */


   // declare histograms
   std::map<std::string, TH3*> promptDCA[20][20];
   std::map<std::string, TH3*> nonPromptDCA[20][20];
   std::map<std::string, TH3*> dataDCA[20][20];
   
   char histName_promptDCA_unswap[128]; 
   char histName_promptDCA_all[128]; 
   char histName_nonPromptDCA_unswap[128]; 
   char histName_nonPromptDCA_all[128]; 
   char histName_data[128];
 
   for(int ipt=0;ipt<ana::nuofpt;ipt++)
     for(int jy=0;jy<ana::nuofY;jy++)
     {  
        sprintf(histName_promptDCA_unswap,"hDcaVsMassAndMvaPD0_%d_%d",ipt+1,jy+1);
        sprintf(histName_promptDCA_all,"hDcaVsMassAndMvaPD0_All_%d_%d",ipt+1,jy+1);
        sprintf(histName_nonPromptDCA_unswap,"hDcaVsMassAndMvaNPD0_%d_%d",ipt+1,jy+1);
        sprintf(histName_nonPromptDCA_all,"hDcaVsMassAndMvaNPD0_All_%d_%d",ipt+1,jy+1);
        sprintf(histName_data,"hDcaVsMassAndMvaDataD0_%d_%d",ipt+1,jy+1);

       (promptDCA[ipt][jy])["h_match_unswap"] = new TH3D(histName_promptDCA_unswap, histName_promptDCA_unswap, 60, 1.7, 2.0, 100, -0.3, 0.7, ana::nDca, ana::dcaMin, ana::dcaMax);
       (promptDCA[ipt][jy])["h_match_all"] = new TH3D(histName_promptDCA_all,histName_promptDCA_all, 60, 1.7, 2.0, 100, -0.3, 0.7, ana::nDca, ana::dcaMin, ana::dcaMax);

       (nonPromptDCA[ipt][jy])["h_match_unswap"] = new TH3D(histName_nonPromptDCA_unswap,histName_nonPromptDCA_unswap, 60, 1.7, 2.0, 100, -0.3, 0.7, ana::nDca, ana::dcaMin, ana::dcaMax);
       (nonPromptDCA[ipt][jy])["h_match_all"] = new TH3D(histName_nonPromptDCA_all,histName_nonPromptDCA_all, 60, 1.7, 2.0, 100, -0.3, 0.7, ana::nDca, ana::dcaMin, ana::dcaMax);

       (dataDCA[ipt][jy])["hdata"] = new TH3D(histName_data,histName_data, 60, 1.7, 2.0, 100, -0.3, 0.7, ana::nDca, ana::dcaMin, ana::dcaMax);
     } 

   // read trees and fill hisotgrams

   std::string mcPD0ChainStr(
                   TString::Format("%s_mc/VertexCompositeNtuple", ana::whichtree[mode].c_str()));
   TChain* mcPD0Chain = new TChain(mcPD0ChainStr.c_str());
   TFileCollection* fcMCPD0 = new TFileCollection(mcPD0List.c_str(), "", mcPD0List.c_str());
   mcPD0Chain->AddFileInfoList(fcMCPD0->GetList());
   d0mc *mcPD0 = new d0mc(mcPD0Chain);

   std::cout << "started filling histograms of mc promptd0" << std::endl;
   loopTree(mcPD0, promptDCA, true);
   std::cout << "ended filling histograms of mc promptd0" << std::endl;

   delete mcPD0;


   std::string mcNPD0ChainStr(
                   TString::Format("%s_mc/VertexCompositeNtuple", ana::whichtree[mode].c_str()));
   TChain* mcNPD0Chain = new TChain(mcNPD0ChainStr.c_str());
   TFileCollection* fcMCNPD0 = new TFileCollection(mcNPD0List.c_str(), "", mcNPD0List.c_str());
   mcNPD0Chain->AddFileInfoList(fcMCNPD0->GetList());
   d0mc *mcNPD0 = new d0mc(mcNPD0Chain);

   std::cout << "started filling histograms of mc non-prompt d0" << std::endl;
   loopTree(mcNPD0, nonPromptDCA, true);
   std::cout << "ended filling histograms of mc non-prompt d0" << std::endl;

   delete mcNPD0;

   std::string dataD0ChainStr = 
                   Form("%s/VertexCompositeNtuple", ana::whichtree[mode].c_str());
   TChain* dataD0Chain= new TChain(dataD0ChainStr.c_str());
   TFileCollection* fcDataD0 = new TFileCollection(dataList.c_str(), "", dataList.c_str());
   dataD0Chain->AddFileInfoList(fcDataD0->GetList());
   d0data *dataD0 = new d0data(dataD0Chain);


   std::cout << "started filling histograms of data d0" << std::endl;
   loopTree(dataD0, dataDCA, false);
   std::cout << "ended filling histograms of data d0" << std::endl;

   delete dataD0;

   
   const std::string dirHists = Form("if [ ! -d \"hists\" ]; then\n"
                              "    mkdir hists \n"
                              "fi");
   gSystem->Exec(dirHists.c_str());

   for(int ipt=0;ipt<ana::nuofpt;ipt++)
     for(int jy=0;jy<ana::nuofY;jy++)
     {
       TFile f4(Form("hists_test/%s_hists_%s/%s_hists_pT%.1f-%.1f_y%.1f-%.1f.root",ana::whichtree[mode].c_str(),ana::whichbk[ana::bk_mode].c_str(),ana::whichtree[mode].c_str(),ana::ptbin[ipt],ana::ptbin[ipt+1],ana::ybin[jy],ana::ybin[jy+1]), "recreate");
       std::cout << Form("ready to write ouput file: %s_hists_%s/%s_hists_pT%.1f-%.1f_y%.1f-%.1f.root",ana::whichtree[mode].c_str(),ana::whichbk[ana::bk_mode].c_str(), ana::whichtree[mode].c_str(),ana::ptbin[ipt],ana::ptbin[ipt+1],ana::ybin[jy],ana::ybin[jy+1]) << std::endl;
       f4.cd();

       for(const auto& h : promptDCA[ipt][jy]) {h.second->Write(); delete h.second;}
       for(const auto& h : nonPromptDCA[ipt][jy]) {h.second->Write(); delete h.second;}
      // for(const auto& h : dataDCA[ipt][jy]) {h.second->Write(); delete h.second;}
     }

   std::cout << "ended writing output file" << std::endl;  
}

void loopTree(d0tree* d0, std::map<std::string, TH3*> d0hists[][20], bool isMC)
{
   int nentries = d0->GetEntries();
   std::cout << "total entries: " << nentries << std::endl;
   //nentries = 10000;
   for(int ientry=0; ientry<nentries; ientry++){
      d0->GetEntry(ientry);
      if(!passKinematicalCut(d0)) continue;
      int ptbin = whichPt(d0->Pt());
      int ybin = whichY(d0->Y());
      if(ptbin==-1 || ybin==-1) continue;

      double dca3D = d0->DecayL3D() * std::sin(d0->PointingAngle3D());
      double mass = d0->Mass();
      double mva = d0->Mva();
      if(isMC){
         if(!d0->MatchGEN()) continue;
         if(!d0->IsSwap()) (d0hists[ptbin][ybin])["h_match_unswap"]->Fill(mass, mva, dca3D);
         (d0hists[ptbin][ybin])["h_match_all"]->Fill(mass, mva, dca3D);
      } else {
         (d0hists[ptbin][ybin])["hdata"]->Fill(mass, mva, dca3D);
      }
   }
}

bool passKinematicalCut(d0tree* d0)
{
   bool passPointingAngle = std::fabs(d0->PointingAngle3D()) < 1;
   //passPointingAngle = true;
   bool passTrkEta = std::fabs(d0->etaD1()) < 2.4 && std::fabs(d0->etaD2()) < 2.4;
   bool passTrkPt = d0->PtD1() > 0.7 && d0->PtD2() > 0.7;
   bool passTrkPtErr = d0->PtErrD1()/d0->PtD1() < 0.1 && d0->PtErrD2()/d0->PtD2() < 0.1;
   bool passTrkPurity = d0->highPurityD1() && d0->highPurityD2();
   bool passTrkNhits = d0->nHitD1() >=11 && d0->nHitD2() >=11; 
   //bool passDeltaEta = std::fabs(d0->etaD1()-d0->etaD2())<1;
   bool passMVA = d0->Mva() > -0.3;
   passMVA = true;

   if(passPointingAngle && passTrkEta 
         && passTrkPt && passTrkPtErr && passTrkPurity && passTrkNhits && passMVA) return true;
   return false;
}
