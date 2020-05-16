#include <iostream>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TBrowser.h>
#include <cstring>
#include <TLatex.h>

using namespace std;


void DistributionDrawer(int EventType, char* CorrelationType){
	TCanvas* c = new TCanvas("c","c",1600,800);

	TString* filename = new TString();
	filename->Form("tpiOutput/EventType%d/%s.root", EventType, CorrelationType);
	TFile* file = new TFile(filename->Data());
	delete filename;

	vector<string> names2D = {
			"cnumepNonIdEPTrue", "cdenepNonIdEPTrue",
			"cnumepNonIdEPTrueNoQS", "cdenepNonIdEPTrueNoQS"
	};
	vector<string> names1D = {"hpt", "heta", "hphiP", "hevmultAll"};
	vector<string> names2F = {"hevmultPID"};

	//Drawing divided histograms (TH2D)
	for(int i=0; i<4; i+=2){
		names2D[i].append(CorrelationType);
		names2D[i + 1].append(CorrelationType);

		double scale_num = ((TH2D*)gDirectory->Get(Form(names2D[i].c_str())))->Integral();
		double scale_den = ((TH2D*)gDirectory->Get(Form(names2D[i+1].c_str())))->Integral();

   		TH2D* num = (TH2D*)gDirectory->Get(Form(names2D[i].c_str()));
	   	TH2D* den = ((TH2D*)gDirectory->Get(Form(names2D[i+1].c_str())));
	   	num->Scale(1./scale_num);
   		den->Scale(1./scale_den);
   		num->Divide(den);
		
		TH2D* hdiv = num;
	
		//Title
		char* title = new char[100];
		if (names2D[i].find("NoQS") != string::npos)
		{
			strcpy(title, "Korelacje katowe");
		}
		else {
			strcpy(title, "Korelacje katowe z uwzglednieniem statystyki kwantowej");
		}
		hdiv->SetTitle(title);
		//hdiv->SetName(names2D[i+1]);


		hdiv->GetXaxis()->SetTitle("#Delta#phi");
		hdiv->GetXaxis()->SetTitleSize(0.07);
		hdiv->GetXaxis()->SetTitleOffset(1.08);

		hdiv->GetYaxis()->SetTitle("#Delta#eta");
		hdiv->GetYaxis()->SetTitleSize(0.07);
		hdiv->GetYaxis()->SetTitleOffset(1.08);

		hdiv->GetZaxis()->SetTitle("C(#Delta#eta, #Delta#phi)");
		hdiv->GetZaxis()->SetTitleSize(0.07);
		hdiv->GetZaxis()->SetTitleOffset(0.60);

		hdiv->SetStats(false);


		hdiv->Draw("surf1");


		char* file_name = new char[40];
		strcpy(file_name, names2D[i].c_str());
		strcat(file_name, ".png");

		c->SaveAs(file_name);
	}


	//Drawing the rest
	//TH1D
	for(string i:names1D){
		i.append(CorrelationType);
		TH1D* h = (TH1D*)file->Get(i.c_str());

		h->SetStats(false);
		h->SetFillColor(kRed);
        h->SetFillStyle(3006);
		if(i.find("heta") != string::npos){
			h->GetXaxis()->SetTitle("#eta");
			h->SetTitle("Rozklad kata #eta");
		}
		else if(i.find("hevmultAll") != string::npos){
			h->GetXaxis()->SetTitle("Krotnosc zderzenia");
			h->GetXaxis()->SetRange(11, 21);
			h->SetTitle("Rozklad krotnosci wyprodukowanych czastek");

		}
		else if(i.find("hpt") != string::npos){
			h->GetXaxis()->SetTitle("p_{T} [GeV/c]");
			h->GetYaxis()->SetTitle("dN/dp_{T} [GeV/c]");
			h->SetTitle("Rozklad pedu poprzecznego czastek");

		}
		else{
			h->GetXaxis()->SetTitle("#varphi");
			h->SetMinimum(0);
			h->SetTitle("Rozklad kata #varphi");

		}
		//char* title = new char[40];
		//strcpy(title, i);

		//h->GetXaxis()->SetTitleSize(0.19);
		//h->GetYaxis()->SetTitleSize(0.19);

		//h->SetTitle(title);
		h->Draw("B");


		char* file_name = new char[25];
		strcpy(file_name, i.c_str());
		strcat(file_name, ".png");
		c->SaveAs(file_name);
	}


	//TH2F
	names2F[0].append(CorrelationType);
	TH2F* h = (TH2F*)file->Get(names2F[0].c_str());

	h->SetTitle("Rozklad wystepowania krotnosci danej czastki wsrod wszystkich zderzen");
	h->SetStats(false);
	h->GetYaxis()->SetRange(0,20);
	h->Draw("colzz");

	char* file_name = new char[25];
	strcpy(file_name, names2F[0].c_str());
	strcat(file_name, ".png");
	c->SaveAs(file_name);

	delete c;
	file->Close();


}
