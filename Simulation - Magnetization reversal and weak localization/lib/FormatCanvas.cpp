void FormatCanvas(TCanvas *cc, int NumOfPads=0) {

  // Assuming there is only 1 pad and 1 histo has been drawn already
  float TickLengthX = 0.01;
  float TickLengthY = 0.01;
  int LabelAndTitleFont = 43;
  int AxisTitleAndLabelSize_InPixels = 13;
  float AxisTitleToAxisLabelRatio = 1.2;
  bool CenterAxisTitle = 1;
  int FrameLineWidth = 1;
  
  
  // SetFrame
//  TFrame *qq = (TFrame*) cc->FindObject("TFrame");
//  qq->SetLineWidth(FrameLineWidth);
  
  
  for (int i=0; i<=NumOfPads; i++){
  
    cc->cd(i);
    float AxisTitleAndLabelSize = AxisTitleAndLabelSize_InPixels;
    bool WithPalette = 0;
	
    // SetTicks and Grids    
    cc->cd(i)->SetTicks(1,1);
    cc->cd(i)->SetGridy(0);
    cc->cd(i)->SetGridx(0);
    cc->cd(i)->SetFrameLineWidth(FrameLineWidth);
    
    // Set Frame margins
	
	if (NumOfPads==1) {
		cc->cd(i)->SetTopMargin(0.05);
		cc->cd(i)->SetBottomMargin(0.15);
		cc->cd(i)->SetLeftMargin(-0.35);
		cc->cd(i)->SetRightMargin(0.05);
		}
	else if (NumOfPads==2) {
		cc->cd(i)->SetBottomMargin(0.15);
		cc->cd(i)->SetLeftMargin(0.13);
		cc->cd(i)->SetRightMargin(0.08);
		}
    else {
		cc->cd(i)->SetBottomMargin(0.2);
		cc->cd(i)->SetLeftMargin(0.15);
		}

	
    // Set Axis properties
    TH1F *htemp = (TH1F*) gPad->FindObject("htemp");
    TGraph *graph = (TGraph*) gPad->FindObject("Graph");
    TMultiGraph *mg = (TMultiGraph*) cc->cd(i)->FindObject("multigraph");
    if (htemp) {
      // retreive Taxis
      TAxis *Ascisse = htemp->GetXaxis();
      TAxis *Ordinate = htemp->GetYaxis();
      TAxis *Third = htemp->GetZaxis();
      // retreive TPalette
      gPad->Update();
      TPaletteAxis *pal = (TPaletteAxis*)htemp->GetListOfFunctions()->FindObject("palette");
      if (pal) {
	WithPalette = 1;
	// pal->SetY2NDC(1- 0.1*(1-cc->cd(i)->GetBottomMargin()));
	// pal->SetY1NDC(cc->cd(i)->GetBottomMargin() + 0.1*(1-cc->cd(i)->GetBottomMargin()));
	// pal->SetX1NDC(pal->GetX1NDC()+0.005);
	// pal->SetX2NDC(pal->GetX2NDC()-0.02);
	TGaxis* Terzo = pal->GetAxis();
	TAxis *Third = htemp->GetZaxis();
      }
    }
    else if (graph) {
      TAxis *Ascisse = graph->GetXaxis();
      TAxis *Ordinate = graph->GetYaxis();
    }
    else if (mg) {
      TAxis *Ascisse = mg->GetXaxis();
      TAxis *Ordinate = mg->GetYaxis();
    }
    if (Ascisse) {
      Ascisse->SetTickLength(TickLengthX);
	  Ascisse->SetNdivisions(308);
	  Ascisse->SetLabelFont(LabelAndTitleFont);
	  Ascisse->SetLabelSize(AxisTitleAndLabelSize);
	  Ascisse->SetTitleFont(LabelAndTitleFont);
	  Ascisse->SetTitleSize(AxisTitleToAxisLabelRatio*AxisTitleAndLabelSize);
	  Ascisse->SetTitleOffset(1./cc->cd(i)->GetHNDC() * 1.1);
	  Ascisse->SetLabelOffset(0.01);
	  Ascisse->CenterTitle(CenterAxisTitle);
      Ordinate->SetTickLength(TickLengthY);
	  Ordinate->SetNdivisions(305);
	  if 	(NumOfPads>3) Ordinate->SetTitleOffset(1.8);
	  else Ordinate->SetTitleOffset(1./cc->cd(i)->GetHNDC() * 1.3);
	  Ordinate->SetLabelFont(LabelAndTitleFont);
	  Ordinate->SetLabelSize(AxisTitleAndLabelSize);
	  Ordinate->SetTitleFont(LabelAndTitleFont);
	  Ordinate->SetTitleSize(AxisTitleToAxisLabelRatio*AxisTitleAndLabelSize);
	  Ordinate->CenterTitle(CenterAxisTitle);
	  Ordinate->SetDecimals(1);
	  
      if (WithPalette) {
	      cout << "nn " << pal << endl;
	Terzo->SetTickSize(1e-3);
	  Terzo->SetLabelFont(LabelAndTitleFont);
	  Terzo->SetLabelSize(AxisTitleAndLabelSize);
	  Terzo->SetTitleFont(LabelAndTitleFont);
	  Terzo->SetTitleSize(AxisTitleToAxisLabelRatio*AxisTitleAndLabelSize);
	  Terzo->SetTitleOffset(1.2);
	  //Third->SetNdivisions(1);
	  Terzo->CenterTitle(CenterAxisTitle);
	  Third->RotateTitle(1);
	  Terzo->SetMaxDigits(3);
	  Third->Draw();
      }
    }  // if Ascisse
    
    
    // Set Legend properties
    TLegend *ll = (TLegend*) gPad->FindObject("TPave");
    if (!ll) ll = (TLegend*) gPad->FindObject("TLegend");
    if (ll) {
      ll->SetTextFont(43);
      ll->SetTextSizePixels(AxisTitleAndLabelSize_InPixels);
      ll->SetFillColor(0);
      ll->SetLineColor(0);
      //ll->SetTextAlign(22);
    }
	
	
  } // end loop on pads	
  
  
  cc->Update();
  return();
}


void FormatGraph(TGraph *s, int MStyle=1, int MCol=1, float MSize=1., int LStyle=1, int LCol=1, float LSize=1.){
	if (MCol==5) MCol+=10;
	if (LCol==5) LCol+=10;
	s->SetMarkerColor(MCol);
	s->SetMarkerStyle(MStyle);
	s->SetMarkerSize(MSize);
	s->SetLineStyle(LStyle);
	s->SetLineColor(LCol);
	s->SetLineWidth(LSize);
	s->SetFillColor(0);
	
	s->GetXaxis()->CenterTitle(1);
	s->GetYaxis()->CenterTitle(1);
	return();
}
//int MStyle=1, MCol=1, LStyle=1, LCol=1; 
//float MSize=1., LSize=1.;
//FormatGraph(gr, MStyle, MCol, MSize, LStyle, LCol, LSize);

TGraph* TGraph_Rebin(TGraph *inGr, float MinX, float MaxX, int NumOfPoints){
    float NewX[100000];
    float NewY[100000];
    float SelectedX;
    int ProcessedPoints=0;
    
	TGraph *input = (TGraph*) inGr->Clone("newGraph");
    input->Sort();
    Double_t a, b, c, d;
    input->GetPoint(0, a, b);  //x, y
    input->GetPoint(input->GetN()-1, c, d);
	if(MinX==MaxX){
		MinX = a;
		MaxX = c;
		}
    float IntervalLength = (MaxX - MinX)/NumOfPoints;
    float MaxXinInput = c;
    float MinXinInput = a;
    
    for (int i=0; i<NumOfPoints; i++){
        SelectedX = MinX + i*IntervalLength;
        
        if ((SelectedX >= MinXinInput)*(SelectedX <= MaxXinInput)){
            NewY[i] = input->Eval(SelectedX);
            NewX[i] = SelectedX;
            ProcessedPoints++;
        }
    }
    
    TGraph *ggg = new TGraph(ProcessedPoints, NewX, NewY);
    return(ggg);
}


TGraph* SubtractBackgroundFormulaFromTGraph (TGraph *InputGraph, TF1* background) {
  bool debug=0;
  int NumOfEntries = InputGraph->GetN();
  Double_t *Values;  Values = InputGraph->GetY();
  Double_t *Xes;     Xes = InputGraph->GetX();
  Double_t Signal;
  
  TNtuple bufferNtuple("buff", "buff", "index:x:value:signal");
  
  for (int jindex=0; jindex< NumOfEntries; jindex++) {
      Signal = (Values[jindex] - background->Eval(Xes[jindex]));
	  bufferNtuple.Fill(jindex, Xes[jindex], Values[jindex], Signal);
	  if (debug) cout << "  --SubtractBackgroundFormulaFromTGraph: " << Xes[jindex] << " " << Signal[jindex] << endl;
  }
  bufferNtuple.Draw("signal:index", "", "goff");
  TGraph *BkGraph = new TGraph(NumOfEntries, Xes, bufferNtuple.GetV1());
  BkGraph->SetNameTitle("SubtractBackgroundFormulaFromTGraph", "SubtractBackgroundFormulaFromTGraph");
  
  
  if (debug) {
	cout << "  --SubtractBackgroundFormulaFromTGraph:NumOfEntries " << NumOfEntries << endl;
	cout << "  --SubtractBackgroundFormulaFromTGraph:BkGraph->GetN() " << BkGraph->GetN() << endl;
	}
  return(BkGraph);
}
  
  

TMultiGraph* ExtractGraphsFromTNtuple (TNtuple *n, string WhatToPlot  ="Stray_h:m", char *WhatToSelect, char *OverallConditions){
    char PlottingOptions[300] = "lp";
    float thresholdForPeakSearch = 0.0001;
    
    TMultiGraph *multigraph = new TMultiGraph("multigraph", "");
    char buffer[300];
    int color=1;
	int MStyle=1, LStyle=1, LCol=1;
	float MSize=2., LSize=2.;
    
    // extract the values of WhatToSelect
    n->Draw(WhatToSelect, OverallConditions, "goff");
    TH1F *h = (TH1F*) gDirectory->Get("htemp");
    TSpectrum *s = new TSpectrum(20);
    Int_t nfound = s->Search(h,1,"new", thresholdForPeakSearch);
    printf("Found %d candidate peaks to fitn \n",nfound);
    Float_t *xpeaks = s->GetPositionX();
    if (nfound <2) {
        xpeaks[0] = h->GetMean();
        nfound = 1;
    }
	
    //Sort the array of Peaks
    int index[300];
    TMath::Sort(nfound, xpeaks, index, false);
	
	// find smallest distance between extracted values
	float precisionForPlot=1e10;
	for (int ll=0; ll+1<nfound; ll++) {
		float a, b;
		a = xpeaks[index[ll]];
		b = xpeaks[index[ll+1]];
		precisionForPlot = ( precisionForPlot < fabs(a-b) ? precisionForPlot : fabs(a-b) );
	}

/*    
    //Create Legend (large one of large num of entries)
    if (nfound>4) TLegend *leg1 = new TLegend(0.75, 0.5, 0.99, 0.99, WhatToSelect);
    else TLegend *leg1 = new TLegend(0.75, 0.7, 0.99, 0.99, WhatToSelect);
*/    
    // Draw Graphs
	//gStyle->SetPalette(51);
	///TColor col;
	float SelectedValues_Max = TMath::MaxElement(nfound, xpeaks), SelectedValues_Min = TMath::MinElement(nfound, xpeaks);
    for (int i=0; i< nfound; i++){
        float mm = xpeaks[index[i]];
        sprintf(buffer, "%s && abs(%s - %.3e)<%e", OverallConditions, WhatToSelect, mm, precisionForPlot/2);
		if (OverallConditions == "") sprintf(buffer, "abs(%s - %.3e)<%e", WhatToSelect, mm, precisionForPlot/2);
        cout << "  Plot Condition " << i << ": " << buffer << endl;
        n->Draw(WhatToPlot.c_str(), buffer, "goff");
            TGraph *gr = new TGraph(n->GetSelectedRows(),n->GetV2(),n->GetV1());
            gr->Sort();
            gr->SetName("Graph");
            if (fabs(log(fabs(mm))/log(10) ) >3) sprintf(buffer, "%.2e", mm);
            else if (fabs(mm) > 5)           sprintf(buffer, "%.1f", mm);
            else                            sprintf(buffer, "%.2f", mm);
            gr->SetTitle(buffer);
				TColor c;
				int ColorNum = (int) c.GetColorPalette(49*(mm-SelectedValues_Min)/(SelectedValues_Max-SelectedValues_Min));
				cout << "col: " << ColorNum << endl;
				//col.GetColorPalette(); //color = i+1;//
				FormatGraph(gr, MStyle, ColorNum, MSize, LStyle, ColorNum, LSize);
            multigraph->Add(gr);
    }
                
    delete s;
    return(multigraph);
}



TGraph* MakeGraphFromFile(char* InputFile = "2_Summary.csv", char* ColumnNames, char* WhatToDraw, char* Selection ="", bool Normalized=0){
	//gr = MakeGraphFromFile(InputFile.c_str(), ColumnNames.c_str(), WhatToDraw.c_str(), Selection.c_str());
  TNtuple n("n", "n", ColumnNames);
  n.ReadFile(InputFile);
  n.Draw(WhatToDraw, Selection, "goff");
  const int NumOfEvents=n.GetSelectedRows();
  
  if (!NumOfEvents) {
	  cout << "  -->MakeGraphFromFile(): no points for graph from " << InputFile << endl;
	  cout << "    Drawing: " << WhatToDraw << endl;
	  cout << "    Selection: " << Selection << endl;
	  return(0);
  }
  
  Double_t YvaluesNormalized[NumOfEvents];
  Double_t *Yvalues  = n.GetV1();
  
  if (Normalized) {
	float MaxValue = TMath::MaxElement(NumOfEvents, Yvalues);
	float MinValue = TMath::MinElement(NumOfEvents, Yvalues);
	for (int i=0; i<NumOfEvents; i++) YvaluesNormalized[i] = (Yvalues[i]-MinValue)/MaxValue;
	TGraph *gr = new TGraph(NumOfEvents, n.GetV2(), YvaluesNormalized);
  }
  else {
	TGraph *gr = new TGraph(n.GetSelectedRows(), n.GetV2(), Yvalues);
	}
  gr->SetNameTitle(InputFile, InputFile);
  
  
  return(gr);
}



FormatTF1(TF1 *l, int color=2){
	l->SetLineColor(color);
	l->SetFillColor(0);
	l->SetLineStyle(1);
	l->SetLineWidth(3);
	return();
}

