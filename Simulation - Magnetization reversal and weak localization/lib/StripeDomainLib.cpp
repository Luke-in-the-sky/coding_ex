// 
// Library of functions for the simulation of the Stripe Domains and their
// effect on the Weak Localization for graphene
//
// The magnetization reversal is computed by
//    	void Looping_m(float, float, int, double)
// which then saves the results to a text file (the name of which is in GlobalConstants::LogFileTitle)

// This allows for storing the simulation results in one place, with its simulation parameters.
// Such database can then be accessed for analysis. Two routines two do this are
// 		TMultiGraph* PlotSomethingSelecting2quantities(string, string*, string)
//   	TMultiGraph* PlotSomethingForAllZ (TCanvas*, int, string, char*, char*)
//



// =========== ===========
//
// Functions analyzing the results of the simulation fromt the stored data
//
// =========== ===========


// ------------
// Read the simulation results, select a given subset and make plots
int ReadFromDataFile( TCanvas *can, int padindex) {
	GlobalConstants g;
    char buffer[300];
	
	// Check File exists
	ifstream ifile(g.LogFileTitle.c_str());
	if (!ifile) {
	  cout << g.LogFileTitle.c_str() << " -- !! File not found" << endl;
	  return;
	}
        
		
	TNtuple *n = new TNtuple("n", "n", g.LogFile_ListOfColumnNames.c_str());
	n->ReadFile(g.LogFileTitle.c_str());
	n->SetMarkerStyle(4);	 
    can->cd(padindex++);
	sprintf(buffer, "AverRho:h*%f", g.Hstray_sat);
        multigraph = PlotSomethingForAllZ (can, padindex, buffer , "l_phi", g.OverallPlottingConditions.c_str());  
		multigraph->GetHistogram()->GetXaxis()->SetTitle("h");
		multigraph->GetHistogram()->SetTitle("<#rhoH(x,z)>");

            
	TNtuple *Ndata = new TNtuple ("Ndata", "Ndata", "H:rho:temp:ang:curr");
	Ndata->ReadFile("Measurements1.dat");
	Ndata->SetMarkerStyle(4);
	Ndata->Draw("rho/4:H", "temp<3 && ang<20", "goff");
		TGraph *gg = new TGraph(Ndata->GetSelectedRows(), Ndata->GetV2(), Ndata->GetV1());
		//gg->GetXaxis()->SetRangeUser(-0.1, 0.8);
		//gg->SetTitle("#rho at 2K for M_{s}=2500 G");
		//gg->GetXaxis()->SetTitle("H/M_{s}");
	multigraph->Add(gg);
        multigraph->GetYaxis()->SetRangeUser(624, 628);
        
        
        can->cd(padindex++);
	sprintf(buffer, "AverRho:h*%f", g.Hstray_sat);
        multigraph = PlotSomethingForAllZ (can, padindex, buffer , "z", g.OverallPlottingConditions.c_str());  
            multigraph->GetHistogram()->GetXaxis()->SetTitle("h");
            multigraph->GetHistogram()->SetTitle("<#rhoH(x,z)>");
        multigraph->Add(gg);
        multigraph->GetYaxis()->SetRangeUser(625, 628);
        
        
        // for different Ms
        can->cd(padindex++);
        string condition = g.OverallPlottingConditions + " && abs(z-8.6e-7)<1e-9";
	sprintf(buffer, "AverRho:h*%f", g.Hstray_sat);
        multigraph = PlotSomethingForAllZ (can, padindex, buffer , "Ms", condition.c_str());  
            multigraph->GetHistogram()->GetXaxis()->SetTitle("h");
            multigraph->GetHistogram()->SetTitle("<#rhoH(x,z)>");
        multigraph->Add(gg);
        multigraph->GetYaxis()->SetRangeUser(625, 628);
        
        
        // for different Hstray_sat
        can->cd(padindex++);
        string condition = g.OverallPlottingConditions + " && abs(z-8.6e-7)<1e-9 && abs(Ms-300)<2";
	sprintf(buffer, "AverRho:h*%f", g.Hstray_sat);
        multigraph = PlotSomethingForAllZ (can, padindex, buffer , "Hstray_sat", condition.c_str());  
            multigraph->GetHistogram()->GetXaxis()->SetTitle("h");
            multigraph->GetHistogram()->SetTitle("<#rhoH(x,z)>");
        multigraph->Add(gg);
        multigraph->GetYaxis()->SetRangeUser(625, 628);
        
        
        
        // for 1 Hstray_sat
        can->cd(padindex++);
        TMultiGraph *mm = new TMultiGraph("multigraph", "");
        TLegend *leg1 = new TLegend(0.75, 0.7, 0.99, 0.99);
        
               
                 
	Ndata->Draw("rho/4:H", "temp<3 && ang<20", "goff");
        TGraph *gg = new TGraph(Ndata->GetSelectedRows(), Ndata->GetV2(), Ndata->GetV1());
            FormatGraph(gg, 1, 4, 1, 1, 1, 1);
            gg->Sort();
            mm->Add(gg, "p");            
            leg1->AddEntry(gg, "data", "p");
        
        string condition = g.OverallPlottingConditions + " && abs(z-8.6e-7)<1e-9 && abs(Ms-300)<2 && abs(Hstray_sat-90)<2";
	sprintf(buffer, "AverRho:h*%f", g.Hstray_sat);
        n->Draw(buffer, condition.c_str(), "goff");
            TGraph *gg = new TGraph(n->GetSelectedRows(), n->GetV2(), n->GetV1());
            //FormatGraph(TGraph *s, int MStyle=1, int MCol=1, float MSize=1., int LStyle=1, int LCol=1, float LSize=1.){
            FormatGraph(gg, 1, 1, 1, 2, 1, 3);
            gg->Sort();
            mm->Add(gg, "l");
            leg1->AddEntry(gg, "simul", "l");
            
            
            mm->Draw("a");
            mm->GetYaxis()->SetRangeUser(625, 627.5);
            mm->GetXaxis()->SetRangeUser(-200, 1500);
            mm->GetXaxis()->SetTitle("H_{ext} [Oe]");
            mm->GetYaxis()->SetTitle("<#rho> [#Omega]");
            mm->SetTitle("");
            TMultiGraph *Newmm = mm->Clone("multigraph");
       
	return(padindex);
}


// ------------
// Read the simulation results, select a given subset and make plots
TMultiGraph* PlotSomethingSelecting2quantities(string WhatToPlot  ="Stray_h:m", string *WhatToSelect, string OverallConditions){
	const int NumOfSelections=2,
    char PlottingOptions[300] = "lp";
    float thresholdForPeakSearch = 0.0001;
    float precisionForPlot[NumOfSelections];
	int nfound[NumOfSelections];
	float ExtractedValues[NumOfSelections][20]; 
	int index[300];
	char buffer[300];
    int color=1;
	
    GlobalConstants g;
	TMultiGraph *multigraph = new TMultiGraph("multigraph", "");
    TNtuple *n = new TNtuple("n", "n", g.LogFile_ListOfColumnNames.c_str());
    n->ReadFile(g.LogFileTitle.c_str());
    n->SetMarkerStyle(g.MStyle);
	n->SetMarkerSize(g.MSize);	
    
	//Create Legend (large one of large num of entries)
    TLegend *leg1 = new TLegend(0.75, 0.5, 0.99, 0.99);
    //else TLegend *leg1 = new TLegend(0.75, 0.7, 0.99, 0.99, WhatToSelect);
 
    // extract the values of WhatToSelect
	for (int sel=0; sel<NumOfSelections; sel++){
		n->Draw(WhatToSelect[sel].c_str(), OverallConditions.c_str(), "goff");
		TH1F *h = (TH1F*) gDirectory->Get("htemp");
		TSpectrum *s = new TSpectrum(20);
		nfound[sel] = s->Search(h,1,"new", thresholdForPeakSearch);
		printf("%s: found %d candidate peaks to fitn \n",WhatToSelect[sel].c_str(),nfound[sel]);
		Float_t *xpeaks = s->GetPositionX();
		if (nfound[sel] <2) {
			xpeaks[0] = h->GetMean();
			nfound[sel] = 1;
		}
		
		//Sort the array of Peaks
		TMath::Sort(nfound[sel], xpeaks, index, false);
		
		// store extracted values
		for (int ll=0; ll<nfound[sel]; ll++) {
			ExtractedValues[sel][ll] = xpeaks[index[ll]];
			cout << ExtractedValues[sel][ll] << ", ";
		}
		cout << endl;
		delete s;
		
		// find smallest distance between extracted values
		precisionForPlot[sel]=1e10;
		for (int ll=0; ll+1<nfound[sel]; ll++) {
			float a, b;
			a = ExtractedValues[sel][ll];
			b = ExtractedValues[sel][ll+1];
			precisionForPlot[sel] = ( precisionForPlot[sel] < fabs(a-b) ? precisionForPlot[sel] : fabs(a-b) );
		}
	} // loop on extracting values of WhatToSelect

        
    // Draw Graphs
	string SelectionForOneGraph="";
	float Msize=1.;
	int marker = 30;
	char buffer2[300];
	// allocate counter for legend entries
	int legcounterI[300];
	for (int i=0; i< nfound[1]; i++) legcounterI[i]=0;
	for (int j=0; j<nfound[0]; j++){
		marker +=j;
		float SelVal = ExtractedValues[0][j];
		sprintf(buffer2, "abs(%s - %.3e)<%e", WhatToSelect[0].c_str(), SelVal, precisionForPlot[0]/2);
		for (int i=0; i<nfound[1]; i++){
			float mm = ExtractedValues[1][i];
			sprintf(buffer, "%s && %s && abs(%s - %.3e)<%e", OverallConditions.c_str(), buffer2, WhatToSelect[1].c_str(), mm, precisionForPlot[1]/2);
			cout << "  Plot Condition " << i << ": " << buffer << endl;
			n->Draw(WhatToPlot.c_str(), buffer, "goff");
			
			// If selection is not null, make a Graph
			if (n->GetSelectedRows()){
				TGraph *gr = new TGraph(n->GetSelectedRows(),n->GetV2(),n->GetV1());
				gr->Sort();
				gr->SetName("Graph");
				gr->SetTitle(buffer);
				
				//FormatGraph(TGraph *s, int MStyle=1, int MCol=1, float MSize=1., int LStyle=1, int LCol=1, float LSize=1.){
				color=i+1;
				FormatGraph(gr, marker, color, Msize, 1, color++, 1);
				if (!legcounterI[i]) {
					if (fabs(log(fabs(mm))/log(10) ) >3) 	sprintf(buffer, "%.2e", mm);
					else if (fabs(mm) > 5)					sprintf(buffer, "%.1f", mm);
					else                            		sprintf(buffer, "%.2f", mm);
					leg1->AddEntry(gr, buffer, "l");
					legcounterI[i]++;
				}
				multigraph->Add(gr);
			}
		}
	} // loop on "sel"

	multigraph->Draw("apl");
	leg1->Draw();
	
	return (multigraph);
}


// ------------
// Read the simulation results, select a given subset and make plots
TMultiGraph* PlotSomethingForAllZ (TCanvas *can, int padindex, string WhatToPlot  ="Stray_h:m", char *WhatToSelect, char *OverallConditions){
    char PlottingOptions[300] = "lp";
    float thresholdForPeakSearch = 0.0001;
    
    GlobalConstants g;
    TNtuple *n = new TNtuple("n", "n", g.LogFile_ListOfColumnNames.c_str());
    n->ReadFile(g.LogFileTitle.c_str());
    n->SetMarkerStyle(g.MStyle);
	n->SetMarkerSize(g.MSize);
    
    TMultiGraph *multigraph = new TMultiGraph("multigraph", "");
    char buffer[300];
    int color=1;
    
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

    
    //Create Legend (large one of large num of entries)
    if (nfound>4) TLegend *leg1 = new TLegend(0.75, 0.5, 0.99, 0.99, WhatToSelect);
    else TLegend *leg1 = new TLegend(0.75, 0.7, 0.99, 0.99, WhatToSelect);
    
    // Draw Graphs
    for (int i=0; i< nfound; i++){
        float mm = xpeaks[index[i]];
        sprintf(buffer, "%s && abs(%s - %.3e)<%e", OverallConditions, WhatToSelect, mm, precisionForPlot/2);
        cout << "  Plot Condition " << i << ": " << buffer << endl;
        n->Draw(WhatToPlot.c_str(), buffer, "goff");
		
		TGraph *gr = new TGraph(n->GetSelectedRows(),n->GetV2(),n->GetV1());
		gr->Sort();
		gr->SetName("Graph");
		if (fabs(log(fabs(mm))/log(10) ) >3) sprintf(buffer, "%.2e", mm);
		else if (fabs(mm) > 5)           sprintf(buffer, "%.1f", mm);
		else                            sprintf(buffer, "%.2f", mm);
		gr->SetTitle(buffer);
			gr->SetMarkerStyle(g.MStyle); 	gr->SetMarkerSize(g.MSize);
		gr->SetLineColor(color); gr->SetMarkerColor(color++);
		leg1->AddEntry(gr, gr->GetTitle(), PlottingOptions);
		multigraph->Add(gr);
    }
            
    sprintf(buffer, "a%s", PlottingOptions);
    multigraph->Draw(buffer);
    leg1->Draw();
    
    delete s;
    delete n;
    return(multigraph);
}




// =========== ===========
//
// Run the simulation for one magnetization reversal
//
// =========== ===========

void Looping_m(float MinM=0., float MaxM=0.9, int NumOfPoints=5, double FixZ = 1.004e-6) {
    bool verbose =1;
    GlobalConstants g;
    ofstream OutFile;
    OutFile.open(g.LogFileTitle.c_str(), ios::app);
    OutFile <<"# " << g.LogFile_ListOfColumnNames << endl;
    
    for (float m = MinM; m <= MaxM; m+= (MaxM-MinM)/NumOfPoints ){
        // Compute d    
        TF1 LoopM__d_m ("LoopM__d_m", d_m, g.Min_m, g.Max_m, 1);
        LoopM__d_m.SetParameter(0, g.tau);
        double periodicity_d = LoopM__d_m.Eval(m);
        if (verbose) cout << "[m=" << m << "] d: " << periodicity_d << " for m=" << m << " and tau=" << g.tau << endl;
        
        // Compute h(m, d) from Eq.6
        TF1 LoopM_fEq6("LoopM_fEq6", h_md, g.Min_d, g.Max_d, 1);
        LoopM_fEq6.SetParameter(0, periodicity_d);
        float externalField = LoopM_fEq6.Eval(m);
        if (verbose) cout << "[m=" << m << "] h: " << externalField << " for m=" << m << " and tau=" << g.tau << endl;

        // Compute limits for Integration along x-axis
        // NB: Npx will not affect quality/speed of integration, only of the plot. NumOfPeriods and precision will adffect integration
        int NumOfPeriods 	    = g.LoopingM_MaxNumOfPeriodsForInt;
        double EpsForIntegration    = g.LoopingM_EpsForIntegration;
        double MinXforInt = g.MinX;
        double MaxXforInt = g.MinX+NumOfPeriods*periodicity_d;
        if (g.GrapheneWidth) {
            MinXforInt = g.MinX + g.GrapheneOffset;
            MaxXforInt = MinXforInt + g.GrapheneWidth;
        }
        
        // Allocate StrayField Functions
        TF2 *f__StrayField_zComp = new TF2("f__StrayField_zComp", StrayField_zComp, MinXforInt, MaxXforInt, g.t/2, g.MaxZ, 2);
        f__StrayField_zComp->SetParameters(m, periodicity_d);
        TF12 *f__StrayField_zComp_x = new TF12("f__StrayField_zComp_x", f__StrayField_zComp, FixZ, "x");
        f__StrayField_zComp_x->SetNpx(g.h_m__Npx);
        f__StrayField_zComp_x->SetParameters(m, periodicity_d);
        
        // Integrate: use Abs(TF1) when computing integral)
        f__StrayField_zComp_x->AbsValue(1);	
        double r1 = f__StrayField_zComp_x->Integral( MinXforInt, MaxXforInt, f__StrayField_zComp_x->GetParameters(), EpsForIntegration) / (MaxXforInt - MinXforInt);
        // Integrate: use TF1 (not Abs(TF1)) when computing integral)
        f__StrayField_zComp_x->AbsValue(0);	
        double r2 = f__StrayField_zComp_x->Integral( MinXforInt, MaxXforInt, f__StrayField_zComp_x->GetParameters(), EpsForIntegration) / (MaxXforInt - MinXforInt);

        
		// Allocate TotalField Functions
        TF2 *f__TotField_zComp = new TF2("f__TotField_zComp", TotField_zComp, MinXforInt, MaxXforInt, g.t/2, g.MaxZ, 2);
        f__TotField_zComp->SetParameters(m, periodicity_d);
        TF12 *f__TotField_zComp_x = new TF12("f__StrayField_zComp_x", f__TotField_zComp, FixZ, "x");
        f__TotField_zComp_x->SetNpx(g.h_m__Npx);
        f__TotField_zComp_x->SetParameters(m, periodicity_d);
        
        // Integrate: use Abs(TF1) when computing integral)
        f__TotField_zComp_x->AbsValue(1);	
        double r3 = f__TotField_zComp_x->Integral( MinXforInt, MaxXforInt, f__TotField_zComp_x->GetParameters(), EpsForIntegration) / (MaxXforInt - MinXforInt);
        // Integrate: use TF1 (not Abs(TF1)) when computing integral)
        f__TotField_zComp_x->AbsValue(0);	
        double r4 = f__TotField_zComp_x->Integral( MinXforInt, MaxXforInt, f__TotField_zComp_x->GetParameters(), EpsForIntegration) / (MaxXforInt - MinXforInt);


        // Compute Average WL-resistivity
        //
        // Allocate LocalWLfunction
        TF1 *f__WL = new TF1("f__WL", LocalWeakLocalization, MinXforInt, MaxXforInt, 4);
        f__WL->SetParameters(m, periodicity_d, externalField, FixZ);
        // Integrate
        double AverageWLresistivity = f__WL->Integral( MinXforInt, MaxXforInt, f__WL->GetParameters(), EpsForIntegration) / (MaxXforInt - MinXforInt);
		cout << "  Integration interval: (" << MinXforInt << "," << MaxXforInt << ")" << endl;
        
        
        // Write to File
        WL_parameters WLpar;
        OutFile << m << "	" << periodicity_d  << "	" <<  externalField  << "	" << FixZ << "	" << g.t << " " << g.tau << "	" << r1 << "	" << r2 << "	" << r3 << "	" << r4<< "	" << g.Convergence__e_d << "    " << g.Convergence__phi << "	" << NumOfPeriods << "	" << EpsForIntegration << "	" << g.GrapheneWidth << "	" << g.GrapheneOffset<< "	" << AverageWLresistivity << "  " << WLpar.rho0 << "  " << WLpar.l_phi << "  " << WLpar.l_i << "  " << WLpar.l_star << "	" << g.Ms << "  " << g.Hstray_sat << endl;
        if (verbose) {
                cout << "[m=" << m << ", z=" << FixZ << "] <|h|>: " << r1 << endl;
                cout << "[m=" << m << ", z=" << FixZ << "] < h >: " << r2 << " (" << EpsForIntegration << ")" <<  endl << endl;
        }
        
        delete f__StrayField_zComp;
        delete f__StrayField_zComp_x;
        delete f__TotField_zComp;
        delete f__TotField_zComp_x;
        delete f__WL;
    } // endl of loop on m

    return();       
}




// =========== ===========
//
// Function calcultating the total magnetic field and its impact on the Weak Localization
//
// =========== ===========


// ------------
// Compute <rho(m)>
// 		retunrs the value of rho (call this function for fitting)
Double_t ComputeAverageRho(Double_t *Variables, Double_t *par){
    float m = Variables[0];
	float z = par[0];
	float rho	= par[1];
	float lphi	=par[2];
	
	bool verbose =1;
    GlobalConstants g;
    
	// Compute d    
	TF1 LoopM__d_m ("LoopM__d_m", d_m, g.Min_m, g.Max_m, 1);
	LoopM__d_m.SetParameter(0, g.tau);
	double periodicity_d = LoopM__d_m.Eval(m);
	if (verbose) cout << "[m=" << m << "] d: " << periodicity_d << " for m=" << m << " and tau=" << g.tau << endl;
	
	// Compute h(m, d) from Eq.6
	TF1 LoopM_fEq6("LoopM_fEq6", h_md, g.Min_d, g.Max_d, 1);
	LoopM_fEq6.SetParameter(0, periodicity_d);
	float externalField = LoopM_fEq6.Eval(m);

	
	// Compute limits for Integration along x-axis
	// NB: Npx will not affect quality/speed of integration, only of the plot. NumOfPeriods and precision will adffect integration
	int NumOfPeriods 	    = g.LoopingM_MaxNumOfPeriodsForInt;
	double EpsForIntegration    = g.LoopingM_EpsForIntegration;
	double MinXforInt = g.MinX;
	double MaxXforInt = g.MinX+NumOfPeriods*periodicity_d;
	if (g.GrapheneWidth) {
		MinXforInt = g.MinX + g.GrapheneOffset;
		MaxXforInt = MinXforInt + g.GrapheneWidth;
	}
	
	// Allocate LocalWLfunction	
	WL_parameters::l_phi = lphi;
	WL_parameters::rho0	= rho;
	TF1 *f__WL = new TF1("f__WL", LocalWeakLocalization, MinXforInt, MaxXforInt, 4);
	f__WL->SetParameters(m, periodicity_d, externalField, z);
	// Integrate
	double AverageWLresistivity = f__WL->Integral( MinXforInt, MaxXforInt, f__WL->GetParameters(), EpsForIntegration) / (MaxXforInt - MinXforInt);
	
	
	// Write to File
	if (verbose) {
			cout << "Fitting [m=" << m << "]: rho=" << AverageWLresistivity <<" (" << EpsForIntegration << ")" << endl;
	}
	
	delete f__WL;
    return(AverageWLresistivity);    
}


TMultiGraph* CheckLocalWL(TCanvas *can, int &padindex, float FixZ) {
    GlobalConstants g;
    char buffer[300];
    int color =1;
    float m, period;
    const int NumOfPlots = 3;
    float m_List[NumOfPlots] = {0.01, 0.5, 0.8};
    float period_List[NumOfPlots];
    float externField[NumOfPlots];
    

    // Compute d from Eq.5 and external field  h(m, d) from Eq.6
    TF1 LocalWL__d_m ("LocalWL__d_m", d_m, g.Min_m, g.Max_m, 1);
    LocalWL__d_m.SetParameter(0, g.tau);
    TF1 WL_fEq6("WL_fEq6", h_md, g.Min_d, g.Max_d, 1);   
    for (int i=0; i<NumOfPlots; i++) {
        period_List[i] = LocalWL__d_m.Eval(m_List[i]);
        WL_fEq6.SetParameter(0, period_List[i]);
        externField[i] = WL_fEq6.Eval(m_List[i]);
    }
    
    // Draw x-Profile of StrayField
    can->cd(padindex);    
    TGraph *lgr;
    TLegend *legg = new TLegend(.75, .7, .99, .99, "m");
    TMultiGraph *mult = new TMultiGraph("multigraph", "");
    for (int i=0; i<NumOfPlots; i++) {
        lgr = CheckStrayField_ProfileX(FixZ, "x", m_List[i], period_List[i]);
        //void FormatGraph(TGraph *s, int MStyle=1, int MCol=1, float MSize=1., int LStyle=1, int LCol=1, float LSize=1.){
        FormatGraph(lgr, 4, i+1, 1, 1, i+1, 1);
        legg->AddEntry(lgr, lgr->GetTitle(), "l");
        mult->Add(lgr);
    }
    can->cd(padindex++)->Clear();
    mult->Draw("apl");
    mult->GetXaxis()->SetTitle("x [um]");
    sprintf(buffer,"Stray Field h (z=%.1e)", FixZ);
    mult->GetHistogram()->SetTitle(buffer);
    //legg->Draw();
    
    
    // Draw the resistivity plots
    can->cd(padindex++);    
    // Compute minimum x and maximum x value
    double    MinXforInt = g.MinX + g.GrapheneOffset;
    double    MaxXforInt = MinXforInt + g.GrapheneWidth;
        
    TLegend *leg = new TLegend(.75, .7, .99, .99, "m");
    TMultiGraph *multi = new TMultiGraph("multigraph", "");
    TF1 *fWL = new TF1 ("fWL", LocalWeakLocalization, MinXforInt, MaxXforInt, 4);
    
    for (int i=0; i<NumOfPlots; i++) {
        m = m_List[i];
        period = period_List[i];
        cout << m << " " << period << " " << externField[i] << " " << MinXforInt << " " << MaxXforInt<< endl;
        fWL->SetParameters(m, period, externField[i], FixZ);
        fWL->SetNpx(g.h_m__Npx);
            TGraph *gr = new TGraph(fWL);
            gr->SetName("Graph");
            sprintf(buffer, "%.2f", m);
            gr->SetTitle(buffer);
             	gr->SetMarkerStyle(g.MStyle); 	gr->SetMarkerSize(g.MSize);
            gr->SetLineColor(color);
            gr->SetMarkerColor(color++);
            
            leg->AddEntry(gr, gr->GetTitle(), "l");
            multi->Add(gr);
    }
    
    multi->Draw("apl");
    multi->GetXaxis()->SetTitle("x [um]");
    multi->GetYaxis()->SetTitle("#rho [#Omega]");
    leg->Draw();
 
    delete fWL;
    return(multi);
}



// ------------
// Compute local resistivity function
//      fix z, use StrayField_zComp to compute the magnetic field, feed it into the weak localization function
Double_t LocalWeakLocalization(Double_t *Variable, Double_t *par) {
    GlobalConstants g;
    WL_parameters wlPar;
    float x     = Variable[0];
    float m         = par[0];
    float period    = par[1];
    float h         = par[2];
    float FixZ      = par[3];
    
    // Allocate StrayField Functions
    TF2 *fWL__StrayField_zComp = new TF2("fWL__StrayField_zComp", StrayField_zComp, g.MinX, g.MaxX, g.t/2, g.MaxZ, 2);
    fWL__StrayField_zComp->SetParameters(m, period);
    TF12 *fWL__StrayField_zComp_x = new TF12("fWL__StrayField_zComp_x", fWL__StrayField_zComp, FixZ, "x");
    fWL__StrayField_zComp_x->SetNpx(g.h_m__Npx);
    fWL__StrayField_zComp_x->SetParameters(m, period);
    
    // Compute magnetic field: stray field + external field
    //double localField = (fWL__StrayField_zComp_x->Eval(x) + h)* g.Ms * 1e-4; //This will be in T
    double localStrayField = (fWL__StrayField_zComp_x->Eval(x))* g.Ms * 1e-4; //This will be in T
    double ExternalField =  h * g.Hstray_sat * 1e-4; //This will be in T
    
    // Allocate WL Functions
    TF1 f__LocalWL("f__LocalWL", WeakLoc, wlPar.minBfield, wlPar.maxBfield, 4);  // Insert fields in Tesla!
    f__LocalWL.SetNpx(g.h_m__Npx);
    f__LocalWL.SetParameters(wlPar.rho0, wlPar.l_phi, wlPar.l_i, wlPar.l_star);
    //double Output = f__LocalWL.Eval(localField);
    double Output = f__LocalWL.Eval(ExternalField+localStrayField);
    
    delete fWL__StrayField_zComp;
    delete fWL__StrayField_zComp_x;
    return( Output );
}


// ------------
// Draw H(x |m, tau, z)
// 	fix z and m, use StrayField_zComp to get the magnetic potential as funct of x
TGraph* CheckStrayField_ProfileX(double FixValue = 1.004e-6, char *variable, float m=0.5, double periodicity_d=1.e-6) {
	GlobalConstants g;
        bool verbose = 1;
	char buffer[300];
        char PlotTitle[300];

        int NumOfPeriods 	= g.LoopingM_MaxNumOfPeriodsForInt;
        double EpsForIntegration= g.LoopingM_EpsForIntegration;
	
        double MaxXforInt = (g.MinX+NumOfPeriods*periodicity_d)/2;
        double MinXforInt = g.MinX;
	if (g.GrapheneWidth) {
		MaxXforInt = g.GrapheneWidth/2;
		MinXforInt = -1*MaxXforInt;
	}
        
        // Allocate functions & Graph
        TF2 *f__StrayField_zComp = new TF2("f__StrayField_zComp", StrayField_zComp, MinXforInt, MaxXforInt, g.t/2, g.MaxZ, 2);
        f__StrayField_zComp->SetParameters(m, periodicity_d);
        TF12 *f__StrayField_zComp_x = new TF12("f__StrayField_zComp_x", f__StrayField_zComp, FixValue, variable);
        f__StrayField_zComp_x->SetNpx(g.h_m__Npx);
	if(variable=="y") f__StrayField_zComp_x->SetNpx(g.h_Stray__Npx2);;
        f__StrayField_zComp_x->SetParameters(m, periodicity_d);
	TGraph *gr = new TGraph(f__StrayField_zComp_x, "");
        
        // Specifics for x- or z- cut
        if (variable == "y") {
            sprintf(PlotTitle, "H (z | x=%.1f, m=%.2f, #tau=%.0e)/M_{s}",FixValue, m, g.tau);
            gr->GetXaxis()->SetTitle("z [m]");
	}
	else 	{
            sprintf(PlotTitle, "H (x | z=%.3e, m=%.2f, #tau=%.0e)/M_{s}",FixValue, m, g.tau);
            gr->GetXaxis()->SetTitle("x [m]");
           
            // Integrate: use Abs(TF1) when computing integral)
            f__StrayField_zComp_x->AbsValue(1);	
            double r1 = f__StrayField_zComp_x->Integral( MinXforInt, MaxXforInt, f__StrayField_zComp_x->GetParameters(), EpsForIntegration) / (MaxXforInt - MinXforInt);
            // Integrate: use TF1 (not Abs(TF1)) when computing integral)
            f__StrayField_zComp_x->AbsValue(0);	
            double r2 = f__StrayField_zComp_x->Integral( MinXforInt, MaxXforInt, f__StrayField_zComp_x->GetParameters(), EpsForIntegration) / (MaxXforInt - MinXforInt);
    
            // Write to Output
            if (verbose) {
                    cout << "[m=" << m << ", z=" << FixValue << "] <|h|>: " << r1 << endl;
                    cout << "[m=" << m << ", z=" << FixValue << "] < h >: " << r2 << " (" << EpsForIntegration << ")" <<  endl;
            }
        }
		
	// Draw TGraph
	gr->SetName("Graph");
	gr->SetTitle(PlotTitle);
	gr->SetMarkerStyle(g.MStyle);
	gr->SetMarkerSize(g.MSize);
	gr->GetYaxis()->SetTitle("");
	gr->Draw("alp");
	
	delete f__StrayField_zComp;
	delete f__StrayField_zComp_x;
        
        return(gr);
}


// ------------
// Draw H(x,z |m, tau)
// 	fix m, use Eq.A6 to get the magnetic potential, derive it
//	NB: in the paper this is a function of Ms, here we renormalize: phi(x,z)-->phi(x,z)/Ms
void Check__StrayField_zComp(float m, float Periodicity_m) { // *variable is either "x" or "y"
	GlobalConstants g;
	char buffer[300];
        int NumOfPeriods 	    = g.LoopingM_MaxNumOfPeriodsForInt;
        double MinXforInt = g.MinX;
        double MaxXforInt = g.MinX+NumOfPeriods*Periodicity_m;
        
	// Stray Field(x,z |m, d): allocate & initialize function
	TF2 *f__StrayField_zComp = new TF2("f__StrayField_zComp", StrayField_zComp, MinXforInt, MaxXforInt, g.t/2, g.MaxZ, 2);
        
        f__StrayField_zComp->SetTitle("H_{str}(x,z)/M_{s}");
	f__StrayField_zComp->SetNpx(g.NpxFor2D);
	f__StrayField_zComp->SetNpy(g.NpyFor2D);
	f__StrayField_zComp->SetParameters(m, Periodicity_m);
    f__StrayField_zComp->SetLineWidth(1);
	f__StrayField_zComp->Draw("surf3z");
	
	return();
}


Double_t TotField_zComp(Double_t *Variable, Double_t *par) { 
	GlobalConstants g;
	double m 			= par[0];
	double Periodicity_m = par[1];
	double x = Variable[0];
	double z = Variable[1];
	
	// Compute H-field from YIG
	TF2 f__TotField("f__TotField", StrayField_zComp, g.MinX, g.MaxX, g.t/2, g.MaxZ, 2);
	f__TotField.SetParameters(m, Periodicity_m);
	double YigField = f__TotField.Eval(x, z);
	
	// Compute external field (needed to reach magnetization m)
	TF1 f__TotField__h_m("f__TotField__h_m", h_md, g.Min_d, g.Max_d, 1);
	f__TotField__h_m.SetParameter(0, Periodicity_m);
	double ExternalField = f__TotField__h_m.Eval(m);

	return( ExternalField + YigField );
}


// ------------
// H(x,z |m, tau)
// 	fix m, use Eq.A6 to get the magnetic potential, derive it vs z
//	NB: in the paper this is a function of Ms, here we renormalize: phi(x,z)-->phi(x,z)/Ms
Double_t StrayField_zComp(Double_t *Variable, Double_t *par) { 
	GlobalConstants g;
	double m = par[0];
	double Periodicity_m = par[1];
	double x = Variable[0];
	double z = Variable[1];
	
	// phi(x,z |m, d): allocate & initialize function
	TF2 *f__phi = new TF2("f__phi", phi, g.MinX, g.MaxX, g.t/2, g.MaxZ, 2);
	f__phi->SetParameters(m, Periodicity_m);

	// phi(x,z=t |m, d): allocate & initialize function
	TF12 *f__phi_z = new TF12("f__phi_x(z=t)", f__phi, x, "y");
	f__phi_z->SetNpx(g.h_m__Npx);
	f__phi_z->SetParameters(m, Periodicity_m);
	
	double Output = -1.*(f__phi_z->Derivative(z));
	delete f__phi;
	delete f__phi_z;
	
	return( Output );
}





// =========== ===========
//
// Function describing the physics of the Magnetization reversal
//
// =========== ===========


// ------------
// Draw phi(x |m, tau, z)
// 	fix z and m, use Eq.A6 to get the magnetic potential as funct of x
//	NB: in the paper this is a function of Ms, here we renormalize: phi(x,z)-->phi(x,z)/Ms
TGraph* Check__phi_x(float FixValue = 1e-6, char *variable, float m, double Periodicity_m, char* DrawOption="alp") { // *variable is either "x" or "y"
	GlobalConstants g;
	char buffer[300];
        int NumOfPeriods 	    = g.LoopingM_MaxNumOfPeriodsForInt;
        double MinX = g.MinX;
        double MaxX = g.MinX+NumOfPeriods*Periodicity_m; 
        
	// phi(x,z |m, d): allocate & initialize function
	TF2 *f__phi = new TF2("f__phi", phi, MinX, MaxX, g.t/2, g.MaxZ, 2);
	f__phi->SetParameters(m, Periodicity_m);

	// phi(x,z=t |m, d): allocate & initialize function
	TF12 *f__phi_x = new TF12("f__phi_x(z=t)", f__phi, FixValue, variable);
	f__phi_x->SetNpx(g.h_m__Npx);
	f__phi_x->SetParameters(m, Periodicity_m);
	
	// phi(x,z=t |m, d): plot
	TGraph *gr = new TGraph(f__phi_x, "");
	gr->SetName("Graph");
	if (variable == "x") {
		sprintf(buffer, "#phi(x | z=%.3e, m=%.2f, #tau=%.0e)/M_{s}",FixValue, m, g.tau);
		gr->GetXaxis()->SetTitle("x [m]");
	}
	else 	{
		sprintf(buffer, "#phi(z | x=%.1f, m=%.2f, #tau=%.0e)/M_{s}",FixValue, m, g.tau);
		gr->GetXaxis()->SetTitle("z [m]");
	}
	gr->SetTitle(buffer);
	gr->SetMarkerStyle(g.MStyle);
	gr->SetMarkerSize(g.MSize);
	gr->GetYaxis()->SetTitle("");
	gr->Draw(DrawOption);
	
	delete f__phi;
	delete f__phi_x;
	
	return(gr);
}


// ------------
// phi(x,z |m)
// 	fix m, use Eq.A6 to get the magnetic potential
//	NB: in the paper this is a function of Ms, here we renormalize: phi(x,z)-->phi(x,z)/Ms
//  NB: there is a sinh term that will diverge for 2kt = d, so keep t small with respect to d
Double_t phi(Double_t *Variable, Double_t *par) { 
    GlobalConstants g;
    double pi = TMath::Pi(), d, Sum=0., SumForConvergence=0, Addendum, AddendumForConvergence;
    int    MaxN = g.MaxNforSumConvergence;
    double convergenceFactor = g.Convergence__phi;
    double x = Variable[0];
    double z = Variable[1];
    double m = par[0];
    double d = par[1];
    double tau = g.tau;
    bool Converged = 0;

        
    if (fabs(m)<1) {
        // Compute the Summation
        //   NB: phi() is an oscillating functions, so it is not safe to simply wait for Addendum/Small < convergence, you are safer by adding a fixed num of times.
        for (int k=1; k<MaxN; k++){
            // The bit to check for convergence (without the oscillating term)        
            AddendumForConvergence = 2.*d/pow(pi*k,2) * TMath::SinH(pi*k*g.t/d)  * TMath::Exp(-2*pi*k*z/d) ;
            // The bit you need
            Addendum = AddendumForConvergence * TMath::Sin( 0.5*pi*k*(m+1) ) * TMath::Cos(2*pi*k*x/d);
            
            Sum += Addendum;
            SumForConvergence += AddendumForConvergence;
            if (fabs(AddendumForConvergence) < fabs(SumForConvergence*convergenceFactor)) Converged = 1;
            if (Converged) break;
        }
        if (!Converged) cout << "   --!! phi() did not converge [m =" << m << ", MaxN=" << MaxN << ", d=" << d << "]" << endl;
    }
	
    return( 0.5*m*g.t + Sum );
}



// ------------
// Plot h(m)
//	you can plot nfunc of them, to check for convergence
void Checks__h_m(const int nfunc=1) {
	cout << "-- Checks: h_m" << endl;
	GlobalConstants g;
	char buffer[300];
	int color =1;
	TF1 **f__h_m = new TF1*[nfunc];
	TGraph **gr__h_m = new TGraph*[nfunc];
	if (nfunc>1) TLegend *leg = new TLegend(0.5, 0.5, 0.9, 0.9);
	
	for (Int_t i=0;i<nfunc;i++) {
		//int MaxNforSumConvergenceForThisLoop = g.MaxNforSumConvergence/(2*i+1);
		//sprintf(buffer,"MaxN:%i",MaxNforSumConvergenceForThisLoop);
		float tau = g.tau * pow(10,i);
		int MaxNforSumConvergenceForThisLoop = g.MaxNforSumConvergence;
		sprintf(buffer,"#tau:%.1e", tau);
		cout << buffer << endl;
		f__h_m[i] = new TF1(buffer,h_m, g.Min_m, g.Max_m, 2);
		f__h_m[i]->SetParameters(MaxNforSumConvergenceForThisLoop, tau);
		f__h_m[i]->SetNpx(g.h_m__Npx);
		gr__h_m[i] = new TGraph(f__h_m[i], "");
		gr__h_m[i]->SetMarkerStyle(g.MStyle);
		gr__h_m[i]->SetMarkerSize(g.MSize);
		gr__h_m[i]->SetLineColor(color);
		gr__h_m[i]->SetMarkerColor(color++);
		gr__h_m[i]->GetXaxis()->SetTitle("m");
		gr__h_m[i]->GetYaxis()->SetTitle("h");
		if (nfunc>1) leg->AddEntry(gr__h_m[i], gr__h_m[i]->GetName());

		if (i == 0) {
			gr__h_m[i]->Draw("alp");
			gr__h_m[i]->SetName("Graph");
		}
		else        gr__h_m[i]->Draw("lp");
	}
	if (nfunc>1) leg->Draw();
	
	for (Int_t i=0;i<nfunc;i++) delete f__h_m[i];
	delete f__h_m[];
}


void Checks__d_m(int CheckD_m){
    cout << "-- Checks__d_m" << endl;
    int color = 1;
    char buffer[300];
    GlobalConstants g;
    TLegend *leg = new TLegend(.36, .55, .70, .86);
    float StartingTau = g.tau;
    
    // Allocate function
    TF1 *f_CheckD_m = new TF1 ("f_CheckD_m", d_m, g.Min_m, g.Max_m, 1);
    TMultiGraph* multiG = new TMultiGraph("multigraph", "");
    
    for (int i=0; i<CheckD_m; i++) {
        // Update Tau value
        GlobalConstants::tau = StartingTau + i*StartingTau;
        sprintf(buffer,"#tau:%.1e", g.tau);
        cout << buffer << endl;
        
        //Allocate functions
        f_CheckD_m->SetParameter(0, g.tau);
        f_CheckD_m->SetNpx(g.h_m__Npx);
        
        //Plot the graph
        TGraph *gr = new TGraph(f_CheckD_m);
        gr->SetName("Graph");
         	gr->SetMarkerStyle(g.MStyle); 	gr->SetMarkerSize(g.MSize);
        sprintf(buffer, "d(m |#tau = %.2e)", g.tau);
        gr->SetTitle(buffer);
        gr->SetMarkerColor(color);
        gr->SetLineColor(color++);
        multiG->Add(gr);
        leg->AddEntry(gr, gr->GetTitle(), "lp");
    }
    
    multiG->Draw("apl");
    multiG->GetXaxis()->SetTitle("m");
    leg->Draw();
    GlobalConstants::tau =StartingTau;
    
    delete f_CheckD_m;
    return();
}


// ------------
// Optional: Check dependence of h_s on tau (should not vary much)
// fix m at some 0.99 and plot h_tau(tau)
void Checks__h_t(float m){
    GlobalConstants g;
    char buffer[300];
    
    
    // Compute d    
    TF1 h_tau__d_m ("h_m__fEq5", d_m, g.Min_m, g.Max_m, 2);
    h_tau__d_m.SetParameter(0, g.tau);
    double Periodicity_m = h_tau__d_m.Eval(m);
    
    TF1 *fh_tau = new TF1("fh_tau", h_tau, g.MinTau, g.MaxTau, 2);
    fh_tau->SetParameters(m, Periodicity_m);
    fh_tau->SetNpx(g.h_m__Npx);
    TGraph *gr = new TGraph(fh_tau, "");
    gr->SetName("Graph");
    gr->GetXaxis()->SetLimits(g.MinTau,g.MaxTau);
    sprintf(buffer, "h_{m=%.2f}(#tau)", m);
    gr->SetTitle(buffer);
     	gr->SetMarkerStyle(g.MStyle); 	gr->SetMarkerSize(g.MSize);
    gr->GetXaxis()->SetTitle("#tau [m]");
    gr->Draw("alp");
    delete fh_tau;
  }


// ------------
// h(tau |m, d)
// fix m, use Eq.5 to find d(tau) and Eq.6 to find h(tau) = h( m(d(tau)) )
Double_t h_tau(Double_t *x, Double_t *par) { // function of tau for a fixed m
    GlobalConstants g;
    double pi = TMath::Pi();
    int    MaxN = g.MaxNforSumConvergence;
    double tau = x[0];
    double m = par[0];
    double d = par[1];

    // Compute h(m)  
    TF1 h_tau__Eq6("h_m__Eq6", h_md, g.Min_d, g.Max_d, 1);
    h_tau__Eq6.SetParameter(0, d);
    
    return( h_tau__Eq6.Eval(m) );
}



// ------------
// h(m |d, tau)
// 	fix tau, use Eq.5 to find d(m) and Eq.6 to find h(m) = h( m, d(m) )
//	MaxN is a parameter here, so you can check for convergence
Double_t h_m(Double_t *x, Double_t *par) { // function of m for a fixed tau
    GlobalConstants g;
    double pi = TMath::Pi(), d;
    double m = x[0];
    int    MaxN = par[0];
    double tau = par[1]; //g.tau;
	
    
        
    // Compute d    
    TF1 h_m__d_m ("h_m__fEq5", d_m, g.Min_m, g.Max_m, 1);
    h_m__d_m.SetParameter(0, tau);
    d = h_m__d_m.Eval(m);
    
    // Compute h(m)  
    TF1 h_m__Eq6("h_m__Eq6", h_md, g.Min_d, g.Max_d, 1);
    h_m__Eq6.SetParameter(0, d);
	double Output = h_m__Eq6.Eval(m);
	
	//cout << "-- m: " << m << ", d: " << d << ", h_m:" << Output << endl;
    return( Output );
}



Double_t d_m(Double_t *x, Double_t *par) {
    GlobalConstants g;
    float m = x[0];
    float tau = par[0];
    
    // Compute d    
    TF1 f__d_m ("f__d_m", Equation5, g.Min_m, g.Max_m, 2);
    f__d_m.SetParameters(m, tau);
    f__d_m.SetNpx(20);	
    
    float Max_dForLocalSearch = g.Max_d;	
    float Min_dForLocalSearch = g.Min_d;	
    
    while ( f__d_m.Eval(Max_dForLocalSearch/2) > 0 ) Max_dForLocalSearch/=2;
    while ( f__d_m.Eval(Min_dForLocalSearch) / f__d_m.Eval(Max_dForLocalSearch) > 1e5 ) Min_dForLocalSearch += g.Min_d;
    double d = f__d_m.GetX(0., Min_dForLocalSearch, Max_dForLocalSearch, 1e-10, 100);

  //cout << "(" << Min_dForLocalSearch << ", " << Max_dForLocalSearch << ")" << endl;

    return( d );
}


// ------------
// h(m |m, d)
// Solve Eq. 6 for a given {m,d} set
Double_t h_md(Double_t *x, Double_t *par) { 
    GlobalConstants g;
    double m = x[0];
    double d = par[0];
    double Sum=0, Addendum, SumForConvergence=0, AddendumForConvergence;    
    double convergenceFactor = g.Convergence__h_md;
    double pi = TMath::Pi();
    int    MaxN = g.MaxNforSumConvergence;
    bool Converged=0;

    if (m<1) {  
        for (int n=1; n<MaxN; n++){
            // The bit to check for convergence (without the oscillating term)        
            AddendumForConvergence  = 1./pow(n*pi,2) * d/g.t * (1- TMath::Exp(-2*pi*n*g.t/d) );
            // The bit you need
            Addendum                = AddendumForConvergence * TMath::Sin(pi*n*(m+1));
            
            Sum += Addendum;
            SumForConvergence += AddendumForConvergence;
            if (fabs(AddendumForConvergence) < fabs(SumForConvergence*convergenceFactor)) Converged = 1;
            if (Converged) break;
        }
        if (!Converged) cout << "   --!! h_md() did not converge [m =" << m << ", MaxN=" << MaxN << ", d=" << d << "]" << endl;
    }
    
    return (m + Sum);    
}


void Check_Equation5(float m){
    GlobalConstants g;
    float OriginalConvergence = GlobalConstants::Convergence__e_d;
    int color = 1;
    char buffer[300];
    
    TLegend *leg = new TLegend(.6, .25, .85, .55);
    
    TF1 *fEq5 = new TF1("fEq5", Equation5, g.Min_d, g.Max_d, 2);
    fEq5->SetParameters(m, g.tau);
	
    TGraph *gr = new TGraph(fEq5, "");
    gr->SetName("Graph");
    gr->SetTitle("#partiale_{d} / #partiald - 2*[0]/(x*x)");
    gr->SetMarkerColor(color);
    gr->SetLineColor(color++);
    gr->GetXaxis()->SetTitle("d [m]");
    gr->Draw("al");
    sprintf(buffer, "Conv.: %.1e", GlobalConstants::Convergence__e_d);
    leg->AddEntry(gr, buffer, "l");
    delete fEq5;
    
    GlobalConstants::Convergence__e_d/=10;
        TF1 *fEq5 = new TF1("fEq5", Equation5, g.Min_d, g.Max_d, 2);
        fEq5->SetParameters(m, g.tau);
        TGraph *gr = new TGraph(fEq5, "");
        gr->SetMarkerColor(color);
        gr->SetLineColor(color++);
        gr->Draw("l");
        sprintf(buffer, "Conv.: %.1e", GlobalConstants::Convergence__e_d);
        leg->AddEntry(gr, buffer, "l");
        delete fEq5;
    GlobalConstants::Convergence__e_d/=10;
        TF1 *fEq5 = new TF1("fEq5", Equation5, g.Min_d, g.Max_d, 2);
        fEq5->SetParameters(m, g.tau);
        TGraph *gr = new TGraph(fEq5, "");
        gr->SetMarkerColor(color);
        gr->SetLineColor(color++);
        gr->Draw("l");
        sprintf(buffer, "Conv.: %.1e", GlobalConstants::Convergence__e_d);
        leg->AddEntry(gr, buffer, "l");
        delete fEq5;
    
    leg->Draw();
    GlobalConstants::Convergence__e_d=OriginalConvergence;
    return();
}


// ------------
// d, m, tau
// Eq. 5; the solution of this (i.e. the -->GetX(0) ) is d for a given {m,tau} set
Double_t Equation5(Double_t *x, Double_t *par) { // function of d for a fixed m \in {-1,1} and tau
    GlobalConstants g;
    int    MaxN = g.MaxNforSumConvergence;
    double m    = par[0];
    double tau  = par[1];
    double d = x[0];
    
    TF1 LocalE_d ("E_d",e_d,g.Min_d, g.Max_d, 1);
    LocalE_d.SetParameter(0, m);
    Double_t Output = LocalE_d.Derivative(d) - 2*tau/(d*d);
    
    return(Output);
}


// ------------
// e_d(d |m)
// demag energy as a function of d for a given m	
Double_t e_d(Double_t *x, Double_t *par) { // function of d for a fixed m \in {-1,1}
    GlobalConstants g;
    double m = par[0];
    float  convergenceFactor = g.Convergence__e_d;
    int    MaxN = g.MaxNforSumConvergence;
    double d = x[0];
    double pi = TMath::Pi();
    double Sum=0, Addendum, SumForConvergence=0, AddendumForConvergence;
    bool Converged=0;
    
    if (fabs(m)>1) return (m*m);
    
    for (int n=1; n<MaxN; n++){
        // The bit to check for convergence (without the oscillating term)        
        AddendumForConvergence = 4./pow(pi*n,3) * d/g.t * (1 - TMath::Exp(-2*pi*n*g.t/d) );
        // The bit you need
	Addendum = AddendumForConvergence * pow( TMath::Sin(0.5*pi*n*(m+1)) ,2);
        
        Sum += Addendum;
        SumForConvergence += AddendumForConvergence;
        if (fabs(AddendumForConvergence) < fabs(SumForConvergence*convergenceFactor)) Converged = 1;
    //    cout << "ed:" << AddendumForConvergence << " " << SumForConvergence << " "<< AddendumForConvergence/SumForConvergence <<" "<< Converged << endl;
        if (Converged) {
    //    cout << " breaking at " << n << "/" << MaxN << endl;
            break;
        }
    }
    if (!Converged) cout << "   --!! e_d() did not converge [m =" << m << ", MaxN=" << MaxN << ", d=" << d << "]" << endl;
    
    return( pow(m,2) + Sum);
}



// =========== ===========
//
// Functions reading from measurements in ourder to compare with simulation
//
// =========== ===========

TNtuple* temp(float rr){
    TCanvas *can = new TCanvas();
    can->Divide(4,2);
    int padindex=1;
    WL_parameters ww;
	GlobalConstants g;
    char buffer[300];
	
	can->cd(padindex++);
	TNtuple *Ndata = new TNtuple ("Ndata", "Ndata", "H:rho:temp:ang:curr");
	Ndata->ReadFile("Measurements1.dat");
	Ndata->SetMarkerStyle(g.MStyle);
	Ndata->SetMarkerSize(g.MSize);
	Ndata->Draw("rho/4:H/2500", "temp<3 && ang<20", "goff");
	// make the TGraph
	TGraph *gg = new TGraph(Ndata->GetSelectedRows(), Ndata->GetV2(), Ndata->GetV1());
	gg->GetXaxis()->SetRangeUser(-0.1, 0.8);
	gg->SetTitle("#rho at 2K for M_{s}=2500 G");
	gg->GetXaxis()->SetTitle("H/M_{s}");
	gg->Draw("apl");
	
	
	string conditions = "abs(rho0-631)<10";
	
	TNtuple *n = new TNtuple ("n", "n", g.LogFile_ListOfColumnNames.c_str());
	n->ReadFile(g.LogFileTitle.c_str());
	
	can->cd(padindex++);
	n->Draw("l_phi", conditions.c_str());
	can->cd(padindex++);
	n->Draw("l_i", conditions.c_str());
	can->cd(padindex++);
	n->Draw("z", conditions.c_str());
	can->cd(padindex++);
	n->Draw("l_phi:z", "l_phi<2.5e-7", "	");
	
	
	float MinRho=618, MaxRho=628;
	sprintf(buffer, "abs(rho0-%f)<1", rr);
	string generalCond = "abs(l_i-3.3e-7)<1e-8 && (abs(l_phi-2e-7)<1e-8 || abs(l_phi-2.6e-7)<1e-8 || abs(l_phi-3.4e-7)<1e-8 || abs(l_phi-4.e-7)<1e-8) && ";
	string IndividualConditions[3] = {buffer, "abs(z-950e-9)<1e-9", "abs(l_phi-4e-7)<1e-8"};
	
	conditions = generalCond + IndividualConditions[0] + " && " + IndividualConditions[1];
	can->cd(padindex++);
	TMultiGraph *mult = PlotSomethingForAllZ(can, padindex, "AverRho:h", "l_phi", conditions.c_str());
		mult->Add(gg);
		mult->GetYaxis()->SetRangeUser(MinRho, MaxRho);
		mult->GetXaxis()->SetRangeUser(-0.05, 0.8);
	can->cd(padindex++);
	conditions = generalCond + IndividualConditions[0]  + " && " +  IndividualConditions[2];
	TMultiGraph *mult = PlotSomethingForAllZ(can, padindex, "AverRho:h", "z", conditions.c_str());
		mult->Add(gg);
		mult->GetYaxis()->SetRangeUser(MinRho, MaxRho);
		mult->GetXaxis()->SetRangeUser(-0.05, 0.8);
	can->cd(padindex++);
	conditions = generalCond + IndividualConditions[1]  + " && " +  IndividualConditions[2];
	TMultiGraph *mult = PlotSomethingForAllZ(can, padindex, "AverRho:h", "rho0", conditions.c_str());
		mult->Add(gg);
		mult->GetYaxis()->SetRangeUser(MinRho, MaxRho);
		mult->GetXaxis()->SetRangeUser(-0.05, 0.8);
		mult->GetHistogram()->SetTitle("#rho");
		mult->GetXaxis()->SetTitle("H/M_{s}");
	
	FormatCanvas(can, 10);
	return(n);
}


void FitTheData(){
	TNtuple *Ndata = new TNtuple ("Ndata", "Ndata", "H:rho:temp:ang:curr");
	Ndata->ReadFile("Measurements2.dat");
	Ndata->SetMarkerStyle(g.MStyle);
	Ndata->SetMarkerSize(g.MSize);
	Ndata->Draw("rho/4:H/1800", "temp<3 && ang<20", "");
		TGraph *gg = new TGraph(Ndata->GetSelectedRows(), Ndata->GetV2(), Ndata->GetV1());
		gg->SetTitle("#rho at 2K");
		gg->GetXaxis()->SetTitle("H/M_{s}");
		gg->SetMarkerStyle(7);
		gg->Draw("apl");

	TF1 *fitting = new TF1("fitting", ComputeAverageRho, -1, 1, 3);
	fitting->SetParameters(1e-6, 628, 2e-7);
	fitting->SetParLimits(0, 860e-9, 1.1e-6);
	fitting->SetParLimits(1, 625, 640);
	fitting->FixParameter(2, 2e-7);
	//fitting->SetParLimits(2, 1e-7, 5e-7);
		
}
