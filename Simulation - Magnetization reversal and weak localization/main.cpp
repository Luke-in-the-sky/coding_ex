#include "./lib/FormatCanvas.cpp"
#include "./lib/Graphene-WeakLocalization.cpp"
#include "./lib/StripeDomainLib.cpp"

// Simulate the magnetization reversal curve for YIG under different Magnetic Field intensities
//
// The core of the simulation is in the function test(float, float, int, double, bool). 
// You can then loop the simulation over for different values of the simulation's parameters
// by using the function ManyTimesTheMain(). Important parameters of the simulation are allocated
// in a dedicated class named GlobalConstants



// ===========
// Define a class to store global constants for the simulation run
// -----------
class GlobalConstants {
  public:
    double t, Min_m;
    Double_t Max_d, Min_d, Max_m, MinTau, MaxTau, MaxX, MinX, MaxZ, Convergence__h_md, Ms, Hstray_sat;
    int h_m__Npx, MaxNforSumConvergence, NpxFor2D, NpyFor2D, LoopingM_MaxNumOfPeriodsForInt;
    string LogFileTitle, LogFile_ListOfColumnNames;
    static string OverallPlottingConditions;
    static double tau;
    static float GrapheneWidth;
    static float GrapheneOffset;
    static float Convergence__e_d;
    static float Convergence__phi;
    static float LoopingM_EpsForIntegration;

};

GlobalConstants::GlobalConstants() {
	
    // physics
    t		= 1.7e-6;
    Min_m 	= -1.;
    Max_m 	= 3.;
    MinTau 	= 1e-10;
    MaxTau 	= 1e-7;
    Ms      = 300; // 4pi Ms in [G]
    Hstray_sat = 90;
	
	// Spacial window
    Min_d = 1e-7;
    Max_d = 5e-5;
    MinX = 0;
    MaxX = 40e-6;
    MaxZ = t/2 + 0.6e-6;// Max_d/100;
	
    // Numerical Resolution
    h_m__Npx = 50;
    NpxFor2D = 50;
    NpyFor2D = 20;
    Convergence__h_md = 1e-5;
    MaxNforSumConvergence = 1e6;
    LoopingM_MaxNumOfPeriodsForInt = 2;
	
    // Others
    LogFileTitle 	= "LogFile-FixWL.txt";
    LogFile_ListOfColumnNames = "m:d:h:z:t:tau:AbsStray_h:Stray_h:AbsTotField:TotField:ConvergeEd:ConvergPhi:NumOfPeriodsInX:IntegrPrecision:GrapheneWidth:GrapheneOffset:AverRho:rho0:l_phi:l_i:l_star:Ms:Hstray_sat";
}
// Initialize global constants
double GlobalConstants::tau             = 1e-7;
float GlobalConstants::Convergence__e_d	= 1e-7;
float GlobalConstants::Convergence__phi	= 1e-5;
float GlobalConstants::LoopingM_EpsForIntegration= 1e-8;
float GlobalConstants::GrapheneWidth = 1e-5 ;      // if GrapheneWidth==0 ==> MinX=0, MaxX=NumOfPeriods*d
float GlobalConstants::GrapheneOffset = 0.;
string GlobalConstants::OverallPlottingConditions = "Hstray_sat>1";



// ===========
// Define a function to run the simulation for a defined number of times,
// each run with a different set of initial parameters
// -----------
void ManyTimesTheMain() {
    float StartingM     = -0.05;
    float MiddleM       = 0.15;
    float FinishingM    = 10.;
    int NumOfpoints     = 10;
    char buffer[300];
        
    // Set values for the Global Constants
    GlobalConstants g;
	GlobalConstants::tau            = 1e-8;
    GlobalConstants::GrapheneWidth  = 10e-6;
    WL_parameters::rho0	            = 626.5;
    WL_parameters::l_phi 	    	= 1.8e-7;
    WL_parameters::l_i 	            = 2.5e-7;
    WL_parameters::l_star	    	= 0.23e-7;
	
    // Run the simulation for a given set of initial parameters
    for (int i=0; i<3; i++) {
		test( -0.1, .95 , 4, g.t/2 + 10*(i+1)*1.e-9, 0);
		test( 1, 10, 3, g.t/2 + 10*(i+1)*1.e-9, 0);
		}
	test(0.1, 10, 0, g.t/2 + 10, 1);
}



// ===========
// Define the core of the simulation run
// -----------
TCanvas* test(float Min_m=0.001, float Max_m=0.01, int NumOf_m_Points=10, double FixZ = 1.001e-6, bool Draw	= 0) {
    // Drawing Options and segments selectors
    float m = (fabs(Min_m) == 1. ? Min_m - 0.01*Min_m/fabs(Min_m) : Min_m);
    bool FormatCan	= 1 *Draw;
    bool DrawE_d	= 0 *Draw;
    bool CheckEq5	= 0 *Draw;
    int CheckD_m        = 0 *Draw;  // Num of graphs with different tau to superimpose
    bool CheckH_tau     = 0 *Draw;
    int CheckH_m	= 0 *Draw;  // Num of graphs with different tau to superimpose
    bool CheckPhi 	= 0 *Draw;
    bool CheckStrayField = 1 *Draw;
    bool CheckLocalResistivity = 0 *Draw;
    bool ReadDataFile   = 1 *Draw;
    // Computing m
    bool ComputeMoreData= (NumOf_m_Points > 0);
	
	// allocate internal objects and initialize Canvas
    GlobalConstants g;
    char buffer[300];
    int padindex=1;
    TGraph* GraphPointer;	
    TCanvas *can = new TCanvas("can", "can", 800, 400);
    can->Divide(3,3);
        
		
    if (DrawE_d) {
        // Compute E_d. Draw it.
        can->cd(padindex++);
        TF1 *DressedE_d= new TF1("DressedE_d",e_d,g.Min_d, g.Max_d, 1);
        DressedE_d->SetParameter(0, m);
        DressedE_d->SetNpx(g.h_m__Npx);
		//  ..make TGraph
		TGraph *gr = new TGraph(DressedE_d);
		gr->SetName("Graph");
		sprintf(buffer, "e_{d} for m=%.2f", m);
		gr->SetTitle(buffer);
		gr->GetXaxis()->SetTitle("d [m]");
		gr->SetMarkerStyle(4);
		gr->Draw("alp");

     
        // Compute derivative of E_d. Draw it.	
        can->cd(padindex++);
        TMultiGraph *mult = new TMultiGraph("multigraph", "");
        mult->SetName("multigraph");
        TLegend *leg = new TLegend(0.6, 0.5, 0.85, 0.85);
		//  ..make TGraph
        TGraph *gr = (TGraph*)DressedE_d->DrawDerivative("lp");
		gr->SetName("Graph");
		gr->SetTitle("#partiale_{d} / #partiald");
		gr->GetXaxis()->SetTitle("d [m]");
		gr->SetMarkerStyle(4);
		leg->AddEntry(gr, gr->GetTitle(), "l");
		mult->Add(gr);

		
		// Compare with 1/x^2
        TF1 *f5 = new TF1("f5", "2*[0]/(x*x)", g.Min_d, g.Max_d);
		f5->SetParameter(0, g.tau);
		f5->SetNpx(g.h_m__Npx);
		TGraph *gr = new TGraph(f5, "");
		delete f5;
		gr->SetName("Graph");
		gr->SetLineColor(4);
		gr->SetMarkerStyle(4);
		leg->AddEntry(gr, gr->GetTitle(), "l");
		mult->Add(gr);
		
		// Draw all TGraphs
		mult->Draw("alp");
        mult->GetHistogram()->GetXaxis()->SetTitle("d [m]");
        leg->Draw();
	}


	// Compute d    
	cout << "compute d" << endl;	
    TF1 *f__d_m = new TF1 ("f__d_m", d_m, g.Min_m, g.Max_m, 1);
    f__d_m->SetParameter(0, g.tau);
    sprintf(buffer, "d(m |#tau = %.1e)", g.tau);
    double periodicity_d = f__d_m->Eval(m);
    cout << "[] d: " << periodicity_d << " for m=" << m << " and tau=" << g.tau << endl;

	
    if (CheckEq5) {
        can->cd(padindex++);
        Check_Equation5(m);
    }  
    
	
    if(CheckD_m){
        can->cd(padindex++);
        Checks__d_m(CheckD_m);
    }

	
    // Compute h(m, d) from Eq.6
    TF1 *fEq6 = new TF1("fEq6", h_md, g.Min_d, g.Max_d, 1);
    fEq6->SetParameter(0, periodicity_d);
    double extField = fEq6->Eval(m);
    cout << "[] h: " << extField << " for m=" << m << " and tau=" << g.tau << endl;
	
    
    // Draw h(tau) for a fixed m from Eq.5 and Eq.6
    if (CheckH_tau) {
		can->cd(padindex++)->SetLogx();
		Checks__h_t(m);
		}

		
    // Draw h(m) combining Eq.5 and Eq.6
    if (CheckH_m) {
		can->cd(padindex++); 
		Checks__h_m(CheckH_m);
	}
	
	
    // Draw phi(x,z) and 2 profiles
    if (CheckPhi) {
        // 2D plot
        can->cd(padindex++);
        float MaxXforPhiPlot	    = g.MinX+g.LoopingM_MaxNumOfPeriodsForInt*periodicity_d; 
        TF2 *f__phi = new TF2("f__phi", phi, g.MinX, MaxXforPhiPlot, g.t/2, g.MaxZ, 2);
        f__phi->SetTitle("#phi(x,z)/M_{s}");
        f__phi->SetParameters(m, periodicity_d);
        f__phi->SetNpx(g.NpxFor2D);
        f__phi->SetNpy(g.NpyFor2D);
        f__phi->SetLineWidth(1);
        f__phi->Draw("surf3z");
        
        TLegend *leg = new TLegend(.58, .44, .83, .74);
        // x-Profile
        can->cd(padindex++);
        float OriginalConvPhi = GlobalConstants::Convergence__phi;
        GraphPointer = Check__phi_x(FixZ, "x", m, periodicity_d, "alp");
		sprintf(buffer, "Conv.: %.1e", GlobalConstants::Convergence__phi);
		leg->AddEntry(GraphPointer, buffer, "l");
        
        // x-Profile, higher accuracy
        GlobalConstants::Convergence__phi=OriginalConvPhi*100;
		GraphPointer = Check__phi_x(FixZ, "x", m, periodicity_d, "l");
		GraphPointer->SetLineColor(3);
		sprintf(buffer, "Conv.: %.1e", GlobalConstants::Convergence__phi);
		leg->AddEntry(GraphPointer, buffer, "l");
        // x-Profile, higher accuracy
        GlobalConstants::Convergence__phi=OriginalConvPhi/10;
		GraphPointer = Check__phi_x(FixZ, "x", m, periodicity_d, "l");
		GraphPointer->SetLineColor(4);
		sprintf(buffer, "Conv.: %.1e", GlobalConstants::Convergence__phi);
		leg->AddEntry(GraphPointer, buffer, "l");
		GlobalConstants::Convergence__phi=OriginalConvPhi;
        leg->Draw();
        
        // z-Profile
        can->cd(padindex++); 
        Check__phi_x(0, "y", m, periodicity_d, "alp");
    }
	
    
    // Draw Stray Field (x,z)
    if (CheckStrayField) {
        can->cd(padindex++);
        CheckStrayField_ProfileX(FixZ, "x", m, periodicity_d);
        can->cd(padindex++);
        CheckStrayField_ProfileX(0., "y", m, periodicity_d);
    }
	
    
    if (ComputeMoreData) {
            Looping_m(Min_m,Max_m,NumOf_m_Points,FixZ);
    }
	
    
    if (CheckLocalResistivity) {
        can->cd(padindex++);
        CheckLocalWL(can, padindex, FixZ);
    }
    
	
    if(ReadDataFile) {	
            ReadFromDataFile(can, padindex);
    }
    
    
    if(FormatCan) FormatCanvas(can,9);
    
    if (Draw) return(can);
    else return();
    
}