class WL_parameters {
  public:
    static float rho0;
    static float l_phi;
    static float l_i;
    static float l_star;
    static float minBfield;
    static float maxBfield;
    static bool FieldIsInOe;
};
WL_parameters::WL_parameters() {
}

float WL_parameters::rho0		= 627.3 ;
float WL_parameters::l_phi 		= 2.19264e-07;
float WL_parameters::l_i 		= 3.59688e-07;
float WL_parameters::l_star 	= 6.56381e-09;
float WL_parameters::minBfield 	= -0.1;
float WL_parameters::maxBfield 	= 0.1;
bool WL_parameters::FieldIsInOe = 0;


// ========
// Weak Localization formula
// 		Weak Localization for graphene, as in 
//
// Input: 
// Output:	creates 4 TF1 objects, which describe resistivity rho(B) under magnetic field (need to fix 4 parameters)
float WeakLocFormula_GammaMin=1e-3, WeakLocFormula_GammaMax=10.;
float WeakLocFormula_MaxField= 2;

Double_t WeakLoc(Double_t *x, Double_t *par) { // x here
	WL_parameters wl;
	float Rho_0 	= par[0];
	float L_phi=par[1], L_i=par[2], L_star=par[3];
        float h = 6.62606957e-34, e = 1.6e-19; // SI units
        float Constants = e*e/(TMath::Pi()*h); // 1/Ohm
	float A=h/(2*TMath::Pi()*4*e);
	float B_phi = A/(L_phi*L_phi), B_i = A/(L_i*L_i), B_star = A/(L_star*L_star);

	float B= fabs(x[0]);
	if (wl.FieldIsInOe) B= fabs(x[0])*1e-4;
	float FirstTerm 	= DressedF->Eval(B/B_phi);
	float SecondTerm 	= -1.*DressedF->Eval(B/(B_phi + 2*B_i));
	float ThirdTerm 	= -2.*DressedF->Eval(B/(B_phi + B_star + B_i));
	return( Rho_0 -1.*pow(Rho_0, 2)*Constants* ( FirstTerm + SecondTerm + ThirdTerm ) );
}
TF1 *DressedWeakLoc = new TF1("DressedWeakLoc", WeakLoc, -1*WeakLocFormula_MaxField, WeakLocFormula_MaxField, 4);
TF1 *DressedWeakLoc1 = new TF1("DressedWeakLoc", WeakLoc, -1*WeakLocFormula_MaxField, WeakLocFormula_MaxField, 4);
TF1 *DressedWeakLoc2 = new TF1("DressedWeakLoc", WeakLoc, -1*WeakLocFormula_MaxField, WeakLocFormula_MaxField, 4);
TF1 *DressedWeakLoc3 = new TF1("DressedWeakLoc", WeakLoc, -1*WeakLocFormula_MaxField, WeakLocFormula_MaxField, 4);
    
	
	

// ========
// Supporting info
// ========

Double_t LnGamma(Double_t *x, Double_t *par) {
	return (TMath::LnGamma(x[0]));
}
TF1 *DressedLnGamma = new TF1("LnGamma", LnGamma, WeakLocFormula_GammaMin, WeakLocFormula_GammaMax, 0);

Double_t DiGamma(Double_t *x, Double_t *par) {
	return( DressedLnGamma->Derivative(x[0]) );
}
TF1 *DressedDiGamma = new TF1("DiGamma", DiGamma, WeakLocFormula_GammaMin, WeakLocFormula_GammaMax, 0);

Double_t F(Double_t *x, Double_t *par) {
	Double_t input = 0.5 + pow( x[0], -1);
	return( TMath::Log( x[0] ) + DressedDiGamma->Eval(input) );
}
TF1 *DressedF = new TF1("F", F, WeakLocFormula_GammaMin, WeakLocFormula_GammaMax, 0);
