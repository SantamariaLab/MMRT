:Reference : Modeled according to kinetics derived from Magistretti & Alonso 1999
:Comment: corrected rates using q10 = 2.3, target temperature 34, orginal 21

NEURON	{
	SUFFIX Nap
	USEION na READ ena WRITE ina
	RANGE gbar, g, ina, dH, Cp, T0, dS
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gbar = 0.00001 (S/cm2)
 	T0 = 298.15
	Cp=-2860 
	dH=9.011e+04
	dS=62 : at 21 C
}

ASSIGNED	{
	v	(mV)
	ena	(mV)
	ina	(mA/cm2)
	g	(S/cm2)
	celsius (degC)
	mInf
	hInf
	hTau
	hAlpha
	hBeta
	
}

STATE	{
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	rates()
	g = gbar*mInf*h
	ina = g*(v-ena)
}

DERIVATIVE states	{
	rates()
	h' = (hInf-h)/hTau
}

INITIAL{
	rates()
	h = hInf
}

PROCEDURE rates(){
  LOCAL qt, T, kb, hplank, R
  T= celsius + 273.15
  kb = 1.38064852e-23 :m2 kg s-2 K-1
  hplank = 6.62607004e-34 :m2 kg / s
  R= 8.314
  qt=((kb*T/hplank)*exp(-(dH+Cp*(T-T0))/(R*T)+(dS+Cp*(log(T)-log(T0)))/R))
  :  	qt =  2.3^((celsius-21)/10)

  UNITSOFF
	mInf = 1.0/(1+exp((v- -52.6)/-4.6)) : assuming instantaneous activation as modeled by Magistretti and Alonso

	hInf = 1.0/(1+exp((v- -48.8)/10))
	hAlpha = 2.88e-6 * vtrap(v + 17, 4.63)
	hBeta = 6.94e-6 * vtrap(-(v + 64.4), 2.63)

	hTau = (1/(hAlpha + hBeta))/qt
  UNITSON
}

FUNCTION vtrap(x, y) { : Traps for 0 in denominator of rate equations
	UNITSOFF
	if (fabs(x / y) < 1e-6) {
		vtrap = y * (1 - x / y / 2)
	} else {
		vtrap = x / (exp(x / y) - 1)
	}
	UNITSON
}
