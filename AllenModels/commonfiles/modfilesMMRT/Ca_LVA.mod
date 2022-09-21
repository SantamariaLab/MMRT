: Comment: LVA ca channel. Note: mtau is an approximation from the plots
: Reference:		Avery and Johnston 1996, tau from Randall 1997
: Comment: shifted by -10 mv to correct for junction potential
: Comment: corrected rates using q10 = 2.3, target temperature 34, orginal 21

NEURON	{
	SUFFIX Ca_LVA
	USEION ca READ eca WRITE ica
	RANGE gbar, g, ica, dH, Cp, T0, dS
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gbar = 0.00001 (S/cm2)
	T0 = 298.15
	Cp=-3607
	dH=9.718e+04
	dS=86 : at 21 C
}

ASSIGNED	{
	v	(mV)
	eca	(mV)
	ica	(mA/cm2)
	g	(S/cm2)
	celsius (degC)
	mInf
	mTau
	hInf
	hTau
}

STATE	{
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	g = gbar*m*m*h
	ica = g*(v-eca)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
	h' = (hInf-h)/hTau
}

INITIAL{
	rates()
	m = mInf
	h = hInf
}

PROCEDURE rates(){
  LOCAL qt, T, kb, hplank, R
  T= celsius + 273.15
  kb = 1.38064852e-23 :m2 kg s-2 K-1
  hplank = 6.62607004e-34 :m2 kg / s
  R= 8.314
  qt=((kb*T/hplank)*exp(-(dH+Cp*(T-T0))/(R*T)+(dS+Cp*(log(T)-log(T0)))/R))
  :qt = 2.3^((celsius-21)/10)

	UNITSOFF
		v = v + 10
		mInf = 1.0000/(1+ exp((v - -30.000)/-6))
		mTau = (5.0000 + 20.0000/(1+exp((v - -25.000)/5)))/qt
		hInf = 1.0000/(1+ exp((v - -80.000)/6.4))
		hTau = (20.0000 + 50.0000/(1+exp((v - -40.000)/7)))/qt
		v = v - 10
	UNITSON
}
