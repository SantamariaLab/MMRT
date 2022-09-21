: Comment: The persistent component of the K current
: Reference:		Voltage-gated K+ channels in layer 5 neocortical pyramidal neurones from young rats:subtypes and gradients,Korngreen and Sakmann, J. Physiology, 2000


NEURON	{
	SUFFIX K_P
	USEION k READ ek WRITE ik
	RANGE gbar, g, ik, dH, Cp, T0, dS
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gbar = 0.00001 (S/cm2)
	vshift = 0 (mV)
	tauF = 1
	T0 = 298.15
	Cp=-4144
        dH=8.898e+04
	dS=58.2 : at 21 C
}

ASSIGNED	{
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
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
	ik = g*(v-ek)
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

PROCEDURE rates() {
  LOCAL qt, T, kb, hplank, R
  T= celsius + 273.15
  kb = 1.38064852e-23 :m2 kg s-2 K-1
  hplank = 6.62607004e-34 :m2 kg / s
  R= 8.314
  qt=((kb*T/hplank)*exp(-(dH+Cp*(T-T0))/(R*T)+(dS+Cp*(log(T)-log(T0)))/R))
  :qt = 2.3^((celsius-21)/10)
	UNITSOFF
		mInf =  1 / (1 + exp(-(v - (-14.3 + vshift)) / 14.6))
    if (v < -50 + vshift){
    	mTau = tauF * (1.25+175.03*exp(-(v - vshift) * -0.026))/qt
    } else {
      mTau = tauF * (1.25+13*exp(-(v - vshift) * 0.026))/qt
    }
		hInf =  1/(1 + exp(-(v - (-54 + vshift))/-11))
		hTau =  (360+(1010+24*(v - (-55 + vshift)))*exp(-((v - (-75 + vshift))/48)^2))/qt
	UNITSON
}
