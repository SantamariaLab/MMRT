: Based on Kd model of Foust et al. (2011)


NEURON	{
	SUFFIX Kd
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
	T0 = 298.15
	Cp=-4144
        dH=8.898e+04
	dS=55.8 : at 23 C
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
	g = gbar * m * h
	ik = g * (v - ek)
}

DERIVATIVE states	{
	rates()
	m' = (mInf - m) / mTau
	h' = (hInf - h) / hTau
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
  :qt = 2.3^((celsius-23)/10)
  mInf = 1 - 1 / (1 + exp((v - (-43)) / 8))
  mTau = (1/qt)*1
  hInf = 1 / (1 + exp((v - (-67)) / 7.3))
  hTau = 1500/qt
}
