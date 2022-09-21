: Kv2-like channel
: Adapted from model implemented in Keren et al. 2005
: Adjusted parameters to be similar to guangxitoxin-sensitive current in mouse CA1 pyramids from Liu and Bean 2014


NEURON	{
	SUFFIX Kv2like
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
	dS=58.2 : at 21 C
}

ASSIGNED	{
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	g	(S/cm2)
	celsius (degC)
	mInf
	mAlpha
	mBeta
	mTau
	hInf
	h1Tau
	h2Tau
}

STATE	{
	m
	h1
	h2
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	g = gbar * m * m * (0.5 * h1 + 0.5 * h2)
	ik = g * (v - ek)
}

DERIVATIVE states	{
	rates()
	m' = (mInf - m) / mTau
	h1' = (hInf - h1) / h1Tau
	h2' = (hInf - h2) / h2Tau
}

INITIAL{
	rates()
	m = mInf
	h1 = hInf
	h2 = hInf
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
		mAlpha = 0.12 * vtrap( -(v - 43), 11.0)
		mBeta = 0.02 * exp(-(v + 1.27) / 120)
		mInf = mAlpha / (mAlpha + mBeta)
		mTau = 2.5 * (1 / (qt * (mAlpha + mBeta)))

		hInf =  1/(1 + exp((v + 58) / 11))
		h1Tau = (360 + (1010 + 23.7 * (v + 54)) * exp(-((v + 75) / 48)^2)) / qt
		h2Tau = (2350 + 1380 * exp(-0.011 * v) - 210 * exp(-0.03 * v)) / qt
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
