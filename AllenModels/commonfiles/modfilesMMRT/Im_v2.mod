: Based on Im model of Vervaeke et al. (2006)

NEURON	{
	SUFFIX Im_v2
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
	dS=49.2 : at 30 C
}

ASSIGNED	{
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	g	(S/cm2)
	celsius (degC)
	mInf
	mTau
	mAlpha
	mBeta
}

STATE	{
	m
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	g = gbar * m
	ik = g * (v - ek)
}

DERIVATIVE states	{
	rates()
	m' = (mInf - m) / mTau
}

INITIAL{
	rates()
	m = mInf
}

PROCEDURE rates() {
  LOCAL qt, T, kb, hplank, R
  T= celsius + 273.15
  kb = 1.38064852e-23 :m2 kg s-2 K-1
  hplank = 6.62607004e-34 :m2 kg / s
  R= 8.314
  qt=((kb*T/hplank)*exp(-(dH+Cp*(T-T0))/(R*T)+(dS+Cp*(log(T)-log(T0)))/R))
  :qt = 2.3^((celsius-30)/10)

  mAlpha = 0.007 * exp( (6 * 0.4 * (v - (-48))) / 26.12 )
  mBeta = 0.007 * exp( (-6 * (1 - 0.4) * (v - (-48))) / 26.12 )

	mInf = mAlpha / (mAlpha + mBeta)
  mTau = (15 + 1 / (mAlpha + mBeta)) / qt
}
