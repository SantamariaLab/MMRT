: Reference: Colbert and Pan 2002

NEURON	{
	SUFFIX NaTs
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

	malphaF = 0.182
	mbetaF = 0.124
	mvhalf = -40 (mV)
	mk = 6 (mV)

	halphaF = 0.015
	hbetaF = 0.015
	hvhalf = -66 (mV)
	hk = 6 (mV)
 	T0 = 298.15
        Cp=-2860 
        dH=9.011e+04
        dS=59.6 : at 23 C

	Cp_h=-4144
        dH_h=8.898e+04
        dS_h=55.8 : at 23 C

}

ASSIGNED	{
	v	(mV)
	ena	(mV)
	ina	(mA/cm2)
	g	(S/cm2)
	celsius (degC)
	mInf
	mTau
	mAlpha
	mBeta
	hInf
	hTau
	hAlpha
	hBeta
}

STATE	{
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	g = gbar*m*m*m*h
	ina = g*(v-ena)
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
  LOCAL qt,qth, T, kb, hplank, R
  T= celsius + 273.15
  kb = 1.38064852e-23 :m2 kg s-2 K-1
  hplank = 6.62607004e-34 :m2 kg / s
  R= 8.314
  qt=((kb*T/hplank)*exp(-(dH+Cp*(T-T0))/(R*T)+(dS+Cp*(log(T)-log(T0)))/R))

  qth=((kb*T/hplank)*exp(-(dH_h+Cp_h*(T-T0))/(R*T)+(dS_h+Cp_h*(log(T)-log(T0)))/R))

  :qt = 2.3^((celsius-23)/10)

	UNITSOFF
		mAlpha = malphaF * vtrap(-(v - mvhalf), mk)
		mBeta = mbetaF * vtrap((v - mvhalf), mk)

		mInf = mAlpha/(mAlpha + mBeta)
		mTau = (1/(mAlpha + mBeta))/qt

		hAlpha = halphaF * vtrap(v - hvhalf, hk)
		hBeta = hbetaF * vtrap(-(v - hvhalf), hk)

		hInf = hAlpha/(hAlpha + hBeta)
		hTau = (1/(hAlpha + hBeta))/qth
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
