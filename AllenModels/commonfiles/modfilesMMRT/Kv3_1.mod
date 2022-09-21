: Comment: Kv3-like potassium current

NEURON	{
	SUFFIX Kv3_1
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

	T0 = 298.15
	Cp=-4144
        dH=8.898e+04
	dS=46.5 : at 34 C
}

ASSIGNED	{
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	g	(S/cm2)
	mInf
	mTau
}

STATE	{ 
	m
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	g = gbar*m
	ik = g*(v-ek)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
}

INITIAL{
	rates()
	m = mInf
}

PROCEDURE rates(){
  LOCAL qt, T, kb, hplank, R
  T= celsius + 273.15
  kb = 1.38064852e-23 :m2 kg s-2 K-1
  hplank = 6.62607004e-34 :m2 kg / s
  R= 8.314
  qt=((kb*T/hplank)*exp(-(dH+Cp*(T-T0))/(R*T)+(dS+Cp*(log(T)-log(T0)))/R))

	UNITSOFF
		mInf =  1/(1+exp(((v -(18.700 + vshift))/(-9.700))))
		mTau =  0.2*20.000/(qt*(1+exp(((v -(-46.560 + vshift))/(-44.140)))))
	UNITSON
}
