: SK-type calcium-activated potassium current
: Reference : Kohler et al. 1996

NEURON {
       SUFFIX SK
       USEION k READ ek WRITE ik
       USEION ca READ cai
       RANGE gbar, g, ik, dH, Cp, T0, dS

}

UNITS {
      (mV) = (millivolt)
      (mA) = (milliamp)
      (mM) = (milli/liter)
}

PARAMETER {
          v            (mV)
          gbar = .000001 (mho/cm2)
          zTau = 1              (ms)
          ek           (mV)
          cai          (mM)
 	T0 = 298.15
        Cp=-4144
        dH=8.898e+04
        dS= 58.2 :46.5 : at 34 C Q10 correction was not present in original files
}

ASSIGNED {
         zInf
         ik            (mA/cm2)
         g	       (S/cm2)
}

STATE {
      z   FROM 0 TO 1
}

BREAKPOINT {
           SOLVE states METHOD cnexp
           g  = gbar * z
           ik   =  g * (v - ek)
}

DERIVATIVE states {
        rates(cai)
        z' = (zInf - z) / zTau
}

PROCEDURE rates(ca(mM)) {
  LOCAL qt, T, kb, hplank, R
  T= celsius + 273.15
  kb = 1.38064852e-23 :m2 kg s-2 K-1
  hplank = 6.62607004e-34 :m2 kg / s
  R= 8.314
  qt=((kb*T/hplank)*exp(-(dH+Cp*(T-T0))/(R*T)+(dS+Cp*(log(T)-log(T0)))/R))

          if(ca < 1e-7){
	              ca = ca + 1e-07
          }
          zInf = 1/(1 + (0.00043 / ca)^4.8)
	  zTau = 1/qt : assuming parametrization was 34
}

INITIAL {
        rates(cai)
        z = zInf
}
