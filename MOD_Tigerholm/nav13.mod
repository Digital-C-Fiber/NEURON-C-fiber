: This HH model is adopted from cummins et al 2001
: URL: https://bbpteam.epfl.ch/svn/analysis/trunk/IonChannel/xmlTomod/CreateMOD.c 
: Revision: 1499
: Date: 2012-01-28 10:45:44 +0100 (Sat, 28 Jan 2012) 
: Author: rajnish 
:Comment :
:Reference :Nav1.3 sodium channels: rapid repriming and slow closed-state inactivation display quantitative differences after expression in a mammalian cell line and in spinal sensory neurons. J. Neurosci., 2001, 21, 5952-61
: it has only two gates m and h implented by Jenny Tigerholm
: stochastisity has been implemented by Joana Roseira


NEURON {
       SUFFIX nav1p3
       USEION na READ ena WRITE ina
       RANGE gbar, ena, ina, celsiusT, Tshift
       RANGE h_eff
    
       POINTER h0_source
}

UNITS {
      (S) = (siemens)
      (mV) = (millivolts)
      (mA) = (milliamp)
}

PARAMETER {
	  gbar = 0.00001 (S/cm2)
          enainit (mV)
          kvot_qt
          celsiusT
	BBiD = 43

}

ASSIGNED {
	 v	(mV) : NEURON provides this
	 ina	(mA/cm2)
	 g	(S/cm2)
	 minf
	 mtau
	 malpha
	 mbeta
	 h_eff
	 h0_source
		ena	(mV)
         
}

STATE { m }

BREAKPOINT {
	   SOLVE states METHOD cnexp

       h_eff = 1 - h0_source
	   g = gbar * m^3 * h_eff
	   ina = g * (v-ena)
}

INITIAL {
	rates(v) : set hInf, mInf
	: assume that equilibrium has been reached
		m = minf
}

DERIVATIVE states {
	   rates(v)
	   m' = (minf - m)/mtau
}

PROCEDURE rates(v (mV)){
	UNITSOFF
		if(v == -26){
			v = v + 0.000001
		}
		malpha = (0.182 * ((v)- -26))/(1-(exp(-((v)- -26)/9)))
		if(v == -26){
			v = v + 0.000001
		}
		mbeta = (0.124 * (-(v) -26))/(1-(exp(-(-(v) -26)/9)))
		
		
		minf = malpha/(malpha + mbeta)
		mtau = 1/(malpha + mbeta)
	UNITSON
    kvot_qt=1/((2.5^((celsiusT-21)/10)))
         mtau=mtau*kvot_qt
}

