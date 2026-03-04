: nav1p9.mod is the NaV1.9 Na+ current from
: Herzog, Cummins, and Waxman 2001 p1353
: This current is also called the ttx-rp current
: (the tetrodotoxin resistent persistant current)
: implemented by Tom Morse version 2/25/07
: to plot the model run file NaV19_herzog.m

: with noise implemented from
: Guler 2013 

NEURON {
       SUFFIX nav1p9_noise1000
       USEION na  READ ena WRITE ina
       RANGE gbar, ena, slow_inact, m, h, s, gate, ina, se, NNa, phiNa
       RANGE tau_m, tau_h, tau_s, celsiusT, NGFshift1p9
       : if slow_inact=1 then ultra-slow inactivation is included
}

UNITS {
      (S) = (siemens)
      (mV) = (millivolts)
      (mA) = (milliamp)
}

PARAMETER {

	  se = -1
	  
	  gbar = 0 (S/cm2)
	  
	  NNa = 1000  
	  tau = 1 (ms)
	  gamNa = 10 
	  wNa2 = 200 
	  TNa = 800

	  A_am9 = 1.032 (/ms): 1.548 in Baker '05 : A for alpha m(9 etc ...)
	  B_am9 = 6.99 (mV) : -11.01 in Baker '05
	  C_am9 = -14.87115 (mV): -14.871 (mV)  in Baker '05

	  A_ah9 = 0.06435 (/ms) : 0.2574 (/ms) in Baker '05, A for alpha h
	  B_ah9 = 73.26415 (mV) : 63.264  in Baker '05
	  C_ah9 = 3.71928 (mV) : 3.7193  in Baker '05

	  A_as9 = 0.00000016 (/ms) : contributes to ultra slowness
	  B_as9 = 0 (mV)
	  gate = 0 (mV)
	  C_as9 = 12 (mV)

	  A_bs9 = 0.0005 (/ms)
	  B_bs9 = 32 (mV)
	  C_bs9 = 23 (mV)

	  A_bm9 = 5.79 (/ms) : 8.685  in Baker '05 : A for beta m
	  B_bm9 = 130.4 (mV) : -112.4  in Baker '05
	  C_bm9 = 22.9 (mV) : 22.9  in Baker '05

	  A_bh9 = 0.13496 (/ms) : 0.53984  in Baker '05   : A for beta h
	  B_bh9 = 10.27853 (mV) : 0.27853  in Baker '05
	  C_bh9 = -9.09334 (mV) : -9.0933  in Baker '05
	  
	  slow_inact = 1 (1) : to turn on ultra slow inactivation

        kvot_qt
        celsiusT
	NGFshift1p9 = 0 (mV)
	
}

ASSIGNED {
	 v	(mV) : NEURON provides this
	:celsius (degC)
    :i	(mA/cm2)
	 ena	(mV)
	 dt (ms)
	 ina	(mA/cm2)
	 g	(S/cm2)
	 
	 am	(/ms)
	 ah	(/ms)
	 as	(/ms)
	 bm	(/ms)
	 bh	(/ms)
	 bs	(/ms) 
	 tau_h	(ms)
	 tau_m	(ms)
	 tau_s	(ms)
	 minf
	 hinf
	 sinf

	 phiNa
	 xiNa :(/ms) dimensionless
	 
}

STATE { 
	m 
	h 
	s
	qNa
	pNa
}

BREAKPOINT {
	   SOLVE states METHOD cnexp
	   g = gbar * (m * h * s + phiNa)

	   ina = g * (v-ena)
	   :ina = g * (v-69)
}

INITIAL {
	rates(v) : set time constants and infinity values
	: assume that equilibrium has been reached
	
	if (se>0) {set_seed(se)}  
	m = minf
	h = hinf
	s = sinf
	
	qNa=0
	pNa=0
}

DERIVATIVE states {
	LOCAL flag, m0, h0, s0
	rates(v)
	   
	flag=0
	m0=m
	: If the variable leaves [0,1] then random numbers are drawn again
	while (flag==0) {
		m' = (minf - m)/tau_m + etam(am, bm)
		if (m<0 || m>1) {
			m=m0
		} else {flag=1}		
	}
	
	flag=0
	h0=h
	while (flag==0) {
		h' = (hinf - h)/tau_h + etah(ah, bh)
		if (h<0 || h>1) {
			h=h0
		} else {flag=1}
	}
	
	flag=0
	s0=s
	while (flag==0) {
		s' = (sinf - s)/tau_s + etas(as, bs)
		if (s<0 || s>1) { 
			s=s0			
		} else {flag=1}
	}
	
	qNa'= pNa/tau
	pNa'= (-gamNa*pNa/tau-wNa2*(am*(1-m)+bm*m)* qNa) + xiNa

}

FUNCTION etam (am, bm) (/ms) {
    UNITSOFF
	etam = normrand(0,1/sqrt(NNa*dt))*sqrt(am*(1-m)+bm*m)
	UNITSON
}

FUNCTION etah (ah, bh) (/ms) {
	UNITSOFF
	etah = normrand(0,1/sqrt(NNa*dt))*sqrt(ah*(1-h)+bh*h)
	UNITSON
}

FUNCTION etas (as, bs) (/ms) {
	UNITSOFF
	etas = normrand(0,1/sqrt(NNa*dt))*sqrt(as*(1-s)+bs*s)
	UNITSON
}

	
PROCEDURE rates(Vm (mV)) {
	UNITSOFF
	
	am=A_am9/(1+exp((Vm+B_am9)/C_am9))
	ah=A_ah9/(1+exp((Vm+B_ah9)/C_ah9))
	as=A_as9*(exp(-(Vm+gate+B_as9)/C_as9))
	bm=A_bm9/(1+exp((Vm+B_bm9)/C_bm9))
	bh=A_bh9/(1+exp((Vm+B_bh9)/C_bh9))
	bs=A_bs9/(1+exp(-(Vm+gate+B_bs9)/C_bs9))
		
	tau_m = 1.0 / (am + bm) :tau_m = 1.0 / (am(Vm+NGFshift1p9) + bm(Vm+NGFshift1p9))
	minf = am * tau_m   :minf = am(Vm+NGFshift1p9) * tau_m
 
	tau_h = 1.0 / (ah+ bh) :tau_h = 1.0 / (ah(Vm+NGFshift1p9) + bh(Vm+NGFshift1p9))
	hinf = ah * tau_h   :hinf = ah(Vm+NGFshift1p9) * tau_h

	if (slow_inact) {
	    tau_s = 1.0 / (as + bs)
	    	  sinf = as * tau_s
		  } else {
		    tau_s = 0.1	: in a tenth of a millisecond we move to within
		    	  sinf = 1.0 : 1/e factor towards s = 1
			  }

        kvot_qt=1/((2.5^((celsiusT-21)/10)))
        tau_h=tau_h*kvot_qt
        tau_m=tau_m*kvot_qt
        tau_s=tau_s*kvot_qt
	
	xiNa = normrand(0,1/sqrt(tau_m*dt)) * sqrt(gamNa*TNa*(am*(1-m)+bm*m))
	phiNa = sqrt((m*(1-m))/NNa) * h * s * qNa
    
	UNITSON	
}











