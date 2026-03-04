: nattxs_noise1000.mod is a transient ttx-sensitive Na+ current from
: Sheets et al 2007
: with noise implemented from
: Guler 2013 


NEURON {
       SUFFIX nattxs_noise1000					
       USEION na READ ena WRITE ina      
       RANGE gbar, ena, ina, celsiusT, Tshift, NNa, se, phiNa :, m, h, s
}

UNITS {
      (S) = (siemens)		
      (mV) = (millivolts)   
      (mA) = (milliamp)     
}

PARAMETER {
	  se = -1
	  
	  gbar = 0 (S/cm2):0.035135 (S/cm2)  
          enainit (mV)                   
          kvot_qt						 
          celsiusT						 

	  NNa = 1000  
	  tau = 1 (ms)
	  gamNa = 10 
	  wNa2 = 200 
	  TNa = 800

: second commented values are those used in Baker '05
  A_am = 15.5 (/ms)  : 17.235 (/ms) : A for alpha m 
  B_am = -5 (mV)    : 7.58 (mV)
  C_am = -12.08 (mV)   : -11.47 (mV)

  A_ah = 0.38685 (/ms) : 0.23688 (/ms) : A for alpha h
  B_ah = 122.35 (mV)     : 115 (mV)
  C_ah = 15.29 (mV)   : 46.33 (mV)

  A_as = 0.00092 (/ms) : 0.23688 (/ms) : A for alpha h
  B_as = 93.9 (mV)     : 115 (mV)
  C_as = 16.6 (mV)   : 46.33 (mV)

  A_bm = 35.2 (/ms)   : 17.235 (/ms) : A for beta m
  B_bm = 72.7 (mV)    : 66.2 (mV)
  C_bm = 16.7 (mV)    : 19.8 (mV)

  A_bh = 2.00283 (/ms)    : 10.8 (/ms)   : A for beta h
  B_bh = 5.5266 (mV)    : -11.8 (mV)
  C_bh = -12.70195 (mV) : -11.998 (mV)

  A_bs = -132.05 (/ms)    : 10.8 (/ms)   : A for beta h
  B_bs = -384.9 (mV)    : -11.8 (mV)
  C_bs = 28.5 (mV) : -11.998 (mV)

  shift=0 (mV) 
  Tshift=0 (mV) 

}

ASSIGNED {
	v	(mV) : NEURON provides this 
	:celsius (degC)
	ena (mV)
	dt (ms)
	ina (mA/cm2) 
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
	xiNa :(/ms) but dimensionless
}

STATE { 
	m 
	h 
	s 
	qNa
	pNa
} 

BREAKPOINT  {
	SOLVE states METHOD cnexp 	
	g = gbar * (m^3 * h *s + phiNa)
	ina = g * (v-ena) 
}

INITIAL {
	rates(v) : set tau_m, tau_h, hinf, minf  
	: assume that equilibrium has been reached

    :second commented expressions are those in Guler 2013 
    if (se>0) {set_seed(se)}     
	m = minf  :m=am/(am+bm)
	h = hinf  :h=ah/(ah+bh)
	s = sinf  :s=as/(as+bs)
	
	qNa=0
	pNa=0
}

DERIVATIVE states { 
	LOCAL flag, m0, h0, s0
	rates(v)  	
	
	flag=0
	m0=m
	: If the variable leaves [0,1] then random numbers are drawn again
	:second commented expressions are those in Guler 2013
	while (flag==0) {
		m' = (minf - m)/tau_m + etam(am, bm) :m'=am*(1-m)-bm*m+etam(am,bm)
		if (m<0 || m>1) {
			m=m0
		} else {flag=1}		
	}
	
	flag=0
	h0=h
	while (flag==0) {
		h' = (hinf - h)/tau_h + etah(ah, bh) :h'=ah*(1-h)-bh*h+etah(ah,bh)
		if (h<0 || h>1) {
			h=h0
		} else {flag=1}
	}
	
	flag=0
	s0=s
	while (flag==0) {
		s' = (sinf - s)/tau_s + etas(as, bs) :s'=...
		if (s<0 || s>1) { 
			s=s0			
		} else {flag=1}
	}
	
	qNa'= pNa/tau
	pNa'= (-gamNa*pNa/tau-wNa2*(am*(1-m)+bm*m)* qNa) + xiNa
} 


FUNCTION etam (am, bm) (/ms) {
    UNITSOFF
	etam = normrand(0,1/sqrt(3*NNa*dt))*sqrt(am*(1-m)+bm*m)
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
    
    am = A_am/(1+exp((Vm+shift+B_am)/C_am))
    ah = A_ah/(1+exp((Vm+shift+B_ah)/C_ah))
    as = 0.00003 + A_as/(1+exp((Vm+shift+B_as+Tshift)/C_as))
    bm = A_bm/(1+exp((Vm+shift+B_bm)/C_bm))
    bh = -0.00283 + A_bh/(1+exp((Vm+shift+B_bh)/C_bh))
    bs = 132.05 + A_bs/(1+exp((Vm+shift+B_bs+Tshift)/C_bs))

    tau_m = 1.0 / (am + bm)
    minf  = am * tau_m

    tau_h = 1.0 / (ah + bh)
    hinf  = ah * tau_h

    tau_s = 1.0 / (as + bs)
    sinf  = as * tau_s

    kvot_qt = 1/((2.5^((celsiusT-21)/10)))
    tau_m = tau_m*kvot_qt
    tau_h = tau_h*kvot_qt
    tau_s = tau_s*kvot_qt

    xiNa = normrand(0,1/sqrt(tau_m*dt)) * sqrt(gamNa*TNa*(am*(1-m)+bm*m))
    phiNa = sqrt((m^3*(1-m^3))/NNa) * h * s * qNa
	
    UNITSON
}


