: nattxs.mod is a transient ttx-sensitive Na+ current from
: Sheets et al 2007

: nattxs.mod is a transient ttx-sensitive Na+ current from
: Sheets et al 2007

NEURON {
       SUFFIX nattxs
       USEION na READ ena WRITE ina
       RANGE gbar, ena, ina, celsiusT, Tshift
       RANGE h_eff, alphah, betah, gammah, h0
      
      POINTER h0_source :
}

UNITS {
      (S) = (siemens)
      (mV) = (millivolts)
      (mA) = (milliamp)
}

PARAMETER {
	  gbar = 0 (S/cm2):0.035135 (S/cm2)
          enainit (mV)
          kvot_qt
          celsiusT

 
  shift=0 (mV) :10 
  Tshift=0 (mV)

  : use original m
  A_am = 15.5 (/ms)
  B_am = -5 (mV)
  C_am = -12.08 (mV)
  A_bm = 35.2 (/ms)
  B_bm = 72.7 (mV)
  C_bm = 16.7 (mV)

  : original s 
  A_as = 0.00092 (/ms)
  B_as = 93.9 (mV)
  C_as = 16.6 (mV)
  A_bs = -132.05 (/ms)
  B_bs = -384.9 (mV)
  C_bs = 28.5 (mV)

}

ASSIGNED {
	 v	(mV) : NEURON provides this
	 ina	(mA/cm2)
	 g	(S/cm2)

	 tau_m	(ms)
	 minf
   tau_s  (ms)
   sinf
         ena	(mV)
   
   h_eff
   h0_source :
   alphah (/ms)
   betah  (/ms)
   gammah (/ms)

}


STATE { 
        m 
        s
}

 
BREAKPOINT {
	   SOLVE mstates METHOD cnexp    
     :SOLVE hstates METHOD sparse
      
	   h_eff = 1 - h0_source : 

     g = gbar * m^3 * h_eff * s      
	   ina = g * (v-ena)             
}


INITIAL {
	rates(v) : set tau_m, minf, tau_s, sinf

	: assume that equilibrium has been reached
  m = minf
  s = sinf

}

DERIVATIVE mstates {
    rates(v)

    m' = (minf-m)/tau_m
    s' = (sinf-s)/tau_s
}




:original m
FUNCTION alpham(Vm (mV)) (/ms) {
    alpham = A_am/(1+exp((Vm+shift+B_am)/C_am))
}

FUNCTION betam(Vm (mV)) (/ms) {
    betam = A_bm/(1+exp((Vm+shift+B_bm)/C_bm))
}

:original s 
FUNCTION alphas(Vm (mV)) (/ms) {
	 alphas=0.00003+A_as/(1+exp((Vm+shift+B_as+Tshift)/C_as))
}

FUNCTION betas(Vm (mV)) (/ms) {
	 betas=132.05+A_bs/(1+exp((Vm+shift+B_bs+Tshift)/C_bs))
}




: m/s time constant
FUNCTION rates(Vm (mV)) (/ms) {
	 tau_m = 1.0 / (alpham(Vm) + betam(Vm))
         minf = alpham(Vm) * tau_m

   tau_s = 1.0 / (alphas(Vm) + betas(Vm))
         sinf = alphas(Vm) * tau_s


        :Temperature scaling- Q10=2.5, Tref=22, Tref origin=21 
         kvot_qt = 1 / ((2.5^((celsiusT - 21)/10)))    
         tau_m = tau_m * kvot_qt
         tau_s = tau_s * 1 / ((2.5^((celsiusT - 21)/10)))
}