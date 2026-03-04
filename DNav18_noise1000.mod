: The m and h are form sheets 2007
:s and u are form Delmas
: run NaV18_delmas.m to plot the model

: with noise implemented from
: Guler 2013 


NEURON {
	SUFFIX nav1p8_noise1000
	USEION na READ ena WRITE ina
 	RANGE gbar, ena, ina, NNa, se, phiNa :, m, h, s, u
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	se = -1
	
	gbar = 0 (S/cm2) : =220e-9/(100e-12*1e8) (S/cm2) : 220(nS)/100(um)^2       
	kvot_qt
        celsiusT	
	shift_act = 0 (mV)
	shift_inact =0 (mV)
	
	NNa = 1000  
	tau = 1 (ms)
	gamNa = 10 
	wNa2 = 200 
	TNa = 800
}

ASSIGNED {
	v	(mV) : NEURON provides this
	:celsius (degC)
	ena    	(mV)
	dt (ms)
	ina	(mA/cm2)
	g	(S/cm2)
	
	am (/ms)
    bm (/ms)
	as (/ms)
	bs (/ms)
	ah (/ms)
	bh (/ms)
	au (/ms)
	bu (/ms)
	tau_h	(ms)
   	tau_m	(ms)
	tau_s	(ms)
	tau_u	(ms)
	minf
	hinf
	sinf
	uinf
    
	phiNa
	xiNa 
}

STATE {
	m
	h
	s 
	u
	qNa
	pNa
}

BREAKPOINT {
	SOLVE states METHOD cnexp	
	g = gbar * (m^3* h * s * u + phiNa)
	ina = g * (v-ena)
}

INITIAL {
	
    rates(v)
	: assume that equilibrium has been reached
	
	if (se>0) {set_seed(se)}
    m=minf
    h=hinf
    s=sinf
    u=uinf
 
	qNa=0
	pNa=0
}

DERIVATIVE states {
	LOCAL flag, m0, h0, s0, u0
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
    while (flag == 0) {
        s' = (sinf - s)/tau_s + etas(as, bs)
        if (s < 0 || s > 1) {
			s = s0 
		} else { flag = 1 }
    }

    flag=0
    u0=u
    while (flag == 0) {
        u' = (uinf - u)/tau_u + etau(au, bu)
        if (u < 0 || u > 1) {
			u = u0
		} else { flag = 1 }
    }

    qNa'= pNa/tau
    pNa'= (-gamNa*pNa/tau-wNa2*(am*(1-m)+bm*m)*qNa) + xiNa
}

FUNCTION etam(am, bm) (/ms) {
    UNITSOFF
    etam = normrand(0, 1/sqrt(3*NNa*dt))*sqrt(am*(1-m)+bm*m)
    UNITSON
}

FUNCTION etah(ah, bh) (/ms) {
    UNITSOFF
    etah = normrand(0, 1/sqrt(NNa*dt))*sqrt(ah*(1-h)+bh*h)
    UNITSON
}

FUNCTION etas(as, bs) (/ms) {
    UNITSOFF
    etas = normrand(0, 1/sqrt(NNa*dt))*sqrt(as*(1-s)+bs*s)
    UNITSON
}

FUNCTION etau(au, bu) (/ms) {
    UNITSOFF
    etau = normrand(0, 1/sqrt(NNa*dt))*sqrt(au*(1-u)+bu*u)
    UNITSON
}


PROCEDURE rates(Vm (mV)) {
	UNITSOFF
	
    am= 2.85-(2.839)/(1+exp((Vm-1.159)/13.95))
    bm= (7.6205)/(1+exp((Vm+46.463)/8.8289))
	tau_m = 1/(am+bm)
	minf = am/(am+bm)
    
    hinf= 1/(1+exp((Vm+32.2)/4))  
	tau_h=(1.218+42.043*exp(-((Vm+38.1)^2)/(2*15.19^2)))
	ah = (1/(1+exp((Vm+32.2)/4))) / ((1.218+42.043*exp(-((Vm+38.1)^2)/(2*15.19^2)))) 
	: ah = hinf/tau_h
	bh = (1- 1/(1+exp((Vm+32.2)/4))) / ((1.218+42.043*exp(-((Vm+38.1)^2)/(2*15.19^2)))) 
	: bh = (1- hinf) / tau_h
	
	tau_s = 1/(as+bs)			
    sinf = 1/(1 + exp((Vm + 45)/8(mV)))
	as=	0.001(/ms)*5.4203/(1 + exp((Vm + 79.816)/16.269(mV))) : 
	bs= 0.001(/ms)*5.0757/(1 + exp(-(Vm + 15.968)/11.542(mV))) : 
	 	
 	tau_u = 1/(au + bu)	
	uinf = 1/(1 + exp((Vm + 51)/8(mV)))
	au=	0.0002 *2.0434/(1 + exp((Vm + 67.499)/19.51))
	bu= 0.0002 *1.9952/(1 + exp(-(Vm + 30.963)/14.792))

	kvot_qt=1/((2.5^((celsiusT-22)/10)))
        tau_m=tau_m*kvot_qt
        tau_h=tau_h*kvot_qt
        tau_s=tau_s*kvot_qt
        tau_u=tau_u*kvot_qt
		
    xiNa = normrand(0,1/sqrt(tau_m*dt)) * sqrt(gamNa*TNa*(am*(1-m)+bm*m))
    phiNa = sqrt((m^3*(1-m^3))/NNa) * h * s * u* qNa
	
    UNITSON
}



