TITLE na-k dynamics



INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX nakdyn
	USEION na READ ina, nai, nao WRITE nai, nao
	USEION k READ ik, ki,ko WRITE ki, ko
        RANGE nais, naos, kis, kos	:

}

UNITS {
	(molar) = (1/liter)			: moles do not appear in units
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)

	F = (faraday) (coulomb)		

}


PARAMETER {

        naiinf  = 11.4  (mM) : 63.4
        naoinf = 154.0  (mM)
        kiinf = 121.7	(mM) : 144.9
        koinf = 5.6   	(mM)

	nai		(mM)
	nao		(mM)
	ki		(mM)
	ko		(mM)

	theta = 14.5e-3 (um)
	D = 0.1e-6 	(m/s)	:Scriven1981
	

}

STATE {

	nais		(mM)
	naos		(mM)
	kis		(mM)
	kos		(mM)
}

INITIAL {

	nai=naiinf
        ki=kiinf
        nao=naoinf
        ko=koinf

	nais = naiinf
	naos = naoinf
	kis  = kiinf
	kos  = koinf

}

ASSIGNED {

	ina		(mA/cm2)
	ik		(mA/cm2)

	diam 		(um)

}
	
BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state { 
	nais' = -ina*4/F/diam*(1e4)
	naos' = (ina/F - D*(0.1)*(nao - naoinf))/theta*(1e4)
:	naos' = (ina*diam/F - D*(0.1)*(nao - naoinf)*(diam+2*theta))/((diam/2+theta)*theta)*(1e4) :newconcdyn
	kis'  = -ik*4/F/diam*(1e4)
	kos'  = (ik/F - D*(0.1)*(ko - koinf))/theta*(1e4)
:	kos' = (ik*diam/F - D*(0.1)*(ko - koinf)*(diam+2*theta))/((diam/2+theta)*theta)*(1e4) :newconcdyn

:	nais' = 0
:	naos' = 0
:	kis'  = 0
:	kos'  = 0


	nai = nais
	nao = naos
	ki  = kis
	ko  = kos


}



