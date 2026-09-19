: nav13_h.mod for Nav1.3 multi-states h updating
: alpha,beta,gamma from Koester2025

NEURON {
    SUFFIX nav1p3_h
    RANGE h0, h_eff
    RANGE alphah, betah, gammah
    RANGE celsiusT
}

UNITS {
    (mV) = (millivolts)
}

PARAMETER {
    celsiusT
}

ASSIGNED {
    v (mV)

    h_eff
    alphah (/ms)
    betah  (/ms)
    gammah (/ms)
}

STATE{
    h0
    h1
    hw1
    hw2
    hw3
    hw4
    hw5
}

BREAKPOINT {
    SOLVE hstates METHOD sparse
    h_eff = 1 - h0
}

INITIAL {
    rates_h(v)
    SOLVE hstates STEADYSTATE sparse

    h_eff = 1 - h0
}

KINETIC hstates {
    rates_h(v)
    ~ h0  <-> h1  (alphah, 0)
    ~ h1  <-> hw1 (betah, 0)
    ~ hw1 <-> hw2 (gammah, 0)
    ~ hw2 <-> hw3 (gammah, 0)
    ~ hw3 <-> hw4 (gammah, 0)
    ~ hw4 <-> hw5 (gammah, 0)
    ~ hw5 <-> h0  (gammah, 0)

    CONSERVE h0+h1+hw1+hw2+hw3+hw4+hw5 = 1
}

PROCEDURE rates_h(v(mV)){    :h transition rate
  LOCAL x, q10fac
  x = v / 1000.0
  q10fac = 2.5^((celsiusT - 22)/10)   

  alphah = (-96.8351 - 140.6073*x) / ( 2.1407 - 4.0967*exp((x - 0.017261)/(-0.069927)))
  betah  = (-36.6121 - 447.9485*x) / (-0.017397 - 0.016165*exp((x + 0.010331)/(-0.01491)))
  gammah = ( 1000 - 999.9749*x) / ( 0.0044704 + (5.4067e-07)*exp((x + 0.062717)/0.00477)) + 26328.3751

  : To avoid artifacts
  if (alphah < 0) {              
      alphah=0
  }
  if (betah < 0) {
      betah=0
  }
  if (gammah < 0) {
      gammah=0
  }
  
  alphah = alphah / 1000.0       : convert rate from 1/s to 1/ms
  betah = betah / 1000.0
  gammah = gammah / 1000.0

  alphah = alphah * q10fac
  betah = betah * q10fac
  gammah = gammah * q10fac
}
