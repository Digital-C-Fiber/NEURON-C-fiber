: DNav18_h.mod for Nav1.8 multi-states h updating
: alpha,beta,gamma from Koester

NEURON {
    SUFFIX nav1p8_h
    RANGE h0, h_eff
    RANGE alphah, betah, gammah, celsiusT
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

STATE{                 :==========
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
  
   
  :Tigerholm2014
  alphah = (1/(1 + exp((v + 32.2)/4))) / (1.218 + 42.043*exp(-((v + 38.1)*(v + 38.1))/(2*15.192*15.192)))
  betah = (1 - 1/(1 + exp((v + 32.2)/4))) / (1.218 + 42.043*exp(-((v + 38.1)*(v + 38.1))/(2*15.192*15.192)))
  
  :Koester2025
  :alphah = (206.0724 - 3.371*(x+0.01)) / (294.4345 - 328.0076*exp(((x+0.01) - 0.50016)/3.6785))
  :betah  = (17.4508 + 292.8407*(x+0.01)) / ( 0.058965 + 0.0057358*exp(((x+0.01) - 0.061377)/(-0.030721)))
  gammah = (-37.8993 - 15.3631*x) / (-0.00061406 - 15.2606*exp((x + 0.18321)/(-0.021847))) + 0.0014255

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
  
  :alphah = alphah / 1000.0       : convert rate from 1/s to 1/ms
  :betah = betah / 1000.0 
  gammah = gammah / 1000.0 

  alphah = alphah * q10fac
  betah = betah * q10fac
  gammah = gammah * q10fac
}



