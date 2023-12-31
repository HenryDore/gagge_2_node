model gagge_2_node
   ////////////////////////
   // Validation Version //
   ////////////////////////
   //
   // This version uses all parameters and equations the same as the calc2node version.
   // 
   // Simulate for 3600s (1 hour)
   
  constant Real clo = 0.5 "clothing insulation level, clo";
  constant Real met = 1 "metabolic rate, met";
  constant Real wme = 0 "external work, met";
  constant Real pb = 760 "barometric pressure, torr or mmHg";
  constant Real ht = 170 "height, cm";
  constant Real wt = 70 "weight, kg";
  constant Real tskn = 33.7 "setpoint value for skin temperature, °C";
  constant Real tcrn = 36.8 "setpoint value for core temperature, °C";
  constant Real tbn = 36.49 "setpoint value for blood temperature (.1*tskn + .9*tcrn), °C";
  constant Real skbfn = 6.3 "neutral value for skin blood flow";
  constant Real sbc = 5.6697 * 10 ^ (-08) "stephan-Boltzmann constant";
  constant Real sa = ((ht * wt) / 3600 ) ^ .5 "surface Area (m2) according to mosteller formula";
  constant Real w = wme * 58.2 "external work, w/m^2";
  constant Real atm = pb / 760 "atmospheric pressure, atm";
  constant Real rcl = 0.155 * clo "thermal resistance of clothing ensemble, °C m^2/W";
  constant Real facl = 1 + 0.15 * clo "Increase in body surface area due to clothing";
  constant Real tu = 40 "turbulence intensity, %";
  constant Real vel = 0.1 "air velocity m/s";
  constant Real rh = 50 "relative humidity %";
  constant Real csw = 170 "driving coefficient for regulatory sweating, 1";
  constant Real cdil = 120 "driving coefficient for vasodilation, 1";
  constant Real cstr = 0.5 "driving coefficient for vasoconstriction, 1";  
  constant Real lr = 2.2 / atm "Lewis Relation is 2.2 at sea level";
  constant Real rmm = 58.2  "metabolic rate, W";

  parameter Real set_temp = 25 "setpoint for ambient temperature";

  Real m (start = met * 58.2) "metabolic rate, W";
  Real ta (start = 20) "air temperature, °C"; 
  Real tr (start = 20) "mean radiant temperature, °C";
  Real pa "partial vapour pressure of water, mmHg";
  Real tsk (start = tskn, fixed = true) "skin temperature, °C";
  Real tcr (start = tcrn, fixed = true) "core temperature, °C";
  Real tcl (start = tskn) "clothing temperature, °C";
  Real tb  (start = 0.1 * tskn + (1 - 0.1) * tcrn);
  Real skbf (start = skbfn) "skin blood flow, kg/hr m^2";
  Real skbf_ "skin blood flow test";
  Real mshiv (start = 0) "rate of energy released by shivering, W";
  Real alpha (start = 0.1) "fractional skin mass, 1";
  Real esk (start = 0.1 * met) "total evaporative heat loss from the skin W";
  Real wcrit "evaporative efficiency, 1";
  Real icl "";
  Real chc (start = 3 * atm ^ (0.53)) "coductive heat transfer coefficient";
  Real chcv (start = 8.600001 * (vel * atm) ^ 0.53) "convective heat transfer coefficient";
  Real chr (start = 4.7) "radiative heat transfer coefficient";
  Real ctc (start = 4.7 + (3 * atm ^ (0.53))) "?combined convection & radiation coefficient";
  Real ra (start = 1 / (facl * 4.7 + (3 * atm ^ (0.53)))) "resistance of air layer to dry heat transfer, °C m^2 / W";
  Real top (start = ((4.7 * 20) + (3 * atm ^ (0.53)) * 20) /  4.7 + (3 * atm ^ (0.53))) "operative temperature, °C";
  Real dry "total sensible heat loss, W";
  Real hfcs "rate of energy transport between core and skin, W";
  Real eres "heat loss through respiratory evaporation, W";
  Real cres "heat loss through respiratory convection, W";
  Real scr "rate of energy storage in the core, W";
  Real ssk "rate of enery storage in the skin, ";
  Real tcsk "alpha*wt*Cb";
  Real tccr "(1-alpha)*wt*Cb";
  Real deltatsk (start = 0) "rate of change of skin temperature, °C" ;
  Real deltatcr (start = 0) "rate of change of core temperature, °C";
  Real sksig "skin signal ,1";
  Real warms "skin warm signal ,1";
  Real colds "skin cold signal ,1";
  Real crsig "core signal ,1";
  Real warmc "core warm signal ,1";
  Real coldc "core cold signal ,1";
  Real bdsig "blood signal ,1";
  Real warmb "blood warm signal ,1";
  Real coldb "blood cold signal ,1";
  Real regsw "regulatory sweating, g/m^2 hr";
  Real ersw "heat loss through sweating, W";
  Real ersw_"also heat loss through sweating, W";
  Real rea "resistance of air layer to dry heat transfer, °C m^2 / W";
  Real recl "resistance of clothing layer to dry heat transfer, °C m^2 / W";
  Real emax "maximum evapourative capacity, W";
  Real pwet "skin wettedness, 1";
  Real prsw "ratio of actual heat loss due to sweating to maximum heat loss due to sweating";
  Real edif "heat loss due to diffusion of water fvapout rhough the skin";
  Real esk_ "also total evaporative heat loss from the skin W";

function tcl_calculate //calculate tcl, chr, ctc, top & ra iteratively
  input Real tcl;
  input Real tr;
  input Real chc;
  input Real facl;
  input Real ta;
  input Real tsk;
  input Real ctc;
  input Real rcl;
  output Real tcl_;
  output Real chr_;
  output Real ctc_;
  output Real top_;  
  output Real ra_;

  protected
  Real tcl_old;
  Integer finished;
    
  algorithm
    tcl_old := tcl;
    finished := 0;
    while finished == 0 loop
    chr_ := 4 * sbc * (((tcl_old + tr) / 2 + 273.15) ^ 3) * .72;
    ctc_ := chr_ + chc;
    ra_  := 1 / (facl * ctc_);
    top_ := (chr_ * tr + chc * ta) / ctc_;
    tcl_ := (ra_ * tsk + rcl * top_) / (ra_ + rcl);
      if abs(tcl_ - tcl_old) > .01 then
        tcl_old := tcl_;
      else
        finished := 1;
      end if;
  end while;
  //print("tcl_ = " + String(tcl));
end tcl_calculate;

function fnp //make negative part of signals = 0
  input Real x;
  output Real y;
  algorithm
    if x < 0 then
      y := 0;
      else
      y := x;
    end if;
end fnp;

function fnsvp //find saturated vapour pressure
  input Real temperature;
  output Real svp;
  algorithm
  svp := exp(18.6686 - 4030.183 / (temperature + 235));
end fnsvp;

equation
  //constant
  ta = set_temp;
  //step
  //if time < 7200 then ta = 20; else ta = 15; end if; 
  //sine  
  //ta = -1*(2*Modelica.Math.cos((2*Modelica.Constants.pi*time/(3600*24))) + 15)+36;
  
  tr = ta; ///find a better formula for tr

  if clo <= 0 then
    wcrit = 0.38 * vel ^(-0.29);
    icl = 1;
  else
    wcrit = 0.59 * vel ^(-0.08);
    icl   = 0.45;
  end if;
  
  (tcl,chr,ctc,top,ra) = tcl_calculate(tcl,tr,chc,facl,ta,tsk,ctc,rcl);

  dry = (tsk - top) / (ra+rcl);
  hfcs = (tcr - tsk) * (5.28 + 1.163 * skbf);
  eres = 0.0023 * m * (44 - pa);
  cres = 0.0014 * m * (34 - ta);
  scr  = m - hfcs - eres - cres - w; //wm and w in w/2
  ssk  = hfcs - dry - esk;
  tcsk = 0.97 * alpha * wt;
  tccr = 0.97 * (1 - alpha) * wt;
  deltatsk = (ssk * sa) / tcsk / 3600; //°C / s
  deltatcr = (scr * sa) / tccr / 3600; //°C / s
  der(tsk) = deltatsk;
  der(tcr) = deltatcr;
  
  tb  = alpha * tsk + (1 - alpha) * tcr;

  chcv = 8.600001 * (vel * atm) ^ 0.53;
  chc = 3 * atm  ^ 0.53;
  pa = rh * exp(18.6686 - (4030.183/(ta + 235))) / 100;
 
  sksig = tsk - tskn;
  crsig = tcr - tcrn;
  bdsig = tb - tbn;
  warms = fnp(sksig);
  colds = fnp(-sksig);
  warmc = fnp(crsig);
  coldc = fnp(-crsig);
  warmb = fnp(bdsig);
  coldb = fnp(-bdsig);

  skbf_ = (skbfn + cdil * warmc) / (1 + cstr * colds); ///SKBF is in the wrong units as always PLS FIX
  if skbf_ > 90 then 
    skbf = 90;
  elseif skbf_ < 0.5 then 
    skbf = 0.5;
  else
    skbf = skbf_;
  end if;

  regsw = csw * warmb * exp(warms / 10.7); //MAKE SURE THE MAX WORKS

  rea = 1 / (lr * facl * chc);
  recl = rcl / (lr * icl);

  esk = esk_;
  emax = (exp(18.6686 - 4030.183 / (tsk + 235)) - pa) / (rea + recl);
  ersw_ = .68 * regsw;
  if (0.06+0.94 * ersw_ / emax) > wcrit then
    pwet = wcrit;
    prsw = wcrit / 0.94;
    ersw = prsw * emax;
    edif = 0.06 * (1-prsw)*emax;
    esk_ = ersw + edif;
  elseif emax <= 0 then
    pwet = wcrit;
    prsw = wcrit;
    ersw = 0;
    edif = 0;
    esk_ = emax;
  else
    pwet = 0.06 + 0.94 * prsw;
    prsw = ersw / emax;
    ersw = 0.68 * regsw;
    edif = pwet * emax - ersw;
    esk_ = ersw + edif;
  end if;
  
  mshiv = 19.4 * colds * coldc;
  m = rmm + mshiv;
  alpha =  0.0417737 + 0.7451833 / (skbf + 0.585417);

end gagge_2_node;
