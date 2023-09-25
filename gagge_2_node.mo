model gagge_2_node
   // 86400 = 1 day in seconds
  parameter Real vel = 0.1 "air velocity m/s";
  parameter Real rh = 50 "realtive humidity %";
  constant Real clo = 0.5 "clothing insulation level, clo";
  constant Real met = 1 "metabolic rate, met";
  constant Real wme = 0 "external work, met";
  constant Real pb = 760 "barometric pressure, torr or mmHg";
  parameter Real ltime = 60 "exposure time, minutes";
  constant Real ht = 170 "height, cm";
  constant Real wt = 70 "weight, kg";
  parameter Real tu = 40 "turbulence intensity, %";
  parameter Real csw = 170 "driving coefficient for regulatory sweating";
  parameter Real cdil = 120 "driving coefficient for vasodilation";
  parameter Real cstr = 0.5 "driving coefficient for vasoconstriction";  
  constant Real tskn = 33.7 "setpoint (neutral) value for tsk";
  constant Real tcrn = 36.8 "setpoint value for tcr";
  constant Real tbn = 36.8 "setpoint value for tb (.1*tskn + .9*tcrn)";
  constant Real skbfn = 6.3 "neutral value for skbf";
  constant Real sbc = 5.6697 * 10 ^ (-08) "stephan-Boltzmann constant";
  constant Real sa = ((ht * wt) / 3600 ) ^ .5 "surface Area (m2) according to mosteller formula";
  constant Real m = met * 58.2 "metabolic rate, w/m^2";
  constant Real w = wme * 58.2 "external work, w/m^2";
  parameter Real mw = m-w "metabolism - work, w/m^2";
  constant Real atm = pb / 760 "atmospheric pressure, atm";
  parameter Real timeh = ltime / 60 "run time in hours";
  constant Real rcl = 0.155 * clo "thermal resistance of clothing ensemble, °C m^2/W";
  constant Real facl = 1 + 0.15 * clo "Increase in body surface area due to clothing";
  parameter Real lr = 2.2 / atm "Lewis Relation is 2.2 at sea level";
  
  Real ta (start = 20) "air temperature, °C"; 
  Real tr (start = 20) "mean radiant temperature, °C";
  Real pa "?";
  Real tsk (start = tskn) "skin temperature, °C";
  Real tcr (start = tcrn) "core temperature, °C";
  Real tcl (start = tcl_estimate(1)) "clothing temperature, °C";
  Real tb  (start = 0.1 * tskn + (1 - 0.1) * tcrn);
  Real skbf (start = skbfn, max = 90, min = 0.5) "skin blood flow, kg/hr m^2";
  Real skbf_ "skin blood flow test";
  Real mshiv (start = 0) "rate of energy released by shivering, W";
  Real alpha (start = 0.1) "fractional skin mass, ND";
  Real rmm (start = m ) "metabolic rate";
  Real esk (start = 0.1 * met) "total evaporative heat loss from the skin W/m^2";
  Real wcrit;
  Real icl;
  Real chc (start = 3 * atm ^ (0.53)) "?coductive heat transfer coefficient";
  Real chcv (start = 8.600001 * (vel * atm) ^ 0.53) "?convective heat transfer coefficient";
  Real chr (start = 4.7) "?radiative heat transfer coefficient";
  Real ctc (start = 4.7 + (3 * atm ^ (0.53))) "?combined convection & radiation coefficient";
  Real ra (start = 1 / (facl * 4.7 + (3 * atm ^ (0.53)))) "resistance of air layer to dry heat transfer";
  Real top (start = ((4.7 * 20) + (3 * atm ^ (0.53)) * 20) /  4.7 + (3 * atm ^ (0.53))) "operative temperature, °C";
  //Real top (start = 29) "operative temperature, °C";
  Real dry;
  Real hfcs;
  Real eres;
  Real cres;
  Real scr;
  Real ssk;
  Real tcsk;
  Real tccr;
  Real deltatsk (start = 0);
  Real deltatcr (start = 0);
  Real sksig;
  Real warms;
  Real colds;
  Real crsig;
  Real warmc;
  Real coldc;
  Real bdsig;
  Real warmb;
  Real coldb;
  Real regsw (max = 500);
  Real ersw;
  Real ersw_;
 // Real st;
  Real rea;
  Real recl;
  Real emax;
  Real pwet;
  Real prsw;
  Real edif;
  Real esk_;
 // Real ssk_;
  //Real scr_;
  
function tcl_estimate
  input Real temp; 
  output Real tcl;
  algorithm
  tcl := ((4.7 * 25) + (3 * atm ^ (0.53)) * 25) /  4.7 + (3 * atm ^ (0.53))
  + (tskn - ((4.7 * 25) + (3 * atm ^ (0.53)) * 25) /  4.7 + (3 * atm ^ (0.53)))
  / ((4.7 + (3 * atm ^ (0.53))) * ((1 / (facl * 4.7 + (3 * atm ^ (0.53)))) + rcl));
end tcl_estimate;


function tcl_calculate
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

function fnp
  input Real x;
  output Real result;
  algorithm
    if x < 0 then
      result := 0;
      else
      result := x;
    end if;
end fnp;

function fnsvp
  input Real temperature;
  output Real svp;
  algorithm
  svp := exp(18.6686 - 4030.183 / (temperature + 235));
end fnsvp;



equation
  //constant
  ta = 25;
  //step
  //time < 7200 then ta = 20; else ta = 10; end if; 
  //sine  
  //ta = -1*(2*Modelica.Math.cos((2*Modelica.Constants.pi*time/(3600*24))) + 15)+36;
  
  tr = ta;

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


  //tsk = tsk+deltatsk;
  //tcr = tcr+deltatcr;
  der(tsk) = deltatsk;
  der(tcr) = deltatcr;
  
  tb  = alpha * tsk + (1 - alpha) * tcr;

  esk = esk_;
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
    skbf = skbf_/3600;
  end if;

  regsw = csw * warmb * exp(warms / 10.7); //MAKE SURE THE MAX WORKS

  ersw_ = .68 * regsw;

  rea = 1 / (lr * facl * chc);
  recl = rcl / (lr * icl);

  emax = (exp(18.6686 - 4030.183 / (tsk + 235)) - pa) / (rea + recl);

  if (0.06+0.94 * ersw_ / emax) > wcrit then
    pwet = wcrit;
    prsw = wcrit / 0.94;
    ersw = prsw * emax;
    edif = 0.06 * (1-prsw)*emax;
    esk_ = ersw + edif;
  elseif emax == 0 then
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
