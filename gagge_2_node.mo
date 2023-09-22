model gagge_2_node
  parameter Real ta = 25 "air temperature, °C"; 
  parameter Real tr = 25 "mean radiant temperature, °C";
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
  
  Real pa (start = rh * exp(18.6686 - (4030.183/(ta + 235))) / 100) "?";
  Real tsk (start = tskn) "skin temperature, °C";
  Real tcr (start = tcrn) "core temperature, °C";
  Real tcl (start = tcl_estimate(1)) "clothing temperature, °C";
  Real skbf (start = skbfn) "skin blood flow, kg/hr m^2";
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
  Real top (start = ((4.7 * tr) + (3 * atm ^ (0.53)) * ta) /  4.7 + (3 * atm ^ (0.53))) "operative temperature, °C";


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
  input Real ta;
  input Real tsk;
  input Real ctc;
  input Real ra;
  output Real tcl_guess;
    protected
    Real chr_;
    Real ctc_;
    Real tcl_old;
    Real top_;
    Real ra_;
    Integer finished;
  algorithm
    tcl_old := tcl;
    finished := 0;
    while finished == 0 loop
    chr_ := 4 * sbc * (((tcl_old + tr) / 2 + 273.15) ^ 3) * .72;
    ctc_ := chr_ + chc;
    ra_  := 1 / (facl * ctc_);
    top_ := (chr_ * tr + chc * ta) / ctc_;
    tcl_guess := (ra * tsk + rcl * top_) / (ra + rcl);
      if abs(tcl_guess - tcl_old) > .01 then
        tcl_old := tcl_guess;
      else
        finished := 1;
      end if;
  end while;
end tcl_calculate;

function fnp
  input Real x;
  output Real result;
  algorithm
    if value < 0 then
      result := 0;
      else
      result := value;
    end if;
end fnp;

function fnsvp
  input Real temperature;
  output Real svp;
  algorithm
  svp := {exp(18.6686 - 4030.183 / (temperature + 235))};
end fnsvp;



equation

  if clo <= 0 then
    wcrit = 0.38 * vel ^(-0.29);
    icl = 1;
  else
    wcrit = 0.59 * vel ^(-0.08);
    icl   = 0.45;
  end if;
  
  tcl = tcl_calculate(tcl,tr,chc,ta,tsk,ctc,ra);

  tsk = 1;
  tcr = time*2;
  skbf = 3;
  mshiv = 4;
  alpha = 4;
  rmm = 5;
  esk = 5;
  chcv = 6;
  chc = 6;
  pa = 6;
  chr = 7;
  ctc = chr + chc;
  ra = 1/facl*ctc;
  top = 8;
  
end gagge_2_node;
