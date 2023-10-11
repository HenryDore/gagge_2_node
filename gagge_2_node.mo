model gagge_2_node
  // 86400 = 1 day in seconds
  // 3600 = 1 hour in seconds
  import Modelica.Units.SI;
  import Modelica.Units.NonSI;
  import Modelica.Constants;

  constant Real clo = 0.5 "clothing insulation level, clo";
  constant Real met = 1 "metabolic rate, met";
  constant Real wme = 0 "external work, met";
  constant Real pb = 760 "barometric pressure, torr or mmHg";
  constant Real ht = 170 "height, cm";
  constant SI.Mass wt = 70 "weight, kg";
  constant NonSI.Temperature_degC tskn = 33.7 "setpoint value for skin temperature";
  constant NonSI.Temperature_degC tcrn = 36.8 "setpoint value for core temperature";
  constant NonSI.Temperature_degC tbn = 36.49 "setpoint value for mean body temperature (.1*tskn + .9*tcrn)";
  constant Real skbfn = 6.3 "neutral value for skin blood flow";
  constant Real sbc = 5.6697*10^(-8) "stephan-Boltzmann constant";
  constant SI.Area sa = 0.203*(ht/100)^(0.725)*wt^(0.425) "Dubois surface Area";
  constant Real w = wme*58.2 "external work, w/m^2";
  constant Real atm = pb/760 "atmospheric pressure, atm";
  constant Real rcl = 0.155*clo "thermal resistance of clothing ensemble, °C m^2/W";
  constant Real facl = 1 + 0.15*clo "Increase in body surface area due to clothing";
  constant Real tu = 40 "turbulence intensity, %";
  constant Real vel = 0.1 "air velocity m/s";
  constant Real rh = 50 "relative humidity %";
  constant Real csw = 170 "driving coefficient for regulatory sweating, 1";
  constant Real cdil = 120 "driving coefficient for vasodilation, 1";
  constant Real cstr = 0.5 "driving coefficient for vasoconstriction, 1";
  constant Real lr = 2.2/atm "Lewis Relation is 2.2 at sea level";
  constant Real rmm = met*58.2 "basic metabolic rate, W";
  constant Real chc = 3*atm^(0.53) "coductive heat transfer coefficient";
  constant Real chcv = 8.600001*(vel*atm)^0.53 "convective heat transfer coefficient";
  constant Real rea = 1/(lr*facl*chc) "resistance of air layer to dry heat transfer, °C m^2 / W";
  parameter Real set_temp = 22.5 "setpoint for ambient temperature";
  SI.HeatFlux m(start = met*58.2) "metabolic rate";
  NonSI.Temperature_degC ta(start = 20) "air temperature";
  NonSI.Temperature_degC tr(start = 20) "mean radiant temperature";
  Real pa "partial vapour pressure of water, mmHg";
  NonSI.Temperature_degC tsk(start = tskn, fixed = true) "skin temperature";
  NonSI.Temperature_degC tcr(start = tcrn, fixed = true) "core temperature";
  NonSI.Temperature_degC tcl(start = tskn) "clothing temperature, °C";
  NonSI.Temperature_degC tb(start = 0.1*tskn + (1 - 0.1)*tcrn) "mean body temperature";
  Real skbf(start = skbfn) "skin blood flow, kg/hr m^2";
  Real skbf_ "skin blood flow test";
  SI.HeatFlux mshiv(start = 0) "rate of energy released by shivering";
  Real alpha(start = 0.1) "fractional skin mass, 1";
  SI.HeatFlux esk(start = 0.1*met) "total evaporative heat loss from the skin";
  Real wcrit "evaporative efficiency, 1";
  Real icl "";
  Real chr(start = 4.7) "radiative heat transfer coefficient";
  Real ctc(start = 4.7 + (3*atm^(0.53))) "?combined convection & radiation coefficient";
  Real ra(start = 1/(facl*4.7 + (3*atm^(0.53)))) "resistance of air layer to dry heat transfer, °C m^2 / W";
  NonSI.Temperature_degC top(start = ((4.7*20) + (3*atm^(0.53))*20)/4.7 + (3*atm^(0.53))) "operative temperature";
  SI.HeatFlux dry "total sensible heat loss";
  SI.HeatFlux hfcs "rate of energy transport between core and skin";
  SI.HeatFlux eres "heat loss through respiratory evaporation";
  SI.HeatFlux cres "heat loss through respiratory convection";
  SI.HeatFlux scr "rate of energy storage in the core";
  SI.HeatFlux ssk "rate of enery storage in the skin, ";
  Real warms "skin warm signal ,1";
  Real colds "skin cold signal ,1";
  Real warmc "core warm signal ,1";
  Real coldc "core cold signal ,1";
  Real warmb "blood warm signal ,1";
  Real coldb "blood cold signal ,1";
  Real regsw "regulatory sweating, g/m^2 hr";
  SI.HeatFlux ersw "heat loss through sweating";
  SI.HeatFlux ersw_ "also heat loss through sweating";
  Real recl "resistance of clothing layer to dry heat transfer, °C m^2 / W";
  SI.HeatFlux emax "maximum evapourative capacity";
  Real pwet "skin wettedness, 1";
  Real prsw "ratio of actual heat loss due to sweating to maximum heat loss due to sweating";
  SI.HeatFlux edif "heat loss due to diffusion of water fvapout rhough the skin";
  SI.HeatFlux esk_ "also total evaporative heat loss from the skin";

  function tcl_calculate
    //calculate tcl, chr, ctc, top & ra iteratively
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
      chr_ := 4*sbc*(((tcl_old + tr)/2 + 273.15)^3)*.72;
//ASHRAE eqn 35, can add A_r/A_d to improve
      ctc_ := chr_ + chc;
      ra_ := 1/(facl*ctc_);
      top_ := (chr_*tr + chc*ta)/ctc_;
//ASHRAE eqn 8
      tcl_ := (ra_*tsk + rcl*top_)/(ra_ + rcl);
      if abs(tcl_ - tcl_old) > .01 then
        tcl_old := tcl_;
      else
        finished := 1;
      end if;
    end while;
//print("tcl_ = " + String(tcl)); //debug
  end tcl_calculate;

  function fnp
    //make negative part of signals = 0
    input Real x;
    output Real y;
  algorithm
    if x < 0 then
      y := 0;
    else
      y := x;
    end if;
  end fnp;

  //connectors
  Modelica.Blocks.Interfaces.RealInput ambient_temperature annotation(
    Placement(visible = true, transformation(origin = {-120, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput metabolic_rate annotation(
    Placement(visible = true, transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
//constant
//ta = set_temp;
//step
//if time < 7200 then ta = 22.5; else ta = 15; end if;
//sine
//ta = -1*(2*Modelica.Math.cos((2*Modelica.Constants.pi*time/(3600*24))) + 15)+40;
//assume tr = ta. This is an assumption but experimental data from
//DOI 10.1007/s00484-010-0375-4 shows that tr is within 1% of ta on average.
//connector - input
  ta = ambient_temperature;
  tr = ta;
//this is an assumption but it is very close. Experimental data from
  if clo <= 0 then
    wcrit = 0.38*vel^(-0.29);
    icl = 1;
  else
    wcrit = 0.59*vel^(-0.08);
    icl = 0.45;
  end if;
//resistance of clothing layer to dry heat transfer
  recl = rcl/(lr*icl);
//calculate tcl (∴ chr,ctc,top&ra) iteratively
  (tcl, chr, ctc, top, ra) = tcl_calculate(tcl, tr, chc, facl, ta, tsk, ctc, rcl);
//partial vapour pressure of water for ta
  pa = rh*exp(18.6686 - (4030.183/(ta + 235)))/100;
//heat balance
  dry = (tsk - top)/(ra + rcl);
  hfcs = (tcr - tsk)*(5.28 + 1.163*skbf);
  eres = 0.0023*m*(44 - pa);
  cres = 0.0014*m*(34 - ta);
  scr = m - hfcs - eres - cres - w;
//wm and w in w/2
  ssk = hfcs - dry - esk;
  der(tsk) = (ssk*sa)/(0.97*alpha*wt)/3600;
//°C / s
  der(tcr) = (scr*sa)/(0.97*(1 - alpha)*wt)/3600;
//mean body temperature
  tb = alpha*tsk + (1 - alpha)*tcr;
//thermoregulatory signals
  warms = fnp(tsk - tskn);
//sksig
  colds = fnp(-(tsk - tskn));
//-sksig
  warmc = fnp(tcr - tcrn);
//crsig
  coldc = fnp(-(tcr - tcrn));
//-crsig
  warmb = fnp(tb - tbn);
//bdsig
  coldb = fnp(-(tb - tbn));
//-bdsig
//skin blood flow
  skbf_ = (skbfn + cdil*warmc)/(1 + cstr*colds);
///SKBF is in the wrong units as always PLS FIX
  if skbf_ > 90 then
    skbf = 90;
  elseif skbf_ < 0.5 then
    skbf = 0.5;
  else
    skbf = skbf_;
  end if;
//regulatory sweating
  regsw = min(500, (csw*warmb*exp(warms/10.7)));
//is CSW correct?
//skin wettedness conditions
  emax = (exp(18.6686 - 4030.183/(tsk + 235)) - pa)/(rea + recl);
  esk = esk_;
  ersw_ = 0.68*regsw;
  if (0.06 + 0.94*ersw_/emax) > wcrit then
//error
    pwet = wcrit;
    prsw = wcrit/0.94;
    ersw = prsw*emax;
    edif = 0.06*(1 - prsw)*emax;
    esk_ = ersw + edif;
  elseif emax <= 0 then
//error skin wettedness too high
    pwet = wcrit;
    prsw = wcrit;
    ersw = 0;
    edif = 0;
    esk_ = emax;
  else
//regulatory sweating
    pwet = 0.06 + 0.94*prsw;
    prsw = ersw/emax;
    ersw = 0.68*regsw;
    edif = pwet*emax - ersw;
    esk_ = ersw + edif;
  end if;
//shivering
  mshiv = 19.4*colds*coldc;
  m = rmm + mshiv;
//update alpaha from skbf
  alpha = 0.0417737 + 0.7451833/(skbf + 0.585417);
//output connectors
  metabolic_rate = m;
  annotation(
    uses(Modelica(version = "4.0.0")),
    Icon(graphics = {Ellipse(origin = {0, -60}, fillColor = {100, 100, 100}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-60, -20}, {60, 20}}), Rectangle(lineColor = {100, 100, 100}, fillColor = {100, 100, 100}, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-60, 60}, {60, -60}}), Line(origin = {-60, 0}, points = {{0, -60}, {0, 60}}, thickness = 0.5), Line(origin = {60, 0}, points = {{0, 60}, {0, -60}}, thickness = 0.5), Ellipse(origin = {0, 60}, fillColor = {100, 100, 100}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-60, 20}, {60, -20}}), Ellipse(origin = {-1, 60}, fillColor = {150, 150, 150}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-49, 14}, {49, -14}}), Text(origin = {0, -1}, extent = {{-40, -33}, {40, 33}}, textString = "Gagge
2-node")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
    experiment(StartTime = 0, StopTime = 86400, Tolerance = 1e-06, Interval = 60));
end gagge_2_node;
