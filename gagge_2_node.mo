model gagge_2_node
  import Modelica.Units.SI;
  import Modelica.Units.NonSI;
  import Modelica.Constants;
  //constants
  constant Real sbc = 5.6697*10^(-8) "stephan-Boltzmann constant";
  constant NonSI.Temperature_degC tskn = 33.7 "setpoint value for skin temperature";
  constant NonSI.Temperature_degC tcrn = 36.8 "setpoint value for core temperature";
  constant NonSI.Temperature_degC tbn = 36.49 "setpoint value for mean body temperature (.1*tskn + .9*tcrn)";
  constant Real skbfn = 6.3 "neutral value for skin blood flow";
  parameter Real age = 40;
  parameter Integer male = 1;
  parameter Real clo_in = 10 "clo_in = clo value * 100 (clo is multiplied by 10 on input as hack for DOS batch file)";
  Real clo "clothing insulation level, clo";
  Real rcl "thermal resistance of clothing ensemble, °C m^2/W";
  Real facl "Increase in body surface area due to clothing";
  Real sim_time = time;
  constant Real csw = 170 "driving coefficient for regulatory sweating, 1";
  constant Real cdil = 120 "driving coefficient for vasodilation, 1";
  constant Real cstr = 0.5 "driving coefficient for vasoconstriction, 1";
  parameter Real chc = 3*atm^(0.53) "conductive heat transfer coefficient";
  parameter Real chcv = 8.600001*(vel*atm)^0.53 "convective heat transfer coefficient";
  Real rea "resistance of air layer to dry heat transfer, °C m^2 / W";
  //environmental parameters
  parameter Real set_temp = 22.5 "setpoint for ambient temperature";
  constant Real tu = 40 "turbulence intensity, %";
  constant Real vel = 0.1 "air velocity m/s";
  parameter Real rh = 50 "relative humidity %";
  constant Real pb = 760 "barometric pressure, mmHg";
  constant Real atm = pb/760 "atmospheric pressure, atm";
  constant Real lr = 2.2/atm "Lewis Relation is 2.2 at sea level";
  Real pa "partial vapour pressure of water, mmHg";
  Real pa_pa "partial vapour pressure of water, pa";
  //human parameters
  parameter Real ht = 170 "height, cm";
  parameter SI.Mass wt = 70 "weight, kg";
  parameter SI.Area sa = 0.203*(ht/100)^(0.725)*wt^(0.425) "Dubois surface Area";
  parameter Real met = 1 "metabolic rate, met";
  parameter Real wme = 0 "external work, met";
  Real BMR "basal metabolic rate";
  Real w "external work, W";
  Real rmm "basal metabolic rate, W";
  //record Summary
  NonSI.Temperature_degC top(start = ((4.7*20) + (3*atm^(0.53))*20)/4.7 + (3*atm^(0.53))) "operative temperature";
  NonSI.Temperature_degC ta(start = 20) "air temperature";
  NonSI.Temperature_degC tr(start = 20) "mean radiant temperature";
  NonSI.Temperature_degC tsk(start = tskn, fixed = true) "skin temperature";
  NonSI.Temperature_degC tcr(start = tcrn, fixed = true) "core temperature";
  NonSI.Temperature_degC tcl(start = tskn) "clothing temperature, °C";
  NonSI.Temperature_degC tb(start = 0.1*tskn + (1 - 0.1)*tcrn) "mean body temperature";
  //end Summary;
  Real skbf(start = skbfn) "skin blood flow, kg/hr m^2";
  Real skbf_ "skin blood flow test";
  Real alpha(start = 0.1) "fractional skin mass, 1";
  Real wcrit "evaporative efficiency, 1";
  Real icl "";
  Real chr(start = 4.7) "radiative heat transfer coefficient";
  Real ctc(start = 4.7 + (3*atm^(0.53))) "?combined convection & radiation coefficient";
  Real ra(start = 1/(4.7 + (3*atm^(0.53)))) "resistance of air layer to dry heat transfer, °C m^2 / W";
  //calculated values
  SI.HeatFlux m(start = met*58.2) "metabolic rate";
  SI.HeatFlux dry "total sensible heat loss";
  SI.HeatFlux hfcs "rate of energy transport between core and skin";
  SI.HeatFlux eres "heat loss through respiratory evaporation";
  SI.HeatFlux cres "heat loss through respiratory convection";
  SI.HeatFlux scr "rate of energy storage in the core";
  SI.HeatFlux ssk "rate of enery storage in the skin, ";
  SI.HeatFlux ersw "heat loss through sweating";
  SI.HeatFlux ersw_ "also heat loss through sweating";
  SI.HeatFlux emax "maximum evapourative capacity";
  SI.HeatFlux edif "heat loss due to diffusion of water fvapout rhough the skin";
  SI.HeatFlux esk_ "also total evaporative heat loss from the skin";
  SI.HeatFlux mshiv(start = 0) "rate of energy released by shivering";
  SI.HeatFlux esk(start = 0.1*met) "total evaporative heat loss from the skin";
  //signals
  Real warms "skin warm signal ,1";
  Real colds "skin cold signal ,1";
  Real warmc "core warm signal ,1";
  Real coldc "core cold signal ,1";
  Real warmb "blood warm signal ,1";
  Real coldb "blood cold signal ,1";
  Real regsw "regulatory sweating, g/m^2 hr";
  Real pwet "skin wettedness, 1";
  Real prsw "ratio of actual heat loss due to sweating to maximum heat loss due to sweating";
  Real recl "resistance of clothing layer to dry heat transfer, °C m^2 / W";
  //comfort indices
  Real PMV_gagge "";
  Real PMV "Predicted mean vote";
  Real PMV_ISO7730;
  Real PMV_t_cl "Clothing surface temperature";
  Real PMV_h_c "convective heat transfer coefficient";
  Real PMV_f_cl "clothing surfacea area factor";
  Real t_b_c "cold setpoint for calculation of TSENS";
  Real t_b_h "hot setpoint for calculation of TSENS";
  Real TSENS "Thermal sensation";

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
    Placement(transformation(origin = {-120, 0}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-80, 40}, extent = {{-20, -20}, {20, 20}})));

  record summary
    NonSI.Temperature_degC t_ambient "air temperature";
    NonSI.Temperature_degC t_skin "skin temperature";
    NonSI.Temperature_degC t_core "core temperature";
    NonSI.Temperature_degC t_clothing "clothing temperature, °C";
    NonSI.Temperature_degC t_mean_body "mean body temperature";
    SI.HeatFlux q_shivering "rate of energy released by shivering";
    SI.HeatFlux q_sweating "heat loss through sweating";
    SI.HeatFlux q_metabolism "heat generated by metabolism";
    SI.HeatFlux q_sensible "heat loss through convection & radiation";
    SI.HeatFlux q_evaporation "heat loss through evaporation";
    SI.HeatFlux q_c_res "heat loss through convection in breathing";
    SI.HeatFlux q_e_res "heat loss through evaporation in breathing";
    Real PMV "Predicted mean vote";
    Real PPD "Percentage of people dissatisfied";
    Real TSENS "Thermal sensetivity";
  end summary;

  summary RPT1;
  Modelica.Blocks.Interfaces.RealOutput human_results annotation(
    Placement(transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {80, 40}, extent = {{-20, -20}, {20, 20}})));
  Modelica.Blocks.Interfaces.RealOutput PMV_out annotation(
    Placement(transformation(origin = {110, 50}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {80, 80}, extent = {{-20, -20}, {20, 20}})));


  //CALCULATE BMR w and rmm
  //Mifflin-St Jeor equation:
  //BMR (kcal/day) = 10 × weight (kg) + 6.25 × height (cm) – 5 × age (y) + s (kcal/day),
  //where s is +5 for males and -161 for females.
algorithm
  if male == 1 then
    BMR := (10*wt + 6.25*ht - 5*age - 5)*0.048;
  elseif male == 0 then
    BMR := (10*wt + 6.25*ht - 5*age - 161)*0.048;
  else
    BMR := 58.2;
  end if;
  w := wme*58.2;
  rmm := met*58.2;
  
equation
  ta = ambient_temperature - 273.15;
//convert from °K
  tr = ta;
//clo derived variables:
  clo = clo_in/10;
  rcl = 0.155*clo;
  facl = 1 + 0.15*clo;
  rea = 1/(lr*facl*chc);
//================================//
// THERMAL SENSATION CALCULATIONS //
//================================//
// PMV
  PMV_t_cl = 35.7 - 0.028*(m - w) - clo*(3.96*10^(-8)*PMV_f_cl*((tcl + 273)^4 - (tr + 273)^4) + PMV_f_cl*PMV_h_c*(tcl - ta));
  if 2.38*abs(tcl - ta)^0.25 > 12.1*sqrt(vel) then
    PMV_h_c = 2.38*abs(tcl - ta)^0.25;
  else
    PMV_h_c = 12.1*sqrt(vel);
  end if;
  if clo < 0.078 then
    PMV_f_cl = 1.00 + 1.290*clo;
  else
    PMV_f_cl = 1.05*0.645*clo;
  end if;
  PMV_ISO7730 = (0.303*exp(-0.036*m) + 0.028)*((m - w) - 3.05*10^(-3)*(5733 - 6.99*(m - w) - pa_pa) - 0.42*((m - w) - 58.15) - 1.7*10^(-5)*m*(5867 - pa_pa) - 0.0014*m*(34 - ta) - 3.96*10^(-8)*PMV_f_cl*((tcl + 273)^4 - (tr + 273)^4) - PMV_f_cl*PMV_h_c*(tcl - ta));
//PMVG as per derivation
  PMV_gagge = (0.303*exp(-0.036*rmm) + .028)*(rmm - eres - cres - dry - regsw + mshiv);
//gagge PMV
  PMV = (0.303*exp(-2.1*rmm) + 0.028)*((rmm - w) - dry - esk - eres - cres + regsw + mshiv);
  PMV_out = PMV;
// TSENS
  t_b_c = (0.194/58.15)*(rmm - w) + 36.301;
  t_b_h = (0.347/58.15)*(rmm - w) + 36.669;
  if tb < t_b_c then
    TSENS = 0.4685*(tb - t_b_c);
  elseif (tb >= t_b_c) and (tb < t_b_h) then
    TSENS = wcrit*4.7*(tb - t_b_c)/(t_b_h - t_b_c);
  elseif tb >= t_b_h then
    TSENS = wcrit*4.7 + 0.4685*(tb - t_b_h);
  else
    TSENS = 0;
  end if;
  //skin temperature output connector
  human_results = tsk;
//===========//
// REPORTING //
//===========//
  RPT1.t_ambient = ta;
  RPT1.t_skin = tsk;
  RPT1.t_core = tcr;
  RPT1.t_clothing = tcl;
  RPT1.t_mean_body = tb;
  RPT1.q_shivering = mshiv;
  RPT1.q_sweating = ersw;
  RPT1.q_metabolism = m;
  RPT1.q_sensible = dry;
  RPT1.q_evaporation = esk;
  RPT1.q_c_res = cres;
  RPT1.q_e_res = eres;
  RPT1.PMV = PMV;
  RPT1.PPD = 100 - 95*exp(-0.03353*PMV^(4) - 0.2179*PMV^(2));
  RPT1.TSENS = TSENS;

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
  pa_pa = pa*133.322;
//heat balance
  dry = (tsk - top)/(ra + rcl);
  hfcs = (tcr - tsk)*(5.28 + 1.163*skbf);
  eres = 0.0023*m*(44 - pa);
  cres = 0.0014*m*(34 - ta);
  scr = m - hfcs - eres - cres - w;
//wm and w in w/2
//ssk = hfcs - dry - esk;
//added rad_heat from heater into body model
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
  if skbf_ > 90 then
    skbf = 90;
  elseif skbf_ < 0.5 then
    skbf = 0.5;
  else
    skbf = skbf_;
  end if;
//regulatory sweating
  regsw = min(500, (csw*warmb*exp(warms/10.7)));
//skin wettedness conditions
  emax = (exp(18.6686 - 4030.183/(tsk + 235)) - pa)/(rea + recl);
  esk = esk_;
  ersw_ = 0.68*regsw;
  if noEvent(0.06 + 0.94*ersw_/emax >= wcrit) then
//error
    pwet = wcrit;
    prsw = wcrit/0.94;
    ersw = prsw*emax;
    edif = 0.06*(1 - prsw)*emax;
    esk_ = ersw + edif;
  elseif noEvent(emax <= 0) then
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
  annotation(
    uses(Modelica(version = "4.0.0")),
    Icon(graphics = {Ellipse(origin = {0, -60}, fillColor = {255, 170, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-60, -20}, {60, 20}}), Rectangle(lineColor = {100, 100, 100}, fillColor = {255, 170, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-60, 60}, {60, -60}}), Line(origin = {-60, 0}, points = {{0, -60}, {0, 60}}, thickness = 0.5), Line(origin = {60, 0}, points = {{0, 60}, {0, -60}}, thickness = 0.5), Ellipse(origin = {0, 60}, fillColor = {255, 0, 127}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-60, 20}, {60, -20}}), Ellipse(origin = {-1, 60}, fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-49, 14}, {49, -14}}), Text(origin = {0, -1}, extent = {{-40, -33}, {40, 33}}, textString = "Gagge
2-Node", fontName = "Consolas", textStyle = {TextStyle.Bold}), Text(origin = {0, 120}, textColor = {0, 0, 255}, extent = {{-100, 20}, {100, -20}}, textString = "%name")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
    experiment(StartTime = 0, StopTime = 86400, Tolerance = 1e-06, Interval = 60),
    Documentation(info = "<html><head></head><body><p><b></b> Henry Dore (hd255@sussex.ac.uk), 2024</p><div>The 2-node model describes heat transfer inside the components of the body and between the body and its environment, assuming no heat storage within the body. This model has been enhanced and expanded in various ways since its inception, such as the addition of multiple segments, multiple nodes, or an improved sweating model. The 2-node model calculates body temperatures by solving the heat balance of the human body using the following equation:</div><div><br></div><div>Q_met+Q_shiv=Q_diff+Q_sw+Q_rad+Q_conv+Q_res</div><div><br></div><div>Shivering heat generation<span class=\"Apple-tab-span\" style=\"white-space:pre\">	</span>Q_shiv=19.4 A (34-T_sk ) &nbsp;(36.6-T_cr )</div><div><br></div><div>Metabolic heat generation<span class=\"Apple-tab-span\" style=\"white-space:pre\">	</span>Q_met=1.2 met</div><div><br></div><div>Diffusion of vapour heat loss<span class=\"Apple-tab-span\" style=\"white-space:pre\">	</span>Q_diff=A fcl LR H_c &nbsp;w (T_cl-T_a)</div><div><br></div><div>Sweating heat loss<span class=\"Apple-tab-span\" style=\"white-space:pre\">	</span>Q_sw=A λ_sweat &nbsp;Sw Lcf</div><div><br></div><div>Radiative heat loss<span class=\"Apple-tab-span\" style=\"white-space:pre\">	</span>Q_rad=A fcl LR H_c &nbsp;w (T_cl-T_a)</div><div><br></div><div>Convective heat loss<span class=\"Apple-tab-span\" style=\"white-space:pre\">	</span>Q_conv=A fcl H_c &nbsp;(T_cl-T_a)</div><div><br></div><div>Respiration heat loss<span class=\"Apple-tab-span\" style=\"white-space:pre\">	</span>Q_res=0.023 Q_(met ) (44-P_a )+0.0014 (34-T_a)</div><div><br><p>
<img src=\"modelica://gagge_2_node/resources/doc_gagge_2_node.png\" width=\"75%\">
</p><p><br></p></div></body></html>"));
end gagge_2_node;
