model PIR_gagge_smoosh
  //If using whole year weather data:
  //      * Change input weather file to contain year data
  //      * ICl (clothing thermal resistance) and phi (room relative humidity) set to sinusoidal function
  //      * Set simulation start time to 0, stop time to 31536000, number of intervals = 8760 (interval = 3600 s = 1h)
  //If using day weather data:
  //      *Change input weather file to contain day data
  //      * ICl (clothing thermal resistance) and phi (room relative humidity) set to specific values for that day. ICl varies between 0.5 in summer to 1.0 in winter. Phi is the same number as outside value.
  //      * Set simulation start time to 0, stop time to 86400, number of intervals = 144 (interval = 600 s = 10 min)
  import Modelica.Units.SI;
  import Modelica.Units.NonSI;
  import Modelica.Constants;
  //______________________________________________________________________________________________________________________________________
  Real mod_day(unit = "s");
  Real mod_week(unit = "s") "week starting on Monday";
  Real mod_winter(unit = "s");
  //______________________________________________________________________________________________________________________________________
  //Outdoor
  Modelica.Blocks.Sources.Constant pAir(k = 101325);
  //CHANGE FOR Judas IF AT HOME, hdore IF AT WORK!
  Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(calTSky = Buildings.BoundaryConditions.Types.SkyTemperatureCalculation.HorizontalRadiation, computeWetBulbTemperature = true, filNam = Modelica.Utilities.Files.loadResource("C:/Users/hdore/Box/Human Thermal Comfort/Code/WeatherFile_Gatwick-year.txt"), relHumSou = Buildings.BoundaryConditions.Types.DataSource.File) "Weather data reader" annotation(
    Placement(visible = true, transformation(origin = {6, 22}, extent = {{-98, 52}, {-78, 72}}, rotation = 0)));
  //filNam = Modelica.Utilities.Files.loadResource("C:/Program Files/OpenModelica1.21.0-64bit/Buildings 9.1.0/Resources/weatherdata/WeatherFile_Gatwick-1day_29April.txt")) "Weather data reader" annotation(Placement(transformation(extent = {{-98, 52}, {-78, 72}})));
  //______________________________________________________________________________________________________________________________________
  //Office (four elements: one floor, one roof, three walls and one wall with window)
  Real phi "relative humidity";
  parameter SI.Velocity vAir = 0.1;
  parameter Boolean use_vAir_in = true "Get the air velocity from the input connector" annotation(
    Evaluate = true,
    HideResult = true,
    Dialog(group = "Conditional inputs"));
  parameter Boolean use_M_in = true "Get the metabolic rate from the input connector" annotation(
    Evaluate = true,
    HideResult = true,
    Dialog(group = "Conditional inputs"));
  parameter Boolean use_ICl_in = true "Get the clothing insulation from the input connector" annotation(
    Evaluate = true,
    HideResult = true,
    Dialog(group = "Conditional inputs"));
  parameter Boolean use_pAir_in = true "Get the air pressure from the input connector" annotation(
    Evaluate = true,
    HideResult = true,
    Dialog(group = "Conditional inputs"));
  Modelica.Blocks.Interfaces.RealInput pAir_in(final quantity = "Pressure", final unit = "Pa", min = 0) if use_pAir_in "Air pressure";
  Modelica.Blocks.Interfaces.RealInput ICl_in if use_ICl_in "Clothing thermal resistance in clo";
  Modelica.Blocks.Interfaces.RealInput vAir_in if use_vAir_in "Air velocity";
  Modelica.Blocks.Interfaces.RealInput M_in(min = 40, max = 600, final quantity = "HeatFlux", final unit = "W/m2") if use_M_in "Metabolic heat generation in W/m2 (not in met)";
  Buildings.BoundaryConditions.SolarIrradiation.DiffusePerez HDifTil[2](each outSkyCon = true, each outGroCon = true, each til = 1.5707963267949, azi = {0, 1.5707963267949}) "Calculates diffuse solar radiation on titled surface for both directions, south-facing window" annotation(
    Placement(visible = true, transformation(origin = {-60, 26}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Buildings.BoundaryConditions.SolarIrradiation.DirectTiltedSurface HDirTil[2](each til = 1.5707963267949, azi = {0, 1.5707963267949}) "Calculates direct solar radiation on titled surface for both directions" annotation(
    Placement(visible = true, transformation(origin = {-60, 46}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Buildings.ThermalZones.ReducedOrder.SolarGain.CorrectionGDoublePane corGDouPan(n = 2, UWin = 2.7) "Correction factor for solar transmission. n is the vector size for input and output and UWin is thermal transmission coefficient of whole window (https://www.bipv.ch/index.php/en/technology-top-en/thermal-aspects/heat-transfer-coefficient#:~:text=Typical%20values%20vary%20between%206W,double%20glazed%20low%20emission%20window.)" annotation(
    Placement(visible = true, transformation(origin = {-32, 52}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Buildings.ThermalZones.ReducedOrder.RC.FourElements thermalZoneFourElements(VAir = 3*5*2.4 "air volume of the zone, 3 m by 4 m room, 2.4 m ceiling height", hConExt = 3 "convective coefficient of heat transfer of exterior walls (indoor), from Awbi (1998)", hConWin = 3 "convective coefficient of heat transfer of windows(indoor), from Wallentén", gWin = 0.78 "total energy transmittance of windows, from https://www.internorm.solutions/index.php/glossary-of-terms-g-value/", ratioWinConRad = 0.1 "ratio for windows between convective and radiative heat emission", nExt = 1 "number of RC elements of exterior walls (min=1)", RExt = {0.004} "https://www.engineeringtoolbox.com/heat-loss-transmission-d_748.html) vector of resistances of exterior walls, from inside to outside", CExt = {1360000} "vector of heat capacities of exterior walls, from inside to outside, calculated using density, thickness, and specific heat capacity of brick wall, from https://www.greenspec.co.uk/building-design/thermal-mass/", hRad = 5 "coefficient of heat transfer for linearised radiation exchange between walls", AInt = 2*(3*2.4 + 5*2.4) - 2*2 "area of interior walls, excluding floor, roof, and window", hConInt = 3 "convective coefficient of heat transfer of interior walls", nInt = 1 "number of RC elements of interior walls (min=1)", RInt = {0.000668895639141} "unchecked value - vector of resistances of interior walls from port to centre", CInt = {12391363.86} "unchecked value - vector of heat capacities of interior walls, from port to centre", RWin = 0.01642857143 "resistor for windows", RExtRem = 0.1265217391 "unchecked value - resistance of remaining resistor RExtRem between capacity n and outside", AFloor = 15 "area of floor plate", hConFloor = 4 "convective coefficient of heat transfer of floor, from Awbi (1998)", nFloor = 1 "number of RC elements of floor plate (min=1)", RFloor = {0.00331421908725} "vector of resistances of floor plate, from inside to outside", RFloorRem = 0.1265217391 "unchecked value - resistance of remaining resistor between capacity n and outside", CFloor = {1674000} "vector of heat capacities of floor plate, from inside to outside, carpet + concrete floor, from https://help.iesve.com/ve2021/table_6_thermal_conductivity__specific_heat_capacity_and_density.htm and https://www.mrsphysics.co.uk/bge/wp-content/uploads/2016/07/thermal-properties-of-building-materials.pdf", ARoof = 15 "area of roof plate", hConRoof = 0.5 "vector of resistances of roof plate, from inside to outside, from Awbi (1998)", nRoof = 1 "number of RC elements of roof plate (min=1)", RRoof = {0.00331421908725} "vector of resistances of roof plate, from inside to outside", RRoofRem = 0.1265217391 "unchecked value - resistance of remaining resistor between capacity n and outside", CRoof = {1520000} "vector of heat capacities of roof plate, from inside to outside, from https://www.mrsphysics.co.uk/bge/wp-content/uploads/2016/07/thermal-properties-of-building-materials.pdf", nOrientations = 2 "number of orientations (min=1)", AWin = {2, 2} "vector of areas of windows by orientations", ATransparent = {2, 2} "vector of areas of transparent (solar radiation transmitted) elements by orientations", AExt = {3*2.4, 5*2.4} "vector of areas of exterior walls by orientations", redeclare replaceable package Medium = Modelica.Media.Air.SimpleAir, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, extWallRC(thermCapExt(each der_T(fixed = true))), intWallRC(thermCapInt(each der_T(fixed = true))), floorRC(thermCapExt(each der_T(fixed = true))), T_start = 295.15, roofRC(thermCapExt(each der_T(fixed = true)))) "Thermal zone" annotation(
    Placement(visible = true, transformation(origin = {67.7299, 25.9102}, extent = {{-11.8701, -12.4898}, {11.8701, 12.4898}}, rotation = 0)));
  Buildings.ThermalZones.ReducedOrder.EquivalentAirTemperature.VDI6007WithWindow eqAirTemp(n = 2, wfGro = 0 "weight factor of the ground (0, if not considered)", wfWall = {0.3043478260869566, 0.6956521739130435} "weight factor of the walls", wfWin = {0.5, 0.5} "weight factor of windows", withLongwave = true, aExt = 0.65 "coefficient of absorption of exterior walls (outdoor)", hConWallOut = 24 "exterior walls convective coefficient of heat transfer (outdoor), from Chávez-Galán (2014)", hRad = 5 "coefficient of heat transfer for linearised radiation", hConWinOut = 24 "window convective coefficient of heat transfer (outdoor), same as external wall", TGro = 285.15 "temperature of ground in contact with floor plate") "Computes equivalent air temperature" annotation(
    Placement(visible = true, transformation(origin = {-36, -2}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Blocks.Math.Add solRad[2] "Sums up solar radiation of both directions" annotation(
    Placement(visible = true, transformation(origin = {-14, 0}, extent = {{-38, 6}, {-28, 16}}, rotation = 0)));
  Buildings.HeatTransfer.Sources.PrescribedTemperature preTem "Prescribed temperature for exterior walls outdoor surface temperature" annotation(
    Placement(visible = true, transformation(origin = {-96, -46}, extent = {{8, -6}, {20, 6}}, rotation = 0)));
  Buildings.HeatTransfer.Sources.PrescribedTemperature preTem1 "Prescribed temperature for windows outdoor surface temperature" annotation(
    Placement(visible = true, transformation(origin = {-10, 12}, extent = {{8, 14}, {20, 26}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.Convection theConWin "Outdoor convective heat transfer of windows" annotation(
    Placement(visible = true, transformation(origin = {-2, 8}, extent = {{38, 16}, {28, 26}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.Convection theConWall "Outdoor convective heat transfer of walls" annotation(
    Placement(visible = true, transformation(origin = {-90, -26}, extent = {{36, 6}, {26, -4}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow perRad "Radiative heat flow of persons" annotation(
    Placement(visible = true, transformation(origin = {59, -63}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow perCon "Convective heat flow of persons" annotation(
    Placement(visible = true, transformation(origin = {59, -77}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const[2](each k = 0) "sets sunblind signal to zero (open)" annotation(
    Placement(visible = true, transformation(origin = {-10, -6}, extent = {{-20, 14}, {-14, 20}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow macConv "convective heat flow of machines" annotation(
    Placement(visible = true, transformation(origin = {59, -91}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant hConWall(k = 25*11.5) "outdoor coefficient of heat transfer for walls" annotation(
    Placement(visible = true, transformation(origin = {-62, -76}, extent = {{-4, -4}, {4, 4}}, rotation = 90)));
  Modelica.Blocks.Sources.Constant hConWin(k = 20*14) "outdoor coefficient of heat transfer for windows" annotation(
    Placement(visible = true, transformation(origin = {30, 46}, extent = {{4, -4}, {-4, 4}}, rotation = 90)));
  Buildings.HeatTransfer.Sources.PrescribedTemperature preTemFloor "prescribed temperature for floor plate outdoor surface temperature" annotation(
    Placement(visible = true, transformation(origin = {67, -8}, extent = {{-6, -6}, {6, 6}}, rotation = 90)));
  Modelica.Blocks.Sources.Constant TSoil(k = 283.15) "Outdoor surface temperature for floor plate" annotation(
    Placement(visible = true, transformation(origin = {88, -16}, extent = {{-4, -4}, {4, 4}}, rotation = 180)));
  Buildings.BoundaryConditions.WeatherData.Bus weaBus annotation(
    Placement(visible = true, transformation(origin = {-81, 17}, extent = {{-9, -11}, {9, 11}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Buildings.ThermalZones.ReducedOrder.EquivalentAirTemperature.VDI6007 eqAirTempVDI(aExt = 0.65, n = 1, wfWall = {1}, wfWin = {0}, wfGro = 0, hConWallOut = 24, hRad = 5, TGro = 285.15) "Computes equivalent air temperature for roof" annotation(
    Placement(visible = true, transformation(origin = {-24, 84}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Buildings.HeatTransfer.Sources.PrescribedTemperature preTemRoof "Prescribed temperature for roof outdoor surface temperature" annotation(
    Placement(visible = true, transformation(origin = {17, 82}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  Modelica.Thermal.HeatTransfer.Components.Convection theConRoof "Outdoor convective heat transfer of roof" annotation(
    Placement(visible = true, transformation(origin = {45, 69}, extent = {{5, -5}, {-5, 5}}, rotation = -90)));
  Modelica.Blocks.Sources.Constant hConRoof(k = 20) "Outdoor coefficient of heat transfer for roof" annotation(
    Placement(visible = true, transformation(origin = {66, 69}, extent = {{4, -4}, {-4, 4}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant const1(k = 0) "Sets sunblind signal to zero (open)" annotation(
    Placement(visible = true, transformation(origin = {-72, -2}, extent = {{68, 90}, {62, 96}}, rotation = 0)));
  //______________________________________________________________________________________________________________________________________
  //Person
  parameter SI.CoefficientOfHeatTransfer hRad(min = 0, max = 10) = epsilon_cl*4.7 "radiative heat transfer coefficient";
  parameter Real C_dil(unit = "kg/(K.h.m2)") = 75 "coefficient of vasodilation";
  parameter Real C_cst(unit = "kg/(K.h.m2)") = 0.5 "coefficient of vasoconstriction";
  parameter SI.SpecificHeatCapacity C_b = 3490 "J/(kg.K), specific heat of body";
  //  parameter SI.ThermalConductance K_sk = 5.28 "W/K minimum heat conductance of skin tissue";
  parameter Real C_bl(unit = "J/(kg.K)") = 4187 "specific heat of blood";
  SI.DimensionlessRatio alpha = 0.1 "fractional skin mass";
  parameter SI.DimensionlessRatio F_eff = 0.72 "ratio of the effective area of the human bodyfor radiant heat exchange to the surface area of the human body clothing (=0.696~0.725)";
  parameter SI.DimensionlessRatio epsilon_cl = 0.9 "emittance of clothing surface, usually higher than 0.9";
  parameter SI.Density rho_sk = 1020 "skin density";
  parameter SI.Density rho = 985 "person density";
  parameter SI.ThermalConductivity K_person = 0.2 "average thermal conductivity of a person";
  parameter SI.SpecificHeatCapacity c_pa = 1005 "J/(kg.K)specific heat capacity of dry air";
  parameter SI.Density rho_w = 1000 "density of liquid water";
  parameter Real H_rb(unit = "W/(m2.K)") = 6 "radiative heat transfer coefficient of a black surface (=5.7~6.3)";
  parameter SI.MolarMass m_a = 28.97e-3 "(g/mol) molar mass of dry air";
  parameter SI.MolarMass m_w = 18.05e-3 "(g/mol) molar mass of water molecules";
  parameter SI.SpecificHeatCapacity c_pv = 1846 "J/(kg.K)specific heat capacity of water vapour";
  parameter SI.SpecificHeatCapacity c_pw = 4186 "J/(kg.K)specific heat capacity of liquid water";
  SI.Pressure p_vr "water-vapour pressure in the room";
  SI.Pressure p_vo "water-vapour pressure of outdoor air";
  SI.Pressure p_vs_To "saturated water-vapour pressure at outdoor air temperature";
  SI.Pressure p_vs_Ta "saturated water-vapour pressure at indoor air temperature";
  //______________________________________________________________________________________________________________________________________
  parameter SI.Length L_1 = 1.75 "person height";
  parameter Real BMI_1 = 23 "body mass index (from https://www.nhs.uk/live-well/healthy-weight/bmi-calculator/)
        if BMI < 18.5, underweight, if 18.5 <= BMI < 25, healthy, if 25 <= BMI < 30, overweight, if BMI > 30, obese, for white heritage people";
  parameter Real Age_1(unit = "years") = 40;
  parameter Boolean female_1 = false "true = female, false = male";
  Real ICl_1;
  SI.HeatFlux M_1 "1.2*58.2 W/m2 is the standard for office work";
  Modelica.Blocks.Interfaces.RealOutput PMV_1 "predicted mean vote, between -3 and 3 (-3 = cold, -2 = cool, 1 = slightly cool, 0 = neutral, +1 = slightly warm, +2 = warm, +3 = hot)";
  Modelica.Blocks.Interfaces.RealOutput PPD_1 "predicted percentage of dissatisfied, between 0 and 1, predicts the number of thermally dissatisfied people among a large group of people";
  parameter SI.HeatFlux W_1(max = 0) = 0 "rate of mechanical work accomplished (must be non-positive, typically equal to 0)";
  SI.Temperature TOpe_1 "operative temperature";
  SI.Temperature TClo_1(max = 273.15 + 36.8) "surface temperature of clothing";
  SI.Temperature TSki_1(min = 273.15 + 33.5, max = 273.15 + 36.8) "skin temperature";
  SI.Temperature T_cr_1(displayUnit = "degC") = 273.15 + 36.8 "core temperature";
  NonSI.Temperature_degC T_comf_1(displayUnit = "degC") "comfortable room temperature";
  SI.CoefficientOfHeatTransfer hCom_1(min = 0, max = 15) "combined heat transfer coefficient";
  SI.CoefficientOfHeatTransfer hCon_1(min = 0, max = 15) "convective heat transfer coefficient";
  SI.HeatFlux L_b_1 "thermal load of the body";
  Real fCl_1(min = 0) "clothing area factor";
  SI.ThermalInsulance RCl_1 "thermal resistance of clothing 1 clo = 0.155 K.m2/W";
  Real V_b_1 "mass blood flow between core and skin (kg/(h m^2)), it should be L/min";
  SI.Mass m_1 "mass of the person";
  Real BMR_1(unit = "kcal") "Harris-Benedict equation";
  SI.Area A_d_1 "DeBois surface area";
  SI.HeatFlux E_max_1 "maximum evaporative potential from the skin surface to the surrounding room space";
  SI.DimensionlessRatio w_1 "wettedness";
  SI.Length R_person_1 "radius of person (cylinder)";
  Real m_rsw_1(unit = "kg/(s.K.m2)") "sweat drive due to the temperature stimuli";
  Real V_in_1(unit = "(m3/s)/m2") "volumetric rate of inhaled air = volumetric rate of exhaled air (V_out)";
  Real V_w_core_1(unit = "(m3/s)/m2") "volumetric rate of liquid water generated in the body core, which turns into water vapour and is exhaled through the nose and the mouth";
  Real V_w_shell_rho_1 "volumetric rate of liquid water generated in the body shell as sweat and density";
  Real Q_core_1(unit = "J/(m2.K)") "heat capacity of the body core";
  Real Q_shell_1(unit = "J/(m2.K)") "heat capacity of the body shell";
  Real stevrate_1(unit = "g/m2.s") "sweat evaporation rate";
  SI.Pressure p_sk_1 "water-vapour pressure of liquid water at skin-surface temperature";
  SI.Pressure p_vs_Tcr_1 "saturatd water-vapour pressure at body-core temperature";
  //Exergy terms
  SI.HeatFlux E_xm_1 "warm exergy generated by metabolism";
  SI.HeatFlux E_inh_1 "warm/cool and wet/dry exergies of the inhaled humid air";
  SI.HeatFlux E_gen_cr_1 "warm and wet exergies of liquid water generated in the core by metabolism";
  SI.HeatFlux E_gen_sh_1 "warm/ cool and wet/dry exergies of the sum of liquid water generated in the shell by metabolism and dry air to let the liquid water disperse";
  SI.HeatFlux E_rad_abs_1 "warm/cool radiant exergy absorbed by the whole of skin and clothing surfaces";
  SI.HeatFlux E_cons_1 "exergy consumption rate, which is only for thermoregulation";
  SI.HeatFlux E_stored_1 "warm exergy stored in the core and the shell";
  SI.HeatFlux E_exh_1 "warm and wet exergies of the exhaled humid air";
  SI.HeatFlux E_sw_1 "warm/cool exergy of the water vapour originating from the sweat and wet/dry exergy of the humid air containing the evaporated sweat";
  SI.HeatFlux E_rad_dis_1 "warm/cool radiant exergy discharged from the whole of skin and clothing surfaces";
  SI.HeatFlux E_conv_1 "warm/cool exergy transferred by convection from the whole of skin and clothing surfaces into the surrounding air";
  //______________________________________________________________________________________________________________________________________
  //Person/room boundary
  SI.HeatFlowRate rad_1 "radiative profile of the person, in W";
  SI.HeatFlowRate conv_1 "convective profile of the person";
  SI.HeatFlowRate conv_mach_1 "convective profile of machine (desktop computer)";
  Modelica.Blocks.Sources.RealExpression realExpression_rad(y = rad_1) annotation(
    Placement(visible = true, transformation(origin = {19, -63}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression realExpression_conv(y = conv_1) annotation(
    Placement(visible = true, transformation(origin = {19, -77}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression realExpression_mach(y = conv_mach_1) annotation(
    Placement(visible = true, transformation(origin = {19, -91}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
  //______________________________________________________________________________________________________________________________________
  //HVAC
  SI.Temperature setT_function;
  Modelica.Blocks.Sources.RealExpression setT(y = setT_function) annotation(
    Placement(visible = true, transformation(origin = {-2, -36}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
  Buildings.Controls.Continuous.LimPID conHeaCoo(Ti(displayUnit = "ks") = 1, controllerType = Modelica.Blocks.Types.SimpleController.P, k = 1, yMax = 1, yMin = 0) "Heating and cooling controller" annotation(
    Placement(visible = true, transformation(origin = {20, -36}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Blocks.Math.Gain gainHeaCoo(k = 0) "Gain for heating and cooling controller, numerically equal to Q_flow nominal, in W" annotation(
    Placement(visible = true, transformation(origin = {38, -36}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow heaCoo "Ideal heater/cooler with limit" annotation(
    Placement(visible = true, transformation(origin = {56, -36}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor "Sensor for ideal heater/cooler" annotation(
    Placement(visible = true, transformation(origin = {0, 0}, extent = {{84, -42}, {72, -30}}, rotation = 0)));
  Modelica.Units.SI.Efficiency eta(max = 1) "efficiency of electrical heater";
  Real power(unit = "W");
  Real energy(unit = "kWh");
  Real price(unit = "£") "assuming £0.34 as average price for kWh in the UK";
  //______________________________________________________________________________________________________________________________________
  //╔════════════╦═════════════════════╗
  //║gagge_2_node║Declarations         ║
  //╚════════════╩═════════════════════╝
  constant Real clo = 0.5 "clothing insulation level, clo";
  constant Real met = 1 "metabolic rate, met";
  constant Real wme = 0 "external work, met";
  constant Real pb = 760 "barometric pressure, torr or mmHg";
  constant Real ht = 170 "height, cm";
  constant Real wt = 70 "weight, kg";
  constant Real tskn = 33.7 "setpoint value for skin temperature, °C";
  constant Real tcrn = 36.8 "setpoint value for core temperature, °C";
  constant Real tbn = 36.49 "setpoint value for mean body temperature (.1*tskn + .9*tcrn), °C";
  constant Real skbfn = 6.3 "neutral value for skin blood flow";
  constant Real sbc = 5.6697*10^(-8) "stephan-Boltzmann constant";
  constant Real sa = 0.203*(ht/100)^(0.725)*wt^(0.425) "Dubois surface Area, m^2";
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
  Real m(start = met*58.2) "metabolic rate, W";
  Real ta(start = 20) "air temperature, °C";
  Real tr(start = 20) "mean radiant temperature, °C";
  Real pa "partial vapour pressure of water, mmHg";
  Real tsk(start = tskn, fixed = true) "skin temperature, °C";
  Real tcr(start = tcrn, fixed = true) "core temperature, °C";
  Real tcl(start = tskn) "clothing temperature, °C";
  Real tb(start = 0.1*tskn + (1 - 0.1)*tcrn) "mean body temperateu,°C";
  Real skbf(start = skbfn) "skin blood flow, kg/hr m^2";
  Real skbf_ "skin blood flow test";
  Real mshiv(start = 0) "rate of energy released by shivering, W";
  Real alpha_(start = 0.1) "fractional skin mass, 1";
  Real esk(start = 0.1*met) "total evaporative heat loss from the skin W";
  Real wcrit "evaporative efficiency, 1";
  Real icl "";
  Real chr(start = 4.7) "radiative heat transfer coefficient";
  Real ctc(start = 4.7 + (3*atm^(0.53))) "?combined convection & radiation coefficient";
  Real ra(start = 1/(facl*4.7 + (3*atm^(0.53)))) "resistance of air layer to dry heat transfer, °C m^2 / W";
  Real top(start = ((4.7*20) + (3*atm^(0.53))*20)/4.7 + (3*atm^(0.53))) "operative temperature, °C";
  Real dry "total sensible heat loss, W";
  Real hfcs "rate of energy transport between core and skin, W";
  Real eres "heat loss through respiratory evaporation, W";
  Real cres "heat loss through respiratory convection, W";
  Real scr "rate of energy storage in the core, W";
  Real ssk "rate of enery storage in the skin, ";
  Real warms "skin warm signal ,1";
  Real colds "skin cold signal ,1";
  Real warmc "core warm signal ,1";
  Real coldc "core cold signal ,1";
  Real warmb "blood warm signal ,1";
  Real coldb "blood cold signal ,1";
  Real regsw "regulatory sweating, g/m^2 hr";
  Real ersw "heat loss through sweating, W";
  Real ersw_ "also heat loss through sweating, W";
  Real recl "resistance of clothing layer to dry heat transfer, °C m^2 / W";
  Real emax "maximum evapourative capacity, W";
  Real pwet "skin wettedness, 1";
  Real prsw "ratio of actual heat loss due to sweating to maximum heat loss due to sweating";
  Real edif "heat loss due to diffusion of water fvapout rhough the skin";
  Real esk_ "also total evaporative heat loss from the skin W";
  Real PMV "predicted mean vote";
  Real PPD "predicted percentage dissatisfied";
  Real L_b "thermal load of body";
  Real E_xm "warm exergy generated by metabolism";
  Real E_inh "warm/cool and wet/dry exergies of the inhaled humid air";
  Real V_in (unit = "(m3/s)/m2") "volumetric rate of inhaled air = volumetric rate of exhaled air (V_out)";
  Real E_gen_cr "warm and wet exergies of the liquid water generated in the core by metabolism";
  Real V_w_core (unit = "(m3/s)/m2") "volumetric rate of liquid water generated in the body core, which turns into water vapour and is exhaled through the nose and the mouth";
  Real E_gen_sh "warm/cool and wet/dry exergies of the sum of liquid water generated in the shell by metabolism and dry air to let the liquid water disperse";
  Real V_w_shell_rho "volumetric rate of liquid water generated in the body shell as sweat and density";
  Real wetted "skin wettedness";
  Real E_rad_abs "Warm/cool radiant exergy absorbed by the whole of skin and clothing surfaces";
  Real E_stored "warm exergy stored in the core and the shell";
  Real Q_core (unit = "J/(m2.K)") "heat capacity of the body core";
  Real E_exh "warm and wet exergies of the exhaled humid air";
  Real p_vs_Tcr "latent pressure something somthing...";
  Real E_sw "warm/cool exergy of the water vapour originating from the sweat and wet/dry exergy of the humid air containing the evaporated sweat";
  Real E_rad_dis "warm/cool radiant exergy discharged from the whole of skin and clothing surfaces";
  Real E_conv "Warm/cool exergy transferred by convection from the whole of skin and clothing surfaces into the surrounding area";
  Real E_cons "exergy conservation equation";
  
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
  //╔════════════╦═════════════════════╗
  //║gagge_2_node║End of declarations  ║
  //╚════════════╩═════════════════════╝
protected
  Modelica.Blocks.Interfaces.RealInput vAir_in_internal "Needed to connect to conditional connector";
  Modelica.Blocks.Interfaces.RealInput M_in_internal(final quantity = "HeatFlux", final unit = "W/m2") "Needed to connect to conditional connector";
  Modelica.Blocks.Interfaces.RealInput ICl_in_internal "Needed to connect to conditional connector";
  Modelica.Blocks.Interfaces.RealInput pAir_in_internal(final quantity = "Pressure", final unit = "Pa", min = 0) "Needed to connect to conditional connector";
  //______________________________________________________________________________________________________________________________________
initial equation
  assert(W_1 <= 0, "Parameter W must be equal to zero or negative.");
//______________________________________________________________________________________________________________________________________
equation
  //╔════════════╦═════════════════════╗
  //║gagge_2_node║Equations            ║
  //╚════════════╩═════════════════════╝
  ta = thermalZoneFourElements.TAir-273.15; //ambient temperature
  tr = thermalZoneFourElements.TRad-273.15; //mean radiant temperature

  icl = 0.25*Modelica.Math.cos((2*Modelica.Constants.pi*time/(3600*24*365))) + 0.75;
  wcrit = 0.59*vel^(-0.08);

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
  der(tsk) = (ssk*sa)/(0.97*alpha_*wt)/3600;
//°C / s
  der(tcr) = (scr*sa)/(0.97*(1 - alpha_)*wt)/3600;
//mean body temperature
  tb = alpha_*tsk + (1 - alpha_)*tcr;
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
  alpha_ = 0.0417737 + 0.7451833/(skbf + 0.585417);
  
//calculate PMV/PPD

L_b = (rmm - w) - 3.05E-3*(5733 - 6.99*(rmm - w) - pa) - 0.42*((rmm - w) - 58.15) - 1.7E-5*rmm*(5867 - (pa+273.15)) - 0.0014*rmm*(307.15 - thermalZoneFourElements.TAir) - 3.96E-8*facl*((tcl+273.15)^4 - thermalZoneFourElements.TRad^4) - facl*chc*((tcl+273.15) - thermalZoneFourElements.TAir);
  PMV = (0.303*Modelica.Math.exp(-0.036*m) + 0.028)*L_b;
  PPD = 1 - 0.95*Modelica.Math.exp(-(0.03353*PMV^4 + 0.2179*PMV^2));


//EXERGY TERMS
//From "Low energy systems for high performance buildings and communities" Working report of IEA ECBS annex 49
  E_xm = m*(1 - (eqAirTemp.TDryBul/T_cr_1)); //"warm exergy generated by metabolism"
  
  V_in = 1.2e-6*m;
  E_inh = V_in*((c_pa*m_a*(pAir_in - p_vr)/(Modelica.Constants.R*thermalZoneFourElements.TAir) + c_pv*m_w*p_vr/(Modelica.Constants.R*thermalZoneFourElements.TAir))*((thermalZoneFourElements.TAir - eqAirTemp.TDryBul) - eqAirTemp.TDryBul*Modelica.Math.log(thermalZoneFourElements.TAir/eqAirTemp.TDryBul)) + (eqAirTemp.TDryBul/thermalZoneFourElements.TAir)*((pAir_in - p_vr)*Modelica.Math.log((pAir_in - p_vr)/(pAir_in - p_vo)) + p_vr*Modelica.Math.log(p_vr/p_vo))); //"warm/cool and wet/dry exergies of the inhaled humid air"


  V_w_core = 1.2e-6*m*(0.029 - 0.049e-4*p_vr);
  E_gen_cr = V_w_core*rho_w*(c_pw*((tcr+273.15) - eqAirTemp.TDryBul - eqAirTemp.TDryBul*Modelica.Math.log((tcr+273.15)/eqAirTemp.TDryBul)) + (Modelica.Constants.R*eqAirTemp.TDryBul/m_w)*Modelica.Math.log(p_vs_To/p_vo)); //"warm and wet exergies of the liquid water generated in the core by metabolism"

//WHY IS THIS ZERO?

// add me
wetted = -(3.05E-3*(5733 - 6.99*(m) - p_vr) + 0.42*(m) - 58.15)/emax; //this is probably wrong

  V_w_shell_rho = (wetted*(emax/sa))/(2450*1000) "2450 J/g, latent heat value of evaporation of liquid water at 30 degC";
  
  E_gen_sh = V_w_shell_rho*(c_pw*((tsk+273.15)- eqAirTemp.TDryBul - eqAirTemp.TDryBul*Modelica.Math.log((tsk+273.15)/eqAirTemp.TDryBul)) + (Modelica.Constants.R*eqAirTemp.TDryBul/m_w)*(Modelica.Math.log(p_vs_To/p_vo) + ((pAir_in - p_vr)/p_vr)*Modelica.Math.log((pAir_in - p_vr)/(pAir_in - p_vo)))); // warm/cool and wet/dry exergies of the sum of liquid water generated in the shell by metabolism and dry air to let the liquid water disperse
  
E_rad_abs = F_eff*facl*eqAirTemp.aExt*epsilon_cl*H_rb*((eqAirTemp.TEqAir - eqAirTemp.TDryBul)^2)/(eqAirTemp.TEqAir + eqAirTemp.TDryBul); //Warm/cool radiant exergy absorbed by the whole of skin and clothing surfaces
//e_stored "warm exergy stored in the core and the shell"

Q_core = (1 - alpha_)*(wt/sa)*C_b;
E_stored = Q_core*(1 - (eqAirTemp.TDryBul/(tcr+273.15)))*der((tcr+273.15)) + Q_shell_1*(1 - (eqAirTemp.TDryBul/(tsk+273.15)))*der((tsk+273.15));

p_vs_Tcr = Modelica.Math.exp(25.89 - 5319/(tcr+273.15));

E_exh = V_in_1*(((c_pa*m_a/(Modelica.Constants.R*(tcr+273.15)))*(pAir_in - p_vs_Tcr) + (c_pv*m_w/(Modelica.Constants.R*(tcr+273.15)))*p_vs_Tcr)*((tcr+273.15) - eqAirTemp.TDryBul - eqAirTemp.TDryBul*Modelica.Math.log((tcr+273.15)/eqAirTemp.TDryBul)) + (eqAirTemp.TDryBul/(tcr+273.15))*((pAir_in - p_vs_Tcr)*Modelica.Math.log((pAir_in - p_vs_Tcr)/(pAir_in - p_vo)) + (p_vs_Tcr*Modelica.Math.log(p_vs_Tcr/p_vo)))); //"warm and wet exergies of the exhaled humid air"

E_sw = V_w_shell_rho*(c_pv*((tcl+273.15) - eqAirTemp.TDryBul - eqAirTemp.TDryBul*Modelica.Math.log((tcl+273.15)/eqAirTemp.TDryBul)) + (Modelica.Constants.R*eqAirTemp.TDryBul/m_w)*(Modelica.Math.log(p_vr/p_vo) + ((pAir_in - p_vr)/p_vr)*Modelica.Math.log((pAir_in - p_vr)/(pAir_in - p_vo)))); //warm/cool exergy of the water vapour originating from the sweat and wet/dry exergy of the humid air containing the evaporated sweat

E_rad_dis = F_eff*facl*epsilon_cl*H_rb*((tcl+273.15) - eqAirTemp.TDryBul)^2/((tcl+273.15) + eqAirTemp.TDryBul); // warm/cool radiant exergy discharged from the whole of skin and clothing surfaces

E_conv = facl*chc*((tcl+273.15) - thermalZoneFourElements.TAir)*(1 - (eqAirTemp.TDryBul/(tcl+273.15)));

//Complete exergy conservation equation
//#####################################
E_cons = E_xm + E_inh + E_gen_cr + E_gen_sh + E_rad_abs - E_stored - E_exh - E_sw - E_rad_dis - E_conv;
//#####################################


  //╔════════════╦═════════════════════╗
  //║gagge_2_node║End of equations     ║
  //╚════════════╩═════════════════════╝
//Common calculations
  phi = 0.05*Modelica.Math.cos((0.4 + 2*Modelica.Constants.pi*time/(3600*24*365))) + 0.8;
// phi = 0.8;
  p_vr = 1000*0.61094*Modelica.Math.exp(17.625*(thermalZoneFourElements.TAir - 273.15)/(thermalZoneFourElements.TAir - 273.15 + 243.04));
  p_vo = 1000*0.61094*Modelica.Math.exp(17.625*(eqAirTemp.TDryBul - 273.15)/(eqAirTemp.TDryBul - 273.15 + 243.04));
  p_vs_To = Modelica.Math.exp(25.89 - 5319/eqAirTemp.TDryBul);
  p_vs_Ta = Modelica.Math.exp(25.89 - 5319/thermalZoneFourElements.TAir);
//______________________________________________________________________________________________________________________________________
//Person #1
//ICl is a sinusoidal function over whole year
  ICl_1 = 0.25*Modelica.Math.cos((2*Modelica.Constants.pi*time/(3600*24*365))) + 0.75;
// ICl set for specific days
// ICl_1 = 0.5;
  m_1 = L_1^2*BMI_1;
  A_d_1 = 0.007184*((100*L_1)^0.725)*(m_1^0.425);
  if female_1 == true then
    BMR_1 = 447.593 + 9.247*m_1 + 3.098*L_1*100 - 4.330*Age_1;
  else
    BMR_1 = 88.362 + 13.397*m_1 + 4.799*L_1*100 - 5.677*Age_1;
  end if;
  M_1 = 1.2*BMR_1*1000*4.18/(86400*A_d_1*0.7);
//skin blood flow
  if (T_cr_1 - 309.95) >= 0 and (TSki_1 - 306.85) > 0 then
//if core temp >= 36.8 °C AND skin temp > 33.7 °C?
    V_b_1 = 6.3 + (C_dil*(T_cr_1 - 36.8 - 273.15)/(1 - C_cst*(TSki_1 - (273.15 + 33.7))));
//vasodilation & vasoconstriction eqn
  else
    V_b_1 = 6.3;
//core temp is ALWAYS 36.8 °C so V_b_1 is always 6.3
  end if;
  E_max_1 = 2.2*hCon_1*(Modelica.Math.exp(25.89 - 5319/TSki_1) - 1000*0.61094*Modelica.Math.exp(17.625*(thermalZoneFourElements.TAir - 273.15)/(thermalZoneFourElements.TAir - 273.15 + 243.04)))/1000 "pressures in kPa";
  w_1 = (3.05E-3*(5733 - 6.99*(M_1 - W_1) - p_vr) + 0.42*((M_1 - W_1) - 58.15))/E_max_1;
// clothing insulation value in SI units
  RCl_1 = 0.155*ICl_1;
// clothing area factor
  if ICl_1 <= 0.5 then
    fCl_1 = 1 + 0.2*ICl_1;
  else
    fCl_1 = 1.05 + 0.1*ICl_1;
  end if;
//convective heat transfer coefficient
  if 2.38*abs(TClo_1 - thermalZoneFourElements.TRad)^0.25 >= 12.1*(vAir^0.5) then
    hCon_1 = 2.38*abs(TClo_1 - thermalZoneFourElements.TRad)^0.25;
  else
    hCon_1 = 12.1*(vAir^0.5);
  end if;
  hCom_1 = hRad + hCon_1;
  TOpe_1 = (hRad*thermalZoneFourElements.TRad + hCon_1*thermalZoneFourElements.TAir)/hCom_1;
//Skin temperature from Butera (1998): TSki = 35.7 + 273.15 - 0.0275*(M - W);
  R_person_1 = sqrt(m_1/(Modelica.Constants.pi*rho*L_1));
//==================================
//T_skin is constant
  TSki_1 = T_cr_1 - (1 - sqrt((rho_sk - alpha*rho)/rho_sk))*M_1*R_person_1/K_person;
//==================================
//TClo from ISO 7730 and Olesen (1982) Olesen, Bjarne W. Thermal comfort. Vol. 2. Bruel & Kjaer, 1982.
  TClo_1 = TSki_1 - RCl_1*(3.96E-8*fCl_1*(TClo_1^4 - thermalZoneFourElements.TRad^4) + fCl_1*hCon_1*(TClo_1 - thermalZoneFourElements.TAir));
  T_comf_1 = -0.67*(thermalZoneFourElements.TAir - 273.15) - 0.048*(eqAirTemp.TDryBul - 273.15) - 0.35*M_1 - 86.46*ICl_1*0.155 + 0.113*100*phi + 66.53;
  m_rsw_1 = (170*(1 - alpha)*(T_cr_1 - (273.15 + 36.8)) + alpha*(TSki_1 - (273.15 + 33.7)))*(Modelica.Math.exp((TSki_1 - (273.15 + 33.7))/10.7))/(60*60*1000) "170 g/(h.K.m2) is C_sw, the driving coefficient for regulatory sweating";
  stevrate_1 = 1000*m_rsw_1*TSki_1;
  V_in_1 = 1.2e-6*M_1;
  V_w_core_1 = 1.2e-6*M_1*(0.029 - 0.049e-4*p_vr);
  p_sk_1 = Modelica.Math.exp(25.89 - 5319/TSki_1);
  p_vs_Tcr_1 = Modelica.Math.exp(25.89 - 5319/T_cr_1);
  V_w_shell_rho_1 = (w_1*E_max_1)/(2450*1000) "2450 J/g, latent heat value of evaporation of liquid water at 30 degC";
  Q_core_1 = (1 - alpha)*(m_1/A_d_1)*C_b;
  Q_shell_1 = alpha*(m_1/A_d_1)*C_b;
  L_b_1 = (M_1 - W_1) - 3.05E-3*(5733 - 6.99*(M_1 - W_1) - p_vr) - 0.42*((M_1 - W_1) - 58.15) - 1.7E-5*M_1*(5867 - p_vr) - 0.0014*M_1*(307.15 - thermalZoneFourElements.TAir) - 3.96E-8*fCl_1*(TClo_1^4 - thermalZoneFourElements.TRad^4) - fCl_1*hCon_1*(TClo_1 - thermalZoneFourElements.TAir);
  PMV_1 = (0.303*Modelica.Math.exp(-0.036*M_1) + 0.028)*L_b_1;
  PPD_1 = 1 - 0.95*Modelica.Math.exp(-(0.03353*PMV_1^4 + 0.2179*PMV_1^2));
//Exergy terms
//From "Low energy systems for high performance buildings and communities" Working report of IEA ECBS annex 49
//Subscript 1 denotes person #1
//E_xm "warm exergy generated by metabolism"
  E_xm_1 = M_1*(1 - (eqAirTemp.TDryBul/T_cr_1));
//E_inh "warm/cool and wet/dry exergies of the inhaled humid air"
  E_inh_1 = V_in_1*((c_pa*m_a*(pAir_in - p_vr)/(Modelica.Constants.R*thermalZoneFourElements.TAir) + c_pv*m_w*p_vr/(Modelica.Constants.R*thermalZoneFourElements.TAir))*((thermalZoneFourElements.TAir - eqAirTemp.TDryBul) - eqAirTemp.TDryBul*Modelica.Math.log(thermalZoneFourElements.TAir/eqAirTemp.TDryBul)) + (eqAirTemp.TDryBul/thermalZoneFourElements.TAir)*((pAir_in - p_vr)*Modelica.Math.log((pAir_in - p_vr)/(pAir_in - p_vo)) + p_vr*Modelica.Math.log(p_vr/p_vo)));
//E_gen_cr "warm and wet exergies of the liquid water generated in the core by metabolism"
  E_gen_cr_1 = V_w_core_1*rho_w*(c_pw*(T_cr_1 - eqAirTemp.TDryBul - eqAirTemp.TDryBul*Modelica.Math.log(T_cr_1/eqAirTemp.TDryBul)) + (Modelica.Constants.R*eqAirTemp.TDryBul/m_w)*Modelica.Math.log(p_vs_To/p_vo));
//E_gen_sh "warm/cool and wet/dry exergies of the sum of liquid water generated in the shell by metabolism and dry air to let the liquid water disperse"
  E_gen_sh_1 = V_w_shell_rho_1*(c_pw*(TSki_1 - eqAirTemp.TDryBul - eqAirTemp.TDryBul*Modelica.Math.log(TSki_1/eqAirTemp.TDryBul)) + (Modelica.Constants.R*eqAirTemp.TDryBul/m_w)*(Modelica.Math.log(p_vs_To/p_vo) + ((pAir_in - p_vr)/p_vr)*Modelica.Math.log((pAir_in - p_vr)/(pAir_in - p_vo))));
//E_rad_abs "Warm/cool radiant exergy absorbed by the whole of skin and clothing surfaces"
  E_rad_abs_1 = F_eff*fCl_1*eqAirTemp.aExt*epsilon_cl*H_rb*((eqAirTemp.TEqAir - eqAirTemp.TDryBul)^2)/(eqAirTemp.TEqAir + eqAirTemp.TDryBul);
//e_stored "warm exergy stored in the core and the shell"
  if der(T_cr_1) == 0 and der(TSki_1) == 0 then
    E_stored_1 = 0;
  else
    E_stored_1 = Q_core_1*(1 - (eqAirTemp.TDryBul/T_cr_1))*der(T_cr_1) + Q_shell_1*(1 - (eqAirTemp.TDryBul/TSki_1))*der(TSki_1);
  end if;
//E_exh "warm and wet exergies of the exhaled humid air"
  E_exh_1 = V_in_1*(((c_pa*m_a/(Modelica.Constants.R*T_cr_1))*(pAir_in - p_vs_Tcr_1) + (c_pv*m_w/(Modelica.Constants.R*T_cr_1))*p_vs_Tcr_1)*(T_cr_1 - eqAirTemp.TDryBul - eqAirTemp.TDryBul*Modelica.Math.log(T_cr_1/eqAirTemp.TDryBul)) + (eqAirTemp.TDryBul/T_cr_1)*((pAir_in - p_vs_Tcr_1)*Modelica.Math.log((pAir_in - p_vs_Tcr_1)/(pAir_in - p_vo)) + (p_vs_Tcr_1*Modelica.Math.log(p_vs_Tcr_1/p_vo))));
//E_sw "warm/cool exergy of the water vapour originating from the sweat and wet/dry exergy of the humid air containing the evaporated sweat"
  E_sw_1 = V_w_shell_rho_1*(c_pv*(TClo_1 - eqAirTemp.TDryBul - eqAirTemp.TDryBul*Modelica.Math.log(TClo_1/eqAirTemp.TDryBul)) + (Modelica.Constants.R*eqAirTemp.TDryBul/m_w)*(Modelica.Math.log(p_vr/p_vo) + ((pAir_in - p_vr)/p_vr)*Modelica.Math.log((pAir_in - p_vr)/(pAir_in - p_vo))));
//E_rad_dis "warm/cool radiant exergy discharged from the whole of skin and clothing surfaces"
  E_rad_dis_1 = F_eff*fCl_1*epsilon_cl*H_rb*(TClo_1 - eqAirTemp.TDryBul)^2/(TClo_1 + eqAirTemp.TDryBul);
//Warm/cool exergy transferred by convection from the whole of skin and clothing surfaces into the surrounding area
  E_conv_1 = fCl_1*hCon_1*(TClo_1 - thermalZoneFourElements.TAir)*(1 - (eqAirTemp.TDryBul/TClo_1));
//Complete exergy conservation equation
//#####################################
  E_cons_1 = E_xm_1 + E_inh_1 + E_gen_cr_1 + E_gen_sh_1 + E_rad_abs_1 - E_stored_1 - E_exh_1 - E_sw_1 - E_rad_dis_1 - E_conv_1;
//#####################################
//if using day data, comment out mod_week and mod_winter, setting them to the right value
  if mod_week < 432000 then
    if mod_day < 21600 then
      conv_1 = 0;
    elseif mod_day >= 21600 and mod_day < 64800 then
      conv_1 = (M_1*A_d_1 - (Q_core_1*A_d_1*(T_cr_1 - TSki_1)/86400 + Q_shell_1*A_d_1*(TClo_1 - thermalZoneFourElements.TAir)/86400))/2.5 "fraction of radiative heat is 1.5x the convective contribution";
    else
      conv_1 = 0;
    end if;
  else
    conv_1 = 0;
  end if;
  rad_1 = 1.5*conv_1;
  if mod_week < 432000 then
    if mod_day < 21600 then
      conv_mach_1 = 0;
    elseif mod_day >= 21600 and mod_day < 64800 then
      conv_mach_1 = 100;
    else
      conv_mach_1 = 0;
    end if;
  else
    conv_mach_1 = 0;
  end if;
//______________________________________________________________________________________________________________________________________
//time pulses
  mod_day = OpenModelica.Internal.realRem(x = time, y = 86400.0);
  mod_week = OpenModelica.Internal.realRem(x = time, y = 604800.0);
  mod_winter = OpenModelica.Internal.realRem(x = time, y = 270*86400.0);
//constant temperature control during winter (October to April)
  if mod_winter < 122*86400 then
    setT_function = 22 + 273.15;
 else
   setT_function = 273.15;
 end if;
  
  
//Temperature control turned off during evenings and weekends
//  if mod_winter < 122*86400 then
//      if mod_week < 432000 then
//      if mod_day < 21600 then
//        setT_function = 273.15;
//      elseif mod_day >= 21600 and mod_day < 64800 then
//        setT_function = 22+273.15;
//      else
//        setT_function = 273.15;
//      end if;
//    else
//      setT_function= 273.15;
//    end if;
//  else
//     setT_function = 273.15;
//  end if;
  if gainHeaCoo.y < 0 then
    eta = 0.1 "cooling";
  elseif gainHeaCoo.y > 0 then
    eta = 0.9 "heating, boiler";
  else
    eta = 0 "off";
  end if;
// comment out for heating?
//  if eta == 0 then
//  power = 0;
//  else
//  power = abs(heatFlowSensor.Q_flow/eta);
//  end if;

// comment out for??
  if mod_day < 21600 then
   power = 0;
  elseif mod_day >= 21600 and mod_day < 64800 and eta > 0 then
    power = abs(heatFlowSensor.Q_flow/eta);
  else
  power = 0;
  end if;

  der(energy) = power/(1000*3600) "conversions to have kWh";
  price = energy*0.34;
//______________________________________________________________________________________________________________________________________
// Conditional connectors
  connect(vAir_in, vAir_in_internal);
  if not use_vAir_in then
    vAir_in_internal = vAir;
  end if;
  connect(M_in, M_in_internal);
  if not use_M_in then
    M_in_internal = M_1;
  end if;
  connect(ICl_in, ICl_in_internal);
  if not use_ICl_in then
    ICl_in_internal = ICl_1;
  end if;
  connect(pAir_in, pAir_in_internal);
  if not use_pAir_in then
    pAir_in_internal = pAir;
  end if;
//Connectors
  connect(thermalZoneFourElements.TAir, conHeaCoo.u_m) annotation(
    Line(points = {{80, 37}, {80, -54.4}, {20, -54.4}, {20, -43}}, color = {0, 0, 127}));
  connect(conHeaCoo.y, gainHeaCoo.u) annotation(
    Line(points = {{26.6, -36}, {32.6, -36}}, color = {0, 0, 127}));
  connect(gainHeaCoo.y, heaCoo.Q_flow) annotation(
    Line(points = {{42.4, -36}, {48.4, -36}}, color = {0, 0, 127}));
  connect(heaCoo.port, heatFlowSensor.port_b) annotation(
    Line(points = {{64, -36}, {72, -36}}, color = {191, 0, 0}));
  connect(heatFlowSensor.port_a, thermalZoneFourElements.intGainsConv) annotation(
    Line(points = {{84, -36}, {96, -36}, {96, 29}, {80, 29}}, color = {191, 0, 0}));
  connect(HDirTil.inc, corGDouPan.inc) annotation(
    Line(points = {{-53, 44}, {-46, 44}, {-46, 48}, {-39, 48}}, color = {0, 0, 127}, thickness = 0.5));
  connect(eqAirTemp.TEqAirWin, preTem1.T) annotation(
    Line(points = {{-29, 0}, {-6.1, 0}, {-6.1, 32}, {-3, 32}}, color = {0, 0, 127}));
  connect(eqAirTemp.TEqAir, preTem.T) annotation(
    Line(points = {{-29, -2}, {-21.1, -2}, {-21.1, -46}, {-89, -46}}, color = {0, 0, 127}));
  connect(weaDat.weaBus, weaBus) annotation(
    Line(points = {{-72, 84}, {-72, 18}, {-84, 18}, {-84, 17}, {-81, 17}}, color = {255, 204, 51}, thickness = 0.5));
  connect(weaBus.TDryBul, eqAirTemp.TDryBul) annotation(
    Line(points = {{-81, 17}, {-10, 17}, {-10, -6}, {-43, -6}}, color = {255, 204, 51}, thickness = 0.5));
  connect(const.y, eqAirTemp.sunblind) annotation(
    Line(points = {{-24, 11}, {-20.35, 11}, {-20.35, 5}, {-36, 5}}, color = {0, 0, 127}));
  connect(HDifTil.HSkyDifTil, corGDouPan.HSkyDifTil) annotation(
    Line(points = {{-53, 30}, {-46, 30}, {-46, 53}, {-39, 53}}, color = {0, 0, 127}));
  connect(HDirTil.H, corGDouPan.HDirTil) annotation(
    Line(points = {{-53, 46}, {-46, 46}, {-46, 56}, {-39, 56}}, color = {0, 0, 127}));
  connect(HDirTil.H, solRad.u1) annotation(
    Line(points = {{-53, 46}, {-42, 46}, {-42, 14}, {-53, 14}}, color = {0, 0, 127}));
  connect(HDifTil.H, solRad.u2) annotation(
    Line(points = {{-53, 26}, {-53, 8}}, color = {0, 0, 127}));
  connect(HDifTil.HGroDifTil, corGDouPan.HGroDifTil) annotation(
    Line(points = {{-53, 22}, {-46, 22}, {-46, 51}, {-39, 51}}, color = {0, 0, 127}));
  connect(solRad.y, eqAirTemp.HSol) annotation(
    Line(points = {{-41.5, 11}, {-38.25, 11}, {-38.25, 2}, {-43, 2}}, color = {0, 0, 127}));
  connect(weaDat.weaBus, HDifTil[1].weaBus) annotation(
    Line(points = {{-72, 84}, {-72, 57}, {-66, 57}, {-66, 26}}, color = {255, 204, 51}, thickness = 0.5));
  connect(weaDat.weaBus, HDifTil[2].weaBus) annotation(
    Line(points = {{-72, 84}, {-72, 57}, {-66, 57}, {-66, 26}}, color = {255, 204, 51}, thickness = 0.5));
  connect(weaDat.weaBus, HDirTil[1].weaBus) annotation(
    Line(points = {{-72, 84}, {-72, 65}, {-66, 65}, {-66, 46}}, color = {255, 204, 51}, thickness = 0.5));
  connect(weaDat.weaBus, HDirTil[2].weaBus) annotation(
    Line(points = {{-72, 84}, {-72, 65}, {-66, 65}, {-66, 46}}, color = {255, 204, 51}, thickness = 0.5));
  connect(perRad.port, thermalZoneFourElements.intGainsRad) annotation(
    Line(points = {{68, -63}, {100, -63}, {100, 31}, {80, 31}}, color = {191, 0, 0}));
  connect(theConWin.solid, thermalZoneFourElements.window) annotation(
    Line(points = {{36, 29}, {56, 29}}, color = {191, 0, 0}));
  connect(preTem1.port, theConWin.fluid) annotation(
    Line(points = {{10, 32}, {26, 32}, {26, 29}}, color = {191, 0, 0}));
  connect(thermalZoneFourElements.extWall, theConWall.solid) annotation(
    Line(points = {{56, 23}, {56, -25}, {-54, -25}}, color = {191, 0, 0}));
  connect(theConWall.fluid, preTem.port) annotation(
    Line(points = {{-64, -25}, {-64, -30.5}, {-76, -30.5}, {-76, -46}}, color = {191, 0, 0}));
  connect(hConWall.y, theConWall.Gc) annotation(
    Line(points = {{-62, -72}, {-62, -30}, {-59, -30}}, color = {0, 0, 127}));
  connect(hConWin.y, theConWin.Gc) annotation(
    Line(points = {{30, 42}, {30, 34}, {31, 34}}, color = {0, 0, 127}));
  connect(weaBus.TBlaSky, eqAirTemp.TBlaSky) annotation(
    Line(points = {{-81, 17}, {-58, 17}, {-58, -2}, {-43, -2}}, color = {255, 204, 51}, thickness = 0.5));
  connect(macConv.port, thermalZoneFourElements.intGainsConv) annotation(
    Line(points = {{68, -91}, {96, -91}, {96, 29}, {80, 29}}, color = {191, 0, 0}));
  connect(perCon.port, thermalZoneFourElements.intGainsConv) annotation(
    Line(points = {{68, -77}, {96, -77}, {96, 29}, {80, 29}}, color = {191, 0, 0}));
  connect(preTemFloor.port, thermalZoneFourElements.floor) annotation(
    Line(points = {{67, -2}, {67, 13}, {68, 13}}, color = {191, 0, 0}));
  connect(TSoil.y, preTemFloor.T) annotation(
    Line(points = {{84, -16}, {69.4, -16}, {69.4, -15}, {67, -15}}, color = {0, 0, 127}));
  connect(preTemRoof.port, theConRoof.fluid) annotation(
    Line(points = {{17, 76}, {31, 76}, {31, 74}, {45, 74}}, color = {191, 0, 0}));
  connect(theConRoof.solid, thermalZoneFourElements.roof) annotation(
    Line(points = {{45, 64}, {45, 38}, {67, 38}}, color = {191, 0, 0}));
  connect(eqAirTempVDI.TEqAir, preTemRoof.T) annotation(
    Line(points = {{-17, 84}, {-1, 84}, {-1, 89}, {17, 89}}, color = {0, 0, 127}));
  connect(theConRoof.Gc, hConRoof.y) annotation(
    Line(points = {{50, 69}, {62, 69}}, color = {0, 0, 127}));
  connect(eqAirTempVDI.TDryBul, eqAirTemp.TDryBul) annotation(
    Line(points = {{-31, 80}, {-96, 80}, {-96, -6}, {-43, -6}}, color = {0, 0, 127}));
  connect(eqAirTempVDI.TBlaSky, eqAirTemp.TBlaSky) annotation(
    Line(points = {{-31, 84}, {-98, 84}, {-98, -2}, {-43, -2}}, color = {0, 0, 127}));
  connect(eqAirTempVDI.HSol[1], weaBus.HGloHor) annotation(
    Line(points = {{-31, 88}, {-100, 88}, {-100, 17}, {-81, 17}}, color = {0, 0, 127}));
  connect(const1.y, eqAirTempVDI.sunblind[1]) annotation(
    Line(points = {{-10, 91}, {-24, 91}}, color = {0, 0, 127}));
  connect(realExpression_rad.y, perRad.Q_flow) annotation(
    Line(points = {{29, -63}, {50, -63}}, color = {0, 0, 127}));
  connect(realExpression_conv.y, perCon.Q_flow) annotation(
    Line(points = {{29, -77}, {50, -77}}, color = {0, 0, 127}));
  connect(realExpression_mach.y, macConv.Q_flow) annotation(
    Line(points = {{28.9, -91}, {49.9, -91}}, color = {0, 0, 127}));
  connect(setT.y, conHeaCoo.u_s) annotation(
    Line(points = {{6.8, -36}, {12.8, -36}}, color = {0, 0, 127}));
  connect(corGDouPan.solarRadWinTrans, thermalZoneFourElements.solRad) annotation(
    Line(points = {{-25, 52}, {40, 52}, {40, 36}, {55, 36}}, color = {0, 0, 127}));
  annotation(
    Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}), graphics = {Rectangle(origin = {-85, 55}, lineColor = {39, 39, 39}, fillColor = {186, 186, 186}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{15, -45}, {-15, 45}}), Rectangle(origin = {16, 80}, fillColor = {186, 186, 186}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{-60, 18}, {60, -18}}), Rectangle(origin = {-13, 22}, fillColor = {186, 186, 186}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{-55, 38}, {55, -38}}), Rectangle(origin = {-70, -54}, fillColor = {186, 186, 186}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{-24, 36}, {24, -36}}), Rectangle(origin = {70, -12}, fillColor = {159, 159, 159}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{-24, 12}, {24, -12}}), Rectangle(origin = {37, -36}, fillColor = {186, 186, 186}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{-51, 10}, {51, -10}}), Rectangle(origin = {36, -73}, fillColor = {186, 186, 186}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{-36, 25}, {36, -25}}), Text(origin = {-91, 9}, extent = {{-3, 3}, {3, -3}}, textString = "Outdoor", fontSize = 16), Text(origin = {-85, 13}, extent = {{-7, 5}, {7, -5}}, textString = "Outdoor", textStyle = {TextStyle.Bold}), Text(origin = {-85, 13}, extent = {{-7, 5}, {7, -5}}, textString = "Outdoor", textStyle = {TextStyle.Bold}), Text(origin = {-73, -87}, extent = {{-5, 3}, {5, -3}}, textString = "Walls", textStyle = {TextStyle.Bold}), Text(origin = {39, -50}, extent = {{-19, 12}, {19, -12}}, textString = "Person/room boundary", textStyle = {TextStyle.Bold}), Text(origin = {-13, -12}, extent = {{-7, 4}, {7, -4}}, textString = "Window", textStyle = {TextStyle.Bold}), Text(origin = {79, 43}, extent = {{-5, 3}, {5, -3}}, textString = "Office", textStyle = {TextStyle.Bold, TextStyle.Bold}), Text(origin = {12, 66}, extent = {{-4, 2}, {4, -2}}, textString = "Roof", textStyle = {TextStyle.Bold}), Text(origin = {55, -21}, extent = {{-5, 3}, {5, -3}}, textString = "Floor", textStyle = {TextStyle.Bold}), Text(origin = {-7, -43}, extent = {{-5, 3}, {5, -3}}, textString = "HVAC", textStyle = {TextStyle.Bold})}),

//day
//    experiment(StartTime = 0, StopTime = 86400, Tolerance = 1e-06, Interval = 60));
//year
    experiment(StartTime = 0, StopTime = 31536000, Tolerance = 1e-06, Interval = 600));
end PIR_gagge_smoosh;
