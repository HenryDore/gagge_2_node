model simple_gagge
  gagge_2_node gagge_2_node1 annotation(
    Placement(transformation(origin = {0, 16}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant const(k = 273.15 + 15)  annotation(
    Placement(transformation(origin = {-42, 20}, extent = {{-10, -10}, {10, 10}})));
  gagge_2_node gagge_2_node11 annotation(
    Placement(transformation(origin = {0, -28}, extent = {{-10, -10}, {10, 10}})));
  gagge_2_node gagge_2_node12 annotation(
    Placement(transformation(origin = {0, -6}, extent = {{-10, -10}, {10, 10}})));
  gagge_2_node gagge_2_node13 annotation(
    Placement(transformation(origin = {0, 42}, extent = {{-10, -10}, {10, 10}})));
equation
  connect(const.y, gagge_2_node1.ambient_temperature) annotation(
    Line(points = {{-30, 20}, {-8, 20}}, color = {0, 0, 127}));
  connect(const.y, gagge_2_node13.ambient_temperature) annotation(
    Line(points = {{-30, 20}, {-20, 20}, {-20, 46}, {-8, 46}}, color = {0, 0, 127}));
  connect(const.y, gagge_2_node12.ambient_temperature) annotation(
    Line(points = {{-30, 20}, {-20, 20}, {-20, -2}, {-8, -2}}, color = {0, 0, 127}));
  connect(const.y, gagge_2_node11.ambient_temperature) annotation(
    Line(points = {{-30, 20}, {-20, 20}, {-20, -24}, {-8, -24}}, color = {0, 0, 127}));

annotation(
    uses(Modelica(version = "4.0.0")),
  experiment(StartTime = 0, StopTime = 86400, Tolerance = 1e-06, Interval = 60));
end simple_gagge;
