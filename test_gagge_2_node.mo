model test_gagge_2_node
  gagge_2_node human annotation(
    Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Cosine room_temp(amplitude = 2.5, f = 1/86400, offset = 20, phase = 3.141592653589793)  annotation(
    Placement(visible = true, transformation(origin = {-56, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(room_temp.y, human.ambient_temperature) annotation(
    Line(points = {{-44, 0}, {-8, 0}}, color = {0, 0, 127}));
  annotation(
    uses(Modelica(version = "4.0.0")),
  experiment(StartTime = 0, StopTime = 86400, Tolerance = 1e-6, Interval = 60));
end test_gagge_2_node;
