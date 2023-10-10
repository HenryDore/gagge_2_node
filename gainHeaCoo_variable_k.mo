block gainHeaCoo_variable_k "Output the product of a gain value with the input signal"
  Modelica.Blocks.Interfaces.RealInput u "Input signal connector" annotation(
    Placement(visible = true, transformation(origin = {0, 80}, extent = {{-140, -20}, {-100, 20}}, rotation = 0), iconTransformation(origin = {0, 80}, extent = {{-140, -20}, {-100, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput k "Gain value multiplied with input signal" annotation(
    Placement(visible = true, transformation(origin = {-120, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput y "Output signal connector" annotation(
    Placement(visible = true, transformation(origin = {0, 0}, extent = {{100, -10}, {120, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{100, -10}, {120, 10}}, rotation = 0)));
equation
  y = k*u;
  annotation(
    Documentation(info = "<html>
<p>
This block computes output <em>y</em> as
<em>product</em> of gain <em>k</em> with the
input <em>u</em>:
</p>
<blockquote><pre>
y = k * u;
</pre></blockquote>

</html>"),
    Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics = {Polygon(lineColor = {0, 0, 127}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, points = {{-100, -100}, {-100, 100}, {100, 0}, {-100, -100}}), Text(extent = {{-150, -140}, {150, -100}}, textString = ""), Text(textColor = {0, 0, 255}, extent = {{-150, 140}, {150, 100}}, textString = "%name"), Text(origin = {-83, -68}, extent = {{-27, 16}, {27, -16}}, textString = "k")}));
end gainHeaCoo_variable_k;
