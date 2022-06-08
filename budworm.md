# budworm model from Strogatz exercise 8.1.10

$$\begin{align}
		S' = r_S S ( 1 - \frac{S}{K_S}\frac{K_E}{E} )
		E' = r_E E ( 1 - \frac{E}{K_E} ) - \frac{B P}{S}
\end{align}$$

## Equilibria

model undefined along $S=0$ and $E=0$.

for simplicity, 
rescale $t$ by $r_S$
and let
$x = S/K_S$,
$y = E/K_E$,
$β = B P / K_E / K_S / r_E$,
$α = r_E/r_S$,
so that

$$\begin{align}
		x' = x ( 1 - \frac{x}{y} )
		y' = α ( y ( 1 - y ) - \frac{β}{x} )
\end{align}$$

Equilibria satisfy 
$x = y$ and 
$y^2(1-y) = β$

Should make a plot of this as β vs y to show equilibria
     _
\   / \
 \_/   \
--0----1--> y
  

Jacobian Matrix at an equilibrium is
$$\begin{pmatrix}
-x/y    &     -x^2    \\
αβ/x^2  &   α(1-2y) 
\end{pmatrix}
=
\begin{pmatrix}
-1    &     -y^2    \\
αβ/y^2  &   α(1-2y) 
\end{pmatrix}
$$

This has trace $-1 + α(1-2y)$
and determinant  $α( β - (1-2y) )$

interesting things happen when determinant passes through zero (β = 1-2y)
and when trace passes through zero and eigenvalue are complex



```julia
using BifurcationKit, Plots, Setfield, Parameters, LinearAlgebra, ForwardDiff, Revise
const BK = BifurcationKit

function budworm!(dz, z, p, t)
    @unpack Rs, Re, Ks, Ke, P, B = p
    S, E = z
    dz[1] = Rs*S*(1-(S*Ke)/(E*Ks))
    dz[2] = Re*E*(1-E/Ke)-(P*B)/S
    dz
end

budworm(z, p) = budworm!(similar(z), z, p, 0)

dbudworm(z,p) = ForwardDiff.jacobian(x -> budworm(x,p), z)
jet = BK.getJet(budworm, dbudworm)

opts = ContinuationPar(pMin = 0., pMax = 1000., detectBifurcation = 3)

params = (Rs = 100., Re = 95., Ks = 1.1, Ke = 1.2, P = 1., B = 100.)
#what are some reasonable initial parameters for this system? Find better ones than these.

br, = continuation(budworm, dbudworm, [-2.,-1.], params, (@lens _.B), opts; printSolution = (x,p) -> (S = x[1], E = x[2]))

plot(br)

```


<?xml version="1.0" encoding="utf-8"?>
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="600" height="400" viewBox="0 0 2400 1600">
<defs>
  <clipPath id="clip940">
    <rect x="0" y="0" width="2400" height="1600"/>
  </clipPath>
</defs>
<path clip-path="url(#clip940)" d="
M0 1600 L2400 1600 L2400 0 L0 0  Z
  " fill="#ffffff" fill-rule="evenodd" fill-opacity="1"/>
<defs>
  <clipPath id="clip941">
    <rect x="480" y="0" width="1681" height="1600"/>
  </clipPath>
</defs>
<path clip-path="url(#clip940)" d="
M247.945 1423.18 L2352.76 1423.18 L2352.76 47.2441 L247.945 47.2441  Z
  " fill="#ffffff" fill-rule="evenodd" fill-opacity="1"/>
<defs>
  <clipPath id="clip942">
    <rect x="247" y="47" width="2106" height="1377"/>
  </clipPath>
</defs>
<polyline clip-path="url(#clip942)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="
  600.526,1423.18 600.526,47.2441 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="
  1040.04,1423.18 1040.04,47.2441 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="
  1479.56,1423.18 1479.56,47.2441 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="
  1919.08,1423.18 1919.08,47.2441 
  "/>
<polyline clip-path="url(#clip940)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  247.945,1423.18 2352.76,1423.18 
  "/>
<polyline clip-path="url(#clip940)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  600.526,1423.18 600.526,1404.28 
  "/>
<polyline clip-path="url(#clip940)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1040.04,1423.18 1040.04,1404.28 
  "/>
<polyline clip-path="url(#clip940)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1479.56,1423.18 1479.56,1404.28 
  "/>
<polyline clip-path="url(#clip940)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1919.08,1423.18 1919.08,1404.28 
  "/>
<path clip-path="url(#clip940)" d="M560.932 1481.64 L568.57 1481.64 L568.57 1455.28 L560.26 1456.95 L560.26 1452.69 L568.524 1451.02 L573.2 1451.02 L573.2 1481.64 L580.839 1481.64 L580.839 1485.58 L560.932 1485.58 L560.932 1481.64 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M600.283 1454.1 Q596.672 1454.1 594.843 1457.66 Q593.038 1461.2 593.038 1468.33 Q593.038 1475.44 594.843 1479.01 Q596.672 1482.55 600.283 1482.55 Q603.918 1482.55 605.723 1479.01 Q607.552 1475.44 607.552 1468.33 Q607.552 1461.2 605.723 1457.66 Q603.918 1454.1 600.283 1454.1 M600.283 1450.39 Q606.093 1450.39 609.149 1455 Q612.228 1459.58 612.228 1468.33 Q612.228 1477.06 609.149 1481.67 Q606.093 1486.25 600.283 1486.25 Q594.473 1486.25 591.394 1481.67 Q588.339 1477.06 588.339 1468.33 Q588.339 1459.58 591.394 1455 Q594.473 1450.39 600.283 1450.39 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M624.473 1481.64 L640.792 1481.64 L640.792 1485.58 L618.848 1485.58 L618.848 1481.64 Q621.51 1478.89 626.093 1474.26 Q630.7 1469.61 631.88 1468.27 Q634.126 1465.74 635.005 1464.01 Q635.908 1462.25 635.908 1460.56 Q635.908 1457.8 633.964 1456.07 Q632.042 1454.33 628.941 1454.33 Q626.741 1454.33 624.288 1455.09 Q621.857 1455.86 619.079 1457.41 L619.079 1452.69 Q621.904 1451.55 624.357 1450.97 Q626.811 1450.39 628.848 1450.39 Q634.218 1450.39 637.413 1453.08 Q640.607 1455.77 640.607 1460.26 Q640.607 1462.39 639.797 1464.31 Q639.01 1466.2 636.903 1468.8 Q636.325 1469.47 633.223 1472.69 Q630.121 1475.88 624.473 1481.64 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M1000.15 1481.64 L1007.79 1481.64 L1007.79 1455.28 L999.476 1456.95 L999.476 1452.69 L1007.74 1451.02 L1012.42 1451.02 L1012.42 1481.64 L1020.06 1481.64 L1020.06 1485.58 L1000.15 1485.58 L1000.15 1481.64 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M1039.5 1454.1 Q1035.89 1454.1 1034.06 1457.66 Q1032.25 1461.2 1032.25 1468.33 Q1032.25 1475.44 1034.06 1479.01 Q1035.89 1482.55 1039.5 1482.55 Q1043.13 1482.55 1044.94 1479.01 Q1046.77 1475.44 1046.77 1468.33 Q1046.77 1461.2 1044.94 1457.66 Q1043.13 1454.1 1039.5 1454.1 M1039.5 1450.39 Q1045.31 1450.39 1048.37 1455 Q1051.44 1459.58 1051.44 1468.33 Q1051.44 1477.06 1048.37 1481.67 Q1045.31 1486.25 1039.5 1486.25 Q1033.69 1486.25 1030.61 1481.67 Q1027.56 1477.06 1027.56 1468.33 Q1027.56 1459.58 1030.61 1455 Q1033.69 1450.39 1039.5 1450.39 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M1059.71 1451.02 L1078.06 1451.02 L1078.06 1454.96 L1063.99 1454.96 L1063.99 1463.43 Q1065.01 1463.08 1066.03 1462.92 Q1067.05 1462.73 1068.06 1462.73 Q1073.85 1462.73 1077.23 1465.9 Q1080.61 1469.08 1080.61 1474.49 Q1080.61 1480.07 1077.14 1483.17 Q1073.67 1486.25 1067.35 1486.25 Q1065.17 1486.25 1062.9 1485.88 Q1060.66 1485.51 1058.25 1484.77 L1058.25 1480.07 Q1060.33 1481.2 1062.55 1481.76 Q1064.78 1482.32 1067.25 1482.32 Q1071.26 1482.32 1073.6 1480.21 Q1075.93 1478.1 1075.93 1474.49 Q1075.93 1470.88 1073.6 1468.77 Q1071.26 1466.67 1067.25 1466.67 Q1065.38 1466.67 1063.5 1467.08 Q1061.65 1467.5 1059.71 1468.38 L1059.71 1451.02 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M1439.21 1481.64 L1446.85 1481.64 L1446.85 1455.28 L1438.54 1456.95 L1438.54 1452.69 L1446.81 1451.02 L1451.48 1451.02 L1451.48 1481.64 L1459.12 1481.64 L1459.12 1485.58 L1439.21 1485.58 L1439.21 1481.64 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M1478.57 1454.1 Q1474.95 1454.1 1473.13 1457.66 Q1471.32 1461.2 1471.32 1468.33 Q1471.32 1475.44 1473.13 1479.01 Q1474.95 1482.55 1478.57 1482.55 Q1482.2 1482.55 1484 1479.01 Q1485.83 1475.44 1485.83 1468.33 Q1485.83 1461.2 1484 1457.66 Q1482.2 1454.1 1478.57 1454.1 M1478.57 1450.39 Q1484.38 1450.39 1487.43 1455 Q1490.51 1459.58 1490.51 1468.33 Q1490.51 1477.06 1487.43 1481.67 Q1484.38 1486.25 1478.57 1486.25 Q1472.75 1486.25 1469.68 1481.67 Q1466.62 1477.06 1466.62 1468.33 Q1466.62 1459.58 1469.68 1455 Q1472.75 1450.39 1478.57 1450.39 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M1508.73 1469.17 Q1505.39 1469.17 1503.47 1470.95 Q1501.57 1472.73 1501.57 1475.86 Q1501.57 1478.98 1503.47 1480.77 Q1505.39 1482.55 1508.73 1482.55 Q1512.06 1482.55 1513.98 1480.77 Q1515.9 1478.96 1515.9 1475.86 Q1515.9 1472.73 1513.98 1470.95 Q1512.08 1469.17 1508.73 1469.17 M1504.05 1467.18 Q1501.04 1466.44 1499.35 1464.38 Q1497.69 1462.32 1497.69 1459.35 Q1497.69 1455.21 1500.63 1452.8 Q1503.59 1450.39 1508.73 1450.39 Q1513.89 1450.39 1516.83 1452.8 Q1519.77 1455.21 1519.77 1459.35 Q1519.77 1462.32 1518.08 1464.38 Q1516.41 1466.44 1513.43 1467.18 Q1516.81 1467.96 1518.68 1470.26 Q1520.58 1472.55 1520.58 1475.86 Q1520.58 1480.88 1517.5 1483.57 Q1514.44 1486.25 1508.73 1486.25 Q1503.01 1486.25 1499.93 1483.57 Q1496.88 1480.88 1496.88 1475.86 Q1496.88 1472.55 1498.77 1470.26 Q1500.67 1467.96 1504.05 1467.18 M1502.34 1459.79 Q1502.34 1462.48 1504 1463.98 Q1505.69 1465.49 1508.73 1465.49 Q1511.74 1465.49 1513.43 1463.98 Q1515.14 1462.48 1515.14 1459.79 Q1515.14 1457.11 1513.43 1455.6 Q1511.74 1454.1 1508.73 1454.1 Q1505.69 1454.1 1504 1455.6 Q1502.34 1457.11 1502.34 1459.79 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M1879.3 1481.64 L1886.94 1481.64 L1886.94 1455.28 L1878.63 1456.95 L1878.63 1452.69 L1886.89 1451.02 L1891.57 1451.02 L1891.57 1481.64 L1899.2 1481.64 L1899.2 1485.58 L1879.3 1485.58 L1879.3 1481.64 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M1909.46 1481.64 L1917.1 1481.64 L1917.1 1455.28 L1908.79 1456.95 L1908.79 1452.69 L1917.05 1451.02 L1921.73 1451.02 L1921.73 1481.64 L1929.37 1481.64 L1929.37 1485.58 L1909.46 1485.58 L1909.46 1481.64 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M1939.62 1481.64 L1947.26 1481.64 L1947.26 1455.28 L1938.95 1456.95 L1938.95 1452.69 L1947.21 1451.02 L1951.89 1451.02 L1951.89 1481.64 L1959.53 1481.64 L1959.53 1485.58 L1939.62 1485.58 L1939.62 1481.64 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M1289.93 1545.35 L1289.93 1562.76 L1300.24 1562.76 Q1305.43 1562.76 1307.91 1560.63 Q1310.42 1558.46 1310.42 1554.04 Q1310.42 1549.58 1307.91 1547.48 Q1305.43 1545.35 1300.24 1545.35 L1289.93 1545.35 M1289.93 1525.81 L1289.93 1540.13 L1299.44 1540.13 Q1304.15 1540.13 1306.45 1538.38 Q1308.77 1536.6 1308.77 1532.97 Q1308.77 1529.37 1306.45 1527.59 Q1304.15 1525.81 1299.44 1525.81 L1289.93 1525.81 M1283.5 1520.52 L1299.92 1520.52 Q1307.27 1520.52 1311.25 1523.58 Q1315.23 1526.63 1315.23 1532.27 Q1315.23 1536.63 1313.19 1539.21 Q1311.16 1541.79 1307.21 1542.42 Q1311.95 1543.44 1314.56 1546.69 Q1317.2 1549.9 1317.2 1554.74 Q1317.2 1561.11 1312.87 1564.57 Q1308.55 1568.04 1300.56 1568.04 L1283.5 1568.04 L1283.5 1520.52 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><polyline clip-path="url(#clip942)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="
  247.945,1364.42 2352.76,1364.42 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="
  247.945,1149.89 2352.76,1149.89 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="
  247.945,935.363 2352.76,935.363 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="
  247.945,720.836 2352.76,720.836 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="
  247.945,506.309 2352.76,506.309 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="
  247.945,291.781 2352.76,291.781 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:2; stroke-opacity:0.1; fill:none" points="
  247.945,77.254 2352.76,77.254 
  "/>
<polyline clip-path="url(#clip940)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  247.945,1423.18 247.945,47.2441 
  "/>
<polyline clip-path="url(#clip940)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  247.945,1364.42 266.842,1364.42 
  "/>
<polyline clip-path="url(#clip940)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  247.945,1149.89 266.842,1149.89 
  "/>
<polyline clip-path="url(#clip940)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  247.945,935.363 266.842,935.363 
  "/>
<polyline clip-path="url(#clip940)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  247.945,720.836 266.842,720.836 
  "/>
<polyline clip-path="url(#clip940)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  247.945,506.309 266.842,506.309 
  "/>
<polyline clip-path="url(#clip940)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  247.945,291.781 266.842,291.781 
  "/>
<polyline clip-path="url(#clip940)" style="stroke:#000000; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  247.945,77.254 266.842,77.254 
  "/>
<path clip-path="url(#clip940)" d="M117.015 1377.76 L124.654 1377.76 L124.654 1351.4 L116.343 1353.06 L116.343 1348.8 L124.607 1347.14 L129.283 1347.14 L129.283 1377.76 L136.922 1377.76 L136.922 1381.7 L117.015 1381.7 L117.015 1377.76 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M146.366 1375.82 L151.251 1375.82 L151.251 1381.7 L146.366 1381.7 L146.366 1375.82 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M162.246 1377.76 L169.885 1377.76 L169.885 1351.4 L161.575 1353.06 L161.575 1348.8 L169.839 1347.14 L174.514 1347.14 L174.514 1377.76 L182.153 1377.76 L182.153 1381.7 L162.246 1381.7 L162.246 1377.76 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M195.625 1377.76 L211.945 1377.76 L211.945 1381.7 L190 1381.7 L190 1377.76 Q192.663 1375.01 197.246 1370.38 Q201.852 1365.73 203.033 1364.38 Q205.278 1361.86 206.158 1360.12 Q207.061 1358.36 207.061 1356.67 Q207.061 1353.92 205.116 1352.18 Q203.195 1350.45 200.093 1350.45 Q197.894 1350.45 195.44 1351.21 Q193.01 1351.98 190.232 1353.53 L190.232 1348.8 Q193.056 1347.67 195.51 1347.09 Q197.963 1346.51 200 1346.51 Q205.371 1346.51 208.565 1349.2 Q211.76 1351.88 211.76 1356.37 Q211.76 1358.5 210.949 1360.42 Q210.162 1362.32 208.056 1364.92 Q207.477 1365.59 204.375 1368.8 Q201.274 1372 195.625 1377.76 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M116.066 1163.24 L123.705 1163.24 L123.705 1136.87 L115.394 1138.54 L115.394 1134.28 L123.658 1132.61 L128.334 1132.61 L128.334 1163.24 L135.973 1163.24 L135.973 1167.17 L116.066 1167.17 L116.066 1163.24 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M145.417 1161.29 L150.302 1161.29 L150.302 1167.17 L145.417 1167.17 L145.417 1161.29 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M161.297 1163.24 L168.936 1163.24 L168.936 1136.87 L160.626 1138.54 L160.626 1134.28 L168.889 1132.61 L173.565 1132.61 L173.565 1163.24 L181.204 1163.24 L181.204 1167.17 L161.297 1167.17 L161.297 1163.24 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M204.815 1148.54 Q208.172 1149.25 210.047 1151.52 Q211.945 1153.79 211.945 1157.12 Q211.945 1162.24 208.426 1165.04 Q204.908 1167.84 198.426 1167.84 Q196.25 1167.84 193.936 1167.4 Q191.644 1166.99 189.19 1166.13 L189.19 1161.61 Q191.135 1162.75 193.45 1163.33 Q195.764 1163.91 198.287 1163.91 Q202.686 1163.91 204.977 1162.17 Q207.292 1160.43 207.292 1157.12 Q207.292 1154.07 205.139 1152.36 Q203.01 1150.62 199.19 1150.62 L195.163 1150.62 L195.163 1146.78 L199.375 1146.78 Q202.824 1146.78 204.653 1145.41 Q206.482 1144.02 206.482 1141.43 Q206.482 1138.77 204.584 1137.36 Q202.709 1135.92 199.19 1135.92 Q197.269 1135.92 195.07 1136.34 Q192.871 1136.75 190.232 1137.63 L190.232 1133.47 Q192.894 1132.73 195.209 1132.36 Q197.547 1131.99 199.607 1131.99 Q204.931 1131.99 208.033 1134.42 Q211.135 1136.82 211.135 1140.94 Q211.135 1143.81 209.491 1145.8 Q207.848 1147.77 204.815 1148.54 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M114.931 948.708 L122.57 948.708 L122.57 922.342 L114.26 924.009 L114.26 919.75 L122.524 918.083 L127.2 918.083 L127.2 948.708 L134.839 948.708 L134.839 952.643 L114.931 952.643 L114.931 948.708 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M144.283 946.764 L149.167 946.764 L149.167 952.643 L144.283 952.643 L144.283 946.764 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M160.163 948.708 L167.802 948.708 L167.802 922.342 L159.491 924.009 L159.491 919.75 L167.755 918.083 L172.431 918.083 L172.431 948.708 L180.07 948.708 L180.07 952.643 L160.163 952.643 L160.163 948.708 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M202.362 922.157 L190.556 940.606 L202.362 940.606 L202.362 922.157 M201.135 918.083 L207.014 918.083 L207.014 940.606 L211.945 940.606 L211.945 944.495 L207.014 944.495 L207.014 952.643 L202.362 952.643 L202.362 944.495 L186.76 944.495 L186.76 939.981 L201.135 918.083 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M116.413 734.181 L124.052 734.181 L124.052 707.815 L115.742 709.482 L115.742 705.223 L124.005 703.556 L128.681 703.556 L128.681 734.181 L136.32 734.181 L136.32 738.116 L116.413 738.116 L116.413 734.181 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M145.765 732.236 L150.649 732.236 L150.649 738.116 L145.765 738.116 L145.765 732.236 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M161.644 734.181 L169.283 734.181 L169.283 707.815 L160.973 709.482 L160.973 705.223 L169.237 703.556 L173.913 703.556 L173.913 734.181 L181.551 734.181 L181.551 738.116 L161.644 738.116 L161.644 734.181 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M191.042 703.556 L209.399 703.556 L209.399 707.491 L195.325 707.491 L195.325 715.963 Q196.343 715.616 197.362 715.454 Q198.38 715.269 199.399 715.269 Q205.186 715.269 208.565 718.44 Q211.945 721.611 211.945 727.028 Q211.945 732.607 208.473 735.709 Q205 738.787 198.681 738.787 Q196.505 738.787 194.237 738.417 Q191.991 738.046 189.584 737.306 L189.584 732.607 Q191.667 733.741 193.889 734.296 Q196.112 734.852 198.588 734.852 Q202.593 734.852 204.931 732.746 Q207.269 730.639 207.269 727.028 Q207.269 723.417 204.931 721.31 Q202.593 719.204 198.588 719.204 Q196.713 719.204 194.838 719.621 Q192.987 720.037 191.042 720.917 L191.042 703.556 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M115.256 519.653 L122.894 519.653 L122.894 493.288 L114.584 494.955 L114.584 490.695 L122.848 489.029 L127.524 489.029 L127.524 519.653 L135.163 519.653 L135.163 523.589 L115.256 523.589 L115.256 519.653 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M144.607 517.709 L149.491 517.709 L149.491 523.589 L144.607 523.589 L144.607 517.709 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M160.487 519.653 L168.126 519.653 L168.126 493.288 L159.815 494.955 L159.815 490.695 L168.079 489.029 L172.755 489.029 L172.755 519.653 L180.394 519.653 L180.394 523.589 L160.487 523.589 L160.487 519.653 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M200.417 504.445 Q197.269 504.445 195.417 506.598 Q193.588 508.751 193.588 512.501 Q193.588 516.228 195.417 518.403 Q197.269 520.556 200.417 520.556 Q203.565 520.556 205.394 518.403 Q207.246 516.228 207.246 512.501 Q207.246 508.751 205.394 506.598 Q203.565 504.445 200.417 504.445 M209.699 489.792 L209.699 494.052 Q207.94 493.218 206.135 492.779 Q204.352 492.339 202.593 492.339 Q197.963 492.339 195.51 495.464 Q193.079 498.589 192.732 504.908 Q194.098 502.894 196.158 501.829 Q198.218 500.742 200.695 500.742 Q205.903 500.742 208.912 503.913 Q211.945 507.061 211.945 512.501 Q211.945 517.825 208.797 521.042 Q205.649 524.26 200.417 524.26 Q194.422 524.26 191.25 519.677 Q188.079 515.07 188.079 506.343 Q188.079 498.149 191.968 493.288 Q195.857 488.404 202.408 488.404 Q204.167 488.404 205.949 488.751 Q207.755 489.098 209.699 489.792 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M116.32 305.126 L123.959 305.126 L123.959 278.761 L115.649 280.427 L115.649 276.168 L123.913 274.501 L128.589 274.501 L128.589 305.126 L136.228 305.126 L136.228 309.061 L116.32 309.061 L116.32 305.126 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M145.672 303.182 L150.556 303.182 L150.556 309.061 L145.672 309.061 L145.672 303.182 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M161.552 305.126 L169.19 305.126 L169.19 278.761 L160.88 280.427 L160.88 276.168 L169.144 274.501 L173.82 274.501 L173.82 305.126 L181.459 305.126 L181.459 309.061 L161.552 309.061 L161.552 305.126 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M189.723 274.501 L211.945 274.501 L211.945 276.492 L199.399 309.061 L194.514 309.061 L206.32 278.436 L189.723 278.436 L189.723 274.501 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M115.51 90.5988 L123.149 90.5988 L123.149 64.2332 L114.839 65.8999 L114.839 61.6407 L123.103 59.974 L127.779 59.974 L127.779 90.5988 L135.417 90.5988 L135.417 94.534 L115.51 94.534 L115.51 90.5988 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M144.862 88.6544 L149.746 88.6544 L149.746 94.534 L144.862 94.534 L144.862 88.6544 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M160.741 90.5988 L168.38 90.5988 L168.38 64.2332 L160.07 65.8999 L160.07 61.6407 L168.334 59.974 L173.01 59.974 L173.01 90.5988 L180.649 90.5988 L180.649 94.534 L160.741 94.534 L160.741 90.5988 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M200.093 78.1221 Q196.76 78.1221 194.838 79.9045 Q192.94 81.6869 192.94 84.8118 Q192.94 87.9368 194.838 89.7192 Q196.76 91.5016 200.093 91.5016 Q203.426 91.5016 205.348 89.7192 Q207.269 87.9137 207.269 84.8118 Q207.269 81.6869 205.348 79.9045 Q203.449 78.1221 200.093 78.1221 M195.417 76.1313 Q192.408 75.3906 190.718 73.3304 Q189.051 71.2702 189.051 68.3073 Q189.051 64.1638 191.991 61.7564 Q194.954 59.349 200.093 59.349 Q205.255 59.349 208.195 61.7564 Q211.135 64.1638 211.135 68.3073 Q211.135 71.2702 209.445 73.3304 Q207.778 75.3906 204.792 76.1313 Q208.172 76.9184 210.047 79.21 Q211.945 81.5017 211.945 84.8118 Q211.945 89.835 208.866 92.5201 Q205.811 95.2053 200.093 95.2053 Q194.375 95.2053 191.297 92.5201 Q188.241 89.835 188.241 84.8118 Q188.241 81.5017 190.139 79.21 Q192.038 76.9184 195.417 76.1313 M193.704 68.7471 Q193.704 71.4323 195.371 72.9369 Q197.061 74.4415 200.093 74.4415 Q203.102 74.4415 204.792 72.9369 Q206.505 71.4323 206.505 68.7471 Q206.505 66.0619 204.792 64.5573 Q203.102 63.0527 200.093 63.0527 Q197.061 63.0527 195.371 64.5573 Q193.704 66.0619 193.704 68.7471 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><path clip-path="url(#clip940)" d="M28.3562 718.597 L45.7028 731.488 L64.0042 717.929 L64.0042 724.836 L49.9996 735.212 L64.0042 745.588 L64.0042 752.495 L45.3526 738.649 L28.3562 751.317 L28.3562 744.41 L41.0558 734.957 L28.3562 725.504 L28.3562 718.597 Z" fill="#000000" fill-rule="evenodd" fill-opacity="1" /><polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  307.515,1384.24 309.587,1382.83 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  309.587,1382.83 312.613,1380.77 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  312.613,1380.77 317.035,1377.76 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  317.035,1377.76 323.493,1373.37 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  323.493,1373.37 332.928,1366.95 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  332.928,1366.95 346.71,1357.58 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  346.71,1357.58 366.325,1344.26 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  366.325,1344.26 387.044,1330.2 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  387.044,1330.2 407.763,1316.15 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  407.763,1316.15 428.482,1302.11 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  428.482,1302.11 449.201,1288.09 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  449.201,1288.09 469.92,1274.07 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  469.92,1274.07 490.639,1260.07 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  490.639,1260.07 511.358,1246.09 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  511.358,1246.09 532.077,1232.11 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  532.077,1232.11 552.796,1218.15 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  552.796,1218.15 573.515,1204.2 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  573.515,1204.2 594.233,1190.26 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  594.233,1190.26 614.952,1176.33 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  614.952,1176.33 635.671,1162.42 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  635.671,1162.42 656.39,1148.52 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  656.39,1148.52 677.109,1134.62 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  677.109,1134.62 697.828,1120.75 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  697.828,1120.75 718.547,1106.88 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  718.547,1106.88 739.266,1093.03 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  739.266,1093.03 759.985,1079.18 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  759.985,1079.18 780.704,1065.35 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  780.704,1065.35 801.423,1051.53 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  801.423,1051.53 822.142,1037.73 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  822.142,1037.73 842.861,1023.93 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  842.861,1023.93 863.579,1010.15 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  863.579,1010.15 884.298,996.375 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  884.298,996.375 905.017,982.615 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  905.017,982.615 925.736,968.867 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  925.736,968.867 946.455,955.131 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  946.455,955.131 967.174,941.406 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  967.174,941.406 987.893,927.693 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  987.893,927.693 1008.61,913.991 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1008.61,913.991 1029.33,900.301 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1029.33,900.301 1050.05,886.623 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1050.05,886.623 1070.77,872.956 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1070.77,872.956 1091.49,859.3 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1091.49,859.3 1112.21,845.657 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1112.21,845.657 1132.93,832.024 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1132.93,832.024 1153.64,818.403 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1153.64,818.403 1174.36,804.794 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1174.36,804.794 1195.08,791.195 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1195.08,791.195 1215.8,777.609 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1215.8,777.609 1236.52,764.033 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1236.52,764.033 1257.24,750.469 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1257.24,750.469 1277.96,736.916 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1277.96,736.916 1298.68,723.375 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1298.68,723.375 1319.4,709.844 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1319.4,709.844 1340.11,696.325 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1340.11,696.325 1360.83,682.817 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1360.83,682.817 1381.55,669.321 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1381.55,669.321 1402.27,655.835 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1402.27,655.835 1422.99,642.361 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1422.99,642.361 1443.71,628.897 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1443.71,628.897 1464.43,615.445 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1464.43,615.445 1485.15,602.004 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1485.15,602.004 1505.87,588.574 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1505.87,588.574 1526.59,575.155 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1526.59,575.155 1547.3,561.746 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1547.3,561.746 1568.02,548.349 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1568.02,548.349 1588.74,534.963 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1588.74,534.963 1609.46,521.588 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1609.46,521.588 1630.18,508.223 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1630.18,508.223 1650.9,494.87 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1650.9,494.87 1671.62,481.527 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1671.62,481.527 1692.34,468.195 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1692.34,468.195 1713.06,454.874 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1713.06,454.874 1733.77,441.564 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1733.77,441.564 1754.49,428.265 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1754.49,428.265 1775.21,414.976 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1775.21,414.976 1795.93,401.698 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1795.93,401.698 1816.65,388.431 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1816.65,388.431 1837.37,375.174 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1837.37,375.174 1858.09,361.928 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1858.09,361.928 1878.81,348.693 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1878.81,348.693 1899.53,335.468 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1899.53,335.468 1920.24,322.254 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1920.24,322.254 1940.96,309.05 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1940.96,309.05 1961.68,295.857 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1961.68,295.857 1982.4,282.675 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  1982.4,282.675 2003.12,269.503 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  2003.12,269.503 2023.84,256.341 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  2023.84,256.341 2044.56,243.19 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  2044.56,243.19 2065.28,230.049 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  2065.28,230.049 2086,216.919 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  2086,216.919 2106.72,203.799 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  2106.72,203.799 2127.43,190.69 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  2127.43,190.69 2148.15,177.591 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  2148.15,177.591 2168.87,164.502 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  2168.87,164.502 2189.59,151.424 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  2189.59,151.424 2210.31,138.356 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  2210.31,138.356 2231.03,125.298 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  2231.03,125.298 2251.75,112.25 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  2251.75,112.25 2272.47,99.2128 
  "/>
<polyline clip-path="url(#clip942)" style="stroke:#009af9; stroke-linecap:butt; stroke-linejoin:round; stroke-width:4; stroke-opacity:1; fill:none" points="
  2272.47,99.2128 2293.19,86.1857 
  "/>
</svg>


