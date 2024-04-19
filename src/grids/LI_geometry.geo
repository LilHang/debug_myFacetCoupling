// length scale
lc = 0.01;
lcf = 0.004;

// ------outer circle points-------
// lower point
Point(1) = {0.215, 0, 0, lc};
// cneter point
Point(2) = {0.215, 0.215, 0, lc};
// upper point
Point(3) = {0.215, 0.43, 0, lc};

// ------inner circle points-------
// lower point
Point(4) = {0.215, 0.17, 0, lcf};
// upper point
Point(5) = {0.215, 0.26, 0, lcf};

// outer circle
Circle(1) = {1, 2, 3};

// inner circle
Circle(2) = {4, 2, 5};

// left boundaries
Line(3) = {3, 5};
Line(4) = {4, 1};


Curve Loop(1) = {1, 3, -2, 4};
// surface
Plane Surface(1) = {1};

// make elements conforming to lowdim
Curve {2} In Surface {1};

// physical groups
Physical Surface(1) = {1};
Physical Curve(1) = {2};


