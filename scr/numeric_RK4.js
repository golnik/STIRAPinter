// Extends numeric.js with a fixed-step classical RK4 integrator that works
// on states made of Complex (Complex.js) values, so the right-hand side
// f(t,y) can be written directly in terms of complex amplitudes.
//
// Requires numeric.js, Complex.js and complex_exteded.js to be loaded first.

// Solution object: fixed grid of points plus piecewise-linear interpolation.
numeric.RK4Complex = function RK4Complex(x, y) {
  this.x = x;
  this.y = y;
};

numeric.RK4Complex.prototype.at = function at(xi) {
  if (typeof xi !== "number") {
    var n = xi.length, ret = new Array(n);
    for (var i = 0; i < n; i++) ret[i] = this.at(xi[i]);
    return ret;
  }
  var xs = this.x, ys = this.y;
  var floor = Math.floor;
  var i = 0, j = xs.length - 1, k;
  while (j - i > 1) {
    k = floor(0.5 * (i + j));
    if (xs[k] <= xi) i = k; else j = k;
  }
  var t = (xi - xs[i]) / (xs[i + 1] - xs[i]);
  return cAdd(cScale(ys[i], 1 - t), cScale(ys[i + 1], t));
};

// Classical fixed-step RK4 integrator for y' = f(t,y) on [x0,x1] using n
// equal steps (default 100). y is a Complex value or a (nested) array of
// Complex values; f(t,y) must return a value of the same shape as y.
numeric.rk4Complex = function rk4Complex(x0, x1, y0, f, n) {
  if (typeof n === "undefined") n = 100;
  var h = (x1 - x0) / n;
  var xs = [x0], ys = [y0];
  var x = x0, y = y0;
  var add = cAdd, scale = cScale;

  for (var step = 0; step < n; step++) {
    var k1 = f(x, y);
    var k2 = f(x + 0.5 * h, add(y, scale(k1, 0.5 * h)));
    var k3 = f(x + 0.5 * h, add(y, scale(k2, 0.5 * h)));
    var k4 = f(x + h, add(y, scale(k3, h)));

    y = add(y, scale(add(add(k1, scale(k2, 2)), add(scale(k3, 2), k4)), h / 6));
    x = x0 + (step + 1) * h;

    xs.push(x);
    ys.push(y);
  }

  return new numeric.RK4Complex(xs, ys);
};
