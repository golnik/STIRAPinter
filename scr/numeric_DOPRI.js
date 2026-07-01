// Extends numeric.js with a Dormand-Prince RK45 integrator that works on
// states made of Complex (Complex.js) values, so the right-hand side f(t,y)
// can be written directly in terms of complex amplitudes instead of having
// to split every equation into real/imaginary components by hand.
//
// Requires numeric.js, Complex.js and complex_exteded.js to be loaded first.

// Solution object with the same dense-output interpolation as numeric.Dopri,
// rewritten in terms of cAdd/cSub/cScale so it works on Complex states.
numeric.DopriComplex = function DopriComplex(x, y, f, ymid, iterations, msg) {
  this.x = x;
  this.y = y;
  this.f = f;
  this.ymid = ymid;
  this.iterations = iterations;
  this.message = msg;
};

numeric.DopriComplex.prototype._at = function _at(xi, j) {
  function sqr(x) { return x * x; }
  var xs = this.x, ys = this.y, k1 = this.f, ymid = this.ymid;
  var x0 = xs[j], x1 = xs[j + 1], y0 = ys[j], y1 = ys[j + 1];
  var h = x1 - x0;
  var xh = x0 + 0.5 * h;
  var yh = ymid[j];
  var p = cSub(k1[j], cScale(y0, 1 / (x0 - xh) + 2 / (x0 - x1)));
  var q = cSub(k1[j + 1], cScale(y1, 1 / (x1 - xh) + 2 / (x1 - x0)));
  var w = [
    sqr(xi - x1) * (xi - xh) / sqr(x0 - x1) / (x0 - xh),
    sqr(xi - x0) * sqr(xi - x1) / sqr(x0 - xh) / sqr(x1 - xh),
    sqr(xi - x0) * (xi - xh) / sqr(x1 - x0) / (x1 - xh),
    (xi - x0) * sqr(xi - x1) * (xi - xh) / sqr(x0 - x1) / (x0 - xh),
    (xi - x1) * sqr(xi - x0) * (xi - xh) / sqr(x0 - x1) / (x1 - xh)
  ];
  return cAdd(
    cAdd(
      cAdd(
        cAdd(cScale(y0, w[0]), cScale(yh, w[1])),
        cScale(y1, w[2])),
      cScale(p, w[3])),
    cScale(q, w[4]));
};

numeric.DopriComplex.prototype.at = function at(x) {
  var floor = Math.floor;
  if (typeof x !== "number") {
    var n = x.length, ret = new Array(n);
    for (var i = 0; i < n; i++) ret[i] = this.at(x[i]);
    return ret;
  }
  var xs = this.x;
  var i = 0, j = xs.length - 1, k;
  while (j - i > 1) {
    k = floor(0.5 * (i + j));
    if (xs[k] <= x) i = k; else j = k;
  }
  return this._at(x, i);
};

// Dormand-Prince 5(4) adaptive integrator for y' = f(t,y), where y is a
// Complex value or a (nested) array of Complex values. f(t,y) must return
// a value of the same shape as y. Step-size control uses the magnitude
// (.abs()) of the embedded error estimate, since complex numbers have no
// natural ordering; event detection (root finding) from numeric.dopri is
// therefore not supported here.
numeric.dopriComplex = function dopriComplex(x0, x1, y0, f, tol, maxit) {
  if (typeof tol === "undefined") tol = 1e-6;
  if (typeof maxit === "undefined") maxit = 1000;

  var xs = [x0], ys = [y0], k1 = [f(x0, y0)], k2, k3, k4, k5, k6, k7, ymid = [];

  var A2 = 1 / 5;
  var A3 = [3 / 40, 9 / 40];
  var A4 = [44 / 45, -56 / 15, 32 / 9];
  var A5 = [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729];
  var A6 = [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656];
  var b = [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84];
  var bm = [0.5 * 6025192743 / 30085553152,
            0,
            0.5 * 51252292925 / 65400821598,
            0.5 * -2691868925 / 45128329728,
            0.5 * 187940372067 / 1594534317056,
            0.5 * -1776094331 / 19743644256,
            0.5 * 11237099 / 235043384];
  var c = [1 / 5, 3 / 10, 4 / 5, 8 / 9, 1, 1];
  var e = [-71 / 57600, 0, 71 / 16695, -71 / 1920, 17253 / 339200, -22 / 525, 1 / 40];

  var i = 0, er, erinf, it = 0;
  var h = (x1 - x0) / 10;
  var y1;
  var add = cAdd, scale = cScale, norminf = cNorminf;
  var max = Math.max, min = Math.min, pow = Math.pow;

  var ret = new numeric.DopriComplex(xs, ys, k1, ymid, -1, "");

  while (x0 < x1 && it < maxit) {
    ++it;
    if (x0 + h > x1) h = x1 - x0;

    k2 = f(x0 + c[0] * h, add(y0, scale(k1[i], A2 * h)));
    k3 = f(x0 + c[1] * h, add(add(y0, scale(k1[i], A3[0] * h)), scale(k2, A3[1] * h)));
    k4 = f(x0 + c[2] * h, add(add(add(y0, scale(k1[i], A4[0] * h)), scale(k2, A4[1] * h)), scale(k3, A4[2] * h)));
    k5 = f(x0 + c[3] * h, add(add(add(add(y0, scale(k1[i], A5[0] * h)), scale(k2, A5[1] * h)), scale(k3, A5[2] * h)), scale(k4, A5[3] * h)));
    k6 = f(x0 + c[4] * h, add(add(add(add(add(y0, scale(k1[i], A6[0] * h)), scale(k2, A6[1] * h)), scale(k3, A6[2] * h)), scale(k4, A6[3] * h)), scale(k5, A6[4] * h)));

    y1 = add(add(add(add(add(y0, scale(k1[i], h * b[0])), scale(k3, h * b[2])), scale(k4, h * b[3])), scale(k5, h * b[4])), scale(k6, h * b[5]));
    k7 = f(x0 + h, y1);

    er = add(add(add(add(add(scale(k1[i], h * e[0]), scale(k3, h * e[2])), scale(k4, h * e[3])), scale(k5, h * e[4])), scale(k6, h * e[5])), scale(k7, h * e[6]));
    erinf = norminf(er);

    if (erinf > tol) { // reject
      h = 0.2 * h * pow(tol / erinf, 0.25);
      if (x0 + h === x0) {
        ret.message = "Step size became too small";
        break;
      }
      continue;
    }

    ymid[i] = add(add(add(add(add(add(y0,
              scale(k1[i], h * bm[0])),
              scale(k3, h * bm[2])),
              scale(k4, h * bm[3])),
              scale(k5, h * bm[4])),
              scale(k6, h * bm[5])),
              scale(k7, h * bm[6]));

    ++i;
    xs[i] = x0 + h;
    ys[i] = y1;
    k1[i] = k7;

    x0 += h;
    y0 = y1;
    h = min(0.8 * h * pow(tol / erinf, 0.25), 4 * h);
  }
  ret.iterations = it;
  return ret;
};
