// Shared vector arithmetic for states made of Complex (Complex.js) values,
// used by both numeric_DOPRI.js and numeric_RK4.js. A "complex vector" here
// is either a single Complex value or an (arbitrarily nested) array of them.
//
// Requires Complex.js to be loaded first.

function isArray(x) {
  return Object.prototype.toString.call(x) === "[object Array]";
}

function cAdd(a, b) {
  if (isArray(a)) {
    var n = a.length, r = new Array(n);
    for (var i = 0; i < n; i++) r[i] = cAdd(a[i], b[i]);
    return r;
  }
  return a.add(b);
}

function cSub(a, b) {
  if (isArray(a)) {
    var n = a.length, r = new Array(n);
    for (var i = 0; i < n; i++) r[i] = cSub(a[i], b[i]);
    return r;
  }
  return a.sub(b);
}

// s is always a plain real (JS number) scalar, e.g. a step size.
function cScale(a, s) {
  if (isArray(a)) {
    var n = a.length, r = new Array(n);
    for (var i = 0; i < n; i++) r[i] = cScale(a[i], s);
    return r;
  }
  return a.mul(s);
}

function cNorminf(a) {
  if (isArray(a)) {
    var m = 0;
    for (var i = 0; i < a.length; i++) m = Math.max(m, cNorminf(a[i]));
    return m;
  }
  return a.abs();
}
