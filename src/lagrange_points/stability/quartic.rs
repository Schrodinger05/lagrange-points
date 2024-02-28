use num::{
    complex::{Complex64, ComplexFloat},
    zero, FromPrimitive,
};

pub fn quadratic_roots(a0: f64, b0: f64, c0: f64) -> (Complex64, Complex64) {
    let a = Complex64::from_f64(b0 / a0).expect("Coulnd't convert to complex");
    let b = Complex64::from_f64(c0 / a0).expect("Coulnd't convert to complex");

    let a0 = -0.5 * a;
    let delta = a0 * a0 - b;
    let sqrt_delta = Complex64::sqrt(delta);

    (a0 - sqrt_delta, a0 + sqrt_delta)
}

pub fn cubic_root_single(a0: f64, b0: f64, c0: f64, d0: f64) -> Complex64 {
    let mut res: Complex64 = Complex64::new(0.0, 0.0);

    let a = Complex64::from_f64(b0 / a0).expect("Coulnd't convert to complex");
    let b = Complex64::from_f64(c0 / a0).expect("Coulnd't convert to complex");
    let c = Complex64::from_f64(d0 / a0).expect("Coulnd't convert to complex");

    let third = 1.0 / 3.0;
    let a13 = a * third;
    let a2 = a13 * a13;

    let f = third * b - a2;
    let g = a13 * (2.0 * a2 - b) + c;
    let h = 0.25 * g * g + f * f * f;

    fn cubic_root(x: Complex64, third: f64) -> Complex64 {
        if x.re() >= 0.0 {
            x.powf(third)
        } else {
            -(-x).powf(third)
        }
    }

    if f == g && g == h && h == zero() {
        res = -cubic_root(c, third);
    } else if h.re() <= zero() {
        let j = Complex64::sqrt(-f);
        let k = Complex64::acos(-0.5 * g / (j * j * j));
        let m = Complex64::cos(third * k);

        res = 2.0 * j * m - a13;
    } else {
        let sqrt_h = Complex64::sqrt(h);
        let s = cubic_root(-0.5 * g + sqrt_h, third);
        let u = cubic_root(-0.5 * g - sqrt_h, third);
        let s_plus_u = s + u;

        res = s_plus_u - a13
    }

    res
}

pub fn quartic_roots(a0: f64, b0: f64, c0: f64, d0: f64, e0: f64) -> [Complex64; 4] {
    let a = Complex64::from_f64(b0 / a0).expect("Coulnd't convert to complex");
    let b = Complex64::from_f64(c0 / a0).expect("Coulnd't convert to complex");
    let c = Complex64::from_f64(d0 / a0).expect("Coulnd't convert to complex");
    let d = Complex64::from_f64(e0 / a0).expect("Coulnd't convert to complex");
    let mut t: Complex64;

    let a0 = 0.25 * a;
    let a02 = a0 * a0;

    let p = 3.0 * a02 - 0.5 * b;
    let q = a * a02 - b * a0 + 0.5 * c;
    let r = 3.0 * a02 * a02 - b * a02 + c * a0 - d;

    let z0 = cubic_root_single(1.0, p.re(), r.re(), p.re() * r.re() - 0.5 * q.re() * q.re());
    let s = Complex64::sqrt(2.0 * p + 2.0 * z0 + 0.0);

    if s == zero() {
        t = z0 * z0 + r;
    } else {
        t = -q / s;
    }
    let (r0, r1) = quadratic_roots(1.0, s.re(), z0.re() + t.re());
    let (r2, r3) = quadratic_roots(1.0, -s.re(), z0.re() - t.re());

    [(r0 - a0), (r1 - a0), (r2 - a0), (r3 - a0)]
}
