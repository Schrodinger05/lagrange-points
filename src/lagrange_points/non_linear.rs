mod f_of_x;
mod g_of_x;

use f_of_x::f;
use g_of_x::g;

use super::calculation_loop::CommonVals;

pub struct Point {
    pub x: f64,
    pub y: f64,
}

impl Point {
    pub fn create_a_point(x: f64, y: f64) -> Self {
        Self { x, y }
    }
}

pub fn non_linear_system(point: &Point, common: &CommonVals) -> Point {
    const H: f64 = 0.01;
    let mut x = point.x;
    let mut y = point.y;
    let mut a: f64;
    let mut b: f64;
    let mut c: f64;
    let mut d: f64;
    let mut t: f64;
    let mut xm: f64;
    let mut xn: f64;
    let mut p: f64 = 1.0;
    let mut q: f64 = 1.0;

    //begin calculation
    while f64::abs(p) + f64::abs(q) > 1.0e-10 {
        a = f(x + H, y, common);
        b = g(x + H, y, common);
        a = (a - f(x - H, y, common)) / 2.0 / H;
        b = (b - g(x - H, y, common)) / 2.0 / H;
        c = f(x, y + H, common);
        d = g(x, y + H, common);
        c = (c - f(x, y - H, common)) / 2.0 / H;
        d = (d - g(x, y - H, common)) / 2.0 / H;
        t = a * d - b * c;
        xm = f(x, y, common);
        xn = g(x, y, common);
        p = (xm * d - xn * c) / t;
        q = (xn * a - xm * b) / t;
        x = x - p;
        y = y - q;
    }

    Point::create_a_point(x, y)
}
