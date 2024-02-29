use super::calculation_loop::Vals;
use num::complex::ComplexFloat;
mod quartic;

pub fn calculate_stability(results: &mut Vec<Vals>) {
    for result in results {
        let uxx = calculate_uxx(result);
        let uyy = calculate_uyy(result);
        let uxy = calculate_uxy(result);
        let roots = quartic::quartic_roots(1.0, 0.0, 4.0 - uxx - uyy, 0.0, uxx * uyy - uxy.powi(2));
        for root in roots {
            if root.re() != 0.0 {
                result.stability = false;
                break;
            } else if root.im() != 0.0 {
                result.stability = true;
            }
        }
    }
}

fn calculate_uxx(values: &Vals) -> f64 {
    let m = values.mu;
    let d1 = values.d1;
    let d2 = values.d2;
    let x = values.x;
    let y = values.y;
    let z: f64 = 0.0;

    let m21 = m;
    let m22 = m21;

    let res = 1.0
        + (3.0 * (1.0 - 2.0 * m) * (-d1 / 2.0 + 2.0 * m + x).powi(2))
            / (2.0 * ((-d1 / 2.0 + 2.0 * m + x).powi(2) + y.powi(2) + z.powi(2)).powf(5.0 / 2.0));
    let res = res
        - (1.0 - 2.0 * m)
            / (2.0 * ((-d1 / 2.0 + 2.0 * m + x).powi(2) + y.powi(2) + z.powi(2)).powf(3.0 / 2.0));
    let res = res
        + (3.0 * (1.0 - 2.0 * m) * (d1 / 2.0 + 2.0 * m + x).powi(2))
            / (2.0 * ((d1 / 2.0 + 2.0 * m + x).powi(2) + y.powi(2) + z.powi(2)).powf(5.0 / 2.0));
    let res = res
        - (1.0 - 2.0 * m)
            / (2.0 * ((d1 / 2.0 + 2.0 * m + x).powi(2) + y.powi(2) + z.powi(2)).powf(3.0 / 2.0));
    let res = res
        + (3.0 * m22 * (-1.0 - d2 / 2.0 + 2.0 * m + x).powi(2))
            / ((-1.0 - d2 / 2.0 + 2.0 * m + x).powi(2) + y.powi(2) + z.powi(2)).powf(5.0 / 2.0);
    let res = res
        - m22 / ((-1.0 - d2 / 2.0 + 2.0 * m + x).powi(2) + y.powi(2) + z.powi(2)).powf(3.0 / 2.0);
    let res = res
        + (3.0 * m21 * (-1.0 + d2 / 2.0 + 2.0 * m + x).powi(2))
            / ((-1.0 + d2 / 2.0 + 2.0 * m + x).powi(2) + y.powi(2) + z.powi(2)).powf(5.0 / 2.0);
    let res = res
        - m21 / ((-1.0 + d2 / 2.0 + 2.0 * m + x).powi(2) + y.powi(2) + z.powi(2)).powf(3.0 / 2.0);

    res
}

fn calculate_uyy(values: &Vals) -> f64 {
    let m = values.mu;
    let d1 = values.d1;
    let d2 = values.d2;
    let x = values.x;
    let y = values.y;
    let z: f64 = 0.0;

    let m21 = m;
    let m22 = m21;

    let res = 1.0
        + (3.0 * (1.0 - 2.0 * m) * y.powi(2))
            / (2.0 * ((-d1 / 2.0 + 2.0 * m + x).powi(2) + y.powi(2) + z.powi(2)).powf(5.0 / 2.0));
    let res = res
        - (1.0 - 2.0 * m)
            / (2.0 * ((-d1 / 2.0 + 2.0 * m + x).powi(2) + y.powi(2) + z.powi(2)).powf(3.0 / 2.0));
    let res = res
        + (3.0 * m22 * y.powi(2))
            / (2.0 * ((d1 / 2.0 + 2.0 * m + x).powi(2) + y.powi(2) + z.powi(2)).powf(5.0 / 2.0));
    let res = res
        - (1.0 - 2.0 * m)
            / (2.0 * ((d1 / 2.0 + 2.0 * m + x).powi(2) + y.powi(2) + z.powi(2)).powf(3.0 / 2.0));
    let res = res
        + (3.0 * m22 * y.powi(2))
            / ((-1.0 - d2 / 2.0 + 2.0 * m + x).powi(2) + y.powi(2) + z.powi(2)).powf(5.0 / 2.0);
    let res = res
        - m22 / ((-1.0 - d2 / 2.0 + 2.0 * m + x).powi(2) + y.powi(2) + z.powi(2)).powf(3.0 / 2.0);
    let res = res
        + (3.0 * m21 * y.powi(2))
            / ((-1.0 + d2 / 2.0 + 2.0 * m + x).powi(2) + y.powi(2) + z.powi(2)).powf(5.0 / 2.0);
    let res = res
        - m21 / ((-1.0 + d2 / 2.0 + 2.0 * m + x).powi(2) + y.powi(2) + z.powi(2)).powf(3.0 / 2.0);

    res
}

fn calculate_uxy(values: &Vals) -> f64 {
    let m = values.mu;
    let d1 = values.d1;
    let d2 = values.d2;
    let x = values.x;
    let y = values.y;
    let z: f64 = 0.0;

    let m21 = m;
    let m22 = m21;

    let res = (3.0 * (1.0 - 2.0 * m) * (-d1 / 2.0 + 2.0 * m + x) * y) / 2.0
        * ((-d1 / 2.0 + 2.0 * m + x).powi(2) + y.powi(2) + z.powi(2)).powf(5.0 / 2.0);
    let res = res
        + (3.0 * (1.0 - 2.0 * m) * (d1 / 2.0 + 2.0 * m + x) * y) / 2.0
            * ((d1 / 2.0 + 2.0 * m + x).powi(2) + y.powi(2) + z.powi(2)).powf(5.0 / 2.0);
    let res = res
        + (3.0 * m22 * (-1.0 - d2 / 2.0 + 2.0 * m + x) * y)
            / ((-1.0 - d2 / 2.0 + 2.0 * m + x).powi(2) + y.powi(2) + z.powi(2)).powf(5.0 / 2.0);
    let res = res
        + (3.0 * m21 * (-1.0 + d2 / 2.0 + 2.0 * m + x) * y)
            / ((-1.0 + d2 / 2.0 + 2.0 * m + x).powi(2) + y.powi(2) + z.powi(2)).powf(5.0 / 2.0);

    res
}
