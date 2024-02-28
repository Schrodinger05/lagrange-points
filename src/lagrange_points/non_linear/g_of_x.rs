use super::CommonVals;

pub fn g(x: f64, y: f64, common: &CommonVals) -> f64 {
    let mu = common.mu;
    let x11 = common.x11;
    let x12 = common.x12;
    let x21 = common.x21;
    let x22 = common.x22;
    let y11 = common.y11;
    let y12 = common.y12;
    let y21 = common.y21;
    let y22 = common.y22;
    let mut res = 0.0;
    if x.is_normal() {
        let r11 = f64::sqrt((x - x11).powi(2) + (y - y11).powi(2));
        let r12 = f64::sqrt((x - x12).powi(2) + (y - y12).powi(2));
        let r21 = f64::sqrt((x - x21).powi(2) + (y - y21).powi(2));
        let r22 = f64::sqrt((x - x22).powi(2) + (y - y22).powi(2));
        let m11 = (1.0 - 2.0 * mu) / 2.0;
        let m12 = m11;
        let m21 = mu;
        let m22 = m21;
        res = y
            - m11 * (y - y11) / r11.powi(3)
            - m12 * (y - y12) / r12.powi(3)
            - m21 * (y - y21) / r21.powi(3)
            - m22 * (y - y22) / r22.powi(3);
    } else {
        res = 0.0;
    }
    res
}
