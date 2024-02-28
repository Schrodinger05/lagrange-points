use super::non_linear::{non_linear_system, Point};

pub struct Vals {
    pub mu: f64,
    pub d1: f64,
    pub d2: f64,
    pub x: f64,
    pub y: f64,
    pub stability: bool,
}

pub struct CommonVals {
    pub x11: f64,
    pub x12: f64,
    pub x21: f64,
    pub x22: f64,
    pub y11: f64,
    pub y12: f64,
    pub y21: f64,
    pub y22: f64,
    pub mu: f64,
}

pub enum Limits {
    X11,
    X12,
    X21,
    X22,
    NoLimit,
    Num(f64),
}

impl Vals {
    fn new_vals(mu: f64, d1: f64, d2: f64, point: Point) -> Self {
        Self {
            mu,
            d1,
            d2,
            x: point.x,
            y: point.y,
            stability: false,
        }
    }
}

impl CommonVals {
    fn new(mu: f64, d1: f64, d2: f64) -> Self {
        Self {
            x11: -2.0 * mu - d1 / 2.0,
            x12: -2.0 * mu + d1 / 2.0,
            x21: -2.0 * mu - d2 / 2.0 + 1.0,
            x22: -2.0 * mu + d2 / 2.0 + 1.0,
            y11: 0.0,
            y12: 0.0,
            y21: 0.0,
            y22: 0.0,
            mu,
        }
    }
}

fn determine_limit(common: &CommonVals, limit: &Limits) -> f64 {
    match limit {
        Limits::X11 => common.x11,
        Limits::X12 => common.x12,
        Limits::X21 => common.x21,
        Limits::X22 => common.x22,
        Limits::NoLimit => f64::NAN,
        Limits::Num(num) => *num,
    }
}

pub fn calculation_loop(
    initial_guess: &Point,
    guess_step_x: f64,
    guess_step_y: f64,
    upper_limit: Limits,
    lower_limit: Limits,
    can_calculate_y: bool,
) -> Vec<Vals> {
    let mut mu = 0.0;
    let mut d1 = 0.0;
    let mut d2 = 0.0;
    let mut guess: Point;
    let mut results: Vec<Vals> = Vec::new();

    while d2 < 0.131440 {
        while d1 < 0.736068 {
            while mu < 0.25 {
                mu += 0.01;
                // determines a new set of common variables
                let common = CommonVals::new(mu, d1, d2);
                guess = Point::create_a_point(initial_guess.x, initial_guess.y);
                //determines what type of upper and lower values we will use
                let upper = determine_limit(&common, &upper_limit);
                let lower = determine_limit(&common, &lower_limit);
                loop {
                    //increment if increment value is normal
                    if !guess_step_x.is_nan() {
                        guess.x = guess.x + guess_step_x;
                    }
                    if !guess_step_y.is_nan() {
                        guess.y = guess.y + guess_step_y;
                    }
                    //calculate the result
                    let result = non_linear_system(&guess, &common);
                    // determine consditions
                    if can_calculate_y {
                        if upper.is_nan() {
                            if result.y > lower {
                                results.push(Vals::new_vals(mu, d1, d2, result));
                                break;
                            }
                        }
                        if lower.is_nan() {
                            if result.y < upper {
                                results.push(Vals::new_vals(mu, d1, d2, result));
                                break;
                            }
                        }
                        if result.y > lower && result.y < upper {
                            results.push(Vals::new_vals(mu, d1, d2, result));
                            break;
                        }
                    } else {
                        if upper.is_nan() {
                            if result.x > lower {
                                results.push(Vals::new_vals(mu, d1, d2, result));
                                break;
                            }
                        }
                        if lower.is_nan() {
                            if result.x < upper {
                                results.push(Vals::new_vals(mu, d1, d2, result));
                                break;
                            }
                        }
                        if result.x > lower && result.x < upper {
                            results.push(Vals::new_vals(mu, d1, d2, result));
                            break;
                        }
                    }
                }
            }
            d1 += 0.01;
            mu = 0.0;
        }
        d2 += 0.01;
        d1 = 0.0;
    }
    results
}
/*
pub fn l1_loop(vals: &mut Vals) {
    let mut l1_file = fs::File::create("l1.txt").expect("failed to create l1 file");
    l1_file
        .write_all(b"This file has all possible l1 points")
        .expect("Couldn't write to file!");
}
*/
