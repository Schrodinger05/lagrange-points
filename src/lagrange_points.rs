mod non_linear;
use std::{fs, io::Write};

use non_linear::Point;

use self::calculation_loop::{Limits, Vals};

mod calculation_loop;
pub mod stability;

const NONE: f64 = f64::NAN;

pub fn calculate_l1() {
    let initial_guess = Point::create_a_point(-0.12, 0.0);
    let results = calculation_loop::calculation_loop(
        &initial_guess,
        0.0001,
        NONE,
        Limits::X21,
        Limits::X12,
        false,
    );

    print_to_file("l1.txt", results);
}

pub fn calculate_l2() {
    let initial_guess = Point::create_a_point(0.9, 0.0);
    let results = calculation_loop::calculation_loop(
        &initial_guess,
        0.0001,
        NONE,
        Limits::NoLimit,
        Limits::X22,
        false,
    );
    print_to_file("l2.txt", results);
}

pub fn calculate_l3() {
    let initial_guess = Point::create_a_point(-0.9, 0.0);
    let results = calculation_loop::calculation_loop(
        &initial_guess,
        -0.0001,
        NONE,
        Limits::X11,
        Limits::NoLimit,
        false,
    );
    print_to_file("l3.txt", results)
}

pub fn calculate_l4() {
    let initial_guess = Point::create_a_point(0.1, 0.0);
    let results = calculation_loop::calculation_loop(
        &initial_guess,
        NONE,
        0.0001,
        Limits::NoLimit,
        Limits::Num(0.1),
        true,
    );
    print_to_file("l4.txt", results)
}

pub fn calculate_l5() {
    let initial_guess = Point::create_a_point(0.1, 0.0);
    let results = calculation_loop::calculation_loop(
        &initial_guess,
        NONE,
        -0.0001,
        Limits::Num(-0.1),
        Limits::NoLimit,
        true,
    );
    print_to_file("l5.txt", results)
}

fn print_to_file(name: &str, results: Vec<Vals>) {
    let mut file = fs::File::create(name).expect("Couldn't Create File");
    let out = format!("{:36}{:36}{:36}{:36}{:36}\n", "x", "y", "mu", "d1", "d2");
    file.write_all(out.as_bytes())
        .expect("Couldnt make initial write");
    for result in results {
        let out = format!(
            "{:.32}  {:.32}  {:.32}  {:.32}  {:.32}\n",
            result.x, result.y, result.mu, result.d1, result.d2
        );
        file.write_all(out.as_bytes())
            .expect("Could Not Write to file");
    }
    println!("Finished writing to file: {}", name);
}
