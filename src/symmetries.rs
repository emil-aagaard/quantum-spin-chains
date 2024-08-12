pub fn get_symmetry_factors(
    symmetry_signs: &Vec<Vec<[bool; 3]>>,
    time_reversal_eigenvalue: bool,
    reflection_eigenvalue: bool,
    translation_eigenvalue: bool,
) -> Vec<f32> {
    let mut symmetry_factors = Vec::with_capacity(symmetry_signs.len());

    for symmetry_signs_ in symmetry_signs {
        let mut symmetry_factor = 0.0;

        for symmetry_sign in symmetry_signs_ {
            let total_sign = symmetry_sign[0]&time_reversal_eigenvalue
                ^ symmetry_sign[1]&reflection_eigenvalue
                ^ symmetry_sign[2]&translation_eigenvalue;

            symmetry_factor += if total_sign {-1.0} else {1.0}
        }

        symmetry_factors.push(symmetry_factor)
    }

    symmetry_factors
}

pub fn get_eq_class(
    digits: Vec<u8>,
    base: u8,
    length: u8,
    base_powers: &Vec<usize>,
) -> Vec<(usize, [bool; 3])> {
    let all_time_reversed_reflected_digits = time_reverse_reflect_digits(
        digits,
        base,
    );

    let mut eq_class = Vec::new();

    for (time_reversed_reflected_digits, sign) in all_time_reversed_reflected_digits {
        let mut translations = get_translations(
            time_reversed_reflected_digits,
            length,
            base_powers,
            sign,
        );

        eq_class.append(&mut translations)
    }

    eq_class
}

fn time_reverse_reflect_digits(
    digits: Vec<u8>,
    base: u8,
) -> [(Vec<u8>, [bool; 3]); 4] {
    let [digits, time_reversed_digits] = time_reverse_digits(digits, base);
    let [digits, reflected_digits] = reflect_digits(digits);
    let [reflected_digits, time_reversed_reflected_digits] = time_reverse_digits(reflected_digits, base);

    [
        (digits, [false, false, false]),
        (time_reversed_digits, [true, false, false]),
        (reflected_digits, [false, true, false]),
        (time_reversed_reflected_digits, [true, true, false]),
    ]
}

fn time_reverse_digits(
    digits: Vec<u8>,
    base: u8,
) -> [Vec<u8>; 2] {
    let time_reversed_digits = digits
        .iter()
        .map(
            |digit| base - 1 - *digit
        )
        .collect();

    [digits, time_reversed_digits]
}

fn reflect_digits(
    digits: Vec<u8>,
) -> [Vec<u8>; 2] {
    let mut revflected_digits = digits.clone();
    revflected_digits.reverse();

    [digits, revflected_digits]
}

fn get_translations(
    digits: Vec<u8>,
    length: u8,
    base_powers: &Vec<usize>,
    mut symmetry_sign: [bool; 3],
) -> Vec<(usize, [bool; 3])> {
    let basis_state = get_translation(&digits, 0, length, base_powers);
    let mut translations = vec![(basis_state, symmetry_sign)];
    let mut translation;
    
    for translate in 1..length {
        symmetry_sign[2] = !symmetry_sign[2];
        translation = get_translation(&digits, translate, length, base_powers);
        translations.push((translation, symmetry_sign));
    }

    translations
}

fn get_translation(
    digits: &Vec<u8>,
    translate: u8,
    length: u8,
    base_powers: &Vec<usize>,
) -> usize {
    let mut basis_state = 0;

    for index in 0..length {
        let digit_index = ((translate + index) % length) as usize;
        basis_state += digits[digit_index] as usize * base_powers[index as usize];
    }

    basis_state
}