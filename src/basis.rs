use crate::symmetries::get_eq_class;

#[derive(PartialEq, Debug)]
pub struct Representer {
    pub value: usize,
    pub digits: Vec<u8>,
    pub sigmas: Vec<f32>,
}

impl Representer {
    fn new(value: usize, digits: &Vec<u8>, spin: f32) -> Self {
        let mut digits = digits.clone();
        digits.push(digits[0]);
        let sigmas = digits
            .iter()
            .map(|digit| *digit as f32 - spin)
            .collect();

        Self {
            value,
            digits,
            sigmas,
        }
    }
}

#[derive(PartialEq, Debug)]
pub struct BasisStates {
    pub length: usize,
    pub representers: Vec<Representer>,
    pub representer_map: Vec<usize>,
    pub symmetry_signs: Vec<Vec<[bool; 3]>>,
}

impl BasisStates {
    pub fn new(
        base: u8,
        spin: f32,
        length: u8,
        total_s_z: u8,
        base_powers: &Vec<usize>,
    ) -> Self {
        let max_representer_value = (base as usize).pow(length as u32);
        let mut representers = Vec::new();
        let mut representer_map = vec![max_representer_value; max_representer_value];
        let mut symmetry_signs = vec![Vec::new(); max_representer_value];
        let allowed_digit_sum = (spin * length as f32) as u8 + total_s_z;
        let mut index = 0;

        for representer_value in 1..max_representer_value-1 {
            if representer_map[representer_value] == max_representer_value {
                let digits = get_digits(representer_value, base, length);
                let digit_sum: u8 = digits.iter().sum();
                
                if digit_sum == allowed_digit_sum {
                    let representer = Representer::new(
                        representer_value,
                        &digits,
                        spin,
                    );
                    representers.push(representer);

                    let eq_class = get_eq_class(
                        digits,
                        base,
                        length,
                        base_powers,
                    );
        
                    for (other_representer_value, symmetry_sign) in eq_class.iter() {
                        representer_map[*other_representer_value] = index;
                        symmetry_signs[*other_representer_value].push(*symmetry_sign);
                    }

                    index += 1;
                }
            }
        }
        
        Self {
            length: representers.len(),
            representers,
            representer_map,
            symmetry_signs,
        }
    }
}

pub fn get_digits(
    value: usize,
    base: u8,
    length: u8,
) -> Vec<u8> {
    let mut value = value;
    let mut digits = Vec::with_capacity(length as usize);

    for _ in 0..length {
        digits.push((value % base as usize) as u8);
        value /= base as usize;
    }

    digits
}