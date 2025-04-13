use crate::prover::LayerIndexProof;

use fft::fft::FFT;
use polynomials::univariate_polynomial::dense_coefficient_form::UnivariatePolynomial;

use ark_ff::{FftField, PrimeField};

pub fn compute_f_x_squared<T: FftField + PrimeField>(
    idx: usize,
    f_evaluations: (T, T),
    r: T,
    nth_root: Option<T>,
) -> T {
    let (f_positive_x, f_negative_x) = f_evaluations;

    let g_of_x = (f_positive_x + f_negative_x) / T::from(2);
    let h_of_x = (f_positive_x - f_negative_x) / (T::from(2) * nth_root.unwrap().pow([idx as u64]));

    g_of_x + (r * h_of_x)
}

pub fn fold_layer<T: FftField + PrimeField>(evaluations: &[T], r: T) -> Vec<T> {
    let half_length = evaluations.len() / 2;
    let mut folded_layer: Vec<T> = Vec::with_capacity(half_length);
    let nth_root = T::get_root_of_unity(evaluations.len() as u64);

    for idx in 0..half_length {
        let negative_idx = idx + half_length;

        folded_layer.push(compute_f_x_squared(
            idx,
            (evaluations[idx], evaluations[negative_idx]),
            r,
            nth_root,
        ));
    }

    folded_layer
}

pub fn get_f_squared_from_folded_layer<T: FftField + PrimeField>(
    idx: usize,
    folded_layer: &Vec<LayerIndexProof<T>>,
) -> T {
    folded_layer
        .iter()
        .find(|proof| proof.index == idx)
        .unwrap()
        .value
}

pub fn pad_polynomial<T: FftField + PrimeField>(
    values: &[T],
    final_length: usize,
    number_to_add: T,
) -> Vec<T> {
    let mut final_pad = values.to_vec();

    for _ in final_pad.len()..final_length {
        final_pad.push(number_to_add);
    }

    final_pad
}

pub fn perform_reed_solomon<T: FftField + PrimeField>(
    polynomial: UnivariatePolynomial<T>,
    blow_up_factor: usize,
) -> Vec<T> {
    let blown_up_length = polynomial.coefficients.len() * blow_up_factor;
    let padded_polynomial_coefficients = pad_polynomial(
        &polynomial.coefficients,
        blown_up_length.next_power_of_two(),
        T::zero(),
    );

    FFT::convert_to_evaluations(&padded_polynomial_coefficients)
}

pub fn get_layer_proof_indexes(n: usize, given_index: usize) -> (usize, usize) {
    let half_length = n / 2;

    if given_index < half_length {
        return (given_index, given_index + half_length);
    }

    (given_index - half_length, given_index)
}
