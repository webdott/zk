use polynomials::multilinear_polynomial::evaluation_form::{
    BlowUpDirection, MultiLinearPolynomial,
};

use ark_ec::pairing::Pairing;
use ark_ec::PrimeGroup;
use ark_ff::PrimeField;

// Given a set of tau values, we want to generate the lagrange basis array over the boolean hypercube for the number of variables
// We use the check 0 and check one principle
// E.g:
// 000 => (1 - a) * (1 - b) * (1 - c)
// 010 => (1 - a) * (b) * (1 - c)
// 111 => (a) * (b) * (c)
pub fn generate_lagrange_basis_for_n_variables<T: PrimeField>(n: usize, taus: &[T]) -> Vec<T> {
    if n != taus.len() {
        panic!("Length of variables does not match the number of Taus given!")
    }

    let length_of_lagrange_basis = 1 << n;
    let mut lagrange_basis = vec![T::one(); length_of_lagrange_basis];

    for i in 0..length_of_lagrange_basis {
        for (tau_idx, tau) in taus.iter().enumerate() {
            let skip_value = 1 << (n - tau_idx - 1);

            if (i / skip_value) % 2 == 0 {
                // if bit at variable index is 0 (turned off), use check-zero (1 - x)
                lagrange_basis[i] *= T::one() - *tau;
            } else {
                // else use check-one (x)
                lagrange_basis[i] *= tau;
            }
        }
    }

    lagrange_basis
}

// Encrypt lagrange basis by raising the generator to the power of each value in the basis
pub fn encrypt_lagrange_basis<T: PrimeField, P: Pairing>(lagrange_basis: &[T]) -> Vec<P::G1> {
    lagrange_basis
        .iter()
        .map(|&basis_value| P::G1::generator().mul_bigint(basis_value.into_bigint()))
        .collect::<Vec<P::G1>>()
}

pub fn blowup<T: PrimeField>(
    no_of_vars: usize,
    variable_idx: usize,
    evaluation_points: &[T],
) -> Vec<T> {
    let mut blown_up_evals = MultiLinearPolynomial::blow_up_n_times(
        BlowUpDirection::Left,
        variable_idx,
        &evaluation_points,
    );
    blown_up_evals = MultiLinearPolynomial::blow_up_n_times(
        BlowUpDirection::Right,
        no_of_vars - variable_idx - 1,
        &blown_up_evals,
    );

    blown_up_evals
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fr;

    #[test]
    pub fn test_lagrange_basis_for_n_variables_with_same_length_of_taus() {
        assert_eq!(
            generate_lagrange_basis_for_n_variables(3, &[Fr::from(5), Fr::from(2), Fr::from(3)]),
            vec![
                Fr::from(-8),
                Fr::from(12),
                Fr::from(16),
                Fr::from(-24),
                Fr::from(10),
                Fr::from(-15),
                Fr::from(-20),
                Fr::from(30),
            ]
        )
    }
}
