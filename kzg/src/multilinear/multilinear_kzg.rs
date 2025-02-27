use ark_ff::PrimeField;

pub fn get_lagrange_basis_for_n_variables<T: PrimeField>(n: usize, taus: &[T]) -> Vec<T> {
    if n != taus.len() {
        panic!("Length of variables does not match the number of Taus given!")
    }

    let length_of_lagrange_basis = 1 << n;
    let mut lagrange_basis = vec![T::one(); length_of_lagrange_basis];

    for i in 0..length_of_lagrange_basis {
        for (tau_idx, tau) in taus.iter().enumerate() {
            let skip_value = 1 << (n - tau_idx - 1);

            if (i / skip_value) % 2 == 0 {
                lagrange_basis[i] *= T::one() - *tau;
            } else {
                lagrange_basis[i] *= tau;
            }
        }
    }

    lagrange_basis
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::Fq;

    #[test]
    pub fn test_lagrange_basis_for_n_variables_with_same_length_of_taus() {
        assert_eq!(
            get_lagrange_basis_for_n_variables(3, &[Fq::from(5), Fq::from(2), Fq::from(3)]),
            vec![
                Fq::from(-8),
                Fq::from(12),
                Fq::from(16),
                Fq::from(-24),
                Fq::from(10),
                Fq::from(-15),
                Fq::from(-20),
                Fq::from(30),
            ]
        )
    }
}
