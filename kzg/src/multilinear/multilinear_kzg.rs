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

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::Fq;

    #[test]
    pub fn test_lagrange_basis_for_n_variables_with_same_length_of_taus() {
        assert_eq!(
            generate_lagrange_basis_for_n_variables(3, &[Fq::from(5), Fq::from(2), Fq::from(3)]),
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
