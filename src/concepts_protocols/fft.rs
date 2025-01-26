use ark_ff::FftField;

pub struct Polynomial<T: FftField> {
    _marker: std::marker::PhantomData<T>,
}

impl<T: FftField> Polynomial<T> {
    fn fft(coefficients_or_values: &[T], is_inverse: bool) -> Vec<T> {
        // n = len of coefficients_or_values
        // ye = a0, a2, a4....an
        // yo = a1, a3, a5....an-1

        // w -> (roots of unity):
        // *   // (FFT) => e^(2 * PI * i)/n
        // *   // (IFFT) => (1/n) * e^-(2 * PI * i)/n

        // P(w^j) = ye[j] + w^j * (yo[j])
        // P(-w^j) = ye[j] - w^j * (yo[j]); -w^j = w^(j + (n/2))

        let n = coefficients_or_values.len();

        // if it gets to base case, return current coefficients_or_values;
        if n == 1 {
            return vec![coefficients_or_values[0]];
        }

        let (mut even_sequence, mut odd_sequence) = (vec![], vec![]);

        coefficients_or_values
            .iter()
            .enumerate()
            .for_each(|(idx, num)| {
                if idx % 2 == 0 {
                    even_sequence.push(*num);
                } else {
                    odd_sequence.push(*num);
                }
            });

        // recurse to find the further ffts for even and odd sequences
        let (ye, yo) = (
            Self::fft(&even_sequence, is_inverse),
            Self::fft(&odd_sequence, is_inverse),
        );

        let root_of_unity = T::get_root_of_unity(n as u64);

        let w = match root_of_unity {
            Some(rou) => {
                if is_inverse {
                    rou.inverse()
                } else {
                    Some(rou)
                }
            }
            None => None,
        };

        let mut y = vec![T::from(0); n];

        (0..n / 2).into_iter().for_each(|j| {
            let wj = w.unwrap().pow(vec![j as u64]);

            y[j] = ye[j] + wj * yo[j];
            y[j + (n / 2)] = ye[j] - wj * yo[j];
        });

        y
    }

    // Perform Fast Fourier Transforms to convert Polynomial to Values (Samples) Representation
    // This can be done in O(nlogn) time to perform a linear O(n) operation in Sample like evaluation that would have originally taken O(n^2) in Coefficients form
    pub fn transform_to_values(coefficients: &[T]) -> Vec<T> {
        Self::fft(coefficients, false)
    }

    // Perform inverse Fast Fourier Transform to convert Sample representation back to Coefficients
    // This can be done in O(nlogn) time as well to perform a linear O(n) operation in Coefficients form like Multiplication that would have originally taken O(n^2) in Sample form
    pub fn transform_to_coefficients(values: &[T]) -> Vec<T> {
        Self::fft(&values, true)
            .iter()
            .map(|x| *x / T::from(values.len() as u64))
            .collect()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bls12_377::Fr;

    #[test]
    pub fn test_fft_and_ifft() {
        let coefficients = vec![Fr::from(5), Fr::from(3), Fr::from(2), Fr::from(1)];
        let values = Polynomial::transform_to_values(&coefficients);
        let result_coefficients = Polynomial::transform_to_coefficients(&values);

        assert_eq!(result_coefficients, coefficients,)
    }
}
