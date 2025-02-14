use crate::concepts_protocols::fiat_shamir::transcript::Transcript;
use crate::polynomials::multilinear_polynomial::MultiLinearPolynomial;
use crate::polynomials::sum_polynomial::SumPolynomial;
use crate::polynomials::univariate_polynomial::UnivariatePolynomial;
use ark_ff::{BigInteger, PrimeField};
use std::iter;
use std::marker::PhantomData;

pub enum ComposedPolynomial<T: PrimeField> {
    SumPolynomial(SumPolynomial<T>),
    MultilinearPolynomial(MultiLinearPolynomial<T>),
}

#[derive(Debug)]
pub struct SumCheckProof<T: PrimeField> {
    initial_claim_sum: T,
    round_polys: Vec<UnivariatePolynomial<T>>,
}

pub struct Prover<T: PrimeField> {
    _marker: PhantomData<T>,
}

impl<T: PrimeField> Prover<T> {
    // this generates a set of points to partially evaluate a polynomial
    fn generate_evaluation_points(
        transcript: &mut Transcript<T>,
        variables_length: usize,
    ) -> Vec<Option<T>> {
        iter::repeat(())
            .enumerate()
            .map(|(idx, _x)| {
                if idx == 0 {
                    return Some(transcript.sample_challenge());
                }

                return None;
            })
            .take(variables_length)
            .collect()
    }

    fn generate_round_polys(
        initial_polynomial: &ComposedPolynomial<T>,
        transcript: &mut Transcript<T>,
    ) -> Vec<UnivariatePolynomial<T>> {
        let (
            mut resulting_multi_polynomial,
            mut resulting_sum_polynomial,
            mut round_polys,
            number_of_variables,
        ) = match initial_polynomial {
            ComposedPolynomial::SumPolynomial(polynomial) => (
                None,
                Some(polynomial.clone()),
                Vec::with_capacity(polynomial.number_of_variables() as usize),
                polynomial.number_of_variables(),
            ),
            ComposedPolynomial::MultilinearPolynomial(polynomial) => (
                Some(polynomial.clone()),
                None,
                Vec::with_capacity(polynomial.number_of_variables() as usize),
                polynomial.number_of_variables(),
            ),
        };

        // keep adding current polynomial step and sum to the round_polys vec
        (0..number_of_variables).for_each(|_i| {
            let (mut claimed_sum, mut evaluated_polynomial_over_boolean_hypercube) = (
                T::from(0),
                UnivariatePolynomial::new(vec![T::from(0), T::from(0)]),
            );

            if let Some(sum_poly) = &resulting_sum_polynomial {
                let degree = sum_poly.degree();
                let mut evaluation_points = vec![T::from(0); degree + 1];

                for i in 0..degree + 1 {
                    let mut _points = vec![None; sum_poly.number_of_variables() as usize];
                    _points[0] = Some(T::from(i as u8));

                    let res = sum_poly.partial_evaluate(&_points).reduce();

                    evaluation_points[i] = res.iter().sum();
                }

                evaluated_polynomial_over_boolean_hypercube = UnivariatePolynomial::interpolate(
                    vec![T::from(0), T::from(1), T::from(2)],
                    evaluation_points,
                );

                claimed_sum = evaluated_polynomial_over_boolean_hypercube
                    .evaluate_sum_over_boolean_hypercube();
            } else if let Some(multi_poly) = &resulting_multi_polynomial {
                let evaluation_points = multi_poly.get_evaluation_points();
                let (first_half, second_half) =
                    evaluation_points.split_at(evaluation_points.len() / 2);

                let (eval_0, eval_1) = (
                    T::from(first_half.iter().sum::<T>()),
                    T::from(second_half.iter().sum::<T>()),
                );

                claimed_sum = eval_0 + eval_1;

                evaluated_polynomial_over_boolean_hypercube = UnivariatePolynomial::interpolate(
                    vec![T::from(0), T::from(1)],
                    vec![eval_0, eval_1],
                )
            }

            transcript.append(&claimed_sum.into_bigint().to_bytes_le());
            transcript.append(&evaluated_polynomial_over_boolean_hypercube.to_bytes());

            if let Some(sum_poly) = &resulting_sum_polynomial {
                let points = Self::generate_evaluation_points(
                    transcript,
                    sum_poly.number_of_variables() as usize,
                );
                resulting_sum_polynomial = Some(sum_poly.partial_evaluate(&points));
            } else if let Some(multi_poly) = &resulting_multi_polynomial {
                let points = Self::generate_evaluation_points(
                    transcript,
                    multi_poly.number_of_variables() as usize,
                );
                resulting_multi_polynomial = Some(multi_poly.evaluate(points));
            }

            round_polys.push(evaluated_polynomial_over_boolean_hypercube);
        });

        round_polys
    }

    // This creates a sum check proof struct that with the round_polys generated and an initial claim sum
    pub fn generate_sumcheck_proof(init_polynomial: &MultiLinearPolynomial<T>) -> SumCheckProof<T> {
        let mut transcript = Transcript::new();

        // append initial polynomial to transcript to initiate process
        transcript.append(&init_polynomial.to_bytes());

        let round_polys = Self::generate_round_polys(
            &ComposedPolynomial::MultilinearPolynomial(init_polynomial.clone()),
            &mut transcript,
        );

        SumCheckProof {
            initial_claim_sum: init_polynomial.evaluation_sum(),
            round_polys,
        }
    }

    pub fn generate_proof_for_partial_verify(
        initial_claim_sum: T,
        init_poly: SumPolynomial<T>,
    ) -> SumCheckProof<T> {
        let round_polys = Self::generate_round_polys(
            &ComposedPolynomial::SumPolynomial(init_poly),
            &mut Transcript::new(),
        );

        SumCheckProof {
            initial_claim_sum,
            round_polys,
        }
    }
}

struct Verifier<T: PrimeField> {
    _marker: PhantomData<T>,
}

impl<T: PrimeField> Verifier<T> {
    // This check ensures that the last univariate evaluated at a variable is equal to the initial polynomial evaluated at all sampled values
    pub fn perform_oracle_check(
        initial_polynomial: &MultiLinearPolynomial<T>,
        evaluation_values: Vec<Option<T>>,
        final_univariate_poly: &UnivariatePolynomial<T>,
    ) -> bool {
        let last_value = *evaluation_values.last().unwrap();

        *initial_polynomial
            .evaluate(evaluation_values)
            .get_evaluation_points()
            .first()
            .unwrap()
            == final_univariate_poly.evaluate(last_value.expect("Last value is empty"))
    }

    pub fn partial_verify(
        proof: &SumCheckProof<T>,
        transcript: &mut Transcript<T>,
    ) -> (bool, T, Vec<Option<T>>) {
        let mut evaluation_values: Vec<Option<T>> = vec![];
        let mut curr_claimed_sum = proof.initial_claim_sum;

        // This is basically generating all the sampled values e.g(a,b,c)
        // This is done using the same hashing method that the prover used to generate them
        for evaluated_polynomial_over_boolean in &proof.round_polys {
            // from each proof step, get the claimed sum and evaluated polynomial over boolean hypercube.
            // if these two sums don't match, there is no point moving forward
            if evaluated_polynomial_over_boolean.evaluate_sum_over_boolean_hypercube()
                != curr_claimed_sum
            {
                return (false, curr_claimed_sum, evaluation_values);
            }

            transcript.append(&curr_claimed_sum.into_bigint().to_bytes_le());
            transcript.append(&evaluated_polynomial_over_boolean.to_bytes());

            let challenge = transcript.sample_challenge();

            // we then push append the byte equivalent of these values to the transcript and get a challenge which we store in evaluation values array
            evaluation_values.push(Some(challenge));

            // new claim sum is gotten by evaluating the current univariate at the new challenge value
            curr_claimed_sum = evaluated_polynomial_over_boolean.evaluate(challenge);
        }

        (true, curr_claimed_sum, evaluation_values)
    }

    pub fn verify_proof(
        initial_polynomial: &MultiLinearPolynomial<T>,
        proof: SumCheckProof<T>,
    ) -> bool {
        let mut transcript = Transcript::new();

        // append initial polynomial to transcript to initiate process
        transcript.append(&initial_polynomial.to_bytes());

        // check that the initial polynomial evaluated and 0 and 1 is equal to initial claim sum
        if initial_polynomial.evaluation_sum() != proof.initial_claim_sum {
            return false;
        }

        let (partial_verify, _, evaluation_values) = Self::partial_verify(&proof, &mut transcript);

        if !partial_verify {
            return false;
        }

        // if we have a last univariate polynomial variable, perform oracle check, else return false automatically
        let is_correct = match proof.round_polys.last() {
            Some(last_univariate_poly) => Self::perform_oracle_check(
                &initial_polynomial,
                evaluation_values,
                last_univariate_poly,
            ),
            None => false,
        };

        is_correct
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::polynomials::product_polynomial::ProductPolynomial;
    use ark_bn254::Fq;

    #[test]
    fn test_full_sumcheck_pass() {
        let polynomial = vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(5),
        ];

        let initial_polynomial = MultiLinearPolynomial::new(polynomial);

        let sum_check_proof = Prover::generate_sumcheck_proof(&initial_polynomial);

        assert!(Verifier::verify_proof(&initial_polynomial, sum_check_proof));
    }

    #[test]
    fn test_full_sumcheck_fail() {
        let polynomial = vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(65),
        ];

        let sum_check_proof = SumCheckProof {
            initial_claim_sum: Fq::from(10),
            round_polys: vec![
                UnivariatePolynomial::new(vec![Fq::from(3), Fq::from(7)]),
                UnivariatePolynomial::new(vec![Fq::from(9), Fq::from(10)]),
                UnivariatePolynomial::new(vec![Fq::from(10), Fq::from(97)]),
            ],
        };

        assert_eq!(
            Verifier::verify_proof(&MultiLinearPolynomial::new(polynomial), sum_check_proof),
            false
        );
    }

    #[test]
    fn test_partial_sumcheck_pass() {
        let (eval_1, eval_2) = (
            vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)],
            vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(3)],
        );

        let initial_polynomial = SumPolynomial::new(vec![
            ProductPolynomial::new(vec![
                MultiLinearPolynomial::new(eval_1.clone()),
                MultiLinearPolynomial::new(eval_2.clone()),
            ]),
            ProductPolynomial::new(vec![
                MultiLinearPolynomial::new(eval_1),
                MultiLinearPolynomial::new(eval_2),
            ]),
        ]);

        let sum_check_proof =
            Prover::generate_proof_for_partial_verify(Fq::from(12), initial_polynomial);

        assert!(Verifier::partial_verify(&sum_check_proof, &mut Transcript::new()).0);
    }
}
