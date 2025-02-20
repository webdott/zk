use crate::concepts_protocols::fiat_shamir::transcript::Transcript;
use crate::concepts_protocols::sumcheck::sumcheck_protocol::SumCheckProof;
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

pub struct SumcheckProver<T: PrimeField> {
    _marker: PhantomData<T>,
}

impl<T: PrimeField> SumcheckProver<T> {
    // this generates a random challenge and returns points to partially evaluate a polynomial at.
    fn generate_evaluation_points(
        transcript: &mut Transcript<T>,
        variables_length: usize,
    ) -> (T, Vec<Option<T>>) {
        let sample_challenge = transcript.sample_challenge();

        (
            sample_challenge,
            iter::repeat(())
                .enumerate()
                .map(|(idx, _x)| {
                    if idx == 0 {
                        return Some(sample_challenge);
                    }

                    return None;
                })
                .take(variables_length)
                .collect(),
        )
    }

    fn generate_round_polys(
        initial_polynomial: &ComposedPolynomial<T>,
        transcript: &mut Transcript<T>,
    ) -> (Vec<UnivariatePolynomial<T>>, Vec<T>) {
        let (
            mut resulting_multi_polynomial,
            mut resulting_sum_polynomial,
            mut round_polys,
            mut random_challenges,
            number_of_variables,
        ) = match initial_polynomial {
            ComposedPolynomial::SumPolynomial(polynomial) => (
                None,
                Some(polynomial.clone()),
                Vec::with_capacity(polynomial.number_of_variables() as usize),
                Vec::with_capacity(polynomial.number_of_variables() as usize),
                polynomial.number_of_variables(),
            ),
            ComposedPolynomial::MultilinearPolynomial(polynomial) => (
                Some(polynomial.clone()),
                None,
                Vec::with_capacity(polynomial.number_of_variables() as usize),
                Vec::with_capacity(polynomial.number_of_variables() as usize),
                polynomial.number_of_variables(),
            ),
        };

        // The steps for generating the univariate round polys differ based on the type of initial polynomial
        // => In the case of a sum polynomial,
        //    - We get the degree of the sum polynomial and partially evaluate the variable of concern at d+1 points.
        //    - We then reduce and sum at each step to get single evaluation points at which we interpolate at to get a univariate.

        // => In the case of a regular Multilinear poly,
        //    - We keep the variable of concern constant, and we then we partially evaluate the other variables over the boolean hypercube.
        //    - Then evaluate the f(variable) at d+1 points and sum up the values.
        //    - The trick while using multilinear polynomials is that you can half the array and sum up evaluation points.

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

                claimed_sum = evaluation_points[0] + evaluation_points[1];

                evaluated_polynomial_over_boolean_hypercube = UnivariatePolynomial::interpolate(
                    &vec![T::from(0), T::from(1), T::from(2)],
                    &evaluation_points,
                );
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
                    &vec![T::from(0), T::from(1)],
                    &vec![eval_0, eval_1],
                )
            }

            transcript.append(&claimed_sum.into_bigint().to_bytes_le());
            transcript.append(&evaluated_polynomial_over_boolean_hypercube.to_bytes());

            if let Some(sum_poly) = &resulting_sum_polynomial {
                let (challenge, points) = Self::generate_evaluation_points(
                    transcript,
                    sum_poly.number_of_variables() as usize,
                );

                random_challenges.push(challenge);

                resulting_sum_polynomial = Some(sum_poly.partial_evaluate(&points));
            } else if let Some(multi_poly) = &resulting_multi_polynomial {
                let (challenge, points) = Self::generate_evaluation_points(
                    transcript,
                    multi_poly.number_of_variables() as usize,
                );

                random_challenges.push(challenge);

                resulting_multi_polynomial = Some(multi_poly.evaluate(&points));
            }

            round_polys.push(evaluated_polynomial_over_boolean_hypercube);
        });

        (round_polys, random_challenges)
    }

    // This creates a sum check proof, with the round_polys generated and an initial claim sum
    pub fn generate_sumcheck_proof(init_polynomial: &MultiLinearPolynomial<T>) -> SumCheckProof<T> {
        let mut transcript = Transcript::new();

        // append initial polynomial to transcript to initiate process
        transcript.append(&init_polynomial.to_bytes());

        let (round_polys, _) = Self::generate_round_polys(
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
        transcript: &mut Transcript<T>,
    ) -> (SumCheckProof<T>, Vec<T>) {
        let (round_polys, random_points) =
            Self::generate_round_polys(&ComposedPolynomial::SumPolynomial(init_poly), transcript);

        (
            SumCheckProof {
                initial_claim_sum,
                round_polys,
            },
            random_points,
        )
    }
}
