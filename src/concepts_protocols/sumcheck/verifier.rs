use crate::concepts_protocols::fiat_shamir::transcript::Transcript;
use crate::concepts_protocols::sumcheck::sumcheck_protocol::SumCheckProof;
use crate::polynomials::multilinear_polynomial::MultiLinearPolynomial;
use crate::polynomials::univariate_polynomial::UnivariatePolynomial;
use ark_ff::{BigInteger, PrimeField};
use std::marker::PhantomData;

pub struct SumcheckVerifier<T: PrimeField> {
    _marker: PhantomData<T>,
}

impl<T: PrimeField> SumcheckVerifier<T> {
    // This check ensures that the last univariate evaluated at a variable is equal to the initial polynomial evaluated at all sampled values
    pub fn perform_oracle_check(
        initial_polynomial: &MultiLinearPolynomial<T>,
        evaluation_values: Vec<Option<T>>,
        final_univariate_poly: &UnivariatePolynomial<T>,
    ) -> bool {
        let last_value = *evaluation_values.last().unwrap();

        *initial_polynomial
            .evaluate(&evaluation_values)
            .get_evaluation_points()
            .first()
            .unwrap()
            == final_univariate_poly.evaluate(last_value.expect("Last value is empty"))
    }

    // This bit does the partial verification for a proof minus the oracle check.
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

        let (partially_verified, _, evaluation_values) =
            Self::partial_verify(&proof, &mut transcript);

        if !partially_verified {
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
