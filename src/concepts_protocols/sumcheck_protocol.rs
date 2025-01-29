use crate::concepts_protocols::fiat_shamier::transcript::Transcript;
use crate::polynomials::multilinear_polynomial::MultiLinearPolynomial;
use ark_ff::{BigInteger, PrimeField};
use std::iter;

pub struct SumCheckProof<T: PrimeField> {
    initial_claim_sum: T,
    round_steps: Vec<(T, MultiLinearPolynomial<T>)>,
}

pub struct Prover<T: PrimeField> {
    initial_polynomial: MultiLinearPolynomial<T>,
}

impl<T: PrimeField> Prover<T> {
    pub fn new(poly: MultiLinearPolynomial<T>) -> Self {
        Prover {
            initial_polynomial: poly,
        }
    }

    // this generates a set of points to partially evaluate a polynomial
    fn generate_evaluation_points(
        &self,
        transcript: &Transcript<T>,
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

    fn generate_steps(&self) -> Vec<(T, MultiLinearPolynomial<T>)> {
        let mut steps: Vec<(T, MultiLinearPolynomial<T>)> =
            Vec::with_capacity(self.initial_polynomial.get_evaluation_points().len());

        let mut transcript = Transcript::new();

        let mut resulting_polynomial = self.initial_polynomial.clone();

        // keep adding current polynomial step and sum to the steps vec
        (0..self.initial_polynomial.number_of_variables()).for_each(|_i| {
            let (mut eval_0, mut eval_1) = (T::from(0), T::from(0));
            let evaluation_points = resulting_polynomial.get_evaluation_points();

            let half = evaluation_points.len() / 2;

            (0..half).for_each(|x| {
                eval_0 += &evaluation_points[x];
                eval_1 += &evaluation_points[half + x];
            });

            let claimed_sum = eval_0 + eval_1;
            let evaluated_polynomial_over_boolean_hypercube =
                MultiLinearPolynomial::new(vec![eval_0, eval_1]);

            transcript.append(&claimed_sum.into_bigint().to_bytes_le());
            transcript.append(&evaluated_polynomial_over_boolean_hypercube.to_bytes());

            let points = self.generate_evaluation_points(
                &transcript,
                resulting_polynomial.number_of_variables() as usize,
            );

            resulting_polynomial = resulting_polynomial.evaluate(points);

            steps.push((claimed_sum, evaluated_polynomial_over_boolean_hypercube));
        });

        steps
    }

    // This creates a sum check proof struct that with the steps generated and an initial claim sum
    pub fn generate_sumcheck_proof(&self) -> SumCheckProof<T> {
        let round_steps = self.generate_steps();

        SumCheckProof {
            initial_claim_sum: self.initial_polynomial.evaluation_sum(),
            round_steps,
        }
    }
}
struct Verifier<T: PrimeField> {
    initial_polynomial: MultiLinearPolynomial<T>,
}

impl<T: PrimeField> Verifier<T> {
    pub fn new(poly: MultiLinearPolynomial<T>) -> Self {
        Verifier {
            initial_polynomial: poly,
        }
    }

    // This check ensures that the last univariate evaluated at a variable is equal to the initial polynomial evaluated at all sampled values
    pub fn perform_oracle_check(
        &self,
        evaluation_values: Vec<Option<T>>,
        final_univariate_poly: &MultiLinearPolynomial<T>,
    ) -> bool {
        let last_value = *evaluation_values.last().unwrap();

        self.initial_polynomial
            .evaluate(evaluation_values)
            .get_evaluation_points()
            == final_univariate_poly
                .evaluate(vec![last_value])
                .get_evaluation_points()
    }

    pub fn verify_proof(&self, proof: SumCheckProof<T>) -> bool {
        let mut evaluation_values: Vec<Option<T>> = vec![];
        let mut transcript = Transcript::new();

        // check that the initial polynomial evaluated and 0 and 1 is equal to initial claim sum
        if self.initial_polynomial.evaluation_sum() != proof.initial_claim_sum {
            return false;
        }

        // This is basically generating all the sampled values e.g(a,b,c)
        // This is done using the same hashing method that the prover used to generate them
        for step in &proof.round_steps {
            let (claimed_sum, evaluated_polynomial_over_boolean) = step;

            // from each proof step, get the claimed sum and evaluated polynomial over boolean hypercube.
            // if these two sums don't match, there is no point moving forward
            if evaluated_polynomial_over_boolean.evaluation_sum() != *claimed_sum {
                return false;
            }

            transcript.append(&claimed_sum.into_bigint().to_bytes_le());
            transcript.append(&evaluated_polynomial_over_boolean.to_bytes());

            // we then push append the byte equivalent of these values to the transcript and get a challenge which we store in evaluation values array
            evaluation_values.push(Some(transcript.sample_challenge()));
        }

        // if we have a last univariate polynomial variable, perform oracle check, else return false automatically
        if let Some((_a, last_univariate_poly)) = proof.round_steps.last() {
            self.perform_oracle_check(evaluation_values, last_univariate_poly)
        } else {
            false
        }
    }
}

#[cfg(test)]
mod test {
    use crate::concepts_protocols::sumcheck_protocol::{Prover, Verifier};
    use crate::polynomials::multilinear_polynomial::MultiLinearPolynomial;
    use ark_bn254::Fq;

    #[test]
    fn test_hash() {
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
        let prover = Prover::new(MultiLinearPolynomial::new(polynomial.clone()));

        let sum_check_proof = prover.generate_sumcheck_proof();

        let verifier = Verifier::new(MultiLinearPolynomial::new(polynomial));

        assert!(verifier.verify_proof(sum_check_proof));
    }
}
