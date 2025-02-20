use polynomials::univariate_polynomial::dense_coefficient_form::UnivariatePolynomial;

use ark_ff::PrimeField;

#[derive(Debug)]
pub struct SumCheckProof<T: PrimeField> {
    pub initial_claim_sum: T,
    pub round_polys: Vec<UnivariatePolynomial<T>>,
}

#[cfg(test)]
mod test {
    use super::*;

    use fiat_shamir::transcript::Transcript;
    use polynomials::multilinear_polynomial::evaluation_form::MultiLinearPolynomial;
    use polynomials::product_polynomial::ProductPolynomial;
    use polynomials::sum_polynomial::SumPolynomial;

    use crate::prover::SumcheckProver;
    use crate::verifier::SumcheckVerifier;

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

        let initial_polynomial = MultiLinearPolynomial::new(&polynomial);

        let sum_check_proof = SumcheckProver::generate_sumcheck_proof(&initial_polynomial);

        assert!(SumcheckVerifier::verify_proof(
            &initial_polynomial,
            sum_check_proof
        ));
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
            SumcheckVerifier::verify_proof(
                &MultiLinearPolynomial::new(&polynomial),
                sum_check_proof
            ),
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
                MultiLinearPolynomial::new(&eval_1),
                MultiLinearPolynomial::new(&eval_2),
            ]),
            ProductPolynomial::new(vec![
                MultiLinearPolynomial::new(&eval_1),
                MultiLinearPolynomial::new(&eval_2),
            ]),
        ]);

        let (sum_check_proof, _) = SumcheckProver::generate_proof_for_partial_verify(
            Fq::from(12),
            initial_polynomial,
            &mut Transcript::new(),
        );

        assert!(SumcheckVerifier::partial_verify(&sum_check_proof, &mut Transcript::new()).0);
    }
}
