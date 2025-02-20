use polynomials::multilinear_polynomial::evaluation_form::MultiLinearPolynomial;
use sumcheck::sumcheck_protocol::SumCheckProof;

use ark_ff::PrimeField;
use std::fmt::Debug;

#[derive(Debug)]
pub struct GKRProof<T: PrimeField> {
    pub output_poly: MultiLinearPolynomial<T>,
    pub w_polys_evals: Vec<(T, T)>,
    pub sumcheck_proofs: Vec<SumCheckProof<T>>,
}

impl<T: PrimeField> GKRProof<T> {
    pub fn new(
        output_poly: MultiLinearPolynomial<T>,
        w_polys_evals: Vec<(T, T)>,
        sumcheck_proofs: Vec<SumCheckProof<T>>,
    ) -> Self {
        Self {
            output_poly,
            w_polys_evals,
            sumcheck_proofs,
        }
    }
}

#[cfg(test)]
mod tests {
    use arithmetic_circuit::circuit::Circuit;
    use arithmetic_circuit::gate::{Gate, Operation};
    use fiat_shamir::transcript::Transcript;

    use crate::prover::GKRProver;
    use crate::verifier::GKRVerifier;

    use ark_bn254::Fq;

    #[test]
    pub fn test_gkr_sum_check() {
        let mut circuit = Circuit::new(vec![
            vec![
                Gate::new(0, 1, Operation::Add),
                Gate::new(2, 3, Operation::Add),
                Gate::new(4, 5, Operation::Add),
                Gate::new(6, 7, Operation::Mul),
            ],
            vec![
                Gate::new(0, 1, Operation::Mul),
                Gate::new(2, 3, Operation::Add),
            ],
            vec![Gate::new(0, 1, Operation::Add)],
        ]);

        let inputs = vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(5),
            Fq::from(6),
            Fq::from(7),
            Fq::from(8),
        ];

        // initialize both GKR prover and verifier with the same circuit
        let gkr_proof = GKRProver::generate_proof(&mut circuit, &mut Transcript::new(), &inputs);

        assert!(GKRVerifier::verify_proof(
            &inputs,
            &mut circuit,
            &mut Transcript::new(),
            gkr_proof
        ))
    }
}
