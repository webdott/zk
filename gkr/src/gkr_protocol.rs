use kzg::multilinear::prover::MultilinearKZGProof;
use polynomials::multilinear_polynomial::evaluation_form::MultiLinearPolynomial;
use sumcheck::sumcheck_protocol::SumCheckProof;

use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;

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

pub struct GKRProofWithKZG<T: PrimeField, P: Pairing> {
    pub commitment: P::G1,
    pub output_poly: MultiLinearPolynomial<T>,
    pub w_polys_evals: Vec<(T, T)>,
    pub sumcheck_proofs: Vec<SumCheckProof<T>>,
    pub kzg_proofs: Vec<MultilinearKZGProof<T, P>>,
}

impl<T: PrimeField, P: Pairing> GKRProofWithKZG<T, P> {
    pub fn new(
        commitment: P::G1,
        output_poly: MultiLinearPolynomial<T>,
        w_polys_evals: Vec<(T, T)>,
        sumcheck_proofs: Vec<SumCheckProof<T>>,
        kzg_proofs: Vec<MultilinearKZGProof<T, P>>,
    ) -> Self {
        Self {
            commitment,
            output_poly,
            w_polys_evals,
            sumcheck_proofs,
            kzg_proofs,
        }
    }
}

#[cfg(test)]
mod tests {
    use arithmetic_circuit::circuit::Circuit;
    use arithmetic_circuit::gate::{Gate, Operation};
    use fiat_shamir::transcript::Transcript;
    use kzg::multilinear::trusted_setup::TrustedSetup;

    use crate::prover::GKRProver;
    use crate::verifier::GKRVerifier;

    use ark_bls12_381::{Bls12_381, Fr};
    use ark_bn254::Fq;

    pub fn get_test_circuit_and_inputs_fq() -> (Circuit<Fq>, Vec<Fq>) {
        let circuit = Circuit::new(vec![
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

        (circuit, inputs)
    }

    pub fn get_test_circuit_and_inputs_fr() -> (Circuit<Fr>, Vec<Fr>) {
        let circuit = Circuit::new(vec![
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
            Fr::from(1),
            Fr::from(2),
            Fr::from(3),
            Fr::from(4),
            Fr::from(5),
            Fr::from(6),
            Fr::from(7),
            Fr::from(8),
        ];

        (circuit, inputs)
    }

    #[test]
    pub fn test_gkr_sum_check() {
        let (mut circuit, inputs) = get_test_circuit_and_inputs_fq();

        // initialize both GKR prover and verifier with the same circuit
        let gkr_proof = GKRProver::<Fq, Bls12_381>::generate_proof(
            &mut circuit,
            &mut Transcript::new(),
            &inputs,
        );

        assert!(GKRVerifier::<Fq, Bls12_381>::verify_proof(
            &inputs,
            &mut circuit,
            &mut Transcript::new(),
            gkr_proof
        ))
    }

    #[test]
    pub fn test_gkr_sumcheck_with_kzg() {
        let (mut circuit, inputs) = get_test_circuit_and_inputs_fr();

        let trusted_setup: TrustedSetup<Fr, Bls12_381> =
            TrustedSetup::new(&[Fr::from(5), Fr::from(2), Fr::from(3)]);

        // initialize both GKR prover and verifier with the same circuit
        let gkr_proof_with_kzg = GKRProver::<Fr, Bls12_381>::generate_proof_with_kzg(
            &mut circuit,
            &mut Transcript::new(),
            &inputs,
            &trusted_setup.encrypted_lagrange_basis,
        );

        assert!(GKRVerifier::verify_proof_with_kzg(
            &mut circuit,
            &mut Transcript::new(),
            gkr_proof_with_kzg,
            &trusted_setup.encrypted_taus
        ))
    }
}
