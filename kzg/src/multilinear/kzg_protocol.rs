#[cfg(test)]
mod tests {
    use crate::multilinear::trusted_setup::TrustedSetup;

    use crate::multilinear::prover::{MultilinearKZGProof, MultilinearKZGProver};
    use crate::multilinear::verifier::MultilinearKZGVerifier;
    use ark_bls12_381::{Bls12_381, Fr};
    use polynomials::multilinear_polynomial::evaluation_form::MultiLinearPolynomial;

    #[test]
    pub fn test_kzg_protocol_pass() {
        let trusted_setup: TrustedSetup<Fr, Bls12_381> =
            TrustedSetup::new(&[Fr::from(5), Fr::from(2), Fr::from(3)]);
        let polynomial = MultiLinearPolynomial::new(&vec![
            Fr::from(0),
            Fr::from(4),
            Fr::from(0),
            Fr::from(4),
            Fr::from(0),
            Fr::from(4),
            Fr::from(3),
            Fr::from(7),
        ]);
        let openings = vec![Fr::from(6), Fr::from(4), Fr::from(0)];

        let proof: MultilinearKZGProof<Fr, Bls12_381> = MultilinearKZGProver::generate_proof(
            &openings,
            &trusted_setup.encrypted_lagrange_basis,
            &polynomial,
        );

        assert!(MultilinearKZGVerifier::verify_proof(
            proof,
            &openings,
            &trusted_setup.encrypted_taus
        ));
    }
}
