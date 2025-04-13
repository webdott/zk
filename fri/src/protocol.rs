#[cfg(test)]
mod tests {
    use crate::prover::FriProver;
    use crate::utils::perform_reed_solomon;
    use crate::verifier::FriVerifier;

    use fiat_shamir::transcript::GenericTranscript;
    use polynomials::univariate_polynomial::dense_coefficient_form::UnivariatePolynomial;

    use ark_bls12_377::Fr;
    use sha3::{Digest, Keccak256};

    #[test]
    pub fn test_fri_protocol() {
        let init_coefficients =
            UnivariatePolynomial::new(vec![Fr::from(5), Fr::from(3), Fr::from(2)]);
        let blown_up_codeword = perform_reed_solomon(init_coefficients, 2);

        let (final_poly, proof) = FriProver::generate_proof(
            &blown_up_codeword,
            &mut GenericTranscript::new(Keccak256::new()),
            &mut GenericTranscript::new(Keccak256::new()),
        );

        assert!(
            FriVerifier::verify(
                proof,
                &final_poly,
                &mut GenericTranscript::new(Keccak256::new()),
                &mut GenericTranscript::new(Keccak256::new()),
            ),
            "Proof verification failed"
        );
    }
}
