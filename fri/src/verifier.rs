use crate::merkle::MerkleTree;
use crate::prover::FriProof;
use crate::utils::{compute_f_x_squared, get_f_squared_from_folded_layer};

use fiat_shamir::transcript::{GenericHashFunctionTrait, GenericTranscript};
use polynomials::univariate_polynomial::dense_coefficient_form::UnivariatePolynomial;

use ark_ff::{FftField, PrimeField};
use std::marker::PhantomData;

pub struct FriVerifier<T: FftField + PrimeField, F: GenericHashFunctionTrait> {
    _marker: PhantomData<T>,
    _trait: PhantomData<F>,
}

impl<T: FftField + PrimeField, F: GenericHashFunctionTrait> FriVerifier<T, F> {
    fn verify_consistency(
        proof: FriProof<T>,
        commit_transcript: &mut GenericTranscript<T, F>,
        merkle_transcript: &mut GenericTranscript<T, F>,
    ) -> bool {
        for (layer_idx, merkle_root) in proof.layer_merkle_roots.iter().enumerate() {
            commit_transcript.append(merkle_root);

            let n = proof.layer_proofs.len();
            let r = commit_transcript.generate_challenge();
            let mut merkle_tree = MerkleTree::new();
            let nth_root = T::get_root_of_unity(1 << (n - layer_idx - 1) as u64);

            let evaluations_part_of_tree =
                proof.layer_proofs[layer_idx]
                    .iter()
                    .fold(true, |a: bool, b| {
                        a || merkle_tree.verify_proof(
                            &b.value,
                            b.index,
                            &b.proof,
                            &merkle_root,
                            merkle_transcript,
                        )
                    });

            if !evaluations_part_of_tree {
                return false;
            }

            let positive_index = proof.layer_proofs[layer_idx][0].index;

            if layer_idx < proof.layer_proofs.len() - 1 {
                let f_x_squared = compute_f_x_squared(
                    positive_index,
                    (
                        proof.layer_proofs[layer_idx][0].value,
                        proof.layer_proofs[layer_idx][1].value,
                    ),
                    r,
                    nth_root,
                );

                if f_x_squared
                    != get_f_squared_from_folded_layer(
                        positive_index,
                        &proof.layer_proofs[layer_idx + 1],
                    )
                {
                    return false;
                }
            };
        }

        true
    }

    fn verify_degree(polynomial: &UnivariatePolynomial<T>) -> bool {
        polynomial.coefficients.len() == 1
    }

    pub fn verify(
        proof: FriProof<T>,
        final_polynomial: &UnivariatePolynomial<T>,
        commit_transcript: &mut GenericTranscript<T, F>,
        merkle_transcript: &mut GenericTranscript<T, F>,
    ) -> bool {
        Self::verify_degree(final_polynomial)
            && Self::verify_consistency(proof, commit_transcript, merkle_transcript)
    }
}
