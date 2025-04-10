use crate::merkle::{MerkleProof, MerkleTree};
use crate::utils::{fold_layer, pad_polynomial};

use ark_bn254::Fq;
use ark_ff::FftField;
use sha3::digest::core_api::CoreWrapper;
use std::marker::PhantomData;

use fft::fft::Polynomial;
use fiat_shamir::transcript::GenericTranscript;
use polynomials::univariate_polynomial::dense_coefficient_form::UnivariatePolynomial;

pub struct LayerIndexProof<T: FftField> {
    value: T,
    index: usize,
    proof: MerkleProof,
}

pub struct FriProof<T: FftField> {
    layer_merkle_roots: Vec<Vec<u8>>,
    layer_proofs: Vec<Vec<LayerIndexProof<T>>>,
}

impl<T: FftField> FriProof<T> {
    pub fn new(layer_merkle_roots: Vec<Vec<u8>>, proofs: Vec<Vec<LayerIndexProof<T>>>) -> Self {
        Self {
            layer_merkle_roots,
            layer_proofs,
        }
    }
}

impl<T: FftField> FriProof<T> {
    pub fn new() -> Self {
        todo!()
    }
}

pub struct FriProver<T: FftField> {
    _marker: PhantomData<T>,
}

impl<T: FftField> FriProver<T> {
    pub fn generate_proof(
        polynomial: UnivariatePolynomial<T>,
        blow_up_factor: usize,
    ) -> (UnivariatePolynomial<T>, FriProof<T>) {
        let mut commit_transcript: GenericTranscript<Fq, CoreWrapper<Keccak256Core>> =
            GenericTranscript::new(Keccak256::new());

        let blown_up_length = polynomial.coefficients.len() * blow_up_factor;
        let padded_polynomial_coefficients = pad_polynomial(
            &polynomial.coefficients,
            blown_up_length.next_power_of_two(),
            T::zero(),
        );
        let mut layer_evaluations =
            Polynomial::convert_to_evaluations(&padded_polynomial_coefficients);

        let num_of_layers = (blown_up_length as i32).ilog2() as usize;
        let layer_proofs: Vec<Vec<LayerIndexProof<T>>> = Vec::new();
        let mut layer_root_hashes: Vec<Vec<u8>> = Vec::new();

        for _ in 0..num_of_layers + 1 {
            let mut merkle_tree: MerkleTree<Fq, CoreWrapper<Keccak256Core>> = MerkleTree::new();
            let root_hash = merkle_tree.build(
                &layer_evaluations,
                &mut GenericTranscript::new(Keccak256::new()),
            );

            commit_transcript.append(&root_hash);
            let r = commit_transcript.generate_challenge();

            layer_root_hashes.push(root_hash);

            layer_evaluations = fold_layer(&layer_evaluations, &r);
        }

        (
            Polynomial::convert_to_coefficents(&layer_evaluations),
            FriProof::new(layer_root_hashes, layer_proofs),
        )
    }
}
