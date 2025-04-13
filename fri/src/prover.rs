use crate::merkle::{MerkleProof, MerkleTree};
use crate::utils::{fold_layer, get_layer_proof_indexes};

use fft::fft::FFT;
use fiat_shamir::transcript::{GenericHashFunctionTrait, GenericTranscript};
use polynomials::univariate_polynomial::dense_coefficient_form::UnivariatePolynomial;

use ark_ff::{FftField, PrimeField};
use std::marker::PhantomData;

#[derive(Debug)]
pub struct LayerIndexProof<T: FftField + PrimeField> {
    pub value: T,
    pub index: usize,
    pub proof: MerkleProof,
}

impl<T: FftField + PrimeField> LayerIndexProof<T> {
    pub fn new(value: T, index: usize, proof: MerkleProof) -> Self {
        Self {
            value,
            index,
            proof,
        }
    }
}

pub struct FriProof<T: FftField + PrimeField> {
    pub layer_merkle_roots: Vec<Vec<u8>>,
    pub layer_proofs: Vec<Vec<LayerIndexProof<T>>>,
}

impl<T: FftField + PrimeField> FriProof<T> {
    pub fn new(
        layer_merkle_roots: Vec<Vec<u8>>,
        layer_proofs: Vec<Vec<LayerIndexProof<T>>>,
    ) -> Self {
        Self {
            layer_merkle_roots,
            layer_proofs,
        }
    }
}

pub struct FriProver<T: FftField + PrimeField, F: GenericHashFunctionTrait> {
    _marker: PhantomData<T>,
    _trait: PhantomData<F>,
}

impl<T: FftField + PrimeField, F: GenericHashFunctionTrait> FriProver<T, F> {
    fn get_layer_proofs(
        initial_index: usize,
        merkle_trees: &Vec<MerkleTree<T, F>>,
        all_layer_evaluations: &[Vec<T>],
    ) -> Vec<Vec<LayerIndexProof<T>>> {
        let mut given_layer_index = initial_index;
        let mut layer_proofs: Vec<Vec<LayerIndexProof<T>>> = vec![];

        for layer_idx in 0..merkle_trees.len() {
            let (idx, negative_idx) =
                get_layer_proof_indexes(all_layer_evaluations[layer_idx].len(), given_layer_index);

            layer_proofs.push(vec![
                LayerIndexProof::new(
                    all_layer_evaluations[layer_idx][idx],
                    idx,
                    merkle_trees[layer_idx].get_proof(idx),
                ),
                LayerIndexProof::new(
                    all_layer_evaluations[layer_idx][negative_idx],
                    negative_idx,
                    merkle_trees[layer_idx].get_proof(negative_idx),
                ),
            ]);

            given_layer_index = idx;
        }

        layer_proofs
    }

    pub fn generate_proof(
        blown_up_coded_word: &[T],
        commit_transcript: &mut GenericTranscript<T, F>,
        merkle_transcript: &mut GenericTranscript<T, F>,
    ) -> (UnivariatePolynomial<T>, FriProof<T>) {
        let blown_up_length = blown_up_coded_word.len();
        let num_of_layers = (blown_up_length as i32).ilog2() as usize;

        let mut layer_root_hashes: Vec<Vec<u8>> = Vec::new();
        let mut layer_evaluations = blown_up_coded_word.to_vec();
        let mut all_layer_evaluations: Vec<Vec<T>> = Vec::with_capacity(num_of_layers);
        let mut merkle_trees: Vec<MerkleTree<T, F>> = Vec::with_capacity(num_of_layers);

        for layer_idx in 0..num_of_layers + 1 {
            let mut merkle_tree: MerkleTree<T, F> = MerkleTree::new();
            let root_hash = merkle_tree.build(&layer_evaluations, merkle_transcript);

            commit_transcript.append(&root_hash);
            let r = commit_transcript.generate_challenge();

            merkle_trees.push(merkle_tree);
            layer_root_hashes.push(root_hash);
            all_layer_evaluations.push(layer_evaluations.to_vec());

            if layer_idx < num_of_layers {
                layer_evaluations = fold_layer(&layer_evaluations, r);
            }
        }

        let initial_random_index = (*commit_transcript
            .generate_challenge()
            .into_bigint()
            .as_ref()
            .first()
            .unwrap() as usize)
            % blown_up_length;

        let layer_proofs =
            Self::get_layer_proofs(initial_random_index, &merkle_trees, &all_layer_evaluations);

        (
            UnivariatePolynomial::new(FFT::convert_to_coefficents(&layer_evaluations)),
            FriProof::new(layer_root_hashes, layer_proofs),
        )
    }
}
