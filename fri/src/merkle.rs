use ark_ff::{BigInteger, PrimeField};
use std::marker::PhantomData;

use fiat_shamir::transcript::{GenericHashFunctionTrait, GenericTranscript};

pub struct MerkleTree<T: PrimeField, F: GenericHashFunctionTrait> {
    _marker1: PhantomData<T>,
    _marker2: PhantomData<F>,
    hash_layers: Vec<Vec<Vec<u8>>>,
}

#[derive(Debug)]
pub struct MerkleProof {
    hash_path: Vec<Vec<u8>>,
}

impl MerkleProof {
    pub fn new(hash_path: Vec<Vec<u8>>) -> Self {
        Self { hash_path }
    }
}

impl<T: PrimeField, F: GenericHashFunctionTrait> MerkleTree<T, F> {
    pub fn new() -> Self {
        Self {
            _marker1: PhantomData,
            _marker2: PhantomData,
            hash_layers: Vec::new(),
        }
    }

    fn get_hash_partner(&self, hash_index: usize, layer_idx: usize) -> Vec<u8> {
        if hash_index % 2 == 0 {
            self.hash_layers[layer_idx][hash_index + 1].to_vec()
        } else {
            self.hash_layers[layer_idx][hash_index - 1].to_vec()
        }
    }

    fn get_layer_indexes_for_proof_partitions(
        &self,
        index_to_prove: usize,
        indexes_length: usize,
    ) -> Vec<usize> {
        let mut running_index = index_to_prove;
        let mut layer_indexes: Vec<usize> = vec![index_to_prove];

        for _i in 0..indexes_length {
            running_index = running_index / 2;
            layer_indexes.push(running_index);
        }

        layer_indexes
    }

    pub fn get_proof(&self, index_to_prove: usize) -> MerkleProof {
        let mut hash_path = Vec::new();
        let hash_path_length = self.hash_layers.len() - 2;
        let proof_partition_indexes =
            self.get_layer_indexes_for_proof_partitions(index_to_prove, hash_path_length);

        for layer_idx in 0..hash_path_length + 1 {
            hash_path.push(self.get_hash_partner(proof_partition_indexes[layer_idx], layer_idx));
        }

        MerkleProof::new(hash_path)
    }

    pub fn build(&mut self, inputs: &[T], transcript: &mut GenericTranscript<T, F>) -> Vec<u8> {
        let mut current_layer = Vec::from(inputs);
        let input_len = current_layer.len();

        if input_len.next_power_of_two() != input_len {
            // pad input layer with 1's if input length is not a power of 2
            let rem_length = input_len.next_power_of_two() - input_len;
            current_layer.append(&mut vec![T::one(); rem_length]);
        }

        // hash initial layer
        let mut current_hashed_layer = current_layer
            .iter()
            .map(|val| transcript.get_hash(&val.into_bigint().to_bytes_le()))
            .collect::<Vec<_>>();

        self.hash_layers.push(current_hashed_layer.clone());

        while current_hashed_layer.len() > 1 {
            let mut next_hashed_layer = Vec::with_capacity(current_hashed_layer.len() / 2);

            let mut i = 0;

            while i < current_hashed_layer.len() {
                let mut combined_data_to_hash = current_hashed_layer[i].to_vec();
                let mut second_data_to_hash = current_hashed_layer[i + 1].to_vec();

                combined_data_to_hash.append(&mut second_data_to_hash);

                next_hashed_layer.push(transcript.get_hash(&combined_data_to_hash));

                i += 2;
            }

            self.hash_layers.push(next_hashed_layer.to_vec());

            current_hashed_layer = next_hashed_layer;
        }

        current_hashed_layer[0].to_vec()
    }

    pub fn verify_proof(
        &mut self,
        value: &T,
        index_of_value: usize,
        proof: MerkleProof,
        transcript: &mut GenericTranscript<T, F>,
    ) -> bool {
        let root_hash = &self.hash_layers[self.hash_layers.len() - 1][0];
        let hashed_value = transcript.get_hash(&value.into_bigint().to_bytes_le());
        let proof_partition_indexes =
            self.get_layer_indexes_for_proof_partitions(index_of_value, proof.hash_path.len());

        let mut running_hash = hashed_value;

        for (hash_idx, hash) in proof.hash_path.iter().enumerate() {
            let partition_hash_idx = proof_partition_indexes[hash_idx];
            let mut hash_1 = running_hash.to_vec();
            let mut hash_2 = hash.to_vec();

            if partition_hash_idx % 2 == 0 {
                hash_1.append(&mut hash_2);
                running_hash = transcript.get_hash(&hash_1);
            } else {
                hash_2.append(&mut hash_1);
                running_hash = transcript.get_hash(&hash_2);
            }
        }

        root_hash == &running_hash
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use ark_bn254::Fq;
    use sha3::digest::core_api::CoreWrapper;
    use sha3::{Digest, Keccak256, Keccak256Core};

    use fiat_shamir::transcript::GenericTranscript;

    fn get_merkle_tree(values: &[Fq]) -> MerkleTree<Fq, CoreWrapper<Keccak256Core>> {
        let mut merkle_tree: MerkleTree<Fq, CoreWrapper<Keccak256Core>> = MerkleTree::new();

        merkle_tree.build(values, &mut GenericTranscript::new(Keccak256::new()));

        merkle_tree
    }

    #[test]
    pub fn test_verify_merkle_tree() {
        let mut merkle_tree = get_merkle_tree(&[
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(5),
            Fq::from(6),
            Fq::from(7),
            Fq::from(8),
        ]);

        let proof_for_5 = merkle_tree.get_proof(4);

        assert!(merkle_tree.verify_proof(
            &Fq::from(5),
            4,
            proof_for_5,
            &mut GenericTranscript::new(Keccak256::new()),
        ))
    }
}
