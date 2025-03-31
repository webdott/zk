use ark_ff::{BigInteger, PrimeField};
use fiat_shamir::transcript::{GenericHashFunctionTrait, GenericTranscript};
use std::marker::PhantomData;

pub struct MerkleTree<T: PrimeField, F: GenericHashFunctionTrait> {
    _marker1: PhantomData<T>,
    _marker2: PhantomData<F>,
    hash_layers: Vec<Vec<Vec<u8>>>,
}

pub struct MerkleProof {}

impl<T: PrimeField, F: GenericHashFunctionTrait> MerkleTree<T, F> {
    pub fn get_proof() -> MerkleProof {
        todo!()
    }

    pub fn prove(&mut self, inputs: &[T], transcript: &mut GenericTranscript<T, F>) {
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

                next_hashed_layer.push(transcript.get_hash(&second_data_to_hash));

                i += 2;
            }

            self.hash_layers.push(next_hashed_layer.to_vec());

            current_hashed_layer = next_hashed_layer;
        }
    }

    pub fn verify_proof(&mut self, root_hash: &T, proof: MerkleProof) -> bool {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use ark_bn254::Fq;
    use sha3::digest::core_api::CoreWrapper;
    use sha3::{Digest, Keccak256, Keccak256Core};

    use fiat_shamir::transcript::GenericTranscript;

    #[test]
    fn test_merkle() {
        let mut merkle_tree: MerkleTree<Fq, CoreWrapper<Keccak256Core>> = MerkleTree {
            _marker1: PhantomData,
            _marker2: PhantomData,
            hash_layers: vec![],
        };

        merkle_tree.prove(
            &[
                Fq::from(1),
                Fq::from(2),
                Fq::from(3),
                Fq::from(4),
                Fq::from(5),
                Fq::from(6),
                Fq::from(7),
                Fq::from(8),
            ],
            &mut GenericTranscript::new(Keccak256::new()),
        );

        assert!(true)
    }
}
