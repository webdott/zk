use ark_ff::{BigInteger, PrimeField};
use std::marker::PhantomData;

use fiat_shamir::transcript::{GenericHashFunctionTrait, GenericTranscript};

pub struct MerkleTree<T: PrimeField, F: GenericHashFunctionTrait> {
    _marker: PhantomData<F>,
    hash_layers: Vec<Vec<T>>,
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
        current_layer = current_layer
            .iter()
            .map(|val| transcript.get_hash(&val.into_bigint().to_bytes_le()))
            .collect::<Vec<_>>();

        while current_layer.len() > 1 {
            let mut next_layer = Vec::with_capacity(current_layer.len() / 2);

            for mut i in 0..current_layer.len() / 2 {
                next_layer.push(
                    transcript.get_hash(
                        &(current_layer[i] * current_layer[i + 1])
                            .into_bigint()
                            .to_bytes_le(),
                    ),
                );

                i += 1;
            }

            self.hash_layers.push(next_layer.to_vec());

            current_layer = next_layer;
        }

        println!("{}", self.hash_layers);
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
            _marker: PhantomData,
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
