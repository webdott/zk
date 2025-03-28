use ark_ff::{BigInteger, PrimeField};
use fiat_shamir::transcript::{GenericHashFunctionTrait, GenericTranscript};

pub struct MerkleTree<T: PrimeField, F: GenericHashFunctionTrait> {
    hash_layers: Vec<Vec<T>>,
}

pub struct MerkleProof {}

impl<T: PrimeField, F: GenericHashFunctionTrait> MerkleTree<T, F> {
    pub fn getProof() -> MerkleProof {
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

        while current_layer.len() > 1 {
            let hashed_layer = current_layer
                .iter()
                .map(|val| {
                    let hashed_val = transcript.get_hash(&val.into_bigint().to_bytes_le());
                    hashed_val
                })
                .collect::<Vec<_>>();

            self.hash_layers.push(hashed_layer.to_vec());

            current_layer = hashed_layer;
        }

        todo!()
    }

    pub fn verifyProof() -> bool {
        todo!()
    }
}
