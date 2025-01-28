use ark_ff::PrimeField;
use sha3::{Digest, Keccak256};
use std::marker::PhantomData;

pub struct Transcript<T: PrimeField> {
    _marker: PhantomData<T>,
    hasher: Keccak256,
}

impl<T: PrimeField> Transcript<T> {
    pub fn new() -> Self {
        Transcript {
            _marker: Default::default(),
            hasher: Keccak256::new(),
        }
    }

    // update current hasher state with new data
    pub fn append(&mut self, data: &[u8]) {
        self.hasher.update(data);
    }

    pub fn sample_challenge(&self) -> T {
        // uses the current hasher and generates a field value from it
        let hash_result = self.hasher.clone().finalize();

        T::from_le_bytes_mod_order(&hash_result)
    }
}
