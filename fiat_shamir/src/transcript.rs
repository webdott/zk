use ark_ff::PrimeField;
use sha3::digest::Update;
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
        Update::update(&mut self.hasher, data);
    }

    pub fn append_n(&mut self, data: &[&[u8]]) {
        data.iter().for_each(|f| self.append(*f));
    }

    pub fn sample_challenge(&mut self) -> T {
        // uses the current hasher and generates a field value from it
        let hash_result = self.hasher.clone().finalize();

        // we append this result back into the hasher to always generate a separate random value subsequently even without calling updating externally
        self.append(&hash_result);

        T::from_le_bytes_mod_order(&hash_result)
    }

    // squeeze multiple number of challenges
    pub fn sample_n_challenges(&mut self, n: usize) -> Vec<T> {
        (0..n).map(|_| self.sample_challenge()).collect()
    }
}

pub struct GenericTranscript<T: PrimeField, F: GenericHashFunctionTrait> {
    _marker: PhantomData<T>,
    hash_function: F,
}

impl<T: PrimeField, F: GenericHashFunctionTrait> GenericTranscript<T, F> {
    pub fn new(hash_function: F) -> Self {
        Self {
            _marker: PhantomData,
            hash_function,
        }
    }

    pub fn append(&mut self, data: &[u8]) {
        self.hash_function.absorb(data);
    }

    pub fn append_n(&mut self, data: &[&[u8]]) {
        data.iter().for_each(|f| self.append(*f));
    }

    pub fn generate_challenge(&mut self) -> T {
        // uses the current hasher and generates a field value from it
        let hash_result = self.hash_function.squeeze();

        // We append this result back into the hasher
        // We do this to always generate a separate random value subsequently even without calling updating externally
        self.append(&hash_result);

        T::from_le_bytes_mod_order(&hash_result)
    }

    // squeeze multiple number of challenges
    pub fn sample_n_challenges(&mut self, n: usize) -> Vec<T> {
        (0..n).map(|_| self.generate_challenge()).collect()
    }

    // This is used for just getting a hash immediately for a given value (use case: for merkle proof)
    pub fn get_hash(&mut self, data: &[u8]) -> Vec<u8> {
        self.append(&data);

        let hash_result = self.hash_function.squeeze();

        // we don't need to retain the value for next update. So we clear
        self.hash_function.empty();

        hash_result
    }
}

pub trait GenericHashFunctionTrait {
    fn absorb(&mut self, data: &[u8]);
    fn squeeze(&self) -> Vec<u8>;
    fn empty(&mut self);
}

impl GenericHashFunctionTrait for Keccak256 {
    fn absorb(&mut self, data: &[u8]) {
        sha3::Digest::update(self, data);
    }

    fn squeeze(&self) -> Vec<u8> {
        self.clone().finalize().to_vec()
    }

    fn empty(&mut self) {
        self.reset();
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bn254::Fq;
    use sha3::digest::core_api::CoreWrapper;
    use sha3::Keccak256Core;

    #[test]
    fn test_hardcoded_transcript() {
        let mut first_transcript: Transcript<Fq> = Transcript::new();

        let mut second_transcript: Transcript<Fq> = Transcript::new();

        first_transcript.append(b"hello");
        first_transcript.append(b"world");

        second_transcript.append(b"hello");
        second_transcript.append(b"world");

        assert_eq!(
            first_transcript.sample_challenge(),
            second_transcript.sample_challenge()
        );
    }

    #[test]
    fn test_generic_transcript() {
        let mut first_transcript: GenericTranscript<Fq, CoreWrapper<Keccak256Core>> =
            GenericTranscript::new(Keccak256::new());

        let mut second_transcript: GenericTranscript<Fq, CoreWrapper<Keccak256Core>> =
            GenericTranscript::new(Keccak256::new());

        first_transcript.append(b"hello");
        first_transcript.append(b"world");

        second_transcript.append(b"hello");
        second_transcript.append(b"world");

        assert_eq!(
            first_transcript.generate_challenge(),
            second_transcript.generate_challenge()
        );
    }
}
