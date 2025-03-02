use std::marker::PhantomData;

use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;

pub struct MultilinearKZGVerifier<T: PrimeField, P: Pairing> {
    _marker: PhantomData<T>,
    _marker_2: PhantomData<P>,
}

impl<T: PrimeField, P: Pairing> MultilinearKZGVerifier<T, P> {}
