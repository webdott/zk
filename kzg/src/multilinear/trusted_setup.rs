use crate::multilinear::utils::{encrypt_lagrange_basis, generate_lagrange_basis_for_n_variables};
use ark_ec::pairing::Pairing;
use ark_ec::PrimeGroup;
use ark_ff::PrimeField;
use std::marker::PhantomData;

pub struct TrustedSetup<T: PrimeField, P: Pairing> {
    _marker: PhantomData<T>,
    pub encrypted_taus: Vec<P::G2>,
    pub encrypted_lagrange_basis: Vec<P::G1>,
}

impl<T: PrimeField, P: Pairing> TrustedSetup<T, P> {
    pub fn new(taus: &[T]) -> Self {
        let encrypted_lagrange_basis =
            encrypt_lagrange_basis::<T, P>(&generate_lagrange_basis_for_n_variables(taus));
        let encrypted_taus = taus
            .iter()
            .map(|tau| P::G2::generator().mul_bigint(tau.into_bigint()))
            .collect::<Vec<_>>();

        Self {
            _marker: PhantomData,
            encrypted_taus,
            encrypted_lagrange_basis,
        }
    }
}
