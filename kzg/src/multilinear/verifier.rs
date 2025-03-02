use crate::multilinear::prover::MultilinearKZGProof;

use ark_ec::pairing::Pairing;
use ark_ec::PrimeGroup;
use ark_ff::PrimeField;
use std::marker::PhantomData;

pub struct MultilinearKZGVerifier<T: PrimeField, P: Pairing> {
    _marker: PhantomData<T>,
    _marker_2: PhantomData<P>,
}

impl<T: PrimeField, P: Pairing> MultilinearKZGVerifier<T, P> {
    pub fn verify_proof(
        proof: MultilinearKZGProof<T, P>,
        openings: &[T],
        encrypted_taus: &[P::G2],
    ) -> bool {
        // The verifier checks that the LHS and RHS (that contains the proof) are equal
        //        L.H.S                     R.H.S
        //[   g1   ]   [g2]    [       g2        ]   [  g1  ]
        //f_tau - v  *  1  === ∑((tau_i - opening_i) * q_tau_i)

        let g1_v = P::G1::generator().mul_bigint(proof.v.into_bigint());
        let f_tau_minus_v: P::G1 = proof.commitment - g1_v;
        let g2_1 = P::G2::generator().mul_bigint(T::one().into_bigint());

        let lhs = P::pairing(f_tau_minus_v, g2_1);

        let rhs_calc = proof
            .q_taus
            .iter()
            .enumerate()
            .map(|(idx, q_tau_i)| {
                let tau_i = encrypted_taus[idx];
                let g2_a = P::G2::generator().mul_bigint(openings[idx].into_bigint());

                P::pairing(*q_tau_i, tau_i - g2_a)
            })
            .collect::<Vec<_>>();

        // ∑((tau_i - opening_i) * q_tau_i)
        let rhs = rhs_calc.iter().sum();

        lhs == rhs
    }
}
