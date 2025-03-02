use ark_ec::pairing::Pairing;
use ark_ec::PrimeGroup;
use ark_ff::PrimeField;
use std::marker::PhantomData;
use std::ops::Sub;

pub struct MultilinearKZGVerifier<T: PrimeField, P: Pairing> {
    _marker: PhantomData<T>,
    _marker_2: PhantomData<P>,
}

impl<T: PrimeField, P: Pairing> MultilinearKZGVerifier<T, P> {
    pub fn verify_proof(commitment: P::G1, openings: &[T], v_proof: (T, Vec<P::G1>)) -> bool {
        //f_tau - v = (tau - x) * Qi(tau)

        let (v, q_taus) = v_proof;
        let g1_v = P::G1::generator().mul_bigint(-v.into_bigint());
        let f_tau_minus_v: P::G1 = commitment.sub(g1_v);
        let g2_1 = P::G2::generator().mul_bigint(T::one().into_bigint());

        let lhs: P::TargetField = f_tau_minus_v * g2_1; //?

        let rhs_calc = q_taus
            .iter()
            .enumerate()
            .map(|(idx, q_tau)| {
                // ?
            })
            .collect::<Vec<_>>();

        let rhs = rhs_calc.iter().sum();

        lhs == rhs
    }
}
