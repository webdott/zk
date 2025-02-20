use arithmetic_circuit::circuit::Circuit;
use fiat_shamir::transcript::Transcript;
use polynomials::multilinear_polynomial::evaluation_form::MultiLinearPolynomial;
use sumcheck::verifier::SumcheckVerifier;

use crate::gkr_protocol::GKRProof;
use crate::utils::{get_evaluated_muli_addi_at_a, get_folded_polys};

use ark_ff::PrimeField;
use std::marker::PhantomData;

pub struct GKRVerifier<T: PrimeField> {
    _marker: PhantomData<T>,
}

impl<T: PrimeField> GKRVerifier<T> {
    pub fn verify_proof(
        initial_inputs: &[T],
        circuit: &mut Circuit<T>,
        transcript: &mut Transcript<T>,
        proof: GKRProof<T>,
    ) -> bool {
        // performs the same step as prover in output poly
        let length_of_rs = proof.output_poly.number_of_variables();

        transcript.append(&proof.output_poly.to_bytes());

        let mut random_values: Vec<Option<T>> = transcript
            .sample_n_challenges(length_of_rs as usize)
            .into_iter()
            .map(|challenge| Some(challenge))
            .collect();

        for layer_idx in 0..circuit.get_layer_count() {
            let muli_a_b_c = circuit.get_mul_i(layer_idx);
            let addi_a_b_c = circuit.get_add_i(layer_idx);

            let (new_muli_b_c, new_addi_b_c) = match layer_idx {
                0 => get_evaluated_muli_addi_at_a(muli_a_b_c, addi_a_b_c, &random_values),
                _ => {
                    let (alpha, beta) =
                        (transcript.sample_challenge(), transcript.sample_challenge());

                    // Get the new addi's and muli's using alpha beta folding.
                    let (new_muli_b_c, new_addi_b_c) = get_folded_polys(
                        &alpha,
                        &beta,
                        muli_a_b_c,
                        addi_a_b_c,
                        &random_values[0..random_values.len() / 2],
                        &random_values[random_values.len() / 2..],
                    );

                    (new_muli_b_c, new_addi_b_c)
                }
            };

            // Partial verifier checks if partial proof is correct and returns final claim sum and next r values in the process
            let (is_verified, final_claim_sum, next_evaluation_values) =
                SumcheckVerifier::partial_verify(&proof.sumcheck_proofs[layer_idx], transcript);

            // Using the next set of rs gotten from partial prover, we evaluate the new addi's and muli's
            let evaluated_addi_b_c = new_addi_b_c.evaluate(&next_evaluation_values);
            let evaluated_muli_b_c = new_muli_b_c.evaluate(&next_evaluation_values);

            let (new_addi_b_c_eval, new_muli_b_c_eval) = (
                evaluated_addi_b_c.get_evaluation_points().first().unwrap(),
                evaluated_muli_b_c.get_evaluation_points().first().unwrap(),
            );

            // Once we get to the layer before the input, we use the input polynomial instead to build the next_w_i evals,
            let (next_w_i_b_eval, next_w_i_c_eval) = if layer_idx + 1 == circuit.get_layer_count() {
                let (r_b, r_c) = (
                    &next_evaluation_values[0..next_evaluation_values.len() / 2],
                    &next_evaluation_values[next_evaluation_values.len() / 2..],
                );

                let next_w_i = MultiLinearPolynomial::new(&Vec::from(initial_inputs));

                (
                    next_w_i
                        .evaluate(r_b)
                        .get_evaluation_points()
                        .first()
                        .unwrap()
                        .clone(),
                    next_w_i
                        .evaluate(r_c)
                        .get_evaluation_points()
                        .first()
                        .unwrap()
                        .clone(),
                )
            // else use the w_poly evals the prover gives us
            } else {
                proof.w_polys_evals[layer_idx]
            };

            let fbc_eval = (*new_addi_b_c_eval * (next_w_i_b_eval + next_w_i_c_eval))
                + (*new_muli_b_c_eval * (next_w_i_b_eval * next_w_i_c_eval));

            // Now the verifier performs the oracle check not being handled by partial verifier
            // We check if the f_b_c polynomial evaluated at b and c values equal the final claim sum
            if !is_verified || (fbc_eval != final_claim_sum) {
                return false;
            }

            random_values = next_evaluation_values;
        }

        true
    }
}
