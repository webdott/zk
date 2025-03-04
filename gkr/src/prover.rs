use arithmetic_circuit::circuit::Circuit;
use fiat_shamir::transcript::Transcript;
use polynomials::multilinear_polynomial::evaluation_form::MultiLinearPolynomial;
use polynomials::product_polynomial::ProductPolynomial;
use polynomials::sum_polynomial::SumPolynomial;
use sumcheck::prover::SumcheckProver;

use crate::gkr_protocol::GKRProof;
use crate::utils::{get_evaluated_muli_addi_at_a, get_folded_claim_sum, get_folded_polys};

use ark_ff::{BigInteger, PrimeField};
use std::marker::PhantomData;

pub struct GKRProver<T: PrimeField> {
    _marker: PhantomData<T>,
}

impl<T: PrimeField> GKRProver<T> {
    pub fn generate_proof(
        circuit: &mut Circuit<T>,
        transcript: &mut Transcript<T>,
        inputs: &[T],
    ) -> GKRProof<T> {
        // Evaluate the polynomials at the inputs to be able to get w_polys on each layer
        let circuit_evaluations = circuit.evaluate_at_input(Vec::from(inputs));

        // Initialize vecs to store the w_poly_evaluations and the sumcheck proofs at each step
        let (mut w_polys_evals, mut sum_check_proofs) = (
            Vec::with_capacity(circuit.get_layer_count()),
            Vec::with_capacity(circuit.get_layer_count()),
        );

        let length_of_rs = circuit
            .get_w_i(0, &circuit_evaluations)
            .number_of_variables();

        // This variable stores the w_poly for each layer
        let mut running_layer_polynomial = circuit.get_w_i(0, &circuit_evaluations);

        // Commit to the output layer poly by appending to the transcript
        transcript.append(&running_layer_polynomial.to_bytes());

        // generate a number of rs for to evaluate the output layer depending on the number of outputs of the circuit.
        let mut random_values: Vec<Option<T>> = transcript
            .sample_n_challenges(length_of_rs as usize)
            .into_iter()
            .map(|challenge| Some(challenge))
            .collect();

        for layer_idx in 0..circuit.get_layer_count() {
            // Prover is sending the verifier the following at each step:
            //   - W_output poly of the first layer, t
            //   - The evaluations of W_poly of the subsequent layers -> Perform alpha beta folding if more than one output to form one output poly
            //   - Intermediate Sumcheck proof

            let (muli_a_b_c, addi_a_b_c) =
                (circuit.get_mul_i(layer_idx), circuit.get_add_i(layer_idx));

            let (claim_sum, new_muli_b_c, new_addi_b_c) = match layer_idx {
                0 => {
                    let (muli_b_c, addi_b_c) =
                        get_evaluated_muli_addi_at_a(muli_a_b_c, addi_a_b_c, &random_values);

                    (
                        running_layer_polynomial
                            .evaluate(&random_values)
                            .get_evaluation_points()
                            .first()
                            .unwrap()
                            .clone(),
                        muli_b_c,
                        addi_b_c,
                    )
                }
                // perform alpha-beta folding on W poly
                _ => {
                    let (r_b, r_c) = (
                        &random_values[0..random_values.len() / 2],
                        &random_values[random_values.len() / 2..],
                    );

                    let evaluated_running_b_poly = running_layer_polynomial.evaluate(r_b);
                    let evaluated_running_c_poly = running_layer_polynomial.evaluate(r_c);

                    let (w_i_b_eval, w_i_c_eval) = (
                        evaluated_running_b_poly
                            .get_evaluation_points()
                            .first()
                            .unwrap(),
                        evaluated_running_c_poly
                            .get_evaluation_points()
                            .first()
                            .unwrap(),
                    );

                    // commit w's evaluated at rb and rc
                    transcript.append_n(&[
                        &w_i_b_eval.into_bigint().to_bytes_le(),
                        &w_i_c_eval.into_bigint().to_bytes_le(),
                    ]);

                    let (alpha, beta) =
                        (transcript.sample_challenge(), transcript.sample_challenge());

                    //  Get new claim sums, addi and muli polys, alongside evaluations of the current layer's W poly at the random challenges
                    let (new_muli_b_c, new_addi_b_c) =
                        get_folded_polys(&alpha, &beta, muli_a_b_c, addi_a_b_c, r_b, r_c);

                    // append the evaluations to send to the verifier
                    w_polys_evals.push((*w_i_b_eval, *w_i_c_eval));

                    (
                        get_folded_claim_sum(&alpha, &beta, w_i_b_eval, w_i_c_eval),
                        new_muli_b_c,
                        new_addi_b_c,
                    )
                }
            };

            let next_w_i = circuit.get_w_i(layer_idx + 1, &circuit_evaluations);

            // Generate f_b_c -> ( add_i(b, c) * W(b) + W(c) ) + ( mul_i(b, c) * W(b) * W(c) )
            let f_b_c = SumPolynomial::new(vec![
                ProductPolynomial::new(vec![
                    new_muli_b_c,
                    MultiLinearPolynomial::w_mul(&next_w_i, &next_w_i),
                ]),
                ProductPolynomial::new(vec![
                    new_addi_b_c,
                    MultiLinearPolynomial::w_add(&next_w_i, &next_w_i),
                ]),
            ]);

            // Get sumcheck proof and new set of rs to evaluate W and partially evaluate add_i and mul_i at.
            let (sumcheck_proof, random_points) =
                SumcheckProver::generate_proof_for_partial_verify(claim_sum, f_b_c, transcript);

            random_values = random_points.iter().map(|point| Some(*point)).collect();
            running_layer_polynomial = circuit.get_w_i(layer_idx + 1, &circuit_evaluations);

            sum_check_proofs.push(sumcheck_proof);
        }

        GKRProof::new(
            circuit.get_w_i(0, &circuit_evaluations),
            w_polys_evals,
            sum_check_proofs,
        )
    }
}
