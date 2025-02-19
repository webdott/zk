use crate::concepts_protocols::arithmetic_gate::circuit::Circuit;
use crate::concepts_protocols::fiat_shamir::transcript::Transcript;
use crate::concepts_protocols::sumcheck::prover::SumcheckProver;
use crate::concepts_protocols::sumcheck::sumcheck_protocol::SumCheckProof;
use crate::concepts_protocols::sumcheck::verifier::SumcheckVerifier;
use crate::polynomials::multilinear_polynomial::MultiLinearPolynomial;
use crate::polynomials::product_polynomial::ProductPolynomial;
use crate::polynomials::sum_polynomial::SumPolynomial;
use ark_ff::PrimeField;
use std::fmt::Debug;
use std::marker::PhantomData;

pub fn get_folded_claim_sum<T: PrimeField>(
    w_i_b_eval: &T,
    w_i_c_eval: &T,
    alpha: &T,
    beta: &T,
) -> T {
    // Follows the alpha-beta formula:
    //    - ( alpha * (W(rb)) ) + ( beta * (W(rc)) ) => New claim sum
    (*w_i_b_eval * *alpha) + (*w_i_c_eval * *beta)
}

pub fn get_folded_polys<T: PrimeField>(
    alpha: &T,
    beta: &T,
    muli_a_b_c: MultiLinearPolynomial<T>,
    addi_a_b_c: MultiLinearPolynomial<T>,
    r_b: &[Option<T>],
    r_c: &[Option<T>],
) -> (MultiLinearPolynomial<T>, MultiLinearPolynomial<T>) {
    // Follows the alpha-beta formula:
    //    - ( alpha * muli(rb,b,c) ) + ( beta * muli(rc, b, c) ) => New mul_i poly
    //    - ( alpha * addi(rb,b,c) ) + ( beta * addi(rc, b, c) ) => New add_i poly

    let mut eval_points_rb = vec![None; muli_a_b_c.number_of_variables() as usize];
    let mut eval_points_rc = vec![None; muli_a_b_c.number_of_variables() as usize];

    (0..r_b.len()).for_each(|idx| {
        eval_points_rb[idx] = r_b[idx];
        eval_points_rc[idx] = r_c[idx];
    });

    let new_muli_b_c = muli_a_b_c
        .evaluate(eval_points_rb.clone())
        .scalar_mul(*alpha)
        .add(
            &muli_a_b_c
                .evaluate(eval_points_rc.clone())
                .scalar_mul(*beta),
        );

    let new_addi_b_c = addi_a_b_c
        .evaluate(eval_points_rb)
        .scalar_mul(*alpha)
        .add(&addi_a_b_c.evaluate(eval_points_rc).scalar_mul(*beta));

    (new_muli_b_c, new_addi_b_c)
}

#[derive(Debug)]
pub struct GKRProof<T: PrimeField> {
    pub output_poly: MultiLinearPolynomial<T>,
    pub w_polys_evals: Vec<(T, T)>,
    pub sumcheck_proofs: Vec<SumCheckProof<T>>,
}

impl<T: PrimeField> GKRProof<T> {
    pub fn new(
        output_poly: MultiLinearPolynomial<T>,
        w_polys_evals: Vec<(T, T)>,
        sumcheck_proofs: Vec<SumCheckProof<T>>,
    ) -> Self {
        Self {
            output_poly,
            w_polys_evals,
            sumcheck_proofs,
        }
    }
}

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
        circuit.evaluate_at_input(Vec::from(inputs));

        // Initialize vecs to store the w_poly_evaluations and the sumcheck proofs at each step
        let (mut w_polys_evals, mut sum_check_proofs) = (
            Vec::with_capacity(circuit.get_layer_count()),
            Vec::with_capacity(circuit.get_layer_count()),
        );

        let length_of_rs = circuit.get_w_i(0).number_of_variables();

        // This variable stores the w_poly for each layer
        let mut running_layer_polynomial = circuit.get_w_i(0);

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
            //   - W_output poly of the first layer, then just the evaluations of W_poly of the subsequent layers -> Perform alpha beta folding if more than one output to form one output poly
            //   - Intermediate Sumcheck proof
            // perform alpha-beta folding on W poly

            let (muli_a_b_c, addi_a_b_c) =
                (circuit.get_mul_i(layer_idx), circuit.get_add_i(layer_idx));

            let (claim_sum, new_muli_b_c, new_addi_b_c) = match layer_idx {
                0 => {
                    let mut evaluation_points = random_values.clone();

                    evaluation_points.extend(vec![
                        None;
                        (muli_a_b_c.number_of_variables() as usize)
                            - random_values.len()
                    ]);

                    (
                        running_layer_polynomial
                            .evaluate(random_values)
                            .get_evaluation_points()
                            .first()
                            .unwrap()
                            .clone(),
                        muli_a_b_c.evaluate(evaluation_points.clone()),
                        addi_a_b_c.evaluate(evaluation_points),
                    )
                }
                _ => {
                    let (r_b, r_c) = (
                        &random_values[0..random_values.len() / 2],
                        &random_values[random_values.len() / 2..],
                    );

                    let evaluated_running_b_poly =
                        running_layer_polynomial.evaluate(Vec::from(r_b));
                    let evaluated_running_c_poly =
                        running_layer_polynomial.evaluate(Vec::from(r_c));

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

                    let (alpha, beta) =
                        (transcript.sample_challenge(), transcript.sample_challenge());

                    //  Get new claim sums, addi and muli polys, alongside evaluations of the current layer's W poly at the random challenges
                    let (new_muli_b_c, new_addi_b_c) =
                        get_folded_polys(&alpha, &beta, muli_a_b_c, addi_a_b_c, r_b, r_c);

                    w_polys_evals.push((*w_i_b_eval, *w_i_c_eval));

                    (
                        get_folded_claim_sum(&alpha, &beta, w_i_b_eval, w_i_c_eval),
                        new_muli_b_c,
                        new_addi_b_c,
                    )
                }
            };

            let next_w_i = circuit.get_w_i(layer_idx + 1);

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
            running_layer_polynomial = circuit.get_w_i(layer_idx + 1);

            sum_check_proofs.push(sumcheck_proof);
        }

        GKRProof::new(circuit.get_w_i(0), w_polys_evals, sum_check_proofs)
    }
}

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
                0 => {
                    let mut evaluation_points = random_values.clone();

                    evaluation_points.extend(vec![
                        None;
                        (muli_a_b_c.number_of_variables() as usize)
                            - random_values.len()
                    ]);

                    (
                        muli_a_b_c.evaluate(evaluation_points.clone()),
                        addi_a_b_c.evaluate(evaluation_points),
                    )
                }
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
            let evaluated_addi_b_c = new_addi_b_c.evaluate(next_evaluation_values.clone());
            let evaluated_muli_b_c = new_muli_b_c.evaluate(next_evaluation_values.clone());

            let (new_addi_b_c_eval, new_muli_b_c_eval) = (
                evaluated_addi_b_c.get_evaluation_points().first().unwrap(),
                evaluated_muli_b_c.get_evaluation_points().first().unwrap(),
            );

            // Once we get to the layer before the input, we use the input polynomial instead to build the next_w_i evals, else use the w_poly evals the prover gives us
            let (next_w_i_b_eval, next_w_i_c_eval) = if layer_idx + 1 == circuit.get_layer_count() {
                let (r_b, r_c) = (
                    &next_evaluation_values[0..next_evaluation_values.len() / 2],
                    &next_evaluation_values[next_evaluation_values.len() / 2..],
                );

                let next_w_i = MultiLinearPolynomial::new(Vec::from(initial_inputs));

                (
                    next_w_i
                        .evaluate(Vec::from(r_b))
                        .get_evaluation_points()
                        .first()
                        .unwrap()
                        .clone(),
                    next_w_i
                        .evaluate(Vec::from(r_c))
                        .get_evaluation_points()
                        .first()
                        .unwrap()
                        .clone(),
                )
            } else {
                proof.w_polys_evals[layer_idx]
            };

            let fbc_eval = (*new_addi_b_c_eval * (next_w_i_b_eval + next_w_i_c_eval))
                + (*new_muli_b_c_eval * (next_w_i_b_eval * next_w_i_c_eval));

            if !is_verified || (fbc_eval != final_claim_sum) {
                return false;
            }

            random_values = next_evaluation_values;
        }

        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::concepts_protocols::arithmetic_gate::gate::{Gate, Operation};
    use ark_bn254::Fq;
    #[test]
    pub fn test_gkr_sum_check() {
        let mut circuit = Circuit::new(vec![
            vec![
                Gate::new(0, 1, Operation::Add),
                Gate::new(2, 3, Operation::Add),
                Gate::new(4, 5, Operation::Add),
                Gate::new(6, 7, Operation::Mul),
            ],
            vec![
                Gate::new(0, 1, Operation::Mul),
                Gate::new(2, 3, Operation::Add),
            ],
            vec![Gate::new(0, 1, Operation::Add)],
        ]);

        let inputs = vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(5),
            Fq::from(6),
            Fq::from(7),
            Fq::from(8),
        ];

        // initialize both GKR prover and verifier with the same circuit
        let gkr_proof =
            GKRProver::generate_proof(&mut circuit.clone(), &mut Transcript::new(), &inputs);

        assert!(GKRVerifier::verify_proof(
            &inputs,
            &mut circuit,
            &mut Transcript::new(),
            gkr_proof
        ))
    }
}
