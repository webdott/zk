use crate::concepts_protocols::arithmetic_gate::circuit::Circuit;
use crate::concepts_protocols::fiat_shamir::transcript::Transcript;
use crate::concepts_protocols::sumcheck::sumcheck_protocol::{
    SumCheckProof, SumcheckProver, SumcheckVerifier,
};
use crate::polynomials::multilinear_polynomial::MultiLinearPolynomial;
use crate::polynomials::product_polynomial::ProductPolynomial;
use crate::polynomials::sum_polynomial::SumPolynomial;
use ark_ff::PrimeField;
use std::marker::PhantomData;

pub fn get_folded_poly_and_claimed_sum<T: PrimeField>(
    w_poly: MultiLinearPolynomial<T>,
    muli_a_b_c: MultiLinearPolynomial<T>,
    addi_a_b_c: MultiLinearPolynomial<T>,
    r_b: &[Option<T>],
    r_c: &[Option<T>],
    transcript: &mut Transcript<T>,
) -> (T, MultiLinearPolynomial<T>, MultiLinearPolynomial<T>) {
    let (alpha, beta) = (transcript.sample_challenge(), transcript.sample_challenge());

    let w_i_b = w_poly.evaluate(Vec::from(r_b));
    let w_i_c = w_poly.evaluate(Vec::from(r_c));

    let claimed_sum = (*w_i_b.get_evaluation_points().first().unwrap() * alpha)
        + (*w_i_c.get_evaluation_points().first().unwrap() * beta);

    let mut eval_points_rb = vec![None; muli_a_b_c.number_of_variables() as usize];
    let mut eval_points_rc = vec![None; muli_a_b_c.number_of_variables() as usize];

    (0..r_b.len()).for_each(|idx| {
        eval_points_rb[idx] = r_b[idx];
        eval_points_rc[idx] = r_c[idx];
    });

    let new_muli_b_c = muli_a_b_c
        .evaluate(eval_points_rb.clone())
        .scalar_mul(alpha)
        .add(&muli_a_b_c.evaluate(eval_points_rc.clone()).scalar_mul(beta));

    let new_addi_b_c = addi_a_b_c
        .evaluate(eval_points_rb)
        .scalar_mul(alpha)
        .add(&addi_a_b_c.evaluate(eval_points_rc).scalar_mul(beta));

    (claimed_sum, new_muli_b_c, new_addi_b_c)
}

#[derive(Debug)]
pub struct GKRProof<T: PrimeField> {
    pub w_polys: Vec<MultiLinearPolynomial<T>>,
    pub sumcheck_proofs: Vec<SumCheckProof<T>>,
}

impl<T: PrimeField> GKRProof<T> {
    pub fn new(
        w_polys: Vec<MultiLinearPolynomial<T>>,
        sumcheck_proofs: Vec<SumCheckProof<T>>,
    ) -> Self {
        Self {
            w_polys,
            sumcheck_proofs,
        }
    }
}

pub struct GKRProver<T: PrimeField> {
    _marker: PhantomData<T>,
}

impl<T: PrimeField> GKRProver<T> {
    pub fn evaluate_at_input(circuit: &mut Circuit<T>, inputs: &[T]) -> MultiLinearPolynomial<T> {
        circuit.evaluate(Vec::from(inputs))
    }

    pub fn generate_proof(
        circuit: &mut Circuit<T>,
        transcript: &mut Transcript<T>,
        inputs: &[T],
    ) -> GKRProof<T> {
        circuit.evaluate(Vec::from(inputs));

        let (mut w_polys, mut sum_check_proofs) = (
            Vec::with_capacity(circuit.get_layer_count()),
            Vec::with_capacity(circuit.get_layer_count()),
        );

        let length_of_rs = circuit
            .get_w_i(0)
            .get_evaluation_points()
            .len()
            .next_power_of_two()
            .ilog2()
            * 2;

        let mut running_polynomial = circuit.get_w_i(0);

        transcript.append(&running_polynomial.to_bytes());

        let mut random_values = vec![];

        (0..length_of_rs).for_each(|_i| {
            random_values.push(Some(transcript.sample_challenge()));
        });

        (0..circuit.get_layer_count()).for_each(|layer_idx| {
            // Things I'm sending the verifier at each step
            //   W (output polys) -> Perform alpha beta folding if more than one output to form one output poly
            //       pass in the output to transcript
            //   Intermediate Sumcheck proof
            //   Mo (Claim sum at each step -> Evaluation of output poly at random sample challenge generated by transcript)

            // perform alpha-beta folding on W poly

            w_polys.push(running_polynomial.clone());

            let muli_a_b_c = circuit.get_mul_i(layer_idx);
            let addi_a_b_c = circuit.get_add_i(layer_idx);

            let (claim_sum, new_muli_b_c, new_addi_b_c) = get_folded_poly_and_claimed_sum(
                running_polynomial.clone(),
                muli_a_b_c,
                addi_a_b_c,
                &random_values[0..random_values.len() / 2],
                &random_values[random_values.len() / 2..],
                transcript,
            );

            let next_w_i = circuit.get_w_i(layer_idx + 1);

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

            let (sumcheck_proof, random_points) = SumcheckProver::generate_proof_for_partial_verify(
                claim_sum,
                f_b_c.clone(),
                transcript,
            );

            random_values = random_points.iter().map(|point| Some(*point)).collect();
            running_polynomial = circuit.get_w_i(layer_idx + 1);

            sum_check_proofs.push(sumcheck_proof);
        });

        GKRProof::new(w_polys, sum_check_proofs)
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
        let length_of_rs = proof
            .w_polys
            .first()
            .unwrap()
            .get_evaluation_points()
            .len()
            .next_power_of_two()
            .ilog2()
            * 2;

        let mut running_polynomial = proof.w_polys.first().unwrap();

        transcript.append(&running_polynomial.to_bytes());

        let mut random_values = vec![];

        (0..length_of_rs).for_each(|_i| {
            random_values.push(Some(transcript.sample_challenge()));
        });

        for layer_idx in 0..circuit.get_layer_count() {
            let muli_a_b_c = circuit.get_mul_i(layer_idx);
            let addi_a_b_c = circuit.get_add_i(layer_idx);

            let (_claim_sum, new_muli_b_c, new_addi_b_c) = get_folded_poly_and_claimed_sum(
                running_polynomial.clone(),
                muli_a_b_c,
                addi_a_b_c,
                &random_values[0..random_values.len() / 2],
                &random_values[random_values.len() / 2..],
                transcript,
            );

            let next_w_i = if layer_idx + 1 >= proof.w_polys.len() {
                &MultiLinearPolynomial::new(Vec::from(initial_inputs))
            } else {
                &proof.w_polys[layer_idx + 1]
            };

            let f_b_c = SumPolynomial::new(vec![
                ProductPolynomial::new(vec![
                    new_muli_b_c,
                    MultiLinearPolynomial::w_mul(next_w_i, next_w_i),
                ]),
                ProductPolynomial::new(vec![
                    new_addi_b_c,
                    MultiLinearPolynomial::w_add(next_w_i, next_w_i),
                ]),
            ]);

            let (is_verified, final_claim_sum, next_evaluation_values) =
                SumcheckVerifier::partial_verify(&proof.sumcheck_proofs[layer_idx], transcript);

            // println!(
            //     "{:?} {:?}",
            //     f_b_c.evaluate(&next_evaluation_values),
            //     final_claim_sum
            // );

            if !is_verified || (f_b_c.evaluate(&next_evaluation_values) != final_claim_sum) {
                return false;
            }

            running_polynomial = if layer_idx + 1 >= proof.w_polys.len() {
                &proof.w_polys[layer_idx]
            } else {
                &proof.w_polys[layer_idx + 1]
            };
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
