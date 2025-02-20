use polynomials::multilinear_polynomial::evaluation_form::MultiLinearPolynomial;

use ark_ff::PrimeField;

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
        .evaluate(&eval_points_rb)
        .scalar_mul(*alpha)
        .add(&muli_a_b_c.evaluate(&eval_points_rc).scalar_mul(*beta));

    let new_addi_b_c = addi_a_b_c
        .evaluate(&eval_points_rb)
        .scalar_mul(*alpha)
        .add(&addi_a_b_c.evaluate(&eval_points_rc).scalar_mul(*beta));

    (new_muli_b_c, new_addi_b_c)
}

pub fn get_evaluated_muli_addi_at_a<T: PrimeField>(
    muli_a_b_c: MultiLinearPolynomial<T>,
    addi_a_b_c: MultiLinearPolynomial<T>,
    random_values: &[Option<T>],
) -> (MultiLinearPolynomial<T>, MultiLinearPolynomial<T>) {
    // Gets the points for "a" to partially evaluate muli and addi at
    let evaluation_points = (0..(muli_a_b_c.number_of_variables() as usize))
        .map(|idx| {
            if idx < random_values.len() {
                return random_values[idx];
            }

            return None;
        })
        .collect::<Vec<_>>();

    (
        muli_a_b_c.evaluate(&evaluation_points),
        addi_a_b_c.evaluate(&evaluation_points),
    )
}
