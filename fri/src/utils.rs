use ark_ff::FftField;

pub fn fold_layer<T: FftField>(evaluations: &[T], r: &T) -> Vec<T> {
    todo!()
}

pub fn pad_polynomial<T: FftField>(values: &[T], final_length: usize, number_to_add: T) -> Vec<T> {
    let mut final_pad = values.to_vec();

    for idx in final_pad.len()..final_length {
        final_pad[idx] = number_to_add;
    }

    final_pad
}
