use ark_bn254::Fq;
use zk_cohort::UnivariatePolynomial;

fn main() {
    let poly2 = UnivariatePolynomial {
        coeffients: vec![Fq::from(2), Fq::from(-3), Fq::from(1)],
    };

    let poly = UnivariatePolynomial {
        coeffients: vec![Fq::from(-1), Fq::from(1)],
    };

    let resulting_poly = poly.mul(&poly2);

    println!("{:?}", resulting_poly.coeffients);

    println!(
        "{:?}",
        UnivariatePolynomial::interpolate(vec![Fq::from(8), Fq::from(10), Fq::from(16)])
    );

    let x = 6;
    let v = vec![1, 2, 3];
    let v_iter = v.into_iter();

    let mut curr_x = x;

    let res = v_iter
        .reduce(|acc, curr| {
            let curr_res = acc + (curr * curr_x);
            curr_x *= x;

            curr_res
        })
        .unwrap();

    println!("{:?}", res);
}
