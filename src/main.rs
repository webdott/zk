mod concepts_protocols;
mod polynomials;

use ark_bn254::Fq;

fn main() {
    let poly2 = polynomials::univariate_polynomial::UnivariatePolynomial {
        coeffients: vec![Fq::from(2), Fq::from(-3), Fq::from(1)],
    };

    let poly = polynomials::univariate_polynomial::UnivariatePolynomial {
        coeffients: vec![Fq::from(-1), Fq::from(1)],
    };

    let resulting_poly = poly.mul(&poly2);

    println!("{:?}", resulting_poly.coeffients);

    println!(
        "{:?}",
        polynomials::univariate_polynomial::UnivariatePolynomial::interpolate(
            vec![Fq::from(0), Fq::from(1), Fq::from(2)],
            vec![Fq::from(8), Fq::from(10), Fq::from(16)]
        )
    );
}
