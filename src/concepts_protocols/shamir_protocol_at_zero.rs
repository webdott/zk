use crate::polynomials::univariate_polynomial;
use crate::polynomials::univariate_polynomial::UnivariatePolynomial;
use ark_ff::PrimeField;
use rand::Rng;
use std::marker::PhantomData;

pub struct ShamirProtocol<T: PrimeField> {
    _marker: PhantomData<T>,
    quorom: u32,
    number_of_shares: u32,
}

impl<T: PrimeField> ShamirProtocol<T> {
    pub fn new(quorom: u32, number_of_shares: u32) -> Self {
        ShamirProtocol {
            _marker: PhantomData,
            quorom,
            number_of_shares,
        }
    }

    // Given a secret and the number of passwords to generate from it, generate the password shares
    pub fn generate_shares(&self, secret: T) -> Vec<(T, T)> {
        let mut evaluation_points = vec![T::from(secret)];

        let mut random = rand::thread_rng();

        (0..self.quorom - 1).for_each(|_i| {
            evaluation_points.push(T::from(random.gen_range(0..100)));
        });

        let polynomial = UnivariatePolynomial::new(evaluation_points);

        // Once we get the polynomial, evaluate the polynomial at random set of x_points of length (number of shares) and return them
        std::iter::repeat(())
            .map(|()| T::rand(&mut random))
            .filter(|x| x != &T::from(0))
            .map(|x| (x.clone(), polynomial.evaluate(x)))
            .take(self.number_of_shares as usize)
            .collect()
    }

    // Verify that the shares given to reconstruct a secret is up to the quorom
    fn verify_shares(&self, shares: &[(T, T)]) -> bool {
        shares.len() >= self.quorom as usize
    }

    // Get back the secret given a list of password shares
    pub fn reconstruct_secret(&self, shares: &[(T, T)]) -> Result<T, &str> {
        if !self.verify_shares(shares) {
            return Err("Not enough shares to reconstruct secret");
        }

        // Get the list of x_points and y_points from the shares (Optimise to take just d + 1 points to get polynomial back)
        let (x_points, y_points) =
            shares[0..self.quorom as usize]
                .iter()
                .fold((vec![], vec![]), |acc, curr| {
                    (
                        {
                            let mut vec_a = vec![curr.0];
                            vec_a.extend(acc.0);
                            vec_a
                        },
                        {
                            let mut vec_b = vec![curr.1];
                            vec_b.extend(acc.1);
                            vec_b
                        },
                    )
                });

        // Get back the polynomial we got while generating the shares
        let original_polynomial =
            univariate_polynomial::UnivariatePolynomial::interpolate(x_points, y_points);

        // Evaluate the polynomial at secret's x_point (0 in this case)
        Ok(original_polynomial.evaluate(T::from(0)))
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bn254::Fq;

    #[test]
    pub fn test_generate_shares() {
        let shamir = ShamirProtocol::new(3, 7);

        let shares = shamir.generate_shares(Fq::from(62));

        assert_eq!(shares.len(), 7);
    }

    #[test]
    pub fn test_reconstruct_secret_not_enough() {
        let shamir = ShamirProtocol::new(3, 7);

        shamir.generate_shares(Fq::from(62));

        let secret = shamir.reconstruct_secret(&vec![
            (Fq::from(0), Fq::from(15)),
            (Fq::from(1), Fq::from(91)),
        ]);

        assert_eq!(secret, Err("Not enough shares to reconstruct secret"));
    }

    #[test]
    pub fn test_reconstruct_secret_enough_but_wrong() {
        let shamir = ShamirProtocol::new(3, 7);

        shamir.generate_shares(Fq::from(62));

        let secret = shamir.reconstruct_secret(&vec![
            (Fq::from(0), Fq::from(15)),
            (Fq::from(1), Fq::from(91)),
            (Fq::from(3), Fq::from(15)),
        ]);

        assert_ne!(secret, Ok(Fq::from(62)));
    }

    #[test]
    pub fn test_reconstruct_secret_enough_and_right() {
        let secret = Fq::from(62);
        let shamir = ShamirProtocol::new(3, 7);

        let shares = shamir.generate_shares(secret.clone());

        let regenerated_secret = shamir.reconstruct_secret(&shares);

        assert_eq!(regenerated_secret, Ok(secret));
    }
}
