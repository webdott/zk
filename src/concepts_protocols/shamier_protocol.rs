use crate::polynomials::univariate_polynomial;
use ark_ff::PrimeField;
use rand::Rng;

pub struct ShamierProtocol<T: PrimeField> {
    quorom: u32,
    number_of_shares: u32,
    secret_x: T,
}

impl<T: PrimeField> ShamierProtocol<T> {
    pub fn new(quorom: u32, number_of_shares: u32, secret_x: T) -> Self {
        ShamierProtocol {
            quorom,
            number_of_shares,
            secret_x,
        }
    }

    // Given a secret and the number of passwords to generate from it, generate the password shares
    pub fn generate_shares(&self, secret: T) -> Vec<(T, T)> {
        // initialize the set of points to interpolate on with the secret's x point
        let (x_points, y_points): (&mut Vec<T>, &mut Vec<T>) =
            (&mut vec![self.secret_x], &mut vec![secret]);

        let mut random = rand::thread_rng();

        let n = if self.secret_x < T::from(self.quorom) {
            self.quorom
        } else {
            self.quorom - 1
        };

        // generate a bunch of points to fill up the quorom and give them random values
        for i in 0..n {
            if T::from(i) == self.secret_x {
                continue;
            }

            let random_y: T = T::from(random.gen_range(0..100));

            x_points.push(T::from(i));
            y_points.push(random_y);
        }

        // Interpolate on the points generated i.e secret and other x_points in the quorom
        let polynomial = univariate_polynomial::UnivariatePolynomial::interpolate(
            x_points.clone(),
            y_points.clone(),
        );

        // Once we get the polynomial, evaluate the polynomial at random set of x_points of length (number of shares) and return them
        std::iter::repeat(())
            .map(|()| T::rand(&mut random))
            .filter(|x| x != &self.secret_x)
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

        // Get the list of x_points and y_points from the shares
        let (x_points, y_points) = shares.iter().fold((vec![], vec![]), |acc, curr| {
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

        // Evaluate the polynomial at secret's x_point
        Ok(original_polynomial.evaluate(self.secret_x))
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bn254::Fq;

    #[test]
    pub fn test_generate_shares() {
        let shamier = ShamierProtocol::new(3, 7, Fq::from(4));

        let shares = shamier.generate_shares(Fq::from(62));

        assert_eq!(shares.len(), 7);
    }

    #[test]
    pub fn test_reconstruct_secret_not_enough() {
        let shamier = ShamierProtocol::new(3, 7, Fq::from(4));

        shamier.generate_shares(Fq::from(62));

        let secret = shamier.reconstruct_secret(&vec![
            (Fq::from(0), Fq::from(15)),
            (Fq::from(1), Fq::from(91)),
        ]);

        assert_eq!(secret, Err("Not enough shares to reconstruct secret"));
    }

    #[test]
    pub fn test_reconstruct_secret_enough_but_wrong() {
        // [(0, 84), (1, 33), (2, 14592161914559516814830937163504850059130874104865215775126025263096817472401), (3, 22), (5, 14592161914559516814830937163504850059130874104865215775126025263096817472521), (6, 233)]
        // Secret -> 62
        let shamier = ShamierProtocol::new(3, 7, Fq::from(4));

        shamier.generate_shares(Fq::from(62));

        let secret = shamier.reconstruct_secret(&vec![
            (Fq::from(0), Fq::from(15)),
            (Fq::from(1), Fq::from(91)),
            (Fq::from(3), Fq::from(15)),
        ]);

        assert_ne!(secret, Ok(Fq::from(62)));
    }

    #[test]
    pub fn test_reconstruct_secret_enough_and_right() {
        // [(0, 84), (1, 33), (2, 14592161914559516814830937163504850059130874104865215775126025263096817472401), (3, 22), (5, 14592161914559516814830937163504850059130874104865215775126025263096817472521), (6, 233)]
        // Secret -> 62
        let secret = Fq::from(62);
        let shamier = ShamierProtocol::new(3, 7, Fq::from(4));

        let shares = shamier.generate_shares(secret.clone());

        let regenerated_secret = shamier.reconstruct_secret(&shares[0..4]);

        assert_eq!(regenerated_secret, Ok(secret));
    }
}
