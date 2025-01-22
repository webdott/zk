use crate::univariate_polynomial;
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

    pub fn generate_shares(&self, secret: T) -> Vec<(T, T)> {
        let (x_points, y_points): (&mut Vec<T>, &mut Vec<T>) =
            (&mut vec![self.secret_x], &mut vec![secret]);

        let mut random = rand::thread_rng();

        let n = if self.secret_x < T::from(self.quorom) {
            self.quorom
        } else {
            self.quorom - 1
        };

        for i in 0..n {
            if T::from(i) == self.secret_x {
                continue;
            }

            let random_y: T = T::from(random.gen_range(0..100));

            x_points.push(T::from(i));
            y_points.push(random_y);
        }

        if self.secret_x >= T::from(self.quorom) {}

        let polynomial = univariate_polynomial::UnivariatePolynomial::interpolate(
            x_points.clone(),
            y_points.clone(),
        );

        let mut shares = vec![];
        let mut shares_num = 0;

        while shares_num < self.number_of_shares {
            if (T::from(shares_num) == self.secret_x) {
                shares_num += 1;
            }

            shares.push((
                T::from(shares_num),
                polynomial.evaluate(T::from(shares_num)),
            ));

            shares_num += 1;
        }

        println!("{:?}", shares);
        shares
    }

    fn verify_shares(&self, shares: &[(T, T)]) -> bool {
        shares.len() >= self.quorom as usize
    }

    pub fn reconstruct_secret(&self, shares: &[(T, T)]) -> Result<T, &str> {
        if !self.verify_shares(shares) {
            return Err("Not enough shares to reconstruct secret");
        }

        let (x_points, y_points) = shares.iter().fold((vec![], vec![]), |acc, curr| {
            (
                {
                    let mut vec_a = vec![curr.0]; // Start with the current value
                    vec_a.extend(acc.0); // Append the previous values
                    vec_a // Return the resulting vector
                },
                {
                    let mut vec_b = vec![curr.1]; // Start with the current value
                    vec_b.extend(acc.1); // Append the previous values
                    vec_b // Return the resulting vector
                },
            )
        });

        let original_polynomial =
            univariate_polynomial::UnivariatePolynomial::interpolate(x_points, y_points);

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
