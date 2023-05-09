use ark_poly::{
    polynomial::{
        multivariate::{SparsePolynomial, SparseTerm, Term},
        DenseMVPolynomial,
    },
    Polynomial,
};
use ark_test_curves::bls12_381::Fq;
use rstest::rstest;
use sumcheck::sumcheck;

// Normal summation of polynomial over the hypercube
fn normal_sum(prover: sumcheck::Prover) -> Fq {
    let n = 2u32.pow(prover.g.num_vars as u32);
    let mut sum = Fq::from(0);
    for i in 0..n {
        sum += prover
            .g
            .evaluate(&sumcheck::bin_repr(i as usize, prover.g.num_vars))
    }
    return sum;
}

#[rstest]
fn sumcheck_test() {
    // Example from Justin Thaler's book: g = 2(x1)^3 + (x1)(x3) + (x2)(x3)
    let g: sumcheck::MPoly = SparsePolynomial::from_coefficients_vec(
        3,
        vec![
            (Fq::from(2), SparseTerm::new(vec![(0, 3)])),
            (Fq::from(1), SparseTerm::new(vec![(0, 1), (2, 1)])),
            (Fq::from(1), SparseTerm::new(vec![(1, 1), (2, 1)])),
        ],
    );
    let g_sum: Fq = normal_sum(sumcheck::Prover::new(&g));
    assert!(g_sum == Fq::from(12));
    assert!(sumcheck::verify(&g, g_sum, 3));

    // some random tests. 10 times, with degree at most 5, number of variables at most 8, and at most 15 non-zero terms.
    for _i in 0..10 {
        let num_var = sumcheck::rand_num(8).unwrap() + 1;
        let deg = sumcheck::rand_num(5).unwrap() + 1;
        let non_zero = sumcheck::rand_num(deg.pow(num_var as u32)).unwrap() % 15 + 1;
        let mut coeff = vec![];
        for _j in 0..non_zero {
            let mut pows = vec![];
            for i in 0..num_var {
                pows.push((i as usize, sumcheck::rand_num(deg).unwrap()));
            }
            coeff.push((
                Fq::from(sumcheck::rand_num(20).unwrap() as u32 + 1),
                SparseTerm::new(pows),
            ));
        }
        let g: sumcheck::MPoly = SparsePolynomial::from_coefficients_vec(num_var, coeff);
        let g_sum: Fq = normal_sum(sumcheck::Prover::new(&g));
        assert!(sumcheck::verify(&g, g_sum, deg));
    }
}
