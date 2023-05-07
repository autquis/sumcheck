#[macro_use]
extern crate lazy_static;

use ark_poly::polynomial::multivariate::{SparsePolynomial, SparseTerm, Term};
use ark_poly::polynomial::DenseMVPolynomial;
use ark_poly::Polynomial;
use ark_test_curves::bls12_381::Fq as ScalarField;
use rstest::rstest;
use sumcheck::temp as sumcheck;

// Sum all evaluations of polynomial `g` over boolean hypercube
fn slow_sum_g(prover: sumcheck::Prover) -> ScalarField {
	let v = prover.g.num_vars;
	let n = 2u32.pow(v as u32);
	(0..n)
		.map(|n| prover.g.evaluate(&sumcheck::n_to_vec(n as usize, v)))
		.sum()
}


lazy_static! {
    // g = 2(x_1)^3 + (x_1)(x_3) + (x_2)(x_3)
    static ref G_0: sumcheck::MultiPoly = SparsePolynomial::from_coefficients_vec(
        3,
        vec![
            (2u32.into(), SparseTerm::new(vec![(0, 3)])),
            (1u32.into(), SparseTerm::new(vec![(0, 1), (2, 1)])),
            (1u32.into(), SparseTerm::new(vec![(1, 1), (2, 1)])),
        ],
    );
    static ref G_0_SUM: ScalarField = slow_sum_g(sumcheck::Prover::new(&G_0));
    // Test with a larger g
    static ref G_1: sumcheck::MultiPoly = SparsePolynomial::from_coefficients_vec(
        4,
        vec![
            (2u32.into(), SparseTerm::new(vec![(0, 3)])),
            (1u32.into(), SparseTerm::new(vec![(0, 1), (2, 1)])),
            (1u32.into(), SparseTerm::new(vec![(1, 1), (2, 1)])),
            (1u32.into(), SparseTerm::new(vec![(3, 1), (2, 1)])),
        ],
    );
    static ref G_1_SUM: ScalarField = slow_sum_g(sumcheck::Prover::new(&G_1));
}

#[rstest]
#[case(&G_0, &G_0_SUM)]
#[case(&G_1, &G_1_SUM)]
fn sumcheck_test(#[case] p: &sumcheck::MultiPoly, #[case] c: &ScalarField) {
    assert!(sumcheck::verify(&p, *c));
}
