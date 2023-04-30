use ark_poly::{
    polynomial::multivariate::{SparsePolynomial, SparseTerm, Term},
    DenseMVPolynomial,
};
use ark_test_curves::bls12_381::Fq;

fn build_poly(
    num_vars: usize,
    coeff: Vec<i128>,
    sparse_terms: Vec<SparseTerm>,
) -> SparsePolynomial<Fq, SparseTerm> {
    assert_eq!(
        coeff.len(),
        sparse_terms.len(),
        "coefficient length: {}, sparse_terms: {} doen't match",
        coeff.len(),
        sparse_terms.len()
    );

    let mut coeff_vec = Vec::<(Fq, SparseTerm)>::new();

    for idx in 0..coeff.len() {
        coeff_vec.push((Fq::from(coeff[idx]), sparse_terms[idx].clone()));
    }

    let poly = SparsePolynomial::from_coefficients_vec(num_vars, coeff_vec);

    return poly;
}

pub struct PState {
    pub vfy_rands: Vec<Fq>,
    pub round: usize
}

fn reduce_poly(poly: &mut SparsePolynomial<Fq, SparseTerm>, pstate: &mut PState) {
    
}

fn main() {
    println!("Hello, world!");
}
