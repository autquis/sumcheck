use ark_ff::Field;
use ark_poly::{
    polynomial::{
        multivariate::{SparsePolynomial as MSP, SparseTerm, Term},
        univariate::SparsePolynomial as USP,
    },
    Polynomial,
};
use ark_test_curves::bls12_381::Fq;
use rand::Rng;

pub type MPoly = MSP<Fq, SparseTerm>;
pub type UPoly = USP<Fq>;

#[derive(Debug, Clone)]
pub struct Prover {
    pub g: MPoly,
    pub vfy_rands: Vec<Fq>,
}

impl Prover {
    pub fn new(g: &MPoly) -> Self {
        Prover {
            g: g.clone(),
            vfy_rands: vec![],
        }
    }

    /* Fixing r_1,...,r_{j-1} in g, get a univariate
    polynomial in Xj by summing over variables after Xj */
    pub fn get_unipoly(&mut self, vfy_msg: Option<Fq>) -> UPoly {
        if let Some(vfy_msg) = vfy_msg {
            self.vfy_rands.push(vfy_msg);
        }
        // number of remaining variables (including Xj)
        let v = self.g.num_vars - self.vfy_rands.len();
        let mut sum = UPoly::from_coefficients_vec(vec![(0, 0u32.into())]);
        let domain = 0..2u32.pow(v as u32 - 1);
        // Sum over `v - 1` variables.
        for point in domain {
            sum = sum + self.evaluate_round(bin_repr(point as usize, v));
        }
        return sum;
    }

    // Evaluates th polynomial of round `j` at a vector `v - 1` remaining variables.
    pub fn evaluate_round(&self, points: Vec<Fq>) -> UPoly {
        let mut sum = UPoly::from_coefficients_vec(vec![]);
        for (coeff, term) in self.g.terms.clone().into_iter() {
            let (coeff_eval, fixed_term) = self.evaluate_term(&term, &points);
            let curr = match fixed_term {
                Some(_) => UPoly::from_coefficients_vec(vec![(
                    fixed_term.unwrap().degree(),
                    coeff * coeff_eval,
                )]),
                None => UPoly::from_coefficients_vec(vec![(0, coeff * coeff_eval)]),
            };
            sum = curr + sum
        }
        return sum;
    }

    // Returns the term by substituting variables and verifier's messages.
    pub fn evaluate_term(&self, term: &SparseTerm, point: &Vec<Fq>) -> (Fq, Option<SparseTerm>) {
        let mut fixed_term: Option<SparseTerm> = None;
        let mut product: Fq = 1u32.into();
        for (var, power) in term.into_iter() {
            // if the variable is Xj, do not touch it.
            if *var == self.vfy_rands.len() {
                fixed_term = Some(SparseTerm::new(vec![(*var, *power)]));
            // if the variable is already fixed by verifier message, substitute the value
            } else if *var < self.vfy_rands.len() {
                product = product * self.vfy_rands[*var].pow(&[*power as u64]);
            // if the variable is after Xj, evaluate it at `point`.
            } else {
                product = product * point[*var - self.vfy_rands.len()].pow(&[*power as u64]);
            }
        }

        return (product, fixed_term);
    }
}

pub fn verify(g: &MPoly, c_1: Fq, deg: usize) -> bool {
    // check the first claim
    let mut p = Prover::new(g);
    let mut gi = p.get_unipoly(None);
    let mut expected_c = gi.evaluate(&0u32.into()) + gi.evaluate(&1u32.into());
    assert!(gi.degree() <= deg);
    assert_eq!(c_1, expected_c);

    for _j in 1..p.g.num_vars {
        let r = rand_elem();
        expected_c = gi.evaluate(&r.unwrap());
        gi = p.get_unipoly(r);
        let new_c = gi.evaluate(&0u32.into()) + gi.evaluate(&1u32.into());
        assert!(gi.degree() <= deg);
        assert_eq!(expected_c, new_c);
    }
    // consulting the polynomial at a random point
    let r = rand_elem();
    expected_c = gi.evaluate(&r.unwrap());
    p.vfy_rands.push(r.unwrap());
    let new_c = p.g.evaluate(&p.vfy_rands);
    assert!(gi.degree() <= deg);
    assert_eq!(expected_c, new_c);
    true
}

pub fn rand_elem() -> Option<Fq> {
    let mut rng = rand::thread_rng();
    let r: Fq = rng.gen();
    Some(r)
}

pub fn rand_num(moduli: usize) -> Option<usize> {
    let mut rng = rand::thread_rng();
    let r: usize = rng.gen();
    Some(r % moduli)
}

// Binary representation of points over hypercube (extending to n bits)
pub fn bin_repr(point: usize, n: usize) -> Vec<Fq> {
    let mut res = vec![];
    let mut num: u32 = point as u32;
    for c in 0u32..n as u32 {
        let base = 2u32.pow(n as u32 - 1 - c);
        if num / base > 0 {
            res.push(1.into());
            num -= base;
        } else {
            res.push(0.into());
        }
    }
    return res;
}

fn main() {}
