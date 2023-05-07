use ark_ff::Field;
use ark_poly::polynomial::multivariate::{SparsePolynomial, SparseTerm, Term};
use ark_poly::polynomial::univariate::SparsePolynomial as UniSparsePolynomial;
use ark_poly::Polynomial;
use ark_test_curves::bls12_381::Fq as ScalarField;
use rand::Rng;

pub type MultiPoly = SparsePolynomial<ScalarField, SparseTerm>;
pub type UniPoly = UniSparsePolynomial<ScalarField>;

// Converts i into an index in {0,1}^v
pub fn n_to_vec(i: usize, n: usize) -> Vec<ScalarField> {
    let mut res = vec!();
    let mut num: u32 = i as u32;
    for c in (0u32..n as u32) {
        let base = 2u32.pow(n as u32  - 1 - c);
        if num / base > 0 {
            res.push(1.into());
            num -= base;
        } else {
            res.push(0.into());
        }
    }
    return res;
}

// Simulates memory of a single prover instance
#[derive(Debug, Clone)]
pub struct Prover {
    pub g: MultiPoly,
    pub r_vec: Vec<ScalarField>,
}

impl Prover {
    pub fn new(g: &MultiPoly) -> Self {
        Prover {
            g: g.clone(),
            r_vec: vec![],
        }
    }

    // Given polynomial g, fix Xj, evaluate over xj+1
    pub fn gen_uni_polynomial(&mut self, r: Option<ScalarField>) -> UniPoly {
        if r.is_some() {
            self.r_vec.push(r.unwrap());
        }
        let v = self.g.num_vars - self.r_vec.len();
        let mut sum = UniPoly::from_coefficients_vec(vec![(0, 0u32.into())]);
        let domain = 0..2u32.pow(v as u32 - 1);
        for point in domain {
            sum = sum + self.evaluate_gj(n_to_vec(point as usize, v));
        }
        return sum;
    }

    // Evaluates gj over a vector permutation of points, folding all evaluated terms together into one univariate polynomial
    pub fn evaluate_gj(&self, points: Vec<ScalarField>) -> UniPoly {
        let mut sum = UniPoly::from_coefficients_vec(vec![]);
        for (coeff, term) in self.g.terms.clone().into_iter() {
            let (coeff_eval, fixed_term) = self.evaluate_term(&term, &points);
            let curr = match fixed_term {
                Some(_) => UniPoly::from_coefficients_vec(vec![(
                    fixed_term.unwrap().degree(),
                    coeff * coeff_eval,
                )]),
                None => UniPoly::from_coefficients_vec(vec![(0, coeff * coeff_eval)]),
            };
            sum = curr + sum
        }
        return sum;
    }

    // Evaluates a term with a fixed univar, returning (new coefficent, fixed term)
    pub fn evaluate_term(
        &self,
        term: &SparseTerm,
        point: &Vec<ScalarField>,
    ) -> (ScalarField, Option<SparseTerm>) {
        let mut fixed_term: Option<SparseTerm> = None;
        let mut product: ScalarField = 1u32.into();
        for (var, power) in term.into_iter() {
            if *var == self.r_vec.len() {
                fixed_term = Some(SparseTerm::new(vec![(*var, *power)]));
            }
            else if *var < self.r_vec.len() {
                product = product * self.r_vec[*var].pow(&[*power as u64]);
            }
            else {
                product = product * point[var - self.r_vec.len()].pow(&[*power as u64]);
            }
        }

        return (product, fixed_term);
    }

    
}

// Verifier procedures
pub fn get_r() -> Option<ScalarField> {
    let mut rng = rand::thread_rng();
    let r: ScalarField = rng.gen();
    Some(r)
}

// A degree look up table for all variables in g
pub fn max_degrees(g: &MultiPoly) -> Vec<usize> {
    let mut lookup: Vec<usize> = vec![0; g.num_vars];
    for (_, term) in &g.terms {
        for (var, power) in term.into_iter() {
            if *power > lookup[*var] {
                lookup[*var] = *power;
            }
        }
    }
    return lookup;
}

// Verify prover's claim c_1
// Presented pedantically:
pub fn verify(g: &MultiPoly, c_1: ScalarField) -> bool {
    // 1st round
    let mut p = Prover::new(g);
    let mut gi = p.gen_uni_polynomial(None);
    let mut expected_c = gi.evaluate(&0u32.into()) + gi.evaluate(&1u32.into());
    assert_eq!(c_1, expected_c);
    let lookup_degree = max_degrees(&g);
    assert!(gi.degree() <= lookup_degree[0]);

    // middle rounds
    for j in 1..p.g.num_vars {
        let r = get_r();
        expected_c = gi.evaluate(&r.unwrap());
        gi = p.gen_uni_polynomial(r);
        let new_c = gi.evaluate(&0u32.into()) + gi.evaluate(&1u32.into());
        assert_eq!(expected_c, new_c);
        assert!(gi.degree() <= lookup_degree[j]);
    }
    // final round
    let r = get_r();
    expected_c = gi.evaluate(&r.unwrap());
    p.r_vec.push(r.unwrap());
    let new_c = p.g.evaluate(&p.r_vec);
    assert_eq!(expected_c, new_c);
    true
}


fn main() {
    println!()
}
