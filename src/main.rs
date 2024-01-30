use ark_bls12_381::{Bls12_381, G1Projective};
use ark_ec::{PairingEngine, ProjectiveCurve};
use ark_ff::{BigInteger, PrimeField, Zero};
use ark_std::rand::rngs::StdRng;
use ark_std::rand::SeedableRng;
use ark_std::UniformRand;

fn main() {
    let (points, scalars) = generate_random_points_and_scalars(10);
    assert_eq!(points.len(), scalars.len());

    let naive_msm_res = naive_msm(&points, &scalars);

    let affine_result = naive_msm_res.into_affine();
    println!("Result of MSM in affine coordinates: {:?}", affine_result);
}
fn naive_msm(
    points: &Vec<G1Projective>,
    scalars: &Vec<<Bls12_381 as PairingEngine>::Fr>,
) -> G1Projective {
    let mut res = G1Projective::zero();
    for (point, scalar) in points.iter().zip(scalars.iter()) {
        res += point.mul(scalar.into_repr());
    }
    res
}
fn generate_random_points_and_scalars(
    num_pairs: usize,
) -> (Vec<G1Projective>, Vec<<Bls12_381 as PairingEngine>::Fr>) {
    // create random points in G1 and random scalars
    let mut rng = StdRng::seed_from_u64(12345);

    let points = (0..num_pairs)
        .map(|_| G1Projective::rand(&mut rng)) // generate random points in G1
        .collect::<Vec<_>>();

    let scalars = (0..num_pairs)
        .map(|_| <Bls12_381 as PairingEngine>::Fr::rand(&mut rng)) // generate random scalars in G1
        .collect::<Vec<_>>();
    (points, scalars)
}
fn c_bit_msm(
    points: &[G1Projective],
    scalars: &[<Bls12_381 as PairingEngine>::Fr],
    c: usize,
) -> G1Projective {
    // use 2^c - 1 buckets
    let num_buckets = (1 << c) - 1;
    let mut buckets = vec![G1Projective::zero(); num_buckets];

    // put points into correct buckets
    for (point, scalar) in points.iter().zip(scalars.iter()) {
        let scalar_value = scalar.into_repr().as_ref()[0] as usize; // scalar is generally quite small.
        if scalar_value != 0 {
            buckets[scalar_value - 1] += point;
        }
    }

    let mut accumulator = G1Projective::zero();

    let mut result = G1Projective::zero();

    //iterate reverselly
    for bucket in buckets.into_iter().rev() {
        accumulator += bucket;

        result += &accumulator;
    }

    result
}
fn b_bit_msm(
    points: &[G1Projective],
    scalars: &[<Bls12_381 as PairingEngine>::Fr],
    b: usize,
    c: usize,
) -> G1Projective {
    //scalars = a1,a2,a3 ... an

    let k = b / c;

    //result of
    let mut t_points: Vec<G1Projective> = Vec::new();

    //compute Tj for j = 0..k
    for j in 0..k {
        let mut scalars_c_bits = vec![<Bls12_381 as PairingEngine>::Fr::zero(); scalars.len()];
        for (i, scalar) in scalars.iter().enumerate() {
            let c_bit_chunk = get_c_bit_chunk(scalar, j, c);
            scalars_c_bits[i] = c_bit_chunk;
        }
        // now scalars_c_bits = vectors of coefficient in Tj
        // use c_bit_msm to compute Tj and push Tj in to t_points array
        t_points.push(c_bit_msm(points, &scalars_c_bits, c));
    }

    let mut result = t_points[0];

    //use T1,T2,...,Tk to compute T
    for j in 1..k {
        for _ in 0..c {
            result = result.double(); // double c times
        }
        result += &t_points[j];
    }

    result
}
fn get_c_bit_chunk(
    scalar: &<Bls12_381 as PairingEngine>::Fr,
    chunk_index: usize,
    chunk_size: usize,
) -> <Bls12_381 as PairingEngine>::Fr {
    // compute the c-bit start bit and end bit index
    let start_bit = chunk_index * chunk_size;
    let end_bit = start_bit + chunk_size;

    // convert the scalar to big endian
    let bits = scalar.into_repr().to_bits_be();

    //get the desired chunk
    let chunk_bits = bits[start_bit..std::cmp::min(end_bit, bits.len())].to_vec();

    //convert chunk_bits to appropraite representation in Fr
    let chunk_repr =
        <<Bls12_381 as PairingEngine>::Fr as PrimeField>::BigInt::from_bits_be(&chunk_bits);
    let res = <Bls12_381 as PairingEngine>::Fr::from_repr(chunk_repr).unwrap();

    res
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ec::PairingEngine;

    use ark_std::test_rng;

    #[test]
    fn test_msm() {
        let (points, scalars) = generate_random_points_and_scalars(50);

        let naive_msm_result = naive_msm(&points, &scalars);
        let b_bit_msm_result = b_bit_msm(&points, &scalars, 256, 4);

        assert_eq!(
            naive_msm_result.into_affine(),
            b_bit_msm_result.into_affine()
        );
    }
}
