use ark_bls12_381::{Bls12_381, G1Projective};
use ark_ec::{PairingEngine, ProjectiveCurve};
use ark_ff::{PrimeField, Zero};

fn main() {
    // 创建一些 G1 群的点
    let points = vec![
        G1Projective::prime_subgroup_generator(),
        G1Projective::prime_subgroup_generator(),
        // ... 这里可以添加更多点
    ];

    // 创建一些标量
    let scalars = vec![
        <Bls12_381 as PairingEngine>::Fr::from(123456789u64),
        <Bls12_381 as PairingEngine>::Fr::from(987654321u64),
        // ... 这里可以添加更多标量
    ];

    assert_eq!(points.len(), scalars.len());
    //
    // // 执行 MSM：逐个计算标量乘法并累加结果
    let mut msm_result = G1Projective::zero();
    for (point, scalar) in points.iter().zip(scalars.iter()) {
        msm_result += point.mul(scalar.into_repr());
    }
    //
    // // 转换为仿射坐标并打印
    let affine_result = msm_result.into_affine();
    println!("Result of MSM in affine coordinates: {:?}", affine_result);
}
