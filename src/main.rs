use ark_bls12_381::{Bls12_381, G1Projective};
use ark_ec::{PairingEngine, ProjectiveCurve};
use ark_ff::{BigInteger, PrimeField, Zero};
use ark_std::rand::rngs::StdRng;
use ark_std::rand::SeedableRng;
use ark_std::UniformRand; // 确保导入正确的 BigInteger 类型

fn main() {
    let (points, scalars) = generate_random_points_and_scalars(10);
    assert_eq!(points.len(), scalars.len());

    //调用函数
    let naive_msm_res = naive_msm(&points, &scalars);

    //
    // // 转换为仿射坐标并打印
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
    // 创建 G1 群的点和对应的标量
    let mut rng = StdRng::seed_from_u64(12345); // 创建随机数生成器

    let points = (0..num_pairs)
        .map(|_| G1Projective::rand(&mut rng)) // 生成随机点
        .collect::<Vec<_>>();

    let scalars = (0..num_pairs)
        .map(|_| <Bls12_381 as PairingEngine>::Fr::rand(&mut rng)) // 生成随机标量
        .collect::<Vec<_>>();
    (points, scalars)
}
fn c_bit_msm(
    points: &[G1Projective],
    scalars: &[<Bls12_381 as PairingEngine>::Fr],
    c: usize,
) -> G1Projective {
    // 使用2^c - 1个桶
    let num_buckets = (1 << c) - 1;
    let mut buckets = vec![G1Projective::zero(); num_buckets];

    // 分配点到对应的桶中
    for (point, scalar) in points.iter().zip(scalars.iter()) {
        let scalar_value = scalar.into_repr().as_ref()[0] as usize; // 由于标量小于2^c，我们可以直接使用它作为索引
        if scalar_value != 0 {
            buckets[scalar_value - 1] += point;
        }
    }

    // 初始化累积变量
    let mut accumulator = G1Projective::zero();

    // 初始化最终结果
    let mut result = G1Projective::zero();

    // 从最后一个桶（即最大的索引）开始迭代
    for bucket in buckets.into_iter().rev() {
        // 累加当前桶的点
        accumulator += bucket;
        // 将当前的累积值添加到最终结果中
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

    //分成 b / c = k组，有k个向量，也就是T1,T2 ... Tk
    //Ti同样有n个元素
    // 初始化向量，用于存储k个c-bit MSM的结果点
    let k = b / c; // 分组数量
    let mut t_points: Vec<G1Projective> = Vec::new();

    // 对每组标量的c位块进行c-bit MSM
    for j in 0..k {
        let mut scalars_c_bits = vec![<Bls12_381 as PairingEngine>::Fr::zero(); scalars.len()];
        for (i, scalar) in scalars.iter().enumerate() {
            let c_bit_chunk = get_c_bit_chunk(scalar, j, c);
            scalars_c_bits[i] = c_bit_chunk;
        }
        // println!("original scalars: {:?}", scalars);
        // println!("scalars_c_bits: {:?}", scalars_c_bits);
        // 计算每组的c-bit MSM并存储结果点
        t_points.push(c_bit_msm(points, &scalars_c_bits, c));
    }

    // 从第一个结果点开始，将其作为初始累加结果
    let mut result = t_points[0];

    // 遍历其余的结果点，每次将累加结果翻倍c次，然后加上下一个结果点
    for j in 1..k {
        for _ in 0..c {
            result = result.double(); // 翻倍c次
        }
        result += &t_points[j]; // 加上下一个结果点
    }

    result // 返回最终的累加结果点
}
fn get_c_bit_chunk(
    scalar: &<Bls12_381 as PairingEngine>::Fr,
    chunk_index: usize,
    chunk_size: usize,
) -> <Bls12_381 as PairingEngine>::Fr {
    // 将标量转换为其 BigInteger 表示
    let scalar_repr = scalar.into_repr();

    // 计算位块的起始和结束位置
    let start_bit = chunk_index * chunk_size;
    let end_bit = start_bit + chunk_size;

    // 获取 BigInteger 的位表示
    let bits = scalar_repr.to_bits_be();
    println!("bits:{:?}", bits);
    // 提取所需的位块
    let chunk_bits = bits[start_bit..std::cmp::min(end_bit, bits.len())].to_vec();

    println!("chunk_bits:{:?}", chunk_bits);
    // 构建一个新的 BigInteger 表示
    let chunk_repr =
        <<Bls12_381 as PairingEngine>::Fr as PrimeField>::BigInt::from_bits_be(&chunk_bits);
    println!("chunk_repr:{:?}", chunk_repr);

    // 将 BigInteger 转换回 Fr 元素
    let res = <Bls12_381 as PairingEngine>::Fr::from_repr(chunk_repr).unwrap();
    println!("Fr chunk_repr:{:?}", res);
    res
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ec::PairingEngine;

    use ark_std::test_rng;

    #[test]
    fn test_msm() {
        // 创建测试数据
        // let points = vec![
        //     G1Projective::prime_subgroup_generator(),
        //     G1Projective::prime_subgroup_generator().double(),
        // ];
        // let scalars = vec![
        //     <Bls12_381 as PairingEngine>::Fr::from(2u64),
        //     <Bls12_381 as PairingEngine>::Fr::from(5u64),
        // ];
        let (points, scalars) = generate_random_points_and_scalars(3);
        // 执行 MSM

        let naive_msm_result = naive_msm(&points, &scalars);
        let b_bit_msm_result = b_bit_msm(&points, &scalars, 256, 4);

        // 转换结果为仿射坐标进行比较
        assert_eq!(
            naive_msm_result.into_affine(),
            b_bit_msm_result.into_affine()
        );
    }
    #[test]
    fn test_get_c_bit_chunk() {
        let mut rng = test_rng();

        // 生成一个随机标量
        let scalar = <Bls12_381 as PairingEngine>::Fr::rand(&mut rng);

        // 定义要提取的位块的大小和索引
        let chunk_index = 0; // 可以选择不同的索引进行测试
        let chunk_size = 16; // 例如，提取16位

        // 使用函数提取位块
        let chunk = get_c_bit_chunk(&scalar, chunk_index, chunk_size);

        // 将原始标量和提取的位块转换为位表示形式
        let scalar_bits = scalar.into_repr().to_bits_be(); //最大在低位
                                                           // let chunk_bits = chunk.into_repr().to_bits_be();
                                                           // // println!("original scalar {:?}", scalar);
                                                           // // println!("scalar_bits:{:?}", scalar_bits);
                                                           // // 计算期望的位块
                                                           // let expected_bits = &scalar_bits[chunk_index * chunk_size..][..chunk_size];

        // println!("expected bits:{:?}", expected_bits);
        // println!("chunk_bits:{:?}", chunk_bits);

        // 比较提取的位块与期望的位块
        // assert_eq!(
        //     chunk_bits, expected_bits,
        //     "Extracted chunk does not match expected bits."
        // );
    }
}
