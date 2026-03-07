#![allow(unused)]

use crate::{Volume, ray};
use approx::AbsDiffEq;
use core::{
    cmp::{PartialEq, PartialOrd},
    mem,
};
use vectral::{
    point::Point,
    utils::num::{ClosedAdd, ClosedDiv, ClosedMul, ClosedNeg, ClosedSub, One, Signed, Sqrt, Zero},
    vector::Vector,
};

#[derive(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub struct Cylinder<T, const N: usize> {
    pub origin: Point<T, N>,
    pub radius: T,
    pub half_height: T,
}

impl<T, const N: usize> Cylinder<T, N> {
    #[must_use]
    #[inline]
    pub const fn new(origin: Point<T, N>, radius: T, half_height: T) -> Self {
        Self {
            origin,
            radius,
            half_height,
        }
    }
}

impl<T, const N: usize> Volume<N> for Cylinder<T, N> {
    type Scalar = T;

    #[inline]
    fn contains(&self, _point: &Point<Self::Scalar, N>) -> bool {
        todo!()
    }
    
    fn origin(&self) -> Point<Self::Scalar, N> {
        todo!()
    }
    
    fn set_origin(&mut self, new_origin: Point<Self::Scalar, N>) {
        todo!()
    }
    
    fn support_point(&self, direction: &Vector<Self::Scalar, N>) -> Point<Self::Scalar, N> {
        todo!()
    }
}

impl<T: ClosedAdd + Copy, const N: usize> Cylinder<T, N> {
    #[must_use]
    #[inline]
    pub fn height(&self) -> T {
        self.half_height + self.half_height
    }
}

impl<T: ClosedAdd + ClosedMul + ClosedSub + Copy + One + Zero, const N: usize> Cylinder<T, N> {
    #[inline]
    fn end_caps(&self) -> (Point<T, N>, Point<T, N>) {
        let y = Vector::Y * self.half_height;
        (self.origin + y, self.origin - y)
    }
}

impl<T, const N: usize> ray::Intersect<N> for Cylinder<T, N>
where
    T: AbsDiffEq
        + ClosedAdd
        + ClosedDiv
        + ClosedMul
        + ClosedSub
        + Copy
        + One
        + PartialOrd
        + PartialOrd<<T as AbsDiffEq>::Epsilon>
        + Signed
        + Sqrt
        + Zero,
    <T as AbsDiffEq>::Epsilon: Copy + ClosedNeg + PartialEq<T>,
{
    type Scalar = T;

    fn intersection_with(
        &self,
        ray: &ray::Ray<Self::Scalar, N>,
    ) -> Option<ray::Intersection<Self::Scalar, N>> {
        // https://davidjcobb.github.io/articles/ray-cylinder-intersection
        let (cap_top, cap_bottom) = self.end_caps();
        let r1 = ray.origin - cap_bottom;
        let spine = cap_top - cap_bottom;
        let height = self.height();

        let axis = Vector::normalized_unchecked(spine);

        let axis_dot_raydir = Vector::dot(axis, ray.direction);
        let axis_dot_r1 = Vector::dot(axis, r1);
        let r1_len_sq = r1.len_squared();

        let a = T::ONE - (axis_dot_raydir * axis_dot_raydir);
        let b =
            (T::ONE + T::ONE) * (Vector::dot(ray.direction, r1) - (axis_dot_raydir * axis_dot_r1));
        let c = r1_len_sq - axis_dot_r1 * axis_dot_r1 - (self.radius * self.radius);

        let solution = match quadratic_roots(a, b, c) {
            None => return None,
            Some(Solution::Multiple { lower: near, upper: far }) => {
                let hit_pos_near = ray.point_at(near);
                let hit_pos_far = ray.point_at(far);

                let height_offset_near = Vector::dot(cap_top - hit_pos_near, spine);
                let height_offset_far = Vector::dot(cap_top - hit_pos_far, spine);

                // let 
                if near < T::ZERO || height_offset_near < T::default_epsilon() || height_offset_near > height {

                }

                todo!()
            }
            Some(Solution::Single(solution)) => {
                let hit_pos = ray.point_at(solution);
                let height_offset = Vector::dot(cap_top - hit_pos, spine);

                if solution < T::ZERO || height_offset < T::default_epsilon() || height_offset > height {
                    None
                } else {
                    Some(solution)
                }
            }
        };

        todo!()
    }
}

fn quadratic_roots<T>(a: T, b: T, c: T) -> Option<Solution<T>>
where
    T: AbsDiffEq
        + ClosedAdd
        + ClosedDiv
        + ClosedMul
        + ClosedSub
        + Copy
        + One
        + PartialOrd<<T as AbsDiffEq>::Epsilon>
        + Signed
        + Zero
        + Sqrt,
    <T as AbsDiffEq>::Epsilon: Copy + ClosedNeg + PartialEq<T>,
{
    let _2: T = T::ONE + T::ONE;
    let _4: T = _2 + _2;
    let epsilon = T::default_epsilon();

    let discriminant = (b * b) - (_4 * a * c);

    if discriminant > epsilon {
        let b_term = if b < epsilon {
            -b - Sqrt::sqrt(discriminant)
        } else {
            -b + Sqrt::sqrt(discriminant)
        };

        let mut lower = b_term / (_2 * a);
        let mut upper = (_2 * c) / b_term;
        if lower > upper {
            mem::swap(&mut lower, &mut upper);
        }

        Some(Solution::Multiple { lower, upper })
    } else if {
        let is_zero_epsilon = epsilon == T::ZERO;
        if is_zero_epsilon {
            discriminant < epsilon
        } else {
            discriminant > -epsilon && discriminant <= epsilon
        }
    } {
        let root = -(b / _2 * a);
        Some(Solution::Single(root))
    } else {
        None
    }
}

enum Solution<T> {
    Single(T),
    Multiple { lower: T, upper: T },
}

impl<T> Solution<T> {
    #[inline]
    const fn num_solutions(&self) -> usize {
        match *self {
            Self::Single(_) => 1,
            Self::Multiple { .. } => 2,
        }
    }
}
