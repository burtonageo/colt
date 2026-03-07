#![cfg_attr(not(feature = "std"), no_std)]
#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use crate::ray::Intersection as RayIntersection;
use core::cmp::PartialOrd;
use core::{
    mem::{self, MaybeUninit},
    ops::Neg,
    ptr,
};
use vectral::{
    point::Point,
    utils::num::{ClosedAdd, ClosedDiv, ClosedMul, ClosedSub, One, Zero},
    vector::Vector,
};

#[cfg(feature = "std")]
extern crate std;

pub mod ball;
pub mod bounding_box;
pub mod capsule;
pub mod cone;
pub mod cylinder;
pub mod oriented;
pub mod plane;
pub mod ray;

pub trait Volume<const DIM: usize> {
    type Scalar;

    #[must_use]
    fn origin(&self) -> Point<Self::Scalar, DIM>;

    #[must_use]
    fn set_origin(&mut self, new_origin: Point<Self::Scalar, DIM>);

    #[must_use]
    fn contains(&self, point: &Point<Self::Scalar, DIM>) -> bool;

    #[must_use]
    fn support_point(&self, direction: &Vector<Self::Scalar, DIM>) -> Point<Self::Scalar, DIM>;
}

pub trait Intersect<Shape = Self, const DIM: usize = 3>: Volume<DIM> {
    #[must_use]
    fn intersection_with(&self, other: &Shape) -> Option<Intersection<Self::Scalar, DIM>>;

    #[must_use]
    fn intersects(&self, other: &Shape) -> bool {
        self.intersection_with(other).is_some()
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Intersection<T, const DIM: usize> {
    /// The unit normal of the collision.
    pub collision_normal: Vector<T, DIM>,
    /// The penetration distance into the other object.
    pub penetration_depth: T,
}

impl<T: Neg<Output = T>, const DIM: usize> Intersection<T, DIM> {
    #[must_use]
    #[inline]
    pub fn reversed(self) -> Self {
        Intersection {
            collision_normal: self.collision_normal.neg(),
            penetration_depth: self.penetration_depth,
        }
    }
}

#[inline]
pub fn closest_points_between<V1, V2, const N: usize>(
    v1: &V1,
    v2: &V2,
) -> ClosestPoints<V1::Scalar, N>
where
    V1: ?Sized + Volume<N>,
    V2: ?Sized + Volume<N, Scalar = V1::Scalar>,
    V1::Scalar: Copy + ClosedSub + ClosedAdd + ClosedDiv + ClosedMul + One + PartialOrd + Zero,
{
    let (v1_origin, v2_origin) = (v1.origin(), v2.origin());
    let (v1_to_v2, v2_to_v1) = (
        v1_origin.vector_to(v2_origin),
        v2_origin.vector_to(v1_origin),
    );
    let (v1_support, v2_support) = (v1.support_point(&v1_to_v2), v2.support_point(&v2_to_v1));

    // If the origin of the first shape is closer to the second shape's support point than its own,
    // then the shapes overlap. In this case, average the support points and return it as a single
    // support point.
    if v1_origin.distance_squared_to(v2_support) < v1_origin.distance_squared_to(v1_support) {
        ClosestPoints::Connected(Point::midpoint(v1_support, v2_support))
    } else {
        ClosestPoints::Disconnected(v1_support, v2_support)
    }
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum ClosestPoints<T, const N: usize> {
    Connected(Point<T, N>),
    Disconnected(Point<T, N>, Point<T, N>),
}

impl<T, const N: usize> From<Point<T, N>> for ClosestPoints<T, N> {
    #[inline]
    fn from(value: Point<T, N>) -> Self {
        Self::Connected(value)
    }
}

impl<T, const N: usize> From<(Point<T, N>, Point<T, N>)> for ClosestPoints<T, N> {
    #[inline]
    fn from((point_1, point_2): (Point<T, N>, Point<T, N>)) -> Self {
        Self::Disconnected(point_1, point_2)
    }
}

#[must_use]
#[inline]
pub(crate) fn line_overlap<T>(mut min1: T, mut max1: T, mut min2: T, mut max2: T) -> T
where
    T: ClosedSub + PartialOrd + Zero,
{
    #[inline(always)]
    fn min<T: PartialOrd>(x: T, y: T) -> T {
        if x < y { x } else { y }
    }

    #[inline(always)]
    fn max<T: PartialOrd>(x: T, y: T) -> T {
        if x > y { x } else { y }
    }

    if min1 > max1 {
        mem::swap(&mut min1, &mut max1);
    }

    if min2 > max2 {
        mem::swap(&mut min2, &mut max2);
    }

    max(T::ZERO, min(max1, max2) - max(min1, min2))
}

#[inline]
pub(crate) unsafe fn array_assume_init<T, const N: usize>(array: [MaybeUninit<T>; N]) -> [T; N] {
    let mut return_array = MaybeUninit::<[T; N]>::uninit();
    unsafe {
        ptr::copy_nonoverlapping(
            array.as_ptr() as *const T,
            return_array.as_mut_ptr() as *mut _,
            N,
        );
        MaybeUninit::assume_init(return_array)
    }
}
