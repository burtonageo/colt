use crate::{ball::Ball, ray};
use approx::AbsDiffEq;
use core::{cmp::PartialOrd, ops::Neg};
use vectral::{
    point::Point,
    utils::num::{
        ClosedAdd, ClosedDiv, ClosedMul, ClosedNeg, ClosedSub, One, Scalar, Signed, Sqrt, Zero,
    },
    vector::Vector,
};

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct Plane<T, const N: usize> {
    pub normal: Vector<T, N>,
    pub distance: T,
}

impl<T, const N: usize> Plane<T, N> {
    #[must_use]
    #[inline]
    pub const fn new_unnormalised(normal: Vector<T, N>, distance: T) -> Self {
        Self { normal, distance }
    }
}

impl<T: Scalar + Sqrt, const N: usize> Plane<T, N> {
    #[must_use]
    #[inline]
    pub fn new(mut normal: Vector<T, N>, distance: T) -> Option<Self> {
        let len = normal.len();
        if len > T::ZERO {
            normal = normal / len;
            Some(Plane {
                normal,
                distance: distance * len,
            })
        } else {
            None
        }
    }
}

impl<T: Copy + Zero + ClosedAdd + ClosedMul, const N: usize> Plane<T, N> {
    #[must_use]
    #[inline]
    pub fn position(&self) -> Point<T, N> {
        Point::<T, _>::origin() + (self.normal * self.distance)
    }
}

impl<T: Copy + Signed, const N: usize> Plane<T, N> {
    #[inline]
    pub fn abs_distance(&self) -> T {
        self.distance.abs()
    }
}

impl<T: ClosedNeg, const N: usize> Plane<T, N> {
    #[inline]
    pub fn flipped(self) -> Self {
        self.neg()
    }
}

impl<T: ClosedNeg + Clone, const N: usize> Plane<T, N> {
    #[inline]
    pub fn flip(&mut self) {
        *self = self.clone().neg();
    }
}

impl<T: Neg, const N: usize> Neg for Plane<T, N> {
    type Output = Plane<T::Output, N>;
    #[inline]
    fn neg(self) -> Self::Output {
        Plane {
            normal: self.normal.neg(),
            distance: self.distance.neg(),
        }
    }
}

impl<T: One + ClosedMul + Signed, const N: usize> Plane<T, N> {
    #[must_use]
    #[inline]
    pub fn is_facing_outwards(self) -> bool {
        let product = self.normal.into_iter().fold(T::ONE, |x, y| x * y);
        product.is_negative() == self.distance.is_negative()
    }
}

impl<T: Copy + Zero + ClosedAdd + ClosedMul + ClosedSub, const N: usize> Plane<T, N> {
    #[must_use]
    #[inline]
    pub fn distance_to(&self, point: &Point<T, N>) -> T {
        Vector::dot(*point - self.position(), self.normal)
    }
}

impl<T, const N: usize> ray::Intersect<N> for Plane<T, N>
where
    T: AbsDiffEq
        + ClosedAdd
        + ClosedMul
        + ClosedSub
        + ClosedDiv
        + ClosedNeg
        + Copy
        + PartialOrd<T::Epsilon>
        + PartialOrd
        + Signed
        + Zero,
{
    type Scalar = T;

    #[inline]
    fn intersection_with(
        &self,
        ray: &ray::Ray<Self::Scalar, N>,
    ) -> Option<ray::Intersection<Self::Scalar, N>> {
        let denom = Vector::dot(self.normal, ray.direction);
        if denom.abs() >= T::default_epsilon() {
            let pos = self.position();
            let t = (pos - ray.origin).dot(self.normal) / denom;
            if t >= T::ZERO {
                let is_exterior_face = denom >= T::ZERO;
                let face_normal = if is_exterior_face {
                    self.normal
                } else {
                    -self.normal
                };

                return Some(ray::Intersection {
                    position: ray.point_at(t),
                    face_normal,
                    distance_to_intersection: t,
                    is_exterior_face,
                });
            }
        }

        None
    }
}

#[must_use]
#[inline]
pub(crate) fn ball_inside_plane<T, const N: usize>(ball: &Ball<T, N>, plane: &Plane<T, N>) -> bool
where
    T: Copy + Zero + ClosedAdd + ClosedMul + ClosedSub + ClosedNeg + PartialOrd,
{
    -Plane::distance_to(plane, &ball.center) > ball.radius
}

#[allow(unused)]
#[must_use]
#[inline]
pub(crate) fn ball_outside_plane<T, const N: usize>(ball: &Ball<T, N>, plane: &Plane<T, N>) -> bool
where
    T: Copy + Zero + ClosedAdd + ClosedMul + ClosedSub + PartialOrd,
{
    Plane::distance_to(plane, &ball.center) > ball.radius
}

#[allow(unused)]
#[must_use]
#[inline]
pub(crate) fn ball_intersects_plane<T, const N: usize>(
    ball: &Ball<T, N>,
    plane: &Plane<T, N>,
) -> bool
where
    T: Copy + Zero + ClosedAdd + ClosedMul + ClosedSub + PartialOrd + Signed,
{
    Plane::distance_to(plane, &ball.center).abs() >= ball.radius
}

#[must_use]
#[inline]
pub(crate) fn ball_intersects_plane_point<T, const N: usize>(
    ball: &Ball<T, N>,
    plane: &Plane<T, N>,
) -> Option<(Point<T, N>, T)>
where
    T: Copy + Zero + ClosedAdd + ClosedMul + ClosedSub + PartialOrd + Signed + Sqrt,
{
    let d = Plane::distance_to(plane, &ball.center);

    if d.abs() <= ball.radius {
        let projected = plane.normal * d;
        let point = ball.center - projected;
        let radius = {
            let mut x = (ball.radius * ball.radius) - (d * d);
            if x < T::ZERO {
                x = T::ZERO
            };
            Sqrt::sqrt(x)
        };

        Some((point, radius))
    } else {
        None
    }
}
