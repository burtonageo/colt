use crate::{
    Volume,
    ray::{Intersect as RayIntersect, Intersection as RayIntersection, Ray},
};
use core::ops::DivAssign;
use vectral::{num::ClosedDiv, point::Point, rotation::rotate_point_around, vector::Vector};
use vectral::{
    rotation::Rotation,
    utils::num::{ClosedAdd, ClosedMul, ClosedSub, One, Sqrt, Zero},
};

#[derive(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub struct Oriented<T, R, const N: usize> {
    pub geometry: T,
    pub rotation: R,
}

impl<T, R, const N: usize> Oriented<T, R, N> {
    #[must_use]
    #[inline]
    pub const fn new(geometry: T, rotation: R) -> Self {
        Self { geometry, rotation }
    }
}

impl<T: Volume<N>, R, const N: usize> Volume<N> for Oriented<T, R, N>
where
    T::Scalar: Copy + Zero + ClosedDiv + ClosedSub + ClosedAdd + ClosedMul,
    R: Rotation<N, Scalar = <T as Volume<N>>::Scalar> + Copy,
{
    type Scalar = T::Scalar;

    #[inline]
    fn origin(&self) -> Point<Self::Scalar, N> {
        self.geometry.origin()
    }

    #[inline]
    fn set_origin(&mut self, new_origin: Point<Self::Scalar, N>) {
        self.geometry.set_origin(new_origin);
    }

    #[inline]
    fn contains(&self, point: &Point<Self::Scalar, N>) -> bool {
        let origin = self.origin();
        let rotated_point = rotate_point_around(*point, origin, self.rotation);
        self.geometry.contains(&rotated_point)
    }

    #[inline]
    fn support_point(&self, direction: &Vector<Self::Scalar, N>) -> Point<Self::Scalar, N> {
        let unrotated_support = self.geometry.support_point(direction);
        rotate_point_around(unrotated_support, self.origin(), self.rotation)
    }
}

impl<T, R, const N: usize> RayIntersect<N> for Oriented<T, R, N>
where
    T: RayIntersect<N> + Volume<N, Scalar = <T as RayIntersect<N>>::Scalar>,
    R: Rotation<N, Scalar = <T as RayIntersect<N>>::Scalar> + Copy,
    <T as RayIntersect<N>>::Scalar:
        ClosedAdd + ClosedMul + ClosedSub + Copy + DivAssign + One + PartialEq + Zero + Sqrt,
{
    type Scalar = <T as RayIntersect<N>>::Scalar;

    fn intersection_with(
        &self,
        ray: &Ray<Self::Scalar, N>,
    ) -> Option<RayIntersection<Self::Scalar, N>> {
        let origin = self.geometry.origin();
        let inv_rotation = self.rotation.inverse();
        let rotated_ray = ray.rotated_around(origin, &inv_rotation);

        if let Some(intersection) = self.geometry.intersection_with(&rotated_ray) {
            let mut n = self.rotation.transform_vector(intersection.face_normal);
            n.normalize();
            let pt = {
                let v = intersection.position.vector_to(origin);
                let vr = self.rotation.transform_vector(v);
                intersection.position + v - vr
            };

            let dist = ray.origin.distance_to(pt);

            Some(RayIntersection {
                position: pt,
                face_normal: n,
                distance_to_intersection: dist,
                is_exterior_face: intersection.is_exterior_face,
            })
        } else {
            None
        }
    }

    #[inline]
    fn intersects(&self, ray: &Ray<Self::Scalar, N>) -> bool {
        let rotated_ray = ray.rotated_around(self.geometry.origin(), &self.rotation.inverse());
        self.geometry.intersects(&rotated_ray)
    }
}
