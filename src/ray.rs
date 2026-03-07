use core::ops::{Add, DivAssign, Mul, RangeBounds};
use vectral::{
    matrix::Matrix,
    point::Point,
    rotation::Rotation,
    transform::{Transform, Translate},
    utils::num::{ClosedAdd, ClosedMul, ClosedSub, One, Zero},
    vector::Vector,
};

#[derive(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub struct Ray<T = f32, const N: usize = 3> {
    pub origin: Point<T, N>,
    pub direction: Vector<T, N>,
}

impl<T, const N: usize> Ray<T, N> {
    #[inline]
    pub const fn new(origin: Point<T, N>, direction: Vector<T, N>) -> Self {
        Self { origin, direction }
    }
}

impl<T: Add<Output = T> + Copy + Mul<Output = T>, const N: usize> Ray<T, N> {
    #[inline]
    pub fn point_at(&self, time: T) -> Point<T, N> {
        self.origin + (self.direction * time)
    }
}

impl<T, const DIM: usize> Translate<DIM> for Ray<T, DIM>
where
    T: Zero + One + PartialEq + Copy + DivAssign + ClosedMul + ClosedAdd,
    Matrix<T, { DIM + 1 }, { DIM + 1 }>: Sized,
{
    type Scalar = T;

    #[inline]
    fn translated<Trans: Transform<DIM, Scalar = Self::Scalar>>(&self, transform: &Trans) -> Self {
        Ray {
            origin: self.origin.translated(transform),
            direction: self.direction.translated(transform),
        }
    }

    #[inline]
    fn translate_by<Trans: Transform<DIM, Scalar = Self::Scalar>>(&mut self, transform: &Trans) {
        self.origin.translate_by(transform);
        self.direction.translate_by(transform);
    }
}

impl<T, const N: usize> Ray<T, N>
where
    T: Zero + One + PartialEq + Copy + DivAssign + ClosedMul + ClosedAdd + ClosedSub,
{
    #[must_use]
    #[inline]
    pub fn rotated<R: Rotation<N, Scalar = T>>(self, rotation: R) -> Self {
        let new_dir = rotation.transform_vector(self.direction);
        Self {
            origin: self.origin,
            direction: new_dir,
        }
    }

    #[must_use]
    #[inline]
    pub fn rotated_around<R: ?Sized + Rotation<N, Scalar = T>>(
        self,
        point: Point<T, N>,
        rotation: &R,
    ) -> Self {
        let v = self.origin.vector_to(point);
        let v_rotated = rotation.transform_vector(v);
        let new_dir = rotation.transform_vector(self.direction);

        Self {
            origin: self.origin + (v - v_rotated),
            direction: new_dir,
        }
    }
}

pub trait Intersect<const DIM: usize> {
    type Scalar;

    #[must_use]
    fn intersection_with(
        &self,
        ray: &Ray<Self::Scalar, DIM>,
    ) -> Option<Intersection<Self::Scalar, DIM>>;

    #[must_use]
    #[inline]
    fn intersects(&self, ray: &Ray<Self::Scalar, DIM>) -> bool {
        self.intersection_with(ray).is_some()
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub struct Intersection<T, const N: usize> {
    pub position: Point<T, N>,
    pub face_normal: Vector<T, N>,
    pub distance_to_intersection: T,
    pub is_exterior_face: bool,
}

impl<T, const N: usize> Intersection<T, N> {
    #[must_use]
    #[inline]
    pub fn is_within<U, R: RangeBounds<U>>(&self, range: R) -> bool
    where
        T: PartialOrd<U>,
        U: PartialOrd<T>,
    {
        range.contains(&self.distance_to_intersection)
    }
}

#[cfg(test)]
mod tests {
    use crate::bounding_box::BoundingBox;
    use crate::{oriented::Oriented, ray::Intersect as RayIntersect, ray::Ray};
    use vectral::{
        point::Point,
        rotation::{angle::Angle, quaternion::Quaternion},
        vector::Vector,
    };

    #[test]
    fn test_hit_oriented() {
        let bbox = {
            let bbox = BoundingBox::<_, 3>::new(Point::origin(), Vector::splat(3.0));
            Oriented::new(
                bbox,
                Quaternion::from_angle_axis(Angle::Degrees(45.0), Vector::Y),
            )
        };

        let ray = Ray::new(Point::new([-0.1, 0.0, -5.0]), Vector::new([0.0, 0.0, 1.0]));

        let result = bbox.intersection_with(&ray).unwrap();

        let angle = Vector::angle_between(result.face_normal, ray.direction);

        assert_eq!(angle.in_degrees(), Angle::Degrees(135.0).in_degrees());
        assert!(result.is_exterior_face);
    }
}
