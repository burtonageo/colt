use crate::{
    Intersect, Intersection, Volume,
    bounding_box::BoundingBox,
    line_overlap,
    plane::{Plane, ball_inside_plane, ball_intersects_plane_point},
    ray::{self, Ray},
};
use vectral::{
    matrix::Matrix,
    point::Point,
    transform::{Transform, Translate},
    utils::num::{Bounded, ClosedAdd, ClosedMul, ClosedSub, One, Scalar, Signed, Sqrt, Zero},
    vector::Vector,
};

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct Ball<T, const N: usize> {
    pub center: Point<T, N>,
    pub radius: T,
}

pub type Circle<T = f32> = Ball<T, 2>;
pub type Sphere<T = f32> = Ball<T, 3>;

impl<T, const N: usize> Ball<T, N> {
    #[must_use]
    #[inline]
    pub const fn new(center: Point<T, N>, radius: T) -> Self {
        Self { radius, center }
    }
}

impl<T: ClosedAdd + ClosedMul + One, const N: usize> Ball<T, N> {
    #[must_use]
    #[inline]
    pub fn diameter(self) -> T {
        let two = T::ONE + T::ONE;
        self.radius * two
    }
}

impl<T: Scalar + Sqrt + Copy, const N: usize> Volume<N> for Ball<T, N> {
    type Scalar = T;

    #[inline]
    fn contains(&self, point: &Point<T, N>) -> bool {
        let dist = (self.center - *point).len_squared();
        dist <= (self.radius * self.radius)
    }

    #[inline]
    fn origin(&self) -> Point<T, N> {
        self.center
    }

    #[inline]
    fn set_origin(&mut self, new_origin: Point<T, N>) {
        self.center = new_origin;
    }

    #[inline]
    fn support_point(&self, direction: &Vector<T, N>) -> Point<T, N> {
        let dir_normalized = Vector::normalized_unchecked(*direction);
        self.center + (dir_normalized * self.radius)
    }
}

impl<T: Scalar + Copy, const DIM: usize> Translate<DIM> for Ball<T, DIM>
where
    Matrix<T, { DIM + 1 }, { DIM + 1 }>: Sized,
{
    type Scalar = T;

    #[inline]
    fn translated<Trans: Transform<DIM, Scalar = Self::Scalar>>(&self, transform: &Trans) -> Self {
        Ball {
            center: self.center.translated(transform),
            radius: self.radius,
        }
    }

    #[inline]
    fn translate_by<Trans: Transform<DIM, Scalar = Self::Scalar>>(&mut self, transform: &Trans) {
        self.center.translate_by(transform);
    }
}

impl<T, const N: usize> Intersect<Self, N> for Ball<T, N>
where
    T: Scalar + Sqrt,
{
    #[inline]
    fn intersects(&self, other: &Self) -> bool {
        let dist = (self.center - other.center).len_squared();
        let rad1_sq = self.radius * self.radius;
        let rad2_sq = other.radius * other.radius;
        dist < (rad1_sq + rad2_sq)
    }

    #[inline]
    fn intersection_with(&self, other: &Self) -> Option<Intersection<T, N>> {
        if !self.intersects(other) {
            return None;
        }

        let collision_normal = self.center - other.center;
        let penetration_depth = collision_normal.clone().len();

        Some(Intersection {
            collision_normal,
            penetration_depth,
        })
    }
}

impl<T, const N: usize> Intersect<BoundingBox<T, N>, N> for Ball<T, N>
where
    T: Bounded + Copy + Clone + ClosedAdd + ClosedSub + ClosedMul + Scalar + Zero + Signed + Sqrt,
{
    #[inline]
    fn intersects(&self, rhs: &BoundingBox<T, N>) -> bool {
        // Check if the sphere is inside the box
        if rhs.face_planes().all(|face| ball_inside_plane(self, &face)) {
            return true;
        }

        // Check if the sphere intersects any of the box faces
        for (axis, face) in rhs.face_planes().enumerate() {
            let axis = axis / 2;
            if let Some((point, radius)) = ball_intersects_plane_point(self, &face) {
                if rhs
                    .face_planes()
                    .enumerate()
                    .filter_map(|(axis2, face)| {
                        let axis2 = axis2 / 2;
                        if axis2 == axis { None } else { Some(face) }
                    })
                    .all(|face| Plane::distance_to(&face, &point) <= radius)
                {
                    return true;
                }
            }
        }

        false
    }

    #[inline]
    fn intersection_with(&self, other: &BoundingBox<T, N>) -> Option<Intersection<T, N>> {
        let (mut dim, mut best_depth) = (0, <T as Bounded>::MAX);

        let (self_min_point, self_max_point) = {
            let n = Vector::new([self.radius; N]);
            (self.center - n, self.center + n)
        };

        let (min, max) = other.to_min_max();

        for (d, (((min1, max1), min2), max2)) in self_min_point
            .iter()
            .zip(self_max_point.iter())
            .zip(min.iter())
            .zip(max.iter())
            .enumerate()
        {
            let overlap = line_overlap(*min1, *max1, *min2, *max2);

            if overlap <= T::ZERO {
                // no overlap - it's not intersecting, so return failure
                return None;
            } else if overlap < best_depth {
                // Found a smaller overlap - this is a better candidate
                // for the intersecting face.
                dim = d;
                best_depth = overlap;
            }
        }

        let collision_normal = {
            let (c0, c1) = (&self.center, &other.center);
            // Ensure that the sign of the axis vector is correct
            let mut arr = [T::ZERO; _];
            arr[dim] = (c1[dim] - c0[dim]).sign();
            Vector::new(arr)
        };

        Some(Intersection {
            collision_normal,
            penetration_depth: best_depth,
        })
    }
}

impl<T: Scalar + Sqrt + Signed, const DIM: usize> Ball<T, DIM> {
    #[inline]
    fn calculate_ray_distance_to_intersection_point(&self, ray: &Ray<T, DIM>) -> Option<T> {
        let oc = self.center - ray.origin;
        let a = ray.direction.len_squared();
        let h = Vector::dot(ray.direction, oc);
        let c = oc.len_squared() - (self.radius * self.radius);
        let discriminant = (h * h) - (a * c);

        if discriminant < T::ZERO {
            return None;
        }

        let sqrtd = Sqrt::sqrt(discriminant);

        let mut root = (h - sqrtd) / a;
        if root < T::ZERO {
            root = (h + sqrtd) / a;
            if root < T::ZERO {
                return None;
            }
        }

        Some(root)
    }
}

impl<T: Scalar + Sqrt + Signed, const DIM: usize> ray::Intersect<DIM> for Ball<T, DIM> {
    type Scalar = T;

    #[inline]
    fn intersection_with(
        &self,
        ray: &Ray<Self::Scalar, DIM>,
    ) -> Option<ray::Intersection<Self::Scalar, DIM>> {
        self.calculate_ray_distance_to_intersection_point(ray)
            .map(|t| {
                let hit_pos = ray.point_at(t);
                let mut hit_norm = (hit_pos - self.center).normalized_unchecked();

                let is_front_face = Vector::dot(ray.direction, hit_norm) < T::ZERO;
                if !is_front_face {
                    hit_norm = -hit_norm;
                }

                ray::Intersection {
                    position: hit_pos,
                    face_normal: hit_norm,
                    distance_to_intersection: t,
                    is_exterior_face: is_front_face,
                }
            })
    }

    #[inline]
    fn intersects(&self, ray: &Ray<Self::Scalar, DIM>) -> bool {
        self.calculate_ray_distance_to_intersection_point(ray)
            .is_some()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_ulps_eq;
    use core::ops::Neg;
    use vectral::{point::Point3, vector::Vector3};

    #[test]
    fn test_sphere_raycast() {
        let ray = Ray::new(
            Point3::new([0.0f32, 0.0, 0.0]),
            Vector3::new([0.0, 1.0, 0.0]),
        );
        let sph = Sphere::new(Point3::new([0.0f32, 6.0, 0.0]), 5.0);

        let intersection = ray::Intersect::intersection_with(&sph, &ray).unwrap();
        assert_ulps_eq!(intersection.position, Point3::new([0.0, 1.0, 0.0]));
        assert_ulps_eq!(intersection.face_normal, Vector3::new([0.0, -1.0, 0.0]));
    }

    #[test]
    fn test_sphere_raycast_miss() {
        let ray = Ray::new(
            Point3::new([0.0f32, 0.0, 0.0]),
            Vector3::new([0.0, 1.0, 0.0]),
        );
        let sph = Sphere::new(Point3::new([30.0f32, 30.0, 30.0]), 1.0);
        assert!(!ray::Intersect::intersects(&sph, &ray));
    }

    #[test]
    fn test_sphere_sphere_intersect() {
        let sph0 = Sphere::new(Point3::<f32>::origin(), 3.0);
        let sph1 = Sphere::new(Point3::new([1.0f32, 0.0, 0.0]), 2.0);
        assert!(sph0.intersects(&sph1));

        let i = sph0.intersection_with(&sph1).unwrap();
        assert_ulps_eq!(i.penetration_depth, 1.0f32);
        assert_ulps_eq!(i.collision_normal, Vector3::<f32>::X.neg());
    }

    #[test]
    fn test_sphere_sphere_intersect_miss() {
        let sph0 = Sphere::new(Point3::<f32>::origin(), 3.0);
        let sph1 = Sphere::new(Point3::new([5.0f32, 5.0, 1.0]), 2.0);
        assert!(!sph0.intersects(&sph1));
    }

    #[test]
    fn test_sphere_box_intersect() {
        let sph = Sphere::new(Point3::origin(), 6.0);
        let bx = BoundingBox::from_min_with_size(
            Point3::<f32>::new([0.0, 4.0, 0.0]),
            Vector3::new([3.0, 3.0, 3.0]),
        );
        assert!(sph.intersects(&bx));

        let i = sph.intersection_with(&bx).unwrap();
        assert_ulps_eq!(i.penetration_depth, 2.0f32);
        assert_ulps_eq!(i.collision_normal, Vector3::Y);

        assert_eq!(i, sph.intersection_with(&bx).unwrap());
    }

    #[test]
    fn test_box_enclosing_sphere_intersect() {
        let sph = Sphere::new(Point3::new([50.0f32, 50.0, 50.0]), 1.0);
        let bx = BoundingBox::from_min_with_size(
            Point3::<f32>::origin(),
            Vector3::new([100.0, 100.0, 100.0]),
        );
        assert!(sph.intersects(&bx));

        let sph = Sphere::new(Point3::origin(), 100.0);
        let bx = BoundingBox::from_point_with_extents(Point3::origin(), Vector::splat(1.0));
        assert!(sph.intersects(&bx));
    }

    #[test]
    fn test_sphere_enclosing_box_intersect() {
        let sph = Sphere::new(Point3::origin(), 100.0);
        let bx = BoundingBox::from_min_with_size(
            Point3::<f32>::origin(),
            Vector3::new([50.0, 50.0, 50.0]),
        );
        assert!(sph.intersects(&bx));
    }

    #[test]
    fn test_sphere_box_intersect_miss() {
        let sph = Sphere::new(Point3::new([-5.0f32, -5.0, -5.0]), 3.0);
        let bx = BoundingBox::from_min_with_size(Point3::origin(), Vector3::new([2.0, 2.0, 2.0]));
        assert!(!sph.intersects(&bx));
    }
}
