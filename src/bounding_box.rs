use crate::{
    Intersect, Intersection, RayIntersection, Volume,
    ball::Ball,
    line_overlap,
    plane::Plane,
    ray::{self, Ray},
    array_assume_init,
};
use approx::AbsDiffEq;
use core::{
    cmp::{Ordering, PartialOrd, max_by, min_by},
    marker::PhantomData,
    mem::MaybeUninit,
    ops::{AddAssign, DivAssign, Sub},
};
use vectral::{
    matrix::Matrix,
    point::Point,
    transform::{Transform, Translate},
    utils::{
        num::{
            Bounded, ClosedAdd, ClosedDiv, ClosedMul, ClosedNeg, ClosedSub, Float, One, Scalar,
            Signed, Sqrt, Zero,
        },
        product, zip_map,
    },
    vector::Vector,
};

#[repr(C)]
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct BoundingBox<T = f32, const N: usize = 3> {
    pub center: Point<T, N>,
    pub half_extents: Vector<T, N>,
}

impl<T, const N: usize> BoundingBox<T, N> {
    #[must_use]
    #[inline]
    pub const fn from_point_with_half_extents(center: Point<T, N>, half_extents: Vector<T, N>) -> Self {
        Self { center, half_extents }
    }

    #[deprecated = "use `BoundingBox::from_point_with_half_extents()` instead"]
    #[must_use]
    #[inline]
    pub const fn new(center: Point<T, N>, half_extents: Vector<T, N>) -> Self {
        Self { center, half_extents }
    }
}

impl<T: Copy + ClosedAdd + ClosedDiv + One, const N: usize> BoundingBox<T, N> {
    #[inline]
    pub fn from_point_with_extents(center: Point<T, N>, size: Vector<T, N>) -> Self {
        let half_size = {
            size / (T::ONE + T::ONE)
        };
        Self::from_point_with_half_extents(center, half_size)
    }
}

impl<T, const N: usize> Volume<N> for BoundingBox<T, N>
where
    T: AbsDiffEq + Scalar + Signed + Sqrt + PartialOrd<T::Epsilon>,
{
    type Scalar = T;
    #[inline]
    fn origin(&self) -> Point<Self::Scalar, N> {
        self.center
    }

    #[inline]
    fn set_origin(&mut self, new_origin: Point<Self::Scalar, N>) {
        self.center = new_origin;
    }

    #[inline]
    fn contains(&self, point: &Point<Self::Scalar, N>) -> bool {
        let (min, max) = self.to_min_max();

        let larger_than_min = min.iter().zip(point).all(|(min, p)| p >= min);
        let smaller_than_max = max.iter().zip(point).all(|(max, p)| p <= max);

        larger_than_min && smaller_than_max
    }

    #[inline]
    fn support_point(&self, direction: &Vector<T, N>) -> Point<T, N> {
        let dir = Vector::normalized_unchecked(*direction);
        let ray = Ray::new(self.center, dir);

        let mut max_intersection = None;
        for face in self.face_planes() {
            if let Some(intersection) = ray::Intersect::intersection_with(&face, &ray) {
                if max_intersection.is_none_or(|old_intersection: ray::Intersection<_, _>| {
                    intersection.distance_to_intersection
                        < old_intersection.distance_to_intersection
                }) {
                    max_intersection = Some(intersection);
                }
            }
        }

        max_intersection
            .map(|intersection| intersection.position)
            .unwrap_or(self.center)
    }
}

impl<T, const N: usize> BoundingBox<T, N>
where
    T: Copy + PartialOrd + ClosedSub + ClosedAdd + One + ClosedDiv,
{
    #[must_use]
    #[inline]
    pub fn corners(&self) -> Corners<'_, T, N> {
        let (min, max) = self.to_min_max();
        Corners {
            min,
            max,
            idx: 0,
            _boo: PhantomData,
        }
    }

    #[must_use]
    #[inline]
    pub fn subdivide(&self) -> Subdivisions<'_, T, N> {
        let (min, max) = self.to_min_max();
        Subdivisions {
            min_point: min,
            max_point: max,
            center: self.center,
            idx: 0,
            _boo: PhantomData,
        }
    }
}

impl<T: ClosedAdd + Copy, const N: usize> BoundingBox<T, N> {
    #[must_use]
    #[inline]
    pub fn size(&self) -> Vector<T, N> {
        self.half_extents + self.half_extents
    }
}

impl<T, const N: usize> BoundingBox<T, N> {
    #[must_use]
    #[inline]
    pub fn face_planes(&self) -> FacePlanes<'_, T, N> {
        FacePlanes {
            bounding_box: self,
            curr_plane: 0,
        }
    }
}

impl<T: Copy + Zero, const N: usize> BoundingBox<T, N> {
    #[must_use]
    #[inline]
    pub const fn enclosing(point: Point<T, N>) -> Self {
        Self::from_point_with_half_extents(point, Vector::ZERO)
    }
}

impl<T, const N: usize> BoundingBox<T, N>
where
    T: Copy + PartialOrd + ClosedSub + ClosedAdd + One + ClosedDiv,
{
    #[must_use]
    #[inline]
    pub fn all_corners(self) -> [Point<T, N>; 2usize.pow(N as u32)] {
        let mut corners_array = [const { MaybeUninit::uninit() }; 2usize.pow(N as u32)];

        for (i, corner) in self.corners().enumerate() {
            debug_assert!(1 < corners_array.len());
            unsafe {
                corners_array.get_unchecked_mut(i).write(corner);
            }
        }

        unsafe { array_assume_init(corners_array) }
    }
}

impl<T, const N: usize> BoundingBox<T, N>
where
    T: Copy + Zero + One + ClosedAdd + ClosedNeg + ClosedDiv,
{
    #[must_use]
    #[inline]
    pub fn all_face_planes(self) -> [Plane<T, N>; N * 2] {
        let mut face_planes_array = [const { MaybeUninit::uninit() }; N * 2];

        for (i, face_plane) in self.face_planes().enumerate() {
            unsafe {
                face_planes_array.get_unchecked_mut(i).write(face_plane);
            }
        }

        unsafe { array_assume_init(face_planes_array) }
    }
}

impl<T, const N: usize> BoundingBox<T, N>
where
    T: Copy + PartialOrd + ClosedSub + ClosedAdd + One + ClosedDiv + DivAssign,
{
    #[must_use]
    #[inline]
    pub fn all_subdivisions(self) -> [BoundingBox<T, N>; 2usize.pow(N as u32)] {
        let mut subdivs_array = [const { MaybeUninit::uninit() }; 2usize.pow(N as u32)];

        for (i, subdiv) in self.subdivide().enumerate() {
            unsafe {
                subdivs_array.get_unchecked_mut(i).write(subdiv);
            }
        }

        unsafe { array_assume_init(subdivs_array) }
    }
}

impl<T, const N: usize> BoundingBox<T, N>
where
    T: Copy + ClosedSub + ClosedAdd + One + ClosedDiv,
{
    #[must_use]
    #[inline]
    pub fn from_min_max_unchecked(min_point: Point<T, N>, max_point: Point<T, N>) -> Self {
        let size = max_point - min_point;
        let center = min_point + (size / (T::ONE + T::ONE));
        Self::from_point_with_extents(center, size)
    }
}

impl<T, const N: usize> BoundingBox<T, N>
where
    T: Copy + PartialOrd + ClosedSub + ClosedAdd + One + ClosedDiv,
{
    #[must_use]
    #[inline]
    pub fn from_min_max(min_point: Point<T, N>, max_point: Point<T, N>) -> Self {
        let min = Point::from(zip_map(min_point.into(), max_point.into(), |l, r| {
            min_by(l, r, Self::partial_cmp)
        }));

        let max = Point::from(zip_map(min_point.into(), max_point.into(), |l, r| {
            max_by(l, r, Self::partial_cmp)
        }));

        Self::from_min_max_unchecked(min, max)
    }

    #[must_use]
    #[inline]
    pub fn from_min_with_size(min_point: Point<T, N>, size: Vector<T, N>) -> Self {
        let half_extents = size / (T::ONE + T::ONE);
        Self::from_point_with_half_extents(min_point + half_extents, half_extents)
    }

    #[must_use]
    #[inline]
    pub fn to_min_max(self) -> (Point<T, N>, Point<T, N>) {
        (self.center - self.half_extents, self.center + self.half_extents)
    }
}

impl<T: Copy + PartialOrd + ClosedSub + ClosedAdd + One + ClosedDiv, const N: usize>
    BoundingBox<T, N>
{
    #[must_use]
    #[inline]
    pub fn union_with_point(&self, point: Point<T, N>) -> Self {
        let mut ret = *self;
        ret.enclose(point);
        ret
    }

    #[inline]
    pub fn enclose(&mut self, point: Point<T, N>) {
        let (mut min, mut max) = self.to_min_max();

        min = Point::from(zip_map(min.into(), point.into(), |l, r| {
            min_by(l, r, Self::partial_cmp)
        }));

        max = Point::from(zip_map(max.into(), point.into(), |l, r| {
            max_by(l, r, Self::partial_cmp)
        }));

        *self = Self::from_min_max(min, max);
    }

    #[must_use]
    #[inline]
    pub fn union_with_box(&self, other: &BoundingBox<T, N>) -> Self {
        let mut ret = *self;
        ret.enclose_box(other);
        ret
    }

    #[inline]
    pub fn enclose_box(&mut self, bbox: &BoundingBox<T, N>) {
        let (min, max) = bbox.to_min_max();
        self.enclose(min);
        self.enclose(max);
    }

    #[inline(always)]
    fn partial_cmp(x: &T, y: &T) -> Ordering {
        x.partial_cmp(y).unwrap_or(Ordering::Equal)
    }
}

impl<T: Copy + AddAssign + One + ClosedDiv, const N: usize> BoundingBox<T, N> {
    #[inline]
    pub fn expand(&mut self, delta: T) {
        let two = {
            let mut x = T::ONE;
            x += T::ONE;
            x
        };

        self.half_extents += Vector::splat(delta / two);
    }
}

impl<T: One + Copy + ClosedAdd + ClosedMul + ClosedSub> BoundingBox<T> {
    #[must_use]
    #[inline]
    pub fn surface_area(&self) -> T {
        let d = self.half_extents + self.half_extents;
        let two = T::ONE + T::ONE;
        two * ((d.x * d.y) + (d.x * d.z) + (d.y * d.z))
    }
}

impl<T: Copy + Sub + PartialOrd, const N: usize> BoundingBox<T, N>
where
    T::Output: PartialOrd + Sized,
{
    #[must_use]
    #[inline]
    pub fn max_extent(&self) -> Option<usize> {
        let size = self.half_extents.as_slice();
        let mut max = 0;

        while max < size.len() {
            let curr = &size[max];
            let rest = &size[..max + 1];
            if rest.iter().all(|elem| elem < curr) {
                return Some(max);
            }

            max += 1;
        }

        None
    }
}

impl<T, const N: usize> BoundingBox<T, N>
where
    T: ClosedAdd
        + Copy
        + DivAssign
        + ClosedDiv
        + ClosedSub
        + One
        + ClosedMul
        + Sqrt
        + Zero
        + PartialOrd,
{
    #[must_use]
    #[inline]
    pub fn bounding_circular(self) -> Ball<T, N> {
        let (_, max) = self.to_min_max();
        let radius = Point::distance_to(self.center, max);
        Ball::new(self.center, radius)
    }
}

impl<T, const N: usize> BoundingBox<T, N>
where
    T: Copy + PartialOrd + ClosedSub + ClosedAdd + One + ClosedDiv + ClosedMul + Zero + PartialOrd,
{
    #[must_use]
    #[inline]
    fn calculate_ray_distance_to_intersection_point(&self, ray: &Ray<T, N>) -> Option<(T, usize)> {
        // Taken from:
        // https://people.csail.mit.edu/amy/papers/box-jgt.pdf
        let (min, max) = self.to_min_max();
        let bounds = [min, max];
        let inv_direction = ray.direction.map(|x| T::ONE / x);
        let sign = inv_direction.map(|x| if x < T::ZERO { 1 } else { 0 });

        let calculate_min =
            |axis| (bounds[sign[axis]][axis] - ray.origin[axis]) * inv_direction[axis];

        let calculate_max =
            |axis| (bounds[1 - sign[axis]][axis] - ray.origin[axis]) * inv_direction[axis];

        let mut tmin = calculate_min(0);
        let mut tmax = calculate_max(0);
        let mut axis_tmin = 0;

        for i in 1..N {
            let nmin = calculate_min(i);
            let nmax = calculate_max(i);

            if (tmin > nmax) || (nmin > tmax) {
                return None;
            }

            if nmin > tmin {
                tmin = nmin;
                axis_tmin = i;
            }

            if nmax < tmax {
                tmax = nmax;
            }
        }

        // The min value *should* be the closest face to the ray's origin, and the
        // max value should be the furthest.
        //
        // If the distance is negative, there was no collision
        if tmin > T::ZERO {
            Some((tmin, axis_tmin))
        } else {
            None
        }
    }
}

impl<T: Copy + ClosedAdd, const N: usize> BoundingBox<T, N>
where
    T::Output: One + ClosedMul,
{
    #[must_use]
    #[inline]
    pub fn volume(&self) -> T::Output {
        let array = self.size().to_array();
        product(array)
    }
}

impl<T, const DIM: usize> Translate<DIM> for BoundingBox<T, DIM>
where
    T: Zero
        + PartialOrd
        + One
        + PartialEq
        + Copy
        + DivAssign
        + ClosedMul
        + ClosedSub
        + ClosedAdd
        + PartialOrd
        + ClosedDiv,
    Matrix<T, { DIM + 1 }, { DIM + 1 }>: Sized,
{
    type Scalar = T;

    #[inline]
    fn translated<Trans: Transform<DIM, Scalar = Self::Scalar>>(&self, transform: &Trans) -> Self {
        let mut ret = Self::from_point_with_half_extents(Zero::ZERO, Zero::ZERO);
        for corner in self.corners() {
            ret.enclose(corner.translated(transform));
        }

        ret
    }

    #[inline]
    fn translate_by<Trans: Transform<DIM, Scalar = Self::Scalar>>(&mut self, transform: &Trans) {
        *self = self.translated(transform);
    }
}

impl<T, const N: usize> Intersect<Self, N> for BoundingBox<T, N>
where
    T: AbsDiffEq
        + ClosedSub
        + ClosedAdd
        + ClosedDiv
        + Copy
        + PartialOrd
        + PartialOrd<T::Epsilon>
        + Zero
        + One
        + Bounded
        + Signed
        + Sqrt,
{
    #[inline]
    fn intersects(&self, rhs: &Self) -> bool {
        let (min0, max0) = self.to_min_max();
        let (min1, max1) = rhs.to_min_max();
        min0.iter()
            .zip(max0.iter())
            .zip(min1.iter())
            .zip(max1.iter())
            .all(|(((&min0, &max0), &min1), &max1)| line_overlap(min0, max0, min1, max1) > T::ZERO)
    }

    #[inline]
    fn intersection_with(&self, other: &Self) -> Option<Intersection<T, N>> {
        // @TODO(George): Do I want corner/edge normal detection?
        let (mut dim, mut best_depth) = (0, Bounded::MAX);

        let (min1, max1) = self.to_min_max();
        let (min2, max2) = other.to_min_max();

        for (d, (((min1, max1), min2), max2)) in min1
            .iter()
            .zip(max1.iter())
            .zip(min2.iter())
            .zip(max2.iter())
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
            let (c0, c1) = (self.center, other.center);
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

impl<T, const N: usize> ray::Intersect<N> for BoundingBox<T, N>
where
    T: AbsDiffEq
        + Copy
        + PartialOrd
        + ClosedSub
        + ClosedAdd
        + One
        + ClosedDiv
        + Zero
        + PartialOrd
        + PartialOrd<T::Epsilon>
        + Float,
{
    type Scalar = T;

    #[inline]
    fn intersection_with(
        &self,
        ray: &Ray<Self::Scalar, N>,
    ) -> Option<RayIntersection<Self::Scalar, N>> {
        self.calculate_ray_distance_to_intersection_point(ray)
            .map(|(distance, axis)| {
                let v = self.center.vector_to(ray.origin);
                let mut normal = {
                    let elem = if v[axis] < T::ZERO {
                        T::ONE.neg()
                    } else {
                        T::ONE
                    };
                    let mut v = Vector::<T, N>::ZERO;
                    v[axis] = elem;
                    v
                };

                let is_interior_face = self.contains(&ray.origin);
                if is_interior_face {
                    normal = -normal;
                }

                RayIntersection {
                    position: ray.point_at(distance),
                    face_normal: normal,
                    distance_to_intersection: distance,
                    is_exterior_face: !is_interior_face,
                }
            })
    }

    #[inline]
    fn intersects(&self, ray: &Ray<Self::Scalar, N>) -> bool {
        self.calculate_ray_distance_to_intersection_point(ray)
            .is_some()
    }
}

#[derive(Debug)]
pub struct Corners<'a, T, const N: usize> {
    min: Point<T, N>,
    max: Point<T, N>,
    idx: usize,
    _boo: PhantomData<&'a BoundingBox<T, N>>,
}

impl<'a, T: Copy, const N: usize> Iterator for Corners<'a, T, N> {
    type Item = Point<T, N>;
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.len() == 0 {
            return None;
        }

        let mut point = Point::uninit();
        for n in 0..N {
            if (self.idx >> n) & 1 == 1 {
                point[n].write(self.max[n]);
            } else {
                point[n].write(self.min[n]);
            }
        }

        let point = unsafe { Point::assume_init(point) };

        self.idx += 1;

        Some(point)
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let n = ExactSizeIterator::len(self);
        (n, Some(n))
    }
}

impl<T: Copy, const N: usize> ExactSizeIterator for Corners<'_, T, N> {
    #[inline]
    fn len(&self) -> usize {
        2usize.pow(N as u32).saturating_sub(self.idx)
    }
}

#[derive(Debug)]
pub struct Subdivisions<'a, T, const N: usize> {
    min_point: Point<T, N>,
    max_point: Point<T, N>,
    center: Point<T, N>,
    idx: usize,
    _boo: PhantomData<&'a BoundingBox<T, N>>,
}

impl<'a, T, const N: usize> Iterator for Subdivisions<'a, T, N>
where
    T: ClosedAdd + Copy + DivAssign + ClosedSub + ClosedDiv + One,
{
    type Item = BoundingBox<T, N>;
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.len() == 0 {
            return None;
        }

        let center = self.center;
        let (mut min, mut max) = (Point::splat(T::ONE), Point::splat(T::ONE));

        for n in 0..N {
            if (self.idx >> n) & 1 == 1 {
                min[n] = center[n];
                max[n] = self.max_point[n];
            } else {
                min[n] = self.min_point[n];
                max[n] = center[n];
            }
        }

        self.idx += 1;
        Some(BoundingBox::from_min_max_unchecked(min, max))
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let n = ExactSizeIterator::len(self);
        (n, Some(n))
    }
}

impl<T, const N: usize> ExactSizeIterator for Subdivisions<'_, T, N>
where
    T: ClosedAdd + Copy + DivAssign + ClosedSub + ClosedDiv + One,
{
    #[inline]
    fn len(&self) -> usize {
        2usize.pow(N as u32).saturating_sub(self.idx)
    }
}

#[derive(Debug)]
pub struct FacePlanes<'a, T, const N: usize> {
    bounding_box: &'a BoundingBox<T, N>,
    curr_plane: usize,
}

impl<'a, T, const N: usize> FacePlanes<'a, T, N> {
    const NUM_PLANES: usize = N * 2;
}

impl<'a, T, const N: usize> Iterator for FacePlanes<'a, T, N>
where
    T: Copy + Zero + One + ClosedAdd + ClosedNeg + ClosedDiv,
{
    type Item = Plane<T, N>;
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.curr_plane >= Self::NUM_PLANES {
            return None;
        }

        let i = self.curr_plane;
        let should_neg = i % 2 != 0;
        let idx = if should_neg { i - 1 } else { i } / 2;
        let mut normal = Vector::ZERO;
        normal[idx] = if should_neg { -T::ONE } else { T::ONE };

        self.curr_plane += 1;

        let half_size = self.bounding_box.size() / (T::ONE + T::ONE);
        let distance = (self.bounding_box.center + half_size)[idx];

        Some(Plane::new_unnormalised(normal, distance))
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.len();
        (len, Some(len))
    }
}

impl<'a, T, const N: usize> ExactSizeIterator for FacePlanes<'a, T, N>
where
    T: Copy + Zero + One + ClosedAdd + ClosedNeg + ClosedDiv,
{
    #[inline]
    fn len(&self) -> usize {
        Self::NUM_PLANES - self.curr_plane
    }
}

impl<'a, T, const N: usize> From<&'a BoundingBox<T, N>> for FacePlanes<'a, T, N> {
    #[inline]
    fn from(value: &'a BoundingBox<T, N>) -> Self {
        value.face_planes()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_ulps_eq;
    use core::ops::Neg;
    use vectral::{
        point::{Point2, Point3},
        vector::Vector3,
    };

    #[test]
    fn test_new() {
        let p1 = Point3::splat(-2.0f32);
        let p2 = Point3::splat(2.0);

        let box_1 = BoundingBox::from_min_max(p1, p2);
        let box_2 = BoundingBox::from_min_max(p2, p1);

        assert_eq!(box_1, box_2);
    }

    #[test]
    fn test_contains() {
        const SIZE: f64 = 4.0;
        const MAX: f64 = SIZE / 2.0;
        const MIN: f64 = -MAX;

        let bbox = BoundingBox::from_point_with_extents(Point3::ZERO, Vector::splat(SIZE + f64::EPSILON));

        let (mut x, mut y, mut z) = (MIN, MIN, MIN);

        const STEP: f64 = if cfg!(miri) { 1.0 } else { 0.25 };

        while x <= MAX {
            while y <= MAX {
                while z <= MAX {
                    let point = Point3::new([x, y, z]);
                    assert!(
                        bbox.contains(&point),
                        "Box {:?} erroneously does not contain point {:?}",
                        bbox,
                        point
                    );
                    z += STEP;
                }

                z = MIN;
                y += STEP;
            }

            z = MIN;
            y = MIN;
            x += STEP;
        }

        let outside = Point3::splat(MAX + STEP);
        assert!(
            !bbox.contains(&outside),
            "Box {:?} erroneously contains point {:?}",
            bbox,
            outside
        );

        let outside = Point3::splat(MIN - STEP);
        assert!(
            !bbox.contains(&outside),
            "Box {:?} erroneously contains point {:?}",
            bbox,
            outside
        );
    }

    #[test]
    fn test_volume() {
        let bbox = BoundingBox::from_min_max(Point3::splat(0.0), Point3::splat(1.0));
        assert_eq!(bbox.volume(), 1.0);

        let bbox = BoundingBox::from_min_max(Point3::splat(0.0), Point3::new([2.0, 4.0, 8.0]));
        assert_eq!(bbox.volume(), 64.0);
    }

    #[test]
    fn test_center() {
        let bbox = BoundingBox::from_min_max(Point3::splat(0.0), Point3::new([4.0, 6.0, 8.0]));
        assert_eq!(bbox.center, Point3::new([2.0, 3.0, 4.0]));

        let bbox = BoundingBox::from_min_max(Point3::splat(-2.0), Point3::splat(2.0));
        assert_eq!(bbox.center, Point3::splat(0.0));
    }

    #[test]
    fn test_corners() {
        let bbox = BoundingBox::from_min_max(Point2::splat(0.0), Point2::splat(1.0));
        let mut corners = bbox.corners();

        assert_eq!(corners.size_hint(), (4, Some(4)));

        assert_eq!(corners.next(), Some(Point2::from([0.0, 0.0])));
        assert_eq!(corners.next(), Some(Point2::from([1.0, 0.0])));

        assert_eq!(corners.size_hint(), (2, Some(2)));

        assert_eq!(corners.next(), Some(Point2::from([0.0, 1.0])));
        assert_eq!(corners.next(), Some(Point2::from([1.0, 1.0])));

        assert_eq!(corners.next(), None);

        assert_eq!(corners.size_hint(), (0, Some(0)));

        let bbox = BoundingBox::from_min_max(Point3::splat(2.0), Point3::splat(18.0));
        let all_corners = bbox.all_corners();
        let corners_iter = bbox.corners();

        let all_corners_iter = all_corners.into_iter();
        let mut n = 0;
        for (p0, p1) in all_corners_iter.zip(corners_iter) {
            assert_eq!(p0, p1);
            n += 1;
        }

        assert_eq!(n, 8);
    }

    #[test]
    fn test_subdivision() {
        let bbox = BoundingBox::from_min_max(Point2::splat(0.0), Point2::splat(4.0));
        let mut subdivs = bbox.subdivide();

        assert_eq!(subdivs.size_hint(), (4, Some(4)));

        assert_eq!(
            subdivs.next(),
            Some(BoundingBox::from_min_max(
                Point2::from([0.0, 0.0]),
                Point2::new([2.0, 2.0])
            ))
        );
        assert_eq!(
            subdivs.next(),
            Some(BoundingBox::from_min_max(
                Point2::from([2.0, 0.0]),
                Point2::new([4.0, 2.0])
            ))
        );

        assert_eq!(subdivs.size_hint(), (2, Some(2)));

        assert_eq!(
            subdivs.next(),
            Some(BoundingBox::from_min_max(
                Point2::from([0.0, 2.0]),
                Point2::new([2.0, 4.0])
            ))
        );
        assert_eq!(
            subdivs.next(),
            Some(BoundingBox::from_min_max(
                Point2::from([2.0, 2.0]),
                Point2::new([4.0, 4.0])
            ))
        );

        assert_eq!(subdivs.next(), None);

        assert_eq!(subdivs.size_hint(), (0, Some(0)));
    }

    #[test]
    fn test_raycast() {
        let box_1d = BoundingBox::from_point_with_extents(Point::new([1.0]), Vector::new([1.0]));
        let ray = Ray::new(Point::origin(), Vector::new([-1.0]));

        assert_eq!(
            box_1d.calculate_ray_distance_to_intersection_point(&ray),
            None
        );

        let box_3d = BoundingBox::from_point_with_extents(Point::origin(), Vector::splat(1.0));
        let ray = Ray::new(Point::new([0.0, 0.0, -5.0]), Vector::Z);

        let result = ray::Intersect::intersection_with(&box_3d, &ray).unwrap();
        assert_eq!(result.distance_to_intersection, 4.5);
        assert_eq!(result.face_normal, -Vector::<f64>::Z);
        assert!(result.is_exterior_face);
    }

    #[test]
    fn test_aabb_intersection() {
        let aabb_0 =
            BoundingBox::from_min_with_size(Point3::origin(), Vector3::new([5.0f32, 5.0, 1.0]));
        let aabb_1 = BoundingBox::from_min_with_size(
            Point3::new([0.0f32, 0.0, 0.7]),
            Vector3::new([1.0, 1.0, 1.0]),
        );

        assert!(aabb_0.intersects(&aabb_1));

        let intersection = aabb_0.intersection_with(&aabb_1).unwrap();
        assert_ulps_eq!(intersection.collision_normal, Vector3::Z);
        assert_ulps_eq!(intersection.penetration_depth, 0.3f32);

        // Check the negative z-axis
        let aabb_2 = BoundingBox::from_min_with_size(
            Point3::new([0.0f32, 0.0, -0.7]),
            Vector3::new([1.0, 1.0, 1.0]),
        );

        assert!(aabb_0.intersects(&aabb_2));

        let intersection = aabb_0.intersection_with(&aabb_2).unwrap();
        assert_ulps_eq!(intersection.collision_normal, Vector3::<f32>::Z.neg());
        assert_ulps_eq!(intersection.penetration_depth, 0.3f32);
    }

    #[test]
    fn test_aabb_intersection_miss() {
        let aabb_0 =
            BoundingBox::from_min_with_size(Point3::origin(), Vector3::new([3.0f32, 3.0, 3.0]));
        let aabb_1 = BoundingBox::from_min_with_size(
            Point3::new([6.0f32, 6.0, 6.0]),
            Vector3::new([1.0, 1.0, 1.0]),
        );

        assert!(!aabb_0.intersects(&aabb_1));
        assert!(aabb_0.intersection_with(&aabb_1).is_none());
    }

    #[test]
    fn test_planes() {
        let bbox = BoundingBox::from_point_with_extents(Point::origin(), Vector::new([4.0, 4.0, 4.0]));
        let faces = bbox.face_planes();
        assert_eq!(faces.len(), 6);

        let expected_normals = [
            Vector::<f64, 3>::X,
            -Vector::<f64, _>::X,
            Vector::Y,
            -Vector::<f64, _>::Y,
            Vector::Z,
            -Vector::<f64, _>::Z,
        ];

        for (normal, face) in expected_normals.iter().zip(faces) {
            assert_eq!(face.normal, *normal);
            assert_eq!(face.distance, 2.0);
        }
    }

    #[test]
    fn test_support_point() {
        use crate::ray::Intersect;
        let cube = BoundingBox::from_point_with_extents(Point::origin(), Vector::new([3.0, 8.0, 2.0]));
        let direction = Vector::new([0.3, 0.9, 0.5]).normalized();

        let support_point = cube.support_point(&direction);

        let ray = Ray::new(cube.center, direction);
        let mut smallest_dist = f64::MAX;
        for face_plane in cube.face_planes() {
            if let Some(intersection) = face_plane.intersection_with(&ray) {
                let t = intersection.distance_to_intersection;
                if t < smallest_dist {
                    smallest_dist = t;
                }
            }
        }

        assert_eq!(support_point, ray.point_at(smallest_dist));
    }
}
