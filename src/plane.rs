use crate::{ball::Ball, ray};
use approx::AbsDiffEq;
use core::{cmp::PartialOrd, ops::Neg};
#[cfg(feature = "serde")]
use serde_core::de;
use vectral::{
    point::Point,
    utils::num::{
        ClosedAdd, ClosedDiv, ClosedMul, ClosedNeg, ClosedSub, One, Scalar, Signed, Sqrt, Zero,
    },
    vector::Vector,
};
#[cfg(feature = "serde")]
use core::{marker::PhantomData};

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

impl<T: Copy + Scalar + Sqrt, const N: usize> Plane<T, N> {
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

#[cfg(feature = "serde")]
impl<T: serde_core::Serialize, const N: usize> serde_core::Serialize for Plane<T, N> {
    #[inline]
    fn serialize<S: serde_core::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        if serializer.is_human_readable() {
            use serde_core::ser::SerializeStruct;

            let mut s = serializer.serialize_struct("Plane", 2)?;
            s.serialize_field("normal", &self.normal)?;
            s.serialize_field("distance", &self.distance)?;
            s.end()
        } else {
            use serde_core::ser::SerializeTuple;

            let mut s = serializer.serialize_tuple(2)?;
            s.serialize_element(&self.normal)?;
            s.serialize_element(&self.distance)?;
            s.end()
        }
    }
}

#[cfg(feature = "serde")]
impl<'de, T: de::Deserialize<'de>, const N: usize> de::Deserialize<'de> for Plane<T, N> {
    #[inline]
    fn deserialize<D: de::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        use core::fmt;

        enum Field {
            Normal,
            Distance,
        }

        impl<'de> de::Deserialize<'de> for Field {
            #[inline]
            fn deserialize<D: de::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
                struct FieldVisitor;

                impl de::Visitor<'_> for FieldVisitor {
                    type Value = Field;

                    #[inline]
                    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                        formatter.write_str("`normal` or `distance`")
                    }

                    #[inline]
                    fn visit_str<E: de::Error>(self, v: &str) -> Result<Self::Value, E> {
                        match v {
                            "normal" => Ok(Field::Normal),
                            "distance" => Ok(Field::Distance),
                            _ => return Err(de::Error::unknown_field(v, FIELDS)),
                        }
                    }
                }

                deserializer.deserialize_identifier(FieldVisitor)
            }
        }

        struct Visitor<T, const N: usize>(PhantomData<Ball<T, N>>);

        impl<'de, T: de::Deserialize<'de>, const N: usize> de::Visitor<'de> for Visitor<T, N> {
            type Value = Plane<T, N>;
            #[inline]
            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("struct Plane")
            }

            #[inline]
            fn visit_seq<A: de::SeqAccess<'de>>(self, mut seq: A) -> Result<Self::Value, A::Error> {
                let normal = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(0, &self))?;
                let distance = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(1, &self))?;

                Ok(Plane::new_unnormalised(normal, distance))
            }

            #[inline]
            fn visit_map<A: de::MapAccess<'de>>(self, mut map: A) -> Result<Self::Value, A::Error> {
                let mut fields = (None, None);

                while let Some(key) = map.next_key()? {
                    match key {
                        Field::Normal => {
                            if fields.0.is_some() {
                                return Err(de::Error::duplicate_field("normal"));
                            }
                            fields.0 = Some(map.next_value::<Vector<T, N>>()?);
                        }
                        Field::Distance => {
                            if fields.1.is_some() {
                                return Err(de::Error::duplicate_field("distance"));
                            }
                            fields.1 = Some(map.next_value::<T>()?);
                        }
                    }
                }

                let (normal, distance) = fields;
                let normal = normal.ok_or_else(|| de::Error::missing_field("normal"))?;
                let distance = distance.ok_or_else(|| de::Error::missing_field("distance"))?;

                Ok(Plane::new_unnormalised(normal, distance))
            }
        }

        const FIELDS: &'static [&'static str] = &["normal", "distance"];

        if deserializer.is_human_readable() {
            deserializer.deserialize_struct("Plane", FIELDS, Visitor::<T, N>(PhantomData))
        } else {
            deserializer.deserialize_seq(Visitor::<T, N>(PhantomData))
        }
    }
}
