use vectral::point::Point;

#[derive(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub struct Capsule<T, const N: usize> {
    pub origin: Point<T, N>,
    pub radius: T,
    pub half_height: T,
}

impl<T, const N: usize> Capsule<T, N> {
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
