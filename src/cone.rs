use vectral::point::Point;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct Cone<T, const N: usize> {
    pub position: Point<T, N>,
    pub height: T,
    pub radius: T,
}