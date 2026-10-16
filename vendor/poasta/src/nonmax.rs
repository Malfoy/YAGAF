//! The subset of nonmax::NonMaxU32 needed by POASTA, implemented using std.
//! Values keep their natural ordering while Option uses the nonzero niche.
use std::cmp::Ordering;
use std::num::NonZeroU32;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NonMaxU32(NonZeroU32);

impl From<NonMaxU32> for u32 {
    fn from(value: NonMaxU32) -> Self {
        value.get()
    }
}
impl std::fmt::Display for NonMaxU32 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.get().fmt(f)
    }
}

impl NonMaxU32 {
    pub fn new(value: u32) -> Option<Self> {
        NonZeroU32::new(u32::MAX - value).map(Self)
    }
    pub fn get(self) -> u32 {
        u32::MAX - self.0.get()
    }
}
impl PartialOrd for NonMaxU32 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for NonMaxU32 {
    fn cmp(&self, other: &Self) -> Ordering {
        self.get().cmp(&other.get())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn bounds_order_and_niche() {
        assert!(NonMaxU32::new(u32::MAX).is_none());
        assert_eq!(NonMaxU32::new(0).unwrap().get(), 0);
        assert_eq!(NonMaxU32::new(u32::MAX - 1).unwrap().get(), u32::MAX - 1);
        assert!(NonMaxU32::new(0).unwrap() < NonMaxU32::new(1).unwrap());
        assert_eq!(std::mem::size_of::<Option<NonMaxU32>>(), 4);
    }
}
