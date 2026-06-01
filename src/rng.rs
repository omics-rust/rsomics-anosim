/// SplitMix64 — a tiny, fast, well-distributed PRNG used only to drive the
/// permutation test. Each permutation seeds its own stream, so the p-value is
/// deterministic for a given `--seed` regardless of thread count.
pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    pub fn new(seed: u64) -> Self {
        SplitMix64 { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform integer in `0..bound` via Lemire's multiply-shift with rejection
    /// of the biased low zone.
    pub fn bounded(&mut self, bound: u64) -> u64 {
        let mut x = self.next_u64();
        let mut m = (x as u128) * (bound as u128);
        let mut lo = m as u64;
        if lo < bound {
            let threshold = bound.wrapping_neg() % bound;
            while lo < threshold {
                x = self.next_u64();
                m = (x as u128) * (bound as u128);
                lo = m as u64;
            }
        }
        (m >> 64) as u64
    }
}
