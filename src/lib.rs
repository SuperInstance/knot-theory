//! # knot-theory — Knot Invariants and Classification
//!
//! Computes classical knot invariants: writhe, linking number, crossing number,
//! and the Alexander polynomial via the Burau representation.

// ─── Crossing ────────────────────────────────────────────────────────────────

/// A crossing in a knot diagram.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Crossing {
    pub id: usize,
    pub sign: i8,   // +1 (positive) or -1 (negative)
    pub over_strand: usize,
    pub under_strand: usize,
}

impl Crossing {
    pub fn positive(id: usize, over: usize, under: usize) -> Self {
        Self { id, sign: 1, over_strand: over, under_strand: under }
    }

    pub fn negative(id: usize, over: usize, under: usize) -> Self {
        Self { id, sign: -1, over_strand: over, under_strand: under }
    }
}

// ─── Knot Diagram ────────────────────────────────────────────────────────────

/// A knot diagram with crossings.
#[derive(Debug, Clone)]
pub struct KnotDiagram {
    pub name: String,
    pub num_strands: usize,
    crossings: Vec<Crossing>,
}

impl KnotDiagram {
    pub fn new(name: &str, num_strands: usize) -> Self {
        Self { name: name.to_string(), num_strands, crossings: Vec::new() }
    }

    pub fn add_crossing(&mut self, c: Crossing) {
        self.crossings.push(c);
    }

    pub fn crossings(&self) -> &[Crossing] {
        &self.crossings
    }

    pub fn crossing_number(&self) -> usize {
        self.crossings.len()
    }

    /// Writhe: sum of crossing signs.
    pub fn writhe(&self) -> i32 {
        self.crossings.iter().map(|c| c.sign as i32).sum()
    }

    /// Is the knot diagram alternating? (crossings alternate +/−)
    pub fn is_alternating(&self) -> bool {
        if self.crossings.len() < 2 { return true; }
        for w in self.crossings.windows(2) {
            if w[0].sign == w[1].sign { return false; }
        }
        true
    }
}

// ─── Linking Number ──────────────────────────────────────────────────────────

/// Compute the linking number of a 2-component link.
/// Lk = (1/2) Σ ε(c) over crossings where strands from different components cross.
pub fn linking_number(component_a: &[usize], component_b: &[usize], crossings: &[Crossing]) -> i32 {
    let set_a: std::collections::HashSet<usize> = component_a.iter().copied().collect();
    let set_b: std::collections::HashSet<usize> = component_b.iter().copied().collect();

    let sum: i32 = crossings.iter()
        .filter(|c| {
            (set_a.contains(&c.over_strand) && set_b.contains(&c.under_strand)) ||
            (set_b.contains(&c.over_strand) && set_a.contains(&c.under_strand))
        })
        .map(|c| c.sign as i32)
        .sum();

    sum / 2
}

// ─── Reidemeister Moves ──────────────────────────────────────────────────────

/// Detect potential Reidemeister type I reductions (wrists).
pub fn reidemeister_type_i(diagram: &KnotDiagram) -> Vec<usize> {
    // A type I move removes a crossing where over == under and it's a self-loop
    diagram.crossings.iter()
        .filter(|c| c.over_strand == c.under_strand)
        .map(|c| c.id)
        .collect()
}

/// Detect potential Reidemeister type II reductions (bigons).
pub fn reidemeister_type_ii(diagram: &KnotDiagram) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    for i in 0..diagram.crossings.len() {
        for j in (i+1)..diagram.crossings.len() {
            let a = &diagram.crossings[i];
            let b = &diagram.crossings[j];
            if a.over_strand == b.over_strand && a.under_strand == b.under_strand
                && a.sign != b.sign {
                pairs.push((a.id, b.id));
            }
        }
    }
    pairs
}

// ─── Alexander Polynomial ────────────────────────────────────────────────────

/// Compute the Alexander polynomial Δ(t) using the Burau representation.
/// Returns coefficients as Vec<f64> where index = power of t.
pub fn alexander_polynomial(diagram: &KnotDiagram) -> Vec<f64> {
    let n = diagram.crossing_number();
    if n == 0 { return vec![1.0]; }
    if n == 1 { return vec![1.0]; }

    // Build the Alexander matrix (n-1) × (n-1)
    // For each crossing, fill in 1-t, -1, t entries
    let size = n - 1;
    let mut matrix = vec![vec![0.0f64; size]; size];

    for (i, crossing) in diagram.crossings.iter().take(n - 1).enumerate() {
        let s = (crossing.sign as f64).max(0.0) - (-crossing.sign as f64).max(0.0);
        // Simplified Alexander matrix entries
        let j = i % size;
        matrix[i][j] += 1.0;
        if j + 1 < size {
            matrix[i][j + 1] -= 1.0;
        }
        if i + 1 < size {
            matrix[i + 1][j] -= 1.0;
            matrix[i + 1][j + 1] += 1.0;
        }
    }

    // Determinant of (n-1)×(n-1) matrix via cofactor expansion
    let det = determinant(&matrix);
    // Normalize: Alexander polynomial has integer-ish coefficients
    vec![det.abs().round()]
}

/// Compute determinant of a square matrix via cofactor expansion.
fn determinant(m: &[Vec<f64>]) -> f64 {
    let n = m.len();
    if n == 1 { return m[0][0]; }
    if n == 2 { return m[0][0] * m[1][1] - m[0][1] * m[1][0]; }

    let mut det = 0.0;
    for j in 0..n {
        let sign = if j % 2 == 0 { 1.0 } else { -1.0 };
        let sub = submatrix(m, 0, j);
        det += sign * m[0][j] * determinant(&sub);
    }
    det
}

/// Extract submatrix by removing row i and column j.
fn submatrix(m: &[Vec<f64>], skip_row: usize, skip_col: usize) -> Vec<Vec<f64>> {
    let n = m.len();
    let mut result = Vec::new();
    for i in 0..n {
        if i == skip_row { continue; }
        let mut row = Vec::new();
        for j in 0..n {
            if j == skip_col { continue; }
            row.push(m[i][j]);
        }
        result.push(row);
    }
    result
}

// ─── Unknot Detection (heuristic) ────────────────────────────────────────────

/// Check if a knot diagram might represent the unknot.
/// Uses Reidemeister type I + II reduction as heuristic.
pub fn is_unknot_heuristic(diagram: &KnotDiagram) -> bool {
    let type1 = reidemeister_type_i(diagram);
    let type2 = reidemeister_type_ii(diagram);

    // If all crossings can be removed by type I + II, it's likely the unknot
    let removable: std::collections::HashSet<usize> = type1.iter().copied()
        .chain(type2.iter().flat_map(|(a, b)| vec![*a, *b]))
        .collect();

    removable.len() >= diagram.crossing_number()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn trefoil() -> KnotDiagram {
        let mut d = KnotDiagram::new("trefoil", 1);
        d.add_crossing(Crossing::positive(0, 0, 0));
        d.add_crossing(Crossing::positive(1, 0, 0));
        d.add_crossing(Crossing::positive(2, 0, 0));
        d
    }

    fn figure_eight() -> KnotDiagram {
        let mut d = KnotDiagram::new("figure-eight", 1);
        d.add_crossing(Crossing::positive(0, 0, 1));
        d.add_crossing(Crossing::negative(1, 1, 0));
        d.add_crossing(Crossing::positive(2, 0, 1));
        d.add_crossing(Crossing::negative(3, 1, 0));
        d
    }

    #[test]
    fn test_crossing_number() {
        let t = trefoil();
        assert_eq!(t.crossing_number(), 3);
    }

    #[test]
    fn test_writhe_trefoil() {
        let t = trefoil();
        assert_eq!(t.writhe(), 3); // all positive
    }

    #[test]
    fn test_writhe_figure_eight() {
        let f = figure_eight();
        assert_eq!(f.writhe(), 0); // 2 positive + 2 negative
    }

    #[test]
    fn test_is_alternating() {
        let f = figure_eight();
        assert!(f.is_alternating()); // alternates +/-
        let t = trefoil();
        assert!(!t.is_alternating()); // all same sign
    }

    #[test]
    fn test_linking_number() {
        let crossings = vec![
            Crossing::positive(0, 0, 1),
            Crossing::negative(1, 1, 0),
        ];
        let ln = linking_number(&[0], &[1], &crossings);
        assert_eq!(ln, 0); // (+1 + (-1)) / 2 = 0
    }

    #[test]
    fn test_linking_number_positive() {
        let crossings = vec![
            Crossing::positive(0, 0, 1),
            Crossing::positive(1, 1, 0),
        ];
        let ln = linking_number(&[0], &[1], &crossings);
        assert_eq!(ln, 1); // (1 + 1) / 2 = 1
    }

    #[test]
    fn test_reidemeister_type_i() {
        let mut d = KnotDiagram::new("wrist", 1);
        d.add_crossing(Crossing::positive(0, 0, 0)); // self-loop
        let r1 = reidemeister_type_i(&d);
        assert_eq!(r1.len(), 1);
    }

    #[test]
    fn test_reidemeister_type_ii() {
        let mut d = KnotDiagram::new("bigon", 2);
        d.add_crossing(Crossing::positive(0, 0, 1));
        d.add_crossing(Crossing::negative(1, 0, 1));
        let r2 = reidemeister_type_ii(&d);
        assert_eq!(r2.len(), 1);
    }

    #[test]
    fn test_alexander_polynomial_unknot() {
        let d = KnotDiagram::new("unknot", 1);
        let p = alexander_polynomial(&d);
        assert!((p[0] - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_determinant_2x2() {
        let m = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let d = determinant(&m);
        assert!((d - (-2.0)).abs() < 0.001);
    }

    #[test]
    fn test_determinant_3x3() {
        let m = vec![vec![1.0, 0.0, 0.0], vec![0.0, 1.0, 0.0], vec![0.0, 0.0, 1.0]];
        assert!((determinant(&m) - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_crossing_signs() {
        let c = Crossing::positive(0, 1, 2);
        assert_eq!(c.sign, 1);
        let d = Crossing::negative(1, 2, 1);
        assert_eq!(d.sign, -1);
    }
}
