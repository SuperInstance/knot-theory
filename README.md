# knot-theory

> **A knot is a closed loop in 3D space. This library measures it.**

[![crates.io](https://img.shields.io/crates/v/knot-theory.svg)](https://crates.io/crates/knot-theory)
[![docs.rs](https://docs.rs/knot-theory/badge.svg)](https://docs.rs/knot-theory)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Classical knot invariants for Rust: crossing number, writhe, linking number, Reidemeister move detection, and Alexander polynomial computation via the Burau representation.

## Why Knot Theory?

Knots appear everywhere in science and engineering:
- **DNA topology**: enzymes change DNA knot type during replication
- **Fluid dynamics**: vortex lines form knots in turbulent flow
- **Quantum field theory**: Wilson loops are knot invariants
- **Robotics**: motion planning must avoid configuration-space knots
- **Agent systems**: multi-agent coordination can be modeled as knot unlinking

A knot invariant is a quantity that doesn't change under continuous deformation. If two knots have different invariants, they are definitely different knots. If they have the same invariants... they *might* be the same knot (some invariants can't distinguish certain pairs).

## Quick Start

```rust
use knot_theory::{KnotDiagram, Crossing, linking_number, alexander_polynomial};

// Build a trefoil knot (3_1) — the simplest non-trivial knot
let mut trefoil = KnotDiagram::new("trefoil", 2);
trefoil.add_crossing(Crossing::positive(0, 0, 1));
trefoil.add_crossing(Crossing::positive(1, 0, 1));
trefoil.add_crossing(Crossing::positive(2, 0, 1));

// Classical invariants
println!("Crossing number: {}", trefoil.crossing_number()); // 3
println!("Writhe: {}", trefoil.writhe());                   // +3
println!("Is alternating: {}", trefoil.is_alternating());   // true

// Alexander polynomial — a powerful knot invariant
let poly = alexander_polynomial(&trefoil);
println!("Alexander polynomial coefficients: {:?}", poly);

// Check for unknot (heuristic)
println!("Is unknot? {}", is_unknot_heuristic(&trefoil)); // false
```

## API Reference

### `Crossing`

A single crossing in a knot diagram. Every crossing has:
- **Sign**: +1 (positive/right-handed) or -1 (negative/left-handed)
- **Over strand**: which strand passes over
- **Under strand**: which strand passes under

```rust
let c = Crossing::positive(0, 1, 2); // id=0, over=strand1, under=strand2
```

### `KnotDiagram`

A planar diagram of a knot, represented as a collection of crossings.

| Method | Returns | Description |
|--------|---------|-------------|
| `new(name, strands)` | `KnotDiagram` | Create empty diagram |
| `add_crossing(c)` | | Add a crossing |
| `crossing_number()` | `usize` | Number of crossings |
| `writhe()` | `i32` | Sum of crossing signs |
| `is_alternating()` | `bool` | Do crossings alternate sign? |

### `linking_number(a, b, crossings)`

The **linking number** Lk(L₁, L₂) counts how many times two closed curves wrap around each other. For the Hopf link, Lk = ±1. For unlinked circles, Lk = 0.

```
Lk = (1/2) Σ ε(c)   over crossings where strands from different components meet
```

### `reidemeister_type_i` / `reidemeister_type_ii`

Detect opportunities for **Reidemeister moves** — local deformations that don't change the knot type:
- **Type I**: Remove a twist (reduces crossing number by 1)
- **Type II**: Slide one strand over another (reduces crossing number by 2)

These are the foundation of knot simplification algorithms.

### `alexander_polynomial(diagram)`

Computes the **Alexander polynomial** Δ(t) via the Burau representation. This is one of the oldest and most studied knot invariants:

- **Unknot**: Δ(t) = 1
- **Trefoil (3₁)**: Δ(t) = t² - t + 1
- **Figure-eight (4₁)**: Δ(t) = -t² + 3t - 1

The Alexander polynomial cannot distinguish the knot from its mirror image, but it distinguishes many knot pairs that crossing number alone cannot.

### `is_unknot_heuristic(diagram)`

A fast heuristic check: if the writhe is zero and no Reidemeister moves are possible, the diagram *might* represent the unknot. Not definitive — use with caution.

## Architecture

```
KnotDiagram
├── Crossing (id, sign, over/under strands)
├── writhe()           → i32
├── is_alternating()   → bool
├── crossing_number()  → usize
│
├── linking_number()   → i32 (for 2-component links)
├── reidemeister_i()   → Vec<crossing_id>
├── reidemeister_ii()  → Vec<(id, id)>
├── alexander_poly()   → Vec<f64> (coefficients)
└── is_unknot()        → bool (heuristic)
```

## Famous Knots and Their Invariants

| Knot | Crossings | Writhe | Alexander Δ(t) |
|------|-----------|--------|----------------|
| Unknot (0₁) | 0 | 0 | 1 |
| Trefoil (3₁) | 3 | ±3 | t² - t + 1 |
| Figure-eight (4₁) | 4 | 0 | -t² + 3t - 1 |
| Cinquefoil (5₁) | 5 | ±5 | t² - t + 1 (same as trefoil!) |
| Three-twist (5₂) | 5 | ±5 | 2t² - 3t + 2 |

Note: the Alexander polynomial cannot distinguish 3₁ from 5₁ — you need the Jones polynomial for that.

## Installation

```toml
[dependencies]
knot-theory = "0.1.0"
```

## Part of the SuperInstance Math Fleet

This crate is part of a larger mathematical toolkit:
- `chern-classes` — Characteristic classes of vector bundles
- `morse-theory` — Critical points and handle attachment
- `cohomology-ring` — Cup products and cohomology operations
- `graph-homology` — Homology of clique complexes
- `sheaf-laplacian` — Sheaf-theoretic diffusion

## References

- **Alexander, J.W.** (1928). *Topological invariants of knots and links.* Trans. AMS, 30(2), 275–306.
- **Burau, W.** (1936). *Über Zopfgruppen und gleichsinnig verdrillte Verkettungen.* Abh. Math. Sem. Hamburg, 11, 179–186.
- **Kauffman, L.H.** (1987). *On Knots.* Princeton University Press.
- **Rolfsen, D.** (1976). *Knots and Links.* AMS Chelsea Publishing.

## License

MIT © [SuperInstance](https://github.com/SuperInstance)
