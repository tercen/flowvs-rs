//! R's one-dimensional optimiser (`optimize`), which is Brent's method as implemented in
//! `Brent_fmin` (R sources, `src/library/stats/src/optimize.c`).
//!
//! It matters that this is R's exact variant rather than any Brent: the peak finder and the
//! cofactor search both stop on its tolerance rule, so a different stopping rule moves answers
//! in the last digits the plan's tolerance is measured in.

/// R's `.Machine$double.eps^0.25`, the default `tol` of `optimize`.
pub fn default_tol() -> f64 {
    f64::EPSILON.powf(0.25)
}

/// Minimise `f` on `[ax, bx]`. Returns `(argmin, min)`.
pub fn brent_fmin<F: FnMut(f64) -> f64>(ax: f64, bx: f64, mut f: F, tol: f64) -> (f64, f64) {
    let c = (3.0 - 5.0f64.sqrt()) * 0.5;
    let mut eps = f64::EPSILON;
    eps = eps.sqrt();

    let (mut a, mut b) = (ax, bx);
    let mut v = a + c * (b - a);
    let (mut w, mut x) = (v, v);
    let (mut d, mut e) = (0.0f64, 0.0f64);

    let mut fx = f(x);
    let (mut fv, mut fw) = (fx, fx);
    let tol3 = tol / 3.0;

    loop {
        let xm = (a + b) * 0.5;
        let tol1 = eps * x.abs() + tol3;
        let t2 = tol1 * 2.0;
        if (x - xm).abs() <= t2 - (b - a) * 0.5 {
            break;
        }
        let (mut p, mut q, mut r) = (0.0f64, 0.0f64, 0.0f64);
        if e.abs() > tol1 {
            // fit a parabola through the three best points so far
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = (q - r) * 2.0;
            if q > 0.0 {
                p = -p;
            } else {
                q = -q;
            }
            r = e;
            e = d;
        }
        // Declared before the branches, as in R's `Brent_fmin`: the golden-section and parabolic
        // arms both fall through to the same clamping step below, and collapsing them into an
        // expression would break the line-for-line correspondence this port is checked against.
        #[allow(clippy::needless_late_init)]
        let u;
        if p.abs() >= (q * 0.5 * r).abs() || p <= q * (a - x) || p >= q * (b - x) {
            // golden section
            e = if x < xm { b - x } else { a - x };
            d = c * e;
        } else {
            // parabolic interpolation
            d = p / q;
            let cand = x + d;
            if cand - a < t2 || b - cand < t2 {
                d = tol1;
                if x >= xm {
                    d = -d;
                }
            }
        }
        if d.abs() >= tol1 {
            u = x + d;
        } else if d > 0.0 {
            u = x + tol1;
        } else {
            u = x - tol1;
        }
        let fu = f(u);
        if fu <= fx {
            if u < x {
                b = x;
            } else {
                a = x;
            }
            v = w;
            w = x;
            x = u;
            fv = fw;
            fw = fx;
            fx = fu;
        } else {
            if u < x {
                a = u;
            } else {
                b = u;
            }
            if fu <= fw || w == x {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            } else if fu <= fv || v == x || v == w {
                v = u;
                fv = fu;
            }
        }
    }
    (x, fx)
}

/// Maximise `f` on `[ax, bx]` — R's `optimize(..., maximum = TRUE)`, which minimises `-f`.
pub fn brent_fmax<F: FnMut(f64) -> f64>(ax: f64, bx: f64, mut f: F, tol: f64) -> (f64, f64) {
    let (x, neg) = brent_fmin(ax, bx, |t| -f(t), tol);
    (x, -neg)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_the_minimum_of_a_parabola() {
        let (x, y) = brent_fmin(-5.0, 5.0, |t| (t - 1.25) * (t - 1.25) + 3.0, default_tol());
        assert!((x - 1.25).abs() < 1e-6, "x = {x}");
        assert!((y - 3.0).abs() < 1e-10);
    }

    #[test]
    fn maximum_matches_r_on_a_known_case() {
        // R: optimize(function(x) sin(x), c(0, pi), maximum = TRUE)$maximum -> 1.570799
        let (x, y) = brent_fmax(0.0, std::f64::consts::PI, |t| t.sin(), default_tol());
        assert!((x - 1.5707990).abs() < 1e-5, "x = {x}");
        assert!((y - 1.0).abs() < 1e-8);
    }
}
