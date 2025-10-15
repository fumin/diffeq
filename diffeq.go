// Package diffeq provides tools for solving differential equations.
package diffeq

import (
	"math"

	"github.com/pkg/errors"
)

var (
	ErrTooSmallStep = errors.Errorf("too small step")
)

// DydxFunc returns the derivative dy/dx given x and y.
type DydxFunc func(dydx []float64, x float64, y []float64)

// DormandPrince performs the [Dormand-Prince] method.
//
// [Dormand-Prince]: https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method
func DormandPrince(dydxFunc DydxFunc, xspan [2]float64, y0 []float64) (xs []float64, ys [][]float64, err error) {
	rk := &dormandPrince{}
	return rungeKuttaIntegrate(rk, dydxFunc, xspan, y0)
}

type rungeKutta interface {
	errorOrder() int
	yPlusH(yPlusH, te []float64, dydxFunc DydxFunc, x float64, y []float64, h float64)
}

func rungeKuttaIntegrate(rk rungeKutta, dydxFunc DydxFunc, xspan [2]float64, y0 []float64) (xs []float64, ys [][]float64, err error) {
	xs = append(xs, xspan[0])
	ys = append(ys, y0)

	tol := tolerance{abs: 1e-6, rel: 1e-3}
	h := getFirstStep(rk, dydxFunc, xspan, y0, tol)

	var x float64 = xspan[0]
	y := make([]float64, len(y0))
	copy(y, y0)
	te := make([]float64, len(y0))
	for x < xspan[1] {
		yPlusH := make([]float64, len(y0))
		xPlusH, newH, err := rungeKuttaStep(yPlusH, te, rk, dydxFunc, x, y, h, xspan[1], tol)
		if err != nil {
			return nil, nil, errors.Wrap(err, "")
		}

		// Add to result.
		xs = append(xs, xPlusH)
		ys = append(ys, yPlusH)

		// Update iteration state.
		x = xPlusH
		copy(y, yPlusH)
		h = newH
	}

	return
}

func rungeKuttaStep(yPlusH, te []float64, rungeKutta rungeKutta, dydxFunc DydxFunc, x float64, y []float64, h, xMax float64, tol tolerance) (xPlusH float64, newH float64, err error) {
	minStep := 10 * math.Abs(math.Nextafter(x, math.Inf(1))-x)
	rejected := false
	for {
		if h < minStep {
			return math.NaN(), math.NaN(), ErrTooSmallStep
		}

		xPlusH = x + h
		if xPlusH > xMax {
			xPlusH = xMax
			h = xPlusH - x
		}

		// Compute yPlusH and te.
		rungeKutta.yPlusH(yPlusH, te, dydxFunc, x, y, h)

		// Compute errorNorm.
		scale := make([]float64, len(y))
		for i := range scale {
			scale[i] = tol.abs + tol.rel*max(math.Abs(y[i]), math.Abs(yPlusH[i]))
		}
		var errorNorm float64
		for i := range te {
			errorNorm += math.Pow(te[i]*h/scale[i], 2)
		}
		errorNorm = math.Sqrt(errorNorm / float64(len(te)))

		// Compute new h.
		const maxFactor = 10
		const minFactor = 0.2
		const safety = 0.9
		exponent := -1 / (float64(rungeKutta.errorOrder()) + 1)
		var factor float64
		if errorNorm < 1 {
			if errorNorm == 0 {
				factor = maxFactor
			} else {
				factor = min(safety*math.Pow(errorNorm, exponent), maxFactor)
			}
			if rejected {
				factor = min(1, factor)
			}

			newH = h * factor
			return
		} else {
			h *= max(safety*math.Pow(errorNorm, exponent), minFactor)
			rejected = true
		}
	}
}

// E. Hairer, S. P. Norsett G. Wanner, "Solving Ordinary Differential Equations I: Nonstiff Problems", Sec. II.4.
func getFirstStep(rungeKutta rungeKutta, dydxFunc func(dydx []float64, x float64, y []float64), xspan [2]float64, y0 []float64, tol tolerance) float64 {
	dydx0 := make([]float64, len(y0))
	dydxFunc(dydx0, xspan[0], y0)

	// Compute d0, d1.
	scale := make([]float64, len(y0))
	for i, y := range y0 {
		scale[i] = tol.abs + math.Abs(y)*tol.rel
	}
	var d0, d1 float64
	for i := range scale {
		d0 += math.Pow(y0[i]/scale[i], 2)
		d1 += math.Pow(dydx0[i]/scale[i], 2)
	}
	d0 = math.Sqrt(d0 / float64(len(scale)))
	d1 = math.Sqrt(d1 / float64(len(scale)))

	// Compute h0.
	var h0 float64
	if d0 < 1e-5 || d1 < 1e-6 {
		h0 = 1e-6
	} else {
		h0 = 0.01 * d0 / d1
	}
	h0 = min(h0, xspan[1]-xspan[0])

	// Compute d2.
	y1 := make([]float64, len(y0))
	for i := range y0 {
		y1[i] = y0[i] + h0*dydx0[i]
	}
	dydx1 := make([]float64, len(y0))
	dydxFunc(dydx1, xspan[0]+h0, y1)
	var d2 float64
	for i := range dydx0 {
		d2 += math.Pow((dydx1[i]-dydx0[i])/scale[i], 2)
	}
	d2 = math.Sqrt(d2/float64(len(scale))) / h0

	// Compute h1.
	var h1 float64
	if d1 <= 1e-15 && d2 <= 1e-15 {
		h1 = max(1e-6, h0*1e-3)
	} else {
		order := float64(rungeKutta.errorOrder())
		h1 = math.Pow(0.01/max(d1, d2), 1./(order+1))
	}

	return min(100*h0, h1, xspan[1]-xspan[0])
}

type tolerance struct {
	abs float64
	rel float64
}

type dormandPrince struct{}

func (dp *dormandPrince) errorOrder() int {
	return 4
}

func (dp *dormandPrince) yPlusH(yPlusH, te []float64, dydxFunc DydxFunc, x float64, y []float64, h float64) {
	c := []float64{0, 1. / 5, 3. / 10, 4. / 5, 8. / 9, 1, 1}
	a := [][]float64{
		{0, 0, 0, 0, 0, 0},
		{1. / 5, 0, 0, 0, 0, 0},
		{3. / 40, 9. / 40, 0, 0, 0, 0},
		{44. / 45, -56. / 15, 32. / 9, 0, 0, 0},
		{19372. / 6561, -25360. / 2187, 64448. / 6561, -212. / 729, 0, 0},
		{9017. / 3168, -355. / 33, 46732. / 5247, 49. / 176, -5103. / 18656, 0},
		{35. / 384, 0, 500. / 1113, 125. / 192, -2187. / 6784, 11. / 84},
	}
	b := []float64{35. / 384, 0, 500. / 1113, 125. / 192, -2187. / 6784, 11. / 84, 0}
	e := []float64{-71. / 57600, 0, 71. / 16695, -71. / 1920, 17253. / 339200, -22. / 525, 1. / 40}

	k := make([][]float64, len(c))
	for i := range k {
		k[i] = make([]float64, len(y))
	}

	yItp := make([]float64, len(y))
	for i := range k {
		xItp := x + c[i]*h
		for j := range len(yItp) {
			yItp[j] = y[j]
			for l := range len(k) - 1 {
				yItp[j] += a[i][l] * k[l][j] * h
			}
		}

		dydxFunc(k[i], xItp, yItp)
	}

	for i := range yPlusH {
		yPlusH[i] = y[i]
		for j := range b {
			yPlusH[i] += b[j] * h * k[j][i]
		}
	}
	for i := range te {
		te[i] = 0
		for j := range e {
			te[i] += e[j] * k[j][i]
		}
	}
}
