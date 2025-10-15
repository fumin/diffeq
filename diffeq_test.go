package diffeq

import (
	"flag"
	"fmt"
	"log"
	"math"
	"testing"
)

func TestDormandPrince(t *testing.T) {
	type testcase struct {
		dydx  func([]float64, float64, []float64)
		xspan [2]float64
		y0    []float64
		y     func(float64) []float64
		tol   float64
	}
	tests := []testcase{
		{
			dydx: func(dydx []float64, x float64, y []float64) {
				dydx[0] = -y[1]
				dydx[1] = y[0]
			},
			xspan: [2]float64{0, 4},
			y0:    []float64{1, 1},
			y: func(x float64) []float64 {
				y := make([]float64, 2)
				y[0] = math.Cos(x) - math.Sin(x)
				y[1] = math.Cos(x) + math.Sin(x)
				return y
			},
			tol: 2e-3,
		},
		{
			dydx: func(dydx []float64, x float64, y []float64) {
				dydx[0] = y[1] + x
				dydx[1] = y[0]
			},
			xspan: [2]float64{0, 4},
			y0:    []float64{1, -1},
			y: func(x float64) []float64 {
				y := make([]float64, 2)
				y[0] = 1./2*math.Exp(x) + 3./2*math.Exp(-x) - 1
				y[1] = 1./2*math.Exp(x) - 3./2*math.Exp(-x) - x
				return y
			},
			tol: 2e-3,
		},
	}
	for i, test := range tests {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			t.Parallel()
			xs, ys, err := DormandPrince(test.dydx, test.xspan, test.y0)
			if err != nil {
				t.Fatalf("%+v", err)
			}
			for i := range xs {
				x := xs[i]
				yi := ys[i]
				yTrue := test.y(x)

				for j := range yi {
					diff := math.Abs(yi[j] - yTrue[j])
					if diff > test.tol {
						t.Errorf("x %f y[%d] %f %f diff %f", x, j, yi[j], yTrue[j], diff)
					}
				}
			}
		})
	}
}

func TestMain(m *testing.M) {
	flag.Parse()
	log.SetFlags(log.Lmicroseconds | log.Llongfile | log.LstdFlags)

	m.Run()
}
