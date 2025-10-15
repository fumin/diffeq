package diffeq_test

import (
	"fmt"
	"log"

	"github.com/fumin/diffeq"
)

func Example() {
	// Differential equation:
	// dy/dx = z + x
	// dz/dx = y
	dydx := func(dydx []float64, x float64, y []float64) {
		dydx[0] = y[1] + x
		dydx[1] = y[0]
	}
	// Solve for x between [0, 4].
	xspan := [2]float64{0, 4}
	// Initial values:
	// y(0) = 1
	// z(0) = -1
	y0 := []float64{1, -1}

	// Solve equation with the Dormand-Prince method.
	xs, ys, err := diffeq.DormandPrince(dydx, xspan, y0)
	if err != nil {
		log.Fatalf("%+v", err)
	}

	// Print results.
	fmt.Printf("x, y(x), z(x)\n")
	for i := range xs {
		x := xs[i]
		y := ys[i]
		fmt.Printf("%.3f, %.3f, %.3f\n", x, y[0], y[1])
	}

	// Output:
	// x, y(x), z(x)
	// 0.000, 1.000, -1.000
	// 0.091, 0.917, -0.913
	// 0.943, 0.868, -0.243
	// 1.752, 2.144, 0.872
	// 2.633, 6.065, 4.216
	// 3.742, 20.124, 17.311
	// 4.000, 26.328, 23.273
}
