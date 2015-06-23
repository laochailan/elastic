// This is a simulation of the following elastic body equation:
//
// Let D ⊂ ℝ³ be a domain describing an unstrained
// elastic body and u(x,t) be a function so that u(G,t)
// describes the body shape at time t.
//
// 	ρ ∂²_t u = ∂_i d ∂_i u (1-1/|∂_i u|) + f
//
// where ρ is the body density (assumed to be constant),
// d is the stress tensor and f is an external force density.
package main

import (
	"fmt"
	"math"
)

type body struct {
	w, h   int
	dx, dy float64

	damp float64

	d    []float64
	ρ    []float64
	u    [][2]float64
	dudt [][2]float64
	dudx [][2]float64
	dudy [][2]float64
}

func newBody(w, h int, dx, dy float64) *body {
	b := &body{w: w, h: h, dx: dx, dy: dy}

	b.d = make([]float64, w*h)
	b.ρ = make([]float64, w*h)
	b.u = make([][2]float64, w*h)
	b.dudt = make([][2]float64, w*h)
	b.dudx = make([][2]float64, w*(h-1))
	b.dudy = make([][2]float64, w*(h-1))

	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			b.u[y*w+x][0] = float64(x) * dx //- math.Cos(math.Pi/float64(w)*float64(x)*dx)*2
			b.u[y*w+x][1] = float64(y) * dy

			b.d[y*w+x] = 1
			b.ρ[y*w+x] = 1
		}
	}

	for y := 0; y < h-1; y++ {
		for x := 0; x < w; x++ {
			b.dudx[y*w+x][0] = 1
			b.dudx[y*w+x][1] = 0
			b.dudy[y*w+x][0] = 0
			b.dudy[y*w+x][1] = 1

		}
	}
	return b
}

func (b *body) step(dt float64) {
	w := b.w
	h := b.h
	for y := 1; y < h-2; y++ {
		for x := 1; x < w-2; x++ {
			b.dudx[y*w+x][0] = (b.u[y*w+x+1][0] - b.u[y*w+x][0]) / b.dx
			b.dudx[y*w+x][1] = (b.u[y*w+x+1][1] - b.u[y*w+x][1]) / b.dx
			b.dudy[y*w+x][0] = (b.u[(y+1)*w+x][0] - b.u[y*w+x][0]) / b.dy
			b.dudy[y*w+x][1] = (b.u[(y+1)*w+x][1] - b.u[y*w+x][1]) / b.dy
		}
	}

	for y := 1; y < h-1; y++ {
		for x := 1; x < w-1; x++ {
			dudx0 := b.dudx[y*w+x]
			dudx1 := b.dudx[y*w+x-1]

			dudy0 := b.dudy[y*w+x]
			dudy1 := b.dudy[(y-1)*w+x]

			fx := (b.d[y*w+x-1]*dudx0[0]*(1-1/math.Sqrt(dudx0[0]*dudx0[0]+dudx0[1]*dudx0[1])) -
				b.d[y*w+x]*dudx1[0]*(1-1/math.Sqrt(dudx1[0]*dudx1[0]+dudx1[1]*dudx1[1]))) / b.dx
			fx += (b.d[(y-1)*w+x]*dudy0[0]*(1-1/math.Sqrt(dudy0[0]*dudy0[0]+dudy0[1]*dudy0[1])) -
				b.d[y*w+x]*dudy1[0]*(1-1/math.Sqrt(dudy1[0]*dudy1[0]+dudy1[1]*dudy1[1]))) / b.dy
			fy := (b.d[y*w+x-1]*dudx0[1]*(1-1/math.Sqrt(dudx0[0]*dudx0[0]+dudx0[1]*dudx0[1])) -
				b.d[y*w+x]*dudx1[1]*(1-1/math.Sqrt(dudx1[0]*dudx1[0]+dudx1[1]*dudx1[1]))) / b.dx
			fy += (b.d[(y-1)*w+x]*dudy0[1]*(1-1/math.Sqrt(dudy0[0]*dudy0[0]+dudy0[1]*dudy0[1])) -
				b.d[y*w+x]*dudy1[1]*(1-1/math.Sqrt(dudy1[0]*dudy1[0]+dudy1[1]*dudy1[1]))) / b.dy

			fx += (b.fFac(x, y)*b.dudy[y*w+x][1] - b.fFac(x-1, y)*b.dudy[y*w+x-1][1]) / b.dx
			fx += (-b.fFac(x, y)*b.dudx[y*w+x][1] + b.fFac(x, y-1)*b.dudx[(y-1)*w+x][1]) / b.dy
			fy += (-b.fFac(x, y)*b.dudy[y*w+x][0] + b.fFac(x-1, y)*b.dudy[y*w+x-1][0]) / b.dx
			fy += (b.fFac(x, y)*b.dudx[y*w+x][0] - b.fFac(x, y-1)*b.dudx[(y-1)*w+x][0]) / b.dy

			f := -0.0005

			b.dudt[y*w+x][0] += dt * (fx - b.damp*b.dudt[y*w+x][0]) / b.ρ[y*w+x]
			b.dudt[y*w+x][1] += dt * (fy + f - b.damp*b.dudt[y*w+x][1]) / b.ρ[y*w+x]
		}
	}
	for y := 1; y < h-1; y++ {
		for x := 1; x < w-1; x++ {
			b.u[y*w+x][0] += b.dudt[y*w+x][0] * dt
			b.u[y*w+x][1] += b.dudt[y*w+x][1] * dt
			//fmt.Println(b.u[y*w+x][1])

			if (x < w/4 || x > w*3/4) && b.u[y*w+x][1] < 0 {
				b.u[y*w+x][1] = 0
			}
		}
	}
}

func sign(x float64) float64 {
	if x > 0 {
		return 1
	}
	return -1
}

func (b *body) fFac(x, y int) float64 {
	if x <= 0 || y <= 0 || x >= b.w-2 || y >= b.h-2 {
		return 0
	}

	detDu := b.dudx[y*b.w+x][0]*b.dudy[y*b.w+x][1] - b.dudy[y*b.w+x][0]*b.dudx[y*b.w+x][1]

	return b.d[y*b.w+x] * (detDu - sign(detDu))
}

func main() {
	const output = "out.mkv"

	b := newBody(70, 10, 1, 1)
	b.damp = 0.001
	for i := 0; i < 300; i++ {
		for j := 0; j < 500; j++ {
			b.step(0.08)
		}

		err := b.render(fmt.Sprintf("frames/out%05d.png", i), 100, 100)
		if err != nil {
			fmt.Println("rendering output:", err)
			return
		}
	}

	err := encodeVideo(output, "frames/out%05d.png")
	if err != nil {
		fmt.Println("encoding video:", err)
		return
	}

}
