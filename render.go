package main

import (
	"image"
	"image/color"
	"image/png"
	"os"
	"os/exec"
)

// render renders current state to a png
func (b *body) render(filename string, w, h int) error {
	pixpm := 1.0

	img := image.NewRGBA(image.Rect(0, 0, w, h))

	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			img.SetRGBA(x, y, color.RGBA{200, 200, 200, 255})
			if (x < w/2-b.w/4 || x > w/2+b.w/4) && y > (h+b.h)/2 {
				img.SetRGBA(x, y, color.RGBA{100, 100, 100, 255})
			}
		}
	}

	for y := 1; y < b.h-1; y++ {
		for x := 1; x < b.w-1; x++ {
			if x == b.w-2 && y == b.h-2 {
				continue
			}
			nx := int(0.5+b.u[b.w*y+x][0]*pixpm) + (w-b.w)/2
			ny := h - (int(0.5+b.u[b.w*y+x][1]*pixpm) + (h-b.h)/2)

			if nx < 0 || nx >= w || ny < 0 || ny >= h {
				continue
			}
			for i := -1; i <= 1; i++ {
				for j := -1; j <= 1; j++ {
					img.SetRGBA(nx+i, ny+j, color.RGBA{100, 120, 150, 255})
				}
			}
		}
	}

	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	err = png.Encode(file, img)
	if err != nil {
		return err
	}

	return nil
}

// encodeVideo calls ffmpeg to encode the frames in fmask
func encodeVideo(outname, fmask string) error {
	cmd := exec.Command("ffmpeg", "-v", "error", "-y", "-f", "image2", "-i", fmask, outname)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	err := cmd.Run()
	return err
}
