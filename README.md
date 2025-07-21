# Simple Ray Tracer

A basic ray tracer implementation based on "Ray Tracing in One Weekend" tutorial with custom scene composition and camera configuration.

## Features

- Ray tracing with multiple material types:
  - Lambertian (diffuse) materials
  - Metal materials with adjustable fuzziness
  - Dielectric (glass) materials with refraction
- Customizable camera position and orientation
- Anti-aliasing through multiple samples per pixel
- Gamma correction for proper color output

## Building and Running

```bash
# Compile the ray tracer
g++ -std=c++17 -O3 -o raytracer main.cpp

# Run the ray tracer (outputs image.ppm)
./raytracer
```

## Scene Description

The custom scene includes:
- Multiple spheres with different materials and sizes
- Varied positioning for visual appeal
- Custom camera angle for an interesting perspective

## Output

The ray tracer generates a PPM image file that can be viewed with most image viewers or converted to other formats using tools like ImageMagick:

```bash
# Convert PPM to PNG
convert image.ppm image.png
```