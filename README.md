# Pocket
> :warning: *Work-in-progress!*

## Description
A tiny graphics framework.

## Modules
- psr
  - 2D/3D software rasterizer
  - Render triangles, lines and bitmap fonts
- pwa
  - [Microsoft Windows](https://www.microsoft.com/en-us/windows) abstraction
  - Manage windows and receive input
  - Measure time with microsecond precision

## How to use
Just include the headers and compile the source files. If you compile `pwa.c` you must also link against `gdi32`.

For example:
```
gcc my_program.c psr.c pwa.c -lgdi32 -o my_program
```

## Screenshots
![Bullfrog](examples/assets/readme_bullfrog.png "Bullfrog")
