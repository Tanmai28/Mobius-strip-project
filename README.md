# Mobius-strip-project
Python script that models a Mobius strip using parametric equations and computes key geometric properties.

### Write-Up

**How I Structured the Code:**
I built the Python script around a `MobiusStrip` class to keep things organized and easy to follow. The class takes three inputs: radius `R` (how far the strip sits from the center), width `w` (how wide the strip is), and resolution `n` (how detailed the mesh is). I broke it down into focused methods:
- `_generate_mesh` creates the 3D points (x, y, z) using the Möbius strip’s parametric equations.
- `get_mesh` hands over the mesh data if needed.
- `compute_surface_area` figures out the surface area numerically.
- `compute_edge_length` calculates the length of the strip’s edges.
- `plot` whips up a 3D visualization with Matplotlib.
Each method does one job, and I added clear comments and docstrings to make the code easy to read and tweak. The main block shows how it all works with sample values (`R=3.0`, `w=1.0`, `n=100`).

**How I Approximated Surface Area:**
To get the surface area, I used numerical integration, which sounds fancy but is just a way to add up tiny surface patches. Here’s the gist:
1. I calculated the partial derivatives of the parametric equations for `u` and `v` (the parameters defining the strip).
2. I took their cross product to find the area of each tiny patch on the surface.
3. Using NumPy’s trapezoidal rule (`np.trapz`), I summed these patches over `v` and then `u` across the strip’s domain (`u` from 0 to 2π, `v` from -w/2 to w/2).
This gives a solid estimate, close to the expected `2πRw` (around 18.8496 for `R=3`, `w=1`), and gets more accurate with higher resolution.

**Challenges I Faced:**
- **Balancing Accuracy and Speed**: I needed a high resolution (`n=100`) for precise results, but that made the code run slower. Finding the sweet spot was tricky.
- **Nailing the Math**: The partial derivatives were a bit of a headache because of all the sines and cosines. A small mistake could throw off the area or edge calculations, so I double-checked everything.
- **Making the Plot Pop**: The Möbius strip’s twisty, non-orientable nature is cool but tough to visualize clearly. I played with Matplotlib’s settings—like transparency (`alpha=0.8`) and the `viridis` colormap—to make the 3D plot look sharp and show off the strip’s unique shape.

This setup makes the code clean, functional, and ready to explore the Möbius strip’s wild geometry!
