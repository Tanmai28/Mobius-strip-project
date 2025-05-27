import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class MobiusStrip:
    """Class to model a Möbius strip and compute its geometric properties."""
    
    def __init__(self, R, w, n):
        """
        Initialize Möbius strip with radius R, width w, and resolution n.
        
        Parameters:
        - R: Distance from center to strip (float)
        - w: Width of the strip (float)
        - n: Number of points in each dimension of the mesh (int)
        """
        self.R = R
        self.w = w
        self.n = n
        self.u = np.linspace(0, 2 * np.pi, n)
        self.v = np.linspace(-w/2, w/2, n)
        self.U, self.V = np.meshgrid(self.u, self.v)
        self.mesh = None
        self._generate_mesh()

    def _generate_mesh(self):
        """Generate 3D mesh points using parametric equations."""
        u, v = self.U, self.V
        x = (self.R + v * np.cos(u / 2)) * np.cos(u)
        y = (self.R + v * np.cos(u / 2)) * np.sin(u)
        z = v * np.sin(u / 2)
        self.mesh = (x, y, z)

    def get_mesh(self):
        """Return the 3D mesh coordinates (x, y, z)."""
        return self.mesh

    def compute_surface_area(self):
        """
        Compute surface area using numerical integration.
        Uses the cross product of partial derivatives to approximate area.
        """
        u, v = self.U, self.V
        # Partial derivatives
        du = 2 * np.pi / (self.n - 1)
        dv = self.w / (self.n - 1)

        # Parametric equations
        x = (self.R + v * np.cos(u / 2)) * np.cos(u)
        y = (self.R + v * np.cos(u / 2)) * np.sin(u)
        z = v * np.sin(u / 2)

        # Partial derivatives with respect to u
        x_u = (-(self.R + v * np.cos(u / 2)) * np.sin(u) - 
               (v / 2) * np.sin(u / 2) * np.cos(u))
        y_u = ((self.R + v * np.cos(u / 2)) * np.cos(u) - 
               (v / 2) * np.sin(u / 2) * np.sin(u))
        z_u = (v / 2) * np.cos(u / 2)

        # Partial derivatives with respect to v
        x_v = np.cos(u / 2) * np.cos(u)
        y_v = np.cos(u / 2) * np.sin(u)
        z_v = np.sin(u / 2)

        # Cross product magnitude for surface element
        integrand = np.sqrt((y_u * z_v - z_u * y_v)**2 + 
                           (z_u * x_v - x_u * z_v)**2 + 
                           (x_u * y_v - y_u * x_v)**2)

        # Numerical integration using trapezoidal rule
        area = np.trapz(np.trapz(integrand, dx=dv, axis=0), dx=du)
        return area

    def compute_edge_length(self):
        """
        Compute the length of the boundary edge (v = ±w/2).
        Uses arc length formula along the boundary curves.
        """
        u = self.u
        du = 2 * np.pi / (self.n - 1)
        
        # Parametric equations for edges at v = w/2 and v = -w/2
        v_edges = [self.w / 2, -self.w / 2]
        total_length = 0

        for v in v_edges:
            x = (self.R + v * np.cos(u / 2)) * np.cos(u)
            y = (self.R + v * np.cos(u / 2)) * np.sin(u)
            z = v * np.sin(u / 2)

            # Derivatives with respect to u
            dx_du = (-(self.R + v * np.cos(u / 2)) * np.sin(u) - 
                     (v / 2) * np.sin(u / 2) * np.cos(u))
            dy_du = ((self.R + v * np.cos(u / 2)) * np.cos(u) - 
                     (v / 2) * np.sin(u / 2) * np.sin(u))
            dz_du = (v / 2) * np.cos(u / 2)

            # Arc length integrand
            integrand = np.sqrt(dx_du**2 + dy_du**2 + dz_du**2)
            length = np.trapz(integrand, dx=du)
            total_length += length

        return total_length

    def plot(self):
        """Visualize the Mobius strip in 3D."""
        x, y, z = self.mesh
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(x, y, z, cmap='viridis', alpha=0.8)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('Mobius Strip')
        plt.savefig('mobius_strip.png')
        plt.close()

# Example usage
if __name__ == "__main__":
    # Parameters
    R = 3.0  # Radius
    w = 1.0  # Width
    n = 100  # Resolution

    # Create Mobius strip
    mobius = MobiusStrip(R, w, n)
    
    # Compute properties
    area = mobius.compute_surface_area()
    edge_length = mobius.compute_edge_length()
    
    print(f"Surface Area: {area:.4f} square units")
    print(f"Edge Length: {edge_length:.4f} units")
    
    # Generate visualization
    mobius.plot()
