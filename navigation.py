import numpy as np

class NavigationEngine:
    def __init__(self, max_iter=10, tol=1e-4):
        self.c = 299_792_458  # Speed of light
        self.max_iter = max_iter
        self.tol = tol

    def solve_position(self, pseudoranges, satellite_positions):
        """
        Solves for receiver position and clock bias using Gauss-Newton.
        Args:
            pseudoranges: list of measured pseudoranges [ρ₁, ρ₂, ..., ρn]
            satellite_positions: list of ECEF satellite positions [(x₁,y₁,z₁), ...]
        Returns:
            receiver_pos (x,y,z), clock_bias (in seconds)
        """
        num_sats = len(pseudoranges)
        assert num_sats >= 4, "Need at least 4 satellites for PVT solution"

        # Initial guess: center of Earth + 0 clock bias
        x = np.array([0.0, 0.0, 6371e3, 0.0])  # x, y, z, clock_bias

        for iteration in range(self.max_iter):
            H = np.zeros((num_sats, 4))  # Jacobian
            delta_rho = np.zeros(num_sats)

            for i in range(num_sats):
                sat_pos = satellite_positions[i]
                rho_i = pseudoranges[i]

                dx = x[0] - sat_pos[0]
                dy = x[1] - sat_pos[1]
                dz = x[2] - sat_pos[2]
                r = np.sqrt(dx**2 + dy**2 + dz**2)

                # Predicted pseudorange
                rho_pred = r + self.c * x[3]

                # Residual
                delta_rho[i] = rho_i - rho_pred

                # Jacobian row
                H[i, 0] = -dx / r
                H[i, 1] = -dy / r
                H[i, 2] = -dz / r
                H[i, 3] = -self.c

            # Gauss-Newton update
            dx = np.linalg.lstsq(H, delta_rho, rcond=None)[0]
            x += dx

            if np.linalg.norm(dx) < self.tol:
                break

        receiver_pos = x[:3]
        clock_bias = x[3]
        return receiver_pos, clock_bias

    def ecef_to_lla(self, x, y, z):
        """
        Convert ECEF (x, y, z) to Latitude, Longitude, Altitude (WGS-84)
        """
        # WGS84 constants
        a = 6378137.0  # Semi-major axis
        e = 8.1819190842622e-2  # First eccentricity

        b = np.sqrt(a**2 * (1 - e**2))
        ep = np.sqrt((a**2 - b**2) / b**2)
        p = np.sqrt(x**2 + y**2)
        th = np.arctan2(a * z, b * p)

        lon = np.arctan2(y, x)
        lat = np.arctan2((z + ep**2 * b * np.sin(th)**3),
                         (p - e**2 * a * np.cos(th)**3))
        N = a / np.sqrt(1 - e**2 * np.sin(lat)**2)
        alt = p / np.cos(lat) - N

        # Convert to degrees
        lat = np.degrees(lat)
        lon = np.degrees(lon)
        return lat, lon, alt
