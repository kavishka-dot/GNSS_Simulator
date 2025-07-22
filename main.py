import numpy as np
from transmitter import Transmitter
from medium import Medium
from receiver import Receiver
from navigation import NavigationEngine

# -----------------------------
# Simulation Parameters
# -----------------------------
fs = 4.092e6  # Sampling frequency (Hz)
duration_ms = 1  # Duration in milliseconds

# -----------------------------
# Create Transmitters (Satellites)
# -----------------------------
satellites = []
base_pos = np.array([15600e3, 7540e3, 20140e3])  # ECEF center ref

for i in range(4):  # 4 satellites
    pos_offset = np.random.uniform(-500e3, 500e3, size=3)
    velocity = np.random.uniform(-1000, 1000, size=3)
    sat = Transmitter(prn=i+1, initial_position_ecef=base_pos + pos_offset, velocity_ecef=velocity)
    satellites.append(sat)

# -----------------------------
# Define Receiver True Position
# -----------------------------
receiver_true_pos = np.array([1917034.0, 6029783.0, -817119.0])  # ECEF in meters

# -----------------------------
# Medium: Combine signals and generate pseudoranges
# -----------------------------
channel = Medium(fs=fs, duration_ms=duration_ms)
composite_signal, measurements = channel.propagate_signals(satellites, receiver_true_pos)

# -----------------------------
# Receiver: Extract pseudoranges from measurements
# -----------------------------
# (In a real system, you'd track code phase — here we use channel's simulated pseudoranges)
pseudoranges = []
sat_positions = []

for m in measurements:
    prn = m["prn"]
    rho = m["pseudorange"]
    pseudoranges.append(rho)
    
    for sat in satellites:
        if sat.prn == prn:
            sat_positions.append(sat.get_position())
            break

# -----------------------------
# Navigation: Solve for Position and Clock Bias
# -----------------------------
nav = NavigationEngine()
est_pos, clock_bias = nav.solve_position(pseudoranges, sat_positions)

lat, lon, alt = nav.ecef_to_lla(*est_pos)

# -----------------------------
# Results
# -----------------------------
true_lat, true_lon, true_alt = nav.ecef_to_lla(*receiver_true_pos)
error = np.linalg.norm(est_pos - receiver_true_pos)

print("\n--- GNSS End-to-End Simulation ---")
print(f"True Position (ECEF):     {receiver_true_pos}")
print(f"Estimated Position (ECEF): {est_pos}")
print(f"Positioning Error:        {error:.3f} m")
print(f"Clock Bias:               {clock_bias * 1e9:.2f} ns")
print(f"\nTrue Position (LLA):      ({true_lat:.6f}, {true_lon:.6f}, {true_alt:.2f} m)")
print(f"Estimated Position (LLA): ({lat:.6f}, {lon:.6f}, {alt:.2f} m)")
