# GNSS\_Simulator

A fully modular, end-to-end simulation of a civilian Global Navigation Satellite System (GNSS) receiver pipeline written in Python. This project implements the essential components of a software GNSS receiver, including satellite signal generation, channel effects, receiver processing, and position computation.

---

## 🚀 Overview

This simulator models the following GNSS receiver pipeline:

1. **Satellite Transmitters**: Simulated PRN-coded signals (C/A code)
2. **Signal Propagation Medium**: Models delays, Doppler, noise, and real-world errors
3. **Receiver**: Tracks signals, extracts pseudoranges
4. **Navigation Engine**: Solves for receiver position and clock bias
5. **Main Script**: Combines all components into a test scenario

Designed to replicate realistic GNSS reception and solve the Position, Velocity, Time (PVT) problem up to the Position level (Level 6).

---

## 📐 GNSS Positioning Fundamentals

### Pseudorange Equation

The pseudorange measured by the receiver from satellite $i$ is:

$$
\rho_i = R_i + c \delta t - c \delta t_i^{\text{sat}} + I_i + T_i + M_i + \epsilon_i
$$

Where:

* $\rho_i$: measured pseudorange
* $R_i$: geometric range between receiver and satellite
* $c \delta t$: receiver clock bias
* $c \delta t_i^{\text{sat}}$: satellite clock bias
* $I_i$, $T_i$: ionospheric and tropospheric delays
* $M_i$: multipath delay
* $\epsilon_i$: measurement noise

To solve for the receiver position and clock bias, the following nonlinear system is constructed:

$$
\rho_i = \sqrt{(x - x_i)^2 + (y - y_i)^2 + (z - z_i)^2} + c \delta t
$$

Where $(x_i, y_i, z_i)$ is the satellite's ECEF position and $(x, y, z)$ is the unknown receiver position.

---

## 🧱 Components

### 1. Transmitter (`transmitter.py`)

* Simulates PRN-coded signal generation
* Computes code phase and Doppler shift from satellite-to-receiver geometry
* Emits signal sampled at user-defined frequency

### 2. Medium (`medium.py`)

* Propagates signals from all transmitters to receiver
* Adds:

  * Satellite clock bias $c\delta t^{\text{sat}}$
  * Ionospheric delay $I_i$
  * Tropospheric delay $T_i$
  * Multipath delay $M_i$
  * AWGN $\epsilon_i$
* Returns composite signal and component-wise pseudoranges

### 3. Receiver (`receiver.py`)

* Tracks the satellite signals (simulation)
* Extracts pseudoranges from signal or directly from medium

### 4. Navigation Engine (`navigation.py`)

* Solves for $(x, y, z, \delta t)$ using Gauss-Newton iteration
* Converts ECEF position to LLA (Latitude, Longitude, Altitude)

---

## ⚙️ Simulation Workflow (`main.py`)

```python
# Pseudocode:
1. Define satellites with positions and velocities
2. Define true receiver position (ECEF)
3. Simulate signal propagation with Medium
4. Receiver extracts pseudoranges
5. NavigationEngine solves for position
6. Compare estimated vs true position
```

---

## 📌 Gauss-Newton PVT Solver

The linearized system:

$$
H \cdot \Delta x = \Delta \rho
$$

Where:

* $H$: Jacobian matrix (partial derivatives wrt position and clock bias)
* $\Delta \rho$: residuals between measured and predicted pseudoranges
* $\Delta x$: update to $[x, y, z, c\delta t]$

Updated iteratively until $\|\Delta x\| < \epsilon$.

---

## 📦 File Structure

```
gnss_simulator/
├── main.py
├── transmitter.py
├── medium.py
├── receiver.py
├── navigation.py
└── README.md
```

---

## 🛰️ Future Extensions

* Ephemeris-based orbital motion
* Real NAV message decoding
* Dual-frequency ionosphere mitigation
* Doppler-based velocity estimation
* Kalman filtering and RTK
* Integration with real RF front ends (SDR)
