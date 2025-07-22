import numpy as np

class Medium:
    def __init__(self, fs, duration_ms):
        self.fs = fs
        self.duration_ms = duration_ms
        self.duration_s = duration_ms * 1e-3
        self.total_samples = int(fs * self.duration_s)
        self.c = 299_792_458  # Speed of light (m/s)

        # Adjustable error model parameters
        self.snr_db = 40                 # Signal-to-noise ratio in dB
        self.sat_clock_bias_std = 10e-9  # Satellite clock bias (sec)
        self.iono_delay_mean = 5.0       # Ionospheric delay (m)
        self.tropo_delay_mean = 2.5      # Tropospheric delay (m)
        self.multipath_std = 3.0         # Multipath (m)
        self.receiver_noise_std = 1.0    # Receiver measurement error (m)

    def add_awgn(self, signal, snr_db=None):
        """
        Add Additive White Gaussian Noise (AWGN) to the input signal.
        """
        if snr_db is None:
            snr_db = self.snr_db
        signal_power = np.mean(signal**2)
        snr_linear = 10**(snr_db / 10)
        noise_power = signal_power / snr_linear
        noise = np.random.normal(0, np.sqrt(noise_power), size=signal.shape)
        return signal + noise

    def simulate_satellite_clock_bias(self, prn):
        """ Simulate satellite clock bias in seconds """
        return np.random.normal(0, self.sat_clock_bias_std)

    def simulate_ionospheric_delay(self, prn, receiver_pos):
        """ Simulate ionospheric delay in meters """
        return np.random.normal(self.iono_delay_mean, 0.5)

    def simulate_tropospheric_delay(self, prn, receiver_pos):
        """ Simulate tropospheric delay in meters """
        return np.random.normal(self.tropo_delay_mean, 0.3)

    def simulate_multipath_delay(self):
        """ Simulate multipath delay in meters """
        return np.random.normal(0, self.multipath_std)

    def simulate_measurement_noise(self):
        """ Simulate measurement noise in meters """
        return np.random.normal(0, self.receiver_noise_std)

    def simulate_pseudorange_components(self, tx, receiver_position):
        """
        Compute detailed pseudorange terms for transmitter `tx`.
        Returns pseudorange and all components.
        """
        sat_pos = tx.get_position()
        R_i = np.linalg.norm(receiver_position - sat_pos)  # Geometric range

        # Simulate all error terms
        c_dt_sat = self.c * self.simulate_satellite_clock_bias(tx.prn)
        I_i = self.simulate_ionospheric_delay(tx.prn, receiver_position)
        T_i = self.simulate_tropospheric_delay(tx.prn, receiver_position)
        M_i = self.simulate_multipath_delay()
        eps_i = self.simulate_measurement_noise()

        # Final simulated pseudorange
        rho_i = R_i - c_dt_sat + I_i + T_i + M_i + eps_i

        return rho_i, {
            "R_i": R_i,
            "c_dt_sat": c_dt_sat,
            "I_i": I_i,
            "T_i": T_i,
            "M_i": M_i,
            "eps_i": eps_i
        }

    def propagate_signals(self, transmitters, receiver_position):
        """
        Main interface to:
        - Generate signal from each transmitter
        - Combine all signals into one received waveform
        - Simulate pseudorange components
        Returns:
            - Noisy baseband signal
            - List of (PRN, pseudorange, components)
        """
        composite_signal = np.zeros(self.total_samples)
        measurements = []

        for tx in transmitters:
            # Generate signal
            tx_signal, _, _ = tx.simulate_transmission(
                receiver_pos_ecef=receiver_position,
                fs=self.fs,
                duration_ms=self.duration_ms
            )
            composite_signal += tx_signal

            # Compute true + error-based pseudorange
            rho_i, components = self.simulate_pseudorange_components(tx, receiver_position)
            measurements.append({
                "prn": tx.prn,
                "pseudorange": rho_i,
                "components": components
            })

        # Add noise to the final combined signal
        noisy_signal = self.add_awgn(composite_signal)

        return noisy_signal, measurements
