import numpy as np

class Transmitter:
    def __init__(self, prn, initial_position_ecef, velocity_ecef=np.zeros(3)):
        self.prn = prn  # Satellite PRN number (1–32)
        self.position = np.array(initial_position_ecef, dtype=float)
        self.velocity = np.array(velocity_ecef, dtype=float)  # m/s
        self.chip_rate = 1.023e6
        self.nav_rate = 50  # 50 bps NAV message
        self.c = 299_792_458  # Speed of light

    def update_position(self, dt):
        """ Update satellite position linearly over dt seconds """
        self.position += self.velocity * dt

    def get_position(self):
        """ Return current ECEF position """
        return self.position

    def generate_ca_code(self):
        g2_delays = [
            5, 6, 7, 8, 17, 18, 139, 140, 141, 251,
            252, 254, 255, 256, 257, 258, 469, 470, 471, 472,
            473, 474, 509, 512, 513, 514, 515, 516, 859, 860,
            861, 862
        ]
        g2_delay = g2_delays[self.prn - 1]
        g1 = [1] * 10
        g2 = [1] * 10
        ca = []

        for _ in range(1023):
            g1_out = g1[-1]
            g2_out = g2[g2_delay - 1]
            ca.append(g1_out ^ g2_out)
            g1 = [g1[i - 1] if i > 0 else g1[2] ^ g1[9] for i in range(10)]
            g2_feedback = g2[1] ^ g2[2] ^ g2[5] ^ g2[7] ^ g2[8] ^ g2[9]
            g2 = [g2_feedback] + g2[:-1]

        return np.array([1 if bit == 0 else -1 for bit in ca])

    def generate_signal(self, fs, duration_ms, code_delay_chips=0, doppler_hz=0):
        """
        Generate baseband GNSS signal for this satellite:
        - fs: sampling frequency in Hz
        - duration_ms: signal length in milliseconds
        - code_delay_chips: fractional delay to apply to PRN
        - doppler_hz: carrier frequency shift due to motion
        """
        ca_code = self.generate_ca_code()
        samples_per_chip = int(fs / self.chip_rate)
        total_samples = int(fs * duration_ms * 1e-3)

        ca_upsampled = np.repeat(ca_code, samples_per_chip)
        ca_repeated = np.tile(ca_upsampled, int(np.ceil(total_samples / len(ca_upsampled))))
        ca_repeated = ca_repeated[:total_samples]

        delay_samples = int(code_delay_chips * samples_per_chip)
        ca_delayed = np.roll(ca_repeated, delay_samples)

        t = np.arange(total_samples) / fs
        carrier = np.cos(2 * np.pi * doppler_hz * t)
        signal = ca_delayed * carrier
        return signal

    def simulate_transmission(self, receiver_pos_ecef, fs, duration_ms):
        """
        Compute:
        - Geometric range
        - Time delay
        - Doppler shift
        Then generate signal.
        """
        r = receiver_pos_ecef - self.position
        range_m = np.linalg.norm(r)
        unit_r = r / range_m

        relative_velocity = -np.dot(unit_r, self.velocity)
        doppler_hz = -relative_velocity / self.c * self.chip_rate

        code_delay_chips = range_m / (self.c / self.chip_rate)

        signal = self.generate_signal(fs, duration_ms, code_delay_chips, doppler_hz)
        return signal, range_m, doppler_hz
