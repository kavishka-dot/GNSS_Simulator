import numpy as np

class Receiver:
    def __init__(self, fs):
        self.fs = fs  # Sampling frequency in Hz
        self.speed_of_light = 299_792_458  # m/s
        self.chip_rate = 1.023e6
        self.samples_per_chip = int(fs / self.chip_rate)
        self.samples_per_ms = int(fs * 1e-3)

    def acquire(self, signal, prn, doppler_range=5000, doppler_step=250):
        ca_code = self.generate_ca_code(prn)
        ca_oversampled = np.repeat(ca_code, self.samples_per_chip)
        ca_oversampled = ca_oversampled[:len(signal)]

        ca_oversampled = ca_oversampled / np.linalg.norm(ca_oversampled)
        signal = signal / np.linalg.norm(signal)

        best_peak, best_doppler, best_phase = 0, 0, 0
        for doppler in np.arange(-doppler_range, doppler_range + 1, doppler_step):
            t = np.arange(len(ca_oversampled)) / self.fs
            carrier = np.exp(-1j * 2 * np.pi * doppler * t)
            mixed = signal[:len(ca_oversampled)] * carrier
            corr = np.abs(np.correlate(mixed, ca_oversampled, mode='valid'))
            peak = np.max(corr)
            if peak > best_peak:
                best_peak = peak
                best_doppler = doppler
                best_phase = np.argmax(corr)
        return best_doppler, best_phase, best_peak

    def track(self, signal, prn, init_phase, init_doppler, steps=5):
        ca_code = self.generate_ca_code(prn)
        ca_upsampled = np.repeat(ca_code, self.samples_per_chip)
        ca_upsampled = np.tile(ca_upsampled, 10)

        phase = 0
        code_phase = init_phase
        doppler = init_doppler

        for step in range(steps):
            start = step * self.samples_per_ms
            end = start + self.samples_per_ms
            if end > len(signal): break
            segment = signal[start:end]

            t = np.arange(len(segment)) / self.fs
            carrier = np.exp(-1j * (2 * np.pi * doppler * t + phase))

            early = np.roll(ca_upsampled, int(code_phase - 0.5 * self.samples_per_chip))[:len(segment)]
            prompt = np.roll(ca_upsampled, int(code_phase))[:len(segment)]
            late = np.roll(ca_upsampled, int(code_phase + 0.5 * self.samples_per_chip))[:len(segment)]

            E = np.sum(segment * early * carrier)
            P = np.sum(segment * prompt * carrier)
            L = np.sum(segment * late * carrier)

            dll_error = np.abs(E)**2 - np.abs(L)**2
            pll_error = np.arctan2(P.imag, P.real)

            code_phase += 0.01 * dll_error
            code_phase = np.clip(code_phase, 0, len(ca_upsampled) - 1)

            doppler += 1.5 * pll_error
            phase += pll_error
            phase %= 2 * np.pi

        return code_phase, doppler

    def compute_pseudorange(self, code_phase_samples):
        chips = code_phase_samples / self.samples_per_chip
        return chips * self.speed_of_light / self.chip_rate

    def extract_nav_bits(self, signal, prn, code_phase, doppler, num_bits=300):
        ca_code = self.generate_ca_code(prn)
        ca_upsampled = np.repeat(ca_code, self.samples_per_chip)
        ca_upsampled = np.tile(ca_upsampled, 25)

        samples_per_bit = 20 * self.samples_per_ms
        bits = []

        for i in range(num_bits):
            start = int(code_phase + i * samples_per_bit)
            end = start + samples_per_bit
            if end > len(signal): break

            segment = signal[start:end]
            t = np.arange(len(segment)) / self.fs
            carrier = np.exp(-1j * 2 * np.pi * doppler * t)
            wiped = segment * carrier

            code = np.tile(ca_upsampled, 20)[:samples_per_bit]
            corr = np.dot(wiped.real, code)
            bits.append(1 if corr > 0 else 0)

        return bits

    def find_preamble(self, bitstream, pattern=[1,0,0,0,1,0,1,1]):
        for i in range(len(bitstream) - len(pattern)):
            if bitstream[i:i+8] == pattern:
                return i
        return -1

    def parse_subframe1(self, bits):
        def bits_to_int(b): return int(''.join(str(x) for x in b), 2)
        def twos_complement(v, w): return v - (1 << w) if v & (1 << (w-1)) else v

        if len(bits) < 300: raise ValueError("Insufficient bits for subframe.")
        week_number = bits_to_int(bits[60:70])
        ura_index = bits_to_int(bits[70:74])
        sv_health = bits_to_int(bits[74:80])
        iodc = bits_to_int(bits[80:82] + bits[210:218])
        toc = bits_to_int(bits[218:234]) * (2**4)
        af2 = twos_complement(bits_to_int(bits[240:248]), 8) * (2**-55)
        af1 = twos_complement(bits_to_int(bits[248:266]), 16) * (2**-43)
        af0 = twos_complement(bits_to_int(bits[266:286]), 22) * (2**-31)
        return {
            "Week Number": week_number,
            "URA Index": ura_index,
            "SV Health": sv_health,
            "IODC": iodc,
            "T_oc (s)": toc,
            "af0": af0,
            "af1": af1,
            "af2": af2
        }

    def generate_ca_code(self, prn):
        g2_delays = [
            5, 6, 7, 8, 17, 18, 139, 140, 141, 251,
            252, 254, 255, 256, 257, 258, 469, 470, 471, 472,
            473, 474, 509, 512, 513, 514, 515, 516, 859, 860,
            861, 862
        ]
        g2_delay = g2_delays[prn - 1]
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

        return [1 if bit == 0 else -1 for bit in ca]
