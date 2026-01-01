import pyvisa
import time
import numpy as np

# Initialize VISA resource manager
rm = pyvisa.ResourceManager('C:/Windows/System32/visa32.dll')
scope = rm.open_resource('USB0::0x05FF::0xEE3A::LCRY2251C01483::INSTR')
scope.timeout = float('+inf')


def gather_waveform():
    """
    Returns raw waveform data from channels 1 and 2
    as 1-byte-per-point integers.
    Prints the current horizontal magnification (HOR_MAGNIFY) before acquisition.
    """
    # WaveAce Format: 16-bit
    scope.write("COMM_HEADER OFF")
    scope.write("COMM_FORMAT WORD,BIN")  # 16-bit
    scope.write("COMM_ORDER LSB")

    # Optional: weniger Punkte
    scope.write("WFSU SP,0,NP,12000,FP,12000")

    # Daten abfragen als signed int16
    samples = scope.query_binary_values(

        "C2:WAVEFORM?",
        datatype='h',   # signed 16-bit
        is_big_endian=False,
        container=np.array
    )

    return samples[200:]


# Main loop to log waveforms
filename = "waveform_data.txt"
try:
    with open(filename, "w") as f:
        print(f"Starting waveform collection. Data will be saved in '{filename}'. Press Ctrl+C to stop.")
        while True:
            waveform = gather_waveform()
            f.write(f"CH1: {','.join(map(str, waveform))}\n")
            f.write("\n")
            f.flush()
except KeyboardInterrupt:
    print("Waveform collection stopped.")
finally:
    scope.close()
