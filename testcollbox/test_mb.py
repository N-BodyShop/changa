import numpy as np
import pynbody as pb
from scipy.stats import maxwell
from scipy.stats import kstest

snap = pb.load('uniform_box.000128')
speeds = np.linalg.norm(snap['vel'], axis=1)

xvals = np.linspace(0, np.max(speeds))
param = maxwell.fit(speeds, floc=0)
D, p_value = kstest(speeds, 'maxwell', args=param)

print("Comparing similarity of speeds to MB distribution")
print(f"P-value: {p_value:.3f}")

if p_value > 0.05:
    print("Test passed")
else:
    print("Test failed")
