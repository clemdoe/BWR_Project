# 1D CHF

This tutorial provides 2 examples of 1-D channels with boiling water and
achievement of critical heat flux conditions, both in the case of imposed power
and in the case of imposed temperature.

In the case of imposed power, the power is gradually increased and subsequently
decreased in order to reproduce the well-known hysteresis caused by the
boiling crisis.

In the case of imposed temperature, the temperature is only increased with
time. It shows that with imposed temperature we pass through the Leidenfrost
temperature.

The resulting heat flux vs wall temperature curves can be plotted using the following scripts that are available in each case folder:
``` bash
python3 printResults.py
```

**N.B.1**: The models for CHF and post-CHF have not been validated yet. Use
with care.

**N.B.2**: For the moment the critical heat flux can only be imposed at a
constant value. No models or lookup tables have been implemented yet for its
prediction.
