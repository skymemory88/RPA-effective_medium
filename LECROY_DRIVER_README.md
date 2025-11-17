# LeCroy WavePro Driver for Labber

## Overview

This is an improved driver for the LeCroy WavePro 254HD-MS oscilloscope for use with Keysight Labber measurement software. This version addresses several critical issues with data acquisition, timing control, and trigger mode handling.

## Issues Fixed

### 1. **Specified Time Length Control**
**Problem:** The original driver had no way to set the time range or sample rate, so it would just read whatever was currently displayed on the scope.

**Solution:**
- Added `Time range` parameter to control the total time window (in seconds)
- Added `Sample rate` parameter to set the sampling rate
- Added `Memory length` parameter to control the number of points acquired
- These parameters allow you to precisely specify how much data to capture

### 2. **Trigger Mode Handling**
**Problem:** Different trigger modes ('single', 'normal', 'auto') behaved inconsistently because the driver didn't synchronize with the acquisition state.

**Solution:**
- Added proper `Trigger mode` parameter with support for:
  - **Auto**: Scope always triggers (free-running)
  - **Normal**: Waits for valid trigger condition
  - **Single**: Captures one trace then stops
  - **Stop**: Acquisition stopped
- Implemented `_waitForAcquisition()` method to properly synchronize with scope
- For Single mode, driver now arms the trigger and waits for completion
- For Normal mode, driver waits for acquisition to finish before reading

### 3. **Data Acquisition and Parsing Issues**
**Problem:**
- Used deprecated `np.fromstring()` function
- Incorrect byte order assumption (big-endian vs little-endian)
- Improper handling of data length
- No error handling

**Solution:**
- Replaced `np.fromstring()` with `np.frombuffer()`
- Corrected byte order to little-endian (`<i2`) which is LeCroy standard
- Proper parsing of waveform descriptor structure
- Added try-except blocks for robust error handling
- Correctly extracts:
  - `FIRST_VALID_PNT` and `LAST_VALID_PNT` for valid data range
  - `VERTICAL_GAIN` and `VERTICAL_OFFSET` for voltage scaling
  - `HORIZ_INTERVAL` for time step (dt)
  - `HORIZ_OFFSET` for time offset (t0)
  - `WAVE_ARRAY_COUNT` for number of points

### 4. **Additional Improvements**
- Added `performOpen()` to initialize scope settings on connection
- Set communication format to binary words for faster transfer
- Disabled headers (`CHDR OFF`) for easier parsing
- Added comprehensive INI configuration file with all parameters
- Added per-channel settings (coupling, voltage range, offset)
- Proper handling of channel enable/disable state

## Usage

### Installation

1. Copy both files to your Labber driver directory:
   - `LeCroy_WavePro.py`
   - `LeCroy_WavePro.ini`

2. Place them in: `[Labber Installation]/Labber/Drivers/LeCroy_WavePro/`

3. Restart Labber

### Configuration

#### Timebase Settings
```
Time range: 1.0E-3 s        # Total time window (10 divisions)
Sample rate: 1.0E9 Sa/s     # Sampling rate (optional)
Memory length: 10000         # Number of points (optional)
```

**Note:** For most applications, just set `Time range` and let the scope determine optimal sample rate and memory depth. Enable custom settings only if you need specific values.

#### Trigger Settings
```
Trigger mode: Normal         # Auto/Normal/Single/Stop
Trig source: Channel 1       # Which channel to trigger on
Trig slope: Positive         # Positive or Negative edge
Trig level: 0.0 V           # Trigger threshold
```

#### Channel Settings (per channel)
```
Ch1 - Enabled: True          # Enable channel
Ch1 - Coupling: DC 1M        # DC/AC, 1M/50 ohm
Ch1 - Voltage range: 1.0 V   # Volts per division
Ch1 - Offset: 0.0 V          # Vertical offset
Ch1 - Data: [Read Only]      # Waveform output
```

### Example: Capturing a Single Trace

1. **Setup:**
   - Trigger mode: Single
   - Time range: 1.0E-6 (1 µs)
   - Trig source: Channel 1
   - Trig level: 0.5 V
   - Ch1 - Enabled: True

2. **Acquire:**
   - Read `Ch1 - Data`
   - Driver will:
     - Arm single trigger
     - Wait for trigger event
     - Wait for acquisition to complete
     - Read waveform data
     - Return voltage vs time array

### Example: Continuous Acquisition

1. **Setup:**
   - Trigger mode: Normal (or Auto)
   - Time range: 1.0E-3 (1 ms)
   - Ch1 - Enabled: True

2. **Acquire:**
   - Repeatedly read `Ch1 - Data`
   - Each read waits for a fresh acquisition

### Example: Long Time Capture

1. **Setup:**
   - Time range: 10.0 (10 seconds)
   - Use custom memory length: True
   - Memory length: 1000000 (1M points)
   - This gives 100 kSa/s effective rate

2. **Acquire:**
   - Read `Ch1 - Data`
   - Returns 1M point trace over 10 seconds

## Technical Details

### Binary Data Format

LeCroy scopes return data in this format:
```
C1:WF DAT1,#9NNNNNNNNN<binary data>
```

Where:
- `#9` indicates 9-digit byte count follows
- `NNNNNNNNN` is the number of bytes
- Binary data is 16-bit signed integers, little-endian

### Waveform Descriptor Structure

The descriptor contains metadata at fixed byte offsets:
- Bytes 116-120: Wave array count (int32)
- Bytes 124-128: First valid point (int32)
- Bytes 128-132: Last valid point (int32)
- Bytes 156-160: Vertical gain (float32)
- Bytes 160-164: Vertical offset (float32)
- Bytes 176-180: Horizontal interval (float32)
- Bytes 180-188: Horizontal offset (float64)

All multi-byte values are little-endian.

### Voltage Calculation

```python
voltage = (raw_data * vertical_gain) + vertical_offset
```

### Time Axis

```python
time[i] = horizontal_offset + (i * horizontal_interval)
```

## Troubleshooting

### Problem: "Acquisition timeout" warning
**Cause:** Scope isn't triggering within timeout period (10 seconds)
**Solution:**
- Check trigger level is within signal range
- Try Auto trigger mode for testing
- Verify signal is present on trigger channel

### Problem: Empty waveform returned
**Cause:** Channel is disabled or no data available
**Solution:**
- Verify channel is enabled in settings
- Check scope display shows waveform
- Try manual trigger on scope first

### Problem: Wrong time range captured
**Cause:** Time range not set properly
**Solution:**
- Explicitly set Time range parameter
- Read back current Time range to verify
- Check scope display matches expected range

### Problem: Data looks wrong in Single mode
**Cause:** Reading stale data from previous acquisition
**Solution:**
- Driver now handles this automatically
- It arms trigger and waits for new acquisition
- If still issues, try adding delay between acquisitions

## Command Reference

### Key SCPI Commands Used

```
CHDR OFF                    # Disable headers
CFMT DEF9,WORD,BIN         # Set binary word format
TDIV <value>               # Set time/div
MSIZ <value>               # Set memory size
SARA <value>               # Set sample rate
TRMD AUTO|NORM|SINGLE|STOP # Set trigger mode
C1:TRA ON|OFF              # Channel on/off
C1:CPL D1M|A1M|D50|A50|GND # Channel coupling
C1:VDIV <value>            # Volts per division
C1:OFST <value>            # Voltage offset
C1:WF? DESC                # Get waveform descriptor
C1:WF? DAT1                # Get waveform data
TRSE?                      # Query trigger setup
C1:TRLV <value>            # Set trigger level
C1:TRSL POS|NEG            # Set trigger slope
INR?                       # Query internal state
WAIT                       # Wait for operations complete
```

## Testing Checklist

- [ ] Connect to scope via TCPIP
- [ ] Enable Ch1
- [ ] Set Time range to 1 ms
- [ ] Set Trigger mode to Auto
- [ ] Read Ch1 - Data (should return waveform)
- [ ] Change to Normal trigger mode
- [ ] Read Ch1 - Data (should wait for trigger)
- [ ] Change to Single trigger mode
- [ ] Read Ch1 - Data (should capture once)
- [ ] Verify time axis matches Time range setting
- [ ] Try different time ranges
- [ ] Try multiple channels simultaneously

## Known Limitations

1. **Sample rate control:** Some LeCroy models may not support direct sample rate setting via `SARA` command. In this case, sample rate is determined by time/div and memory depth.

2. **Acquisition wait:** The `_waitForAcquisition()` method uses polling with 50ms interval. For very fast repetitive acquisitions, this may add latency.

3. **Sequence mode:** This driver does not support sequence acquisition mode (multiple segments).

4. **Math/Memory traces:** Only supports physical channels (C1-C4), not math or memory traces.

## Version History

**v0.1.0** (2025-11-17)
- Initial release with fixed time control and trigger handling
- Proper binary data parsing
- Acquisition synchronization
- Comprehensive parameter set

## References

- LeCroy WavePro Remote Control Manual
- Labber Instrument Driver Development Guide
- LeCroy Automation Command Reference

## Support

For issues specific to this driver, check:
1. Scope firmware is up to date
2. TCPIP connection is stable
3. Scope responds to manual SCPI commands via terminal
4. Labber log files for detailed error messages
