#!/usr/bin/env python
"""
Test script for LeCroy WavePro driver

This script demonstrates how to use the driver and validates basic functionality.
Update the SCOPE_ADDRESS with your oscilloscope's IP address.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

# Configuration
SCOPE_ADDRESS = '192.168.1.100'  # Update with your scope's IP address
TEST_CHANNEL = 1
TIME_RANGE = 1.0E-3  # 1 ms

def test_connection():
    """Test basic connection to the scope"""
    print("=" * 60)
    print("Test 1: Connection Test")
    print("=" * 60)

    try:
        from LeCroy_WavePro import Driver

        # Create driver instance
        driver = Driver()

        # Open connection
        print(f"Connecting to scope at {SCOPE_ADDRESS}...")
        driver.performOpen(options={'address': SCOPE_ADDRESS})

        # Query identification
        idn = driver.askAndLog('*IDN?')
        print(f"Connected to: {idn}")

        driver.performClose()
        print("✓ Connection test PASSED\n")
        return True

    except Exception as e:
        print(f"✗ Connection test FAILED: {e}\n")
        return False


def test_timebase_control():
    """Test timebase configuration"""
    print("=" * 60)
    print("Test 2: Timebase Control")
    print("=" * 60)

    try:
        from LeCroy_WavePro import Driver
        from InstrumentConfig import InstrumentQuantity

        driver = Driver()
        driver.performOpen(options={'address': SCOPE_ADDRESS})

        # Create a mock quantity for time range
        class MockQuant:
            def __init__(self, name):
                self.name = name
                self.set_cmd = ''

        # Test setting time range
        quant = MockQuant('Time range')
        print(f"Setting time range to {TIME_RANGE} s...")
        driver.performSetValue(quant, TIME_RANGE)

        # Read back time range
        time_range = driver.performGetValue(quant)
        print(f"Read back time range: {time_range} s")

        # Check if values match (within 10%)
        if abs(time_range - TIME_RANGE) / TIME_RANGE < 0.1:
            print("✓ Timebase control test PASSED\n")
            result = True
        else:
            print(f"✗ Timebase control test FAILED: Set {TIME_RANGE}, got {time_range}\n")
            result = False

        driver.performClose()
        return result

    except Exception as e:
        print(f"✗ Timebase control test FAILED: {e}\n")
        return False


def test_trigger_modes():
    """Test different trigger modes"""
    print("=" * 60)
    print("Test 3: Trigger Modes")
    print("=" * 60)

    try:
        from LeCroy_WavePro import Driver

        driver = Driver()
        driver.performOpen(options={'address': SCOPE_ADDRESS})

        class MockQuant:
            def __init__(self, name):
                self.name = name
                self.set_cmd = ''

        modes = ['Auto', 'Normal', 'Single']

        for mode in modes:
            print(f"Testing {mode} mode...")
            quant = MockQuant('Trigger mode')
            driver.performSetValue(quant, mode)

            # Read back
            read_mode = driver.performGetValue(quant)
            print(f"  Set: {mode}, Read: {read_mode}")

            if mode.upper() in read_mode.upper():
                print(f"  ✓ {mode} mode OK")
            else:
                print(f"  ✗ {mode} mode FAILED")
                driver.performClose()
                return False

        driver.performClose()
        print("✓ Trigger modes test PASSED\n")
        return True

    except Exception as e:
        print(f"✗ Trigger modes test FAILED: {e}\n")
        return False


def test_data_acquisition():
    """Test acquiring waveform data"""
    print("=" * 60)
    print("Test 4: Data Acquisition")
    print("=" * 60)

    try:
        from LeCroy_WavePro import Driver
        from InstrumentConfig import InstrumentQuantity

        driver = Driver()
        driver.performOpen(options={'address': SCOPE_ADDRESS})

        class MockQuant:
            def __init__(self, name):
                self.name = name
                self.set_cmd = ''

            @staticmethod
            def getValueFromCmdString(s):
                return s

            @staticmethod
            def getCmdStringFromValue(v):
                return 'C1'

        # Enable channel
        print(f"Enabling channel {TEST_CHANNEL}...")
        ch_enable = MockQuant(f'Ch{TEST_CHANNEL} - Enabled')
        driver.performSetValue(ch_enable, True)

        # Set time range
        time_quant = MockQuant('Time range')
        driver.performSetValue(time_quant, TIME_RANGE)

        # Set trigger to auto for reliable acquisition
        trig_mode = MockQuant('Trigger mode')
        driver.performSetValue(trig_mode, 'Auto')

        # Add getValue methods to driver for this test
        driver._getValue = driver.getValue
        driver.getValue = lambda name: {
            f'Ch{TEST_CHANNEL} - Enabled': True,
            'Trigger mode': 'Auto'
        }.get(name, driver._getValue(name))

        driver.getCmdStringFromValue = MockQuant.getCmdStringFromValue

        # Acquire data
        print(f"Acquiring data from channel {TEST_CHANNEL}...")
        data_quant = MockQuant(f'Ch{TEST_CHANNEL} - Data')
        data = driver.performGetValue(data_quant)

        # Check data
        if 'y' in data and len(data['y']) > 0:
            print(f"✓ Acquired {len(data['y'])} points")
            print(f"  Voltage range: {np.min(data['y']):.3f} to {np.max(data['y']):.3f} V")
            if 't0' in data and 'dt' in data:
                total_time = len(data['y']) * data['dt']
                print(f"  Time span: {total_time*1e3:.3f} ms")
            print("✓ Data acquisition test PASSED\n")
            result = True
        else:
            print("✗ Data acquisition test FAILED: No data returned\n")
            result = False

        driver.performClose()
        return result

    except Exception as e:
        print(f"✗ Data acquisition test FAILED: {e}\n")
        import traceback
        traceback.print_exc()
        return False


def test_data_plot():
    """Acquire and plot waveform data"""
    print("=" * 60)
    print("Test 5: Data Visualization")
    print("=" * 60)

    try:
        from LeCroy_WavePro import Driver

        driver = Driver()
        driver.performOpen(options={'address': SCOPE_ADDRESS})

        class MockQuant:
            def __init__(self, name):
                self.name = name
                self.set_cmd = ''

            @staticmethod
            def getValueFromCmdString(s):
                return s

            @staticmethod
            def getCmdStringFromValue(v):
                return 'C1'

        # Setup
        ch_enable = MockQuant(f'Ch{TEST_CHANNEL} - Enabled')
        driver.performSetValue(ch_enable, True)

        time_quant = MockQuant('Time range')
        driver.performSetValue(time_quant, TIME_RANGE)

        trig_mode = MockQuant('Trigger mode')
        driver.performSetValue(trig_mode, 'Auto')

        # Helper methods
        driver.getValue = lambda name: {
            f'Ch{TEST_CHANNEL} - Enabled': True,
            'Trigger mode': 'Auto'
        }.get(name, True)
        driver.getCmdStringFromValue = MockQuant.getCmdStringFromValue

        # Acquire
        print(f"Acquiring waveform...")
        data_quant = MockQuant(f'Ch{TEST_CHANNEL} - Data')
        data = driver.performGetValue(data_quant)

        if 'y' in data and len(data['y']) > 0:
            # Create time axis
            if 't0' in data and 'dt' in data:
                time = data['t0'] + np.arange(len(data['y'])) * data['dt']
            else:
                time = np.arange(len(data['y']))

            # Plot
            plt.figure(figsize=(12, 6))
            plt.plot(time * 1e3, data['y'])  # Convert to ms
            plt.xlabel('Time (ms)')
            plt.ylabel('Voltage (V)')
            plt.title(f'LeCroy WavePro - Channel {TEST_CHANNEL}')
            plt.grid(True, alpha=0.3)

            # Save plot
            filename = 'lecroy_test_waveform.png'
            plt.savefig(filename, dpi=150, bbox_inches='tight')
            print(f"✓ Plot saved to {filename}")
            print("✓ Data visualization test PASSED\n")
            result = True
        else:
            print("✗ Data visualization test FAILED: No data\n")
            result = False

        driver.performClose()
        return result

    except ImportError:
        print("⚠ Matplotlib not available, skipping plot test\n")
        return True
    except Exception as e:
        print(f"✗ Data visualization test FAILED: {e}\n")
        return False


def main():
    """Run all tests"""
    print("\n" + "=" * 60)
    print("LeCroy WavePro Driver Test Suite")
    print("=" * 60)
    print(f"Scope address: {SCOPE_ADDRESS}")
    print(f"Test channel: {TEST_CHANNEL}")
    print(f"Time range: {TIME_RANGE} s")
    print()

    # Check if we can import the driver
    try:
        from LeCroy_WavePro import Driver
    except ImportError as e:
        print(f"ERROR: Cannot import driver: {e}")
        print("Make sure LeCroy_WavePro.py is in the same directory or in PYTHONPATH")
        return

    # Run tests
    results = []

    results.append(("Connection", test_connection()))
    results.append(("Timebase Control", test_timebase_control()))
    results.append(("Trigger Modes", test_trigger_modes()))
    results.append(("Data Acquisition", test_data_acquisition()))
    results.append(("Data Visualization", test_data_plot()))

    # Summary
    print("=" * 60)
    print("Test Summary")
    print("=" * 60)

    passed = sum(1 for _, result in results if result)
    total = len(results)

    for name, result in results:
        status = "✓ PASSED" if result else "✗ FAILED"
        print(f"{name:25s} {status}")

    print(f"\nTotal: {passed}/{total} tests passed")

    if passed == total:
        print("\n🎉 All tests passed! Driver is working correctly.")
    else:
        print(f"\n⚠ {total - passed} test(s) failed. Check the output above for details.")


if __name__ == '__main__':
    main()
