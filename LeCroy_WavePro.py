#!/usr/bin/env python

import InstrumentDriver
from VISA_Driver import VISA_Driver
from InstrumentConfig import InstrumentQuantity
import numpy as np
import struct
import time

__version__ = "0.1.0"

class Error(Exception):
    pass

class Driver(VISA_Driver):
    """ This class implements the LeCroy WavePro scope driver"""

    def performOpen(self, options={}):
        """Perform the operation of opening the instrument connection"""
        # Open connection
        VISA_Driver.performOpen(self, options)
        # Set up the scope for remote operation
        self.writeAndLog('CHDR OFF')  # Turn off headers for easier parsing
        # Set communication format
        self.writeAndLog('CFMT DEF9,WORD,BIN')  # 16-bit words, binary

    def performSetValue(self, quant, value, sweepRate=0.0, options={}):
        """Perform the Set Value instrument operation. This function should
        return the actual value set by the instrument"""

        # Handle trigger-related commands
        if quant.name == 'Trig slope':
            sTrig = self.getCmdStringFromValue('Trig source')
            quant.set_cmd = '%s:TRSL' % sTrig
        elif quant.name == 'Trig level':
            sTrig = self.getCmdStringFromValue('Trig source')
            quant.set_cmd = '%s:TRLV' % sTrig

        # Handle timebase settings
        elif quant.name == 'Time range':
            # Time range in seconds, convert to time/div (10 divisions)
            time_div = value / 10.0
            self.writeAndLog('TDIV %.6E' % time_div)
            return value

        elif quant.name == 'Sample rate':
            # Set sample rate (samples per second)
            # Note: LeCroy uses SARA (sample rate) or MSIZ (memory size)
            # Setting sample rate directly if supported
            try:
                self.writeAndLog('SARA %.6E' % value)
            except:
                self.log('Warning: Could not set sample rate directly')
            return value

        elif quant.name == 'Memory length':
            # Set memory depth (number of points)
            self.writeAndLog('MSIZ %d' % int(value))
            return value

        elif quant.name == 'Trigger mode':
            # Handle different trigger modes
            if value == 'Auto':
                self.writeAndLog('TRMD AUTO')
            elif value == 'Normal':
                self.writeAndLog('TRMD NORM')
            elif value == 'Single':
                self.writeAndLog('TRMD SINGLE')
            elif value == 'Stop':
                self.writeAndLog('TRMD STOP')
            return value

        # Run standard VISA case with updated commands
        value = VISA_Driver.performSetValue(self, quant, value, sweepRate, options)
        return value


    def performGetValue(self, quant, options={}):
        """Perform the Get Value instrument operation"""

        # Check for trace data requests
        if quant.name in ('Ch1 - Data', 'Ch2 - Data', 'Ch3 - Data', 'Ch4 - Data'):
            # Extract channel number
            channel = int(quant.name[2])

            # Check if channel is enabled
            if self.getValue('Ch%d - Enabled' % channel):
                # Get current trigger mode to handle acquisition properly
                try:
                    trig_mode = self.getValue('Trigger mode')
                except:
                    trig_mode = 'Normal'  # Default

                # For Single mode, arm the trigger before reading
                if trig_mode == 'Single':
                    self.writeAndLog('TRMD SINGLE')
                    # Wait for trigger and acquisition to complete
                    self._waitForAcquisition(timeout=10.0)
                elif trig_mode == 'Normal':
                    # For Normal mode, ensure we have a fresh acquisition
                    # Wait for acquisition to complete
                    self._waitForAcquisition(timeout=10.0)

                # Get waveform data
                value = self._readTrace(channel)
            else:
                # Channel not enabled, return empty array
                value = InstrumentQuantity.getTraceDict([])

        elif quant.name == 'Trig source':
            # Trigger source, parse from TRSE? response
            sAns = self.askAndLog('TRSE?').strip()
            i1 = sAns.find(',SR,') + 4
            i2 = sAns.find(',HT')
            # Convert response to a number
            value = quant.getValueFromCmdString(sAns[i1:i2])

        elif quant.name == 'Trig level':
            # Trigger level for the selected source
            sTrig = self.getCmdStringFromValue('Trig source')
            sAns = self.askAndLog('%s:TRLV?' % sTrig).strip()
            i1 = sAns.find('V')
            value = float(sAns[:i1])

        elif quant.name == 'Trig slope':
            # Trigger slope for the selected source
            sTrig = self.getCmdStringFromValue('Trig source')
            quant.set_cmd = '%s:TRSL' % sTrig
            value = VISA_Driver.performGetValue(self, quant, options)

        elif quant.name == 'Time range':
            # Get time/div and multiply by 10 for total range
            sAns = self.askAndLog('TDIV?').strip()
            time_div = float(sAns.split()[-1].rstrip('S'))
            value = time_div * 10.0

        elif quant.name == 'Sample rate':
            # Get actual sample rate
            sAns = self.askAndLog('SARA?').strip()
            # Parse sample rate from response
            value = float(sAns.split()[-1].rstrip('Sa/s'))

        elif quant.name == 'Memory length':
            # Get memory depth
            sAns = self.askAndLog('MSIZ?').strip()
            value = float(sAns.split()[-1])

        elif quant.name == 'Trigger mode':
            # Get trigger mode
            sAns = self.askAndLog('TRMD?').strip()
            if 'AUTO' in sAns.upper():
                value = 'Auto'
            elif 'NORM' in sAns.upper():
                value = 'Normal'
            elif 'SINGLE' in sAns.upper():
                value = 'Single'
            elif 'STOP' in sAns.upper():
                value = 'Stop'
            else:
                value = 'Normal'

        else:
            # For all other cases, call VISA driver
            value = VISA_Driver.performGetValue(self, quant, options)

        return value


    def _waitForAcquisition(self, timeout=10.0):
        """Wait for the oscilloscope to complete acquisition"""
        start_time = time.time()
        while (time.time() - start_time) < timeout:
            # Query the trigger status
            try:
                # INR? returns Internal State Change Register
                # Bit 0 is set when a new waveform is available
                status = self.askAndLog('INR?')
                # Check if acquisition is complete
                # For LeCroy, we can also check TRMD? or wait state
                trig_state = self.askAndLog('TRMD?')
                if 'STOP' in trig_state.upper():
                    # Acquisition complete
                    return True

                # Alternative: check using *OPC? or WAIT command
                # Some LeCroy scopes support WAIT command
                try:
                    self.writeAndLog('WAIT')
                    return True
                except:
                    pass

            except:
                pass

            time.sleep(0.05)  # Wait 50ms before checking again

        # Timeout occurred
        self.log('Warning: Acquisition timeout')
        return False


    def _readTrace(self, channel):
        """Read trace data from the specified channel

        Args:
            channel: Channel number (1-4)

        Returns:
            Dictionary with trace data and timing information
        """
        try:
            # Request waveform descriptor
            self.write('C%d:WF? DESC;' % channel, bCheckError=False)
            sDesc = self.read(ignore_termination=True)

            # Parse descriptor header
            # Find the start of data: #9NNNNNNNNN where N is byte count
            indx = sDesc.find(b'#9')
            if indx < 0:
                # Try alternate format
                indx = sDesc.find(b'#8')
                if indx < 0:
                    self.log('Error: Could not find data header')
                    return InstrumentQuantity.getTraceDict([])
                byte_count_len = 8
            else:
                byte_count_len = 9

            # Skip header and byte count
            sDesc = sDesc[indx+2+byte_count_len:]

            # Extract waveform parameters from descriptor
            # See LeCroy Remote Control Manual for WAVEDESC structure
            iFirst = struct.unpack('<i', sDesc[124:128])[0]  # FIRST_VALID_PNT
            iLast = struct.unpack('<i', sDesc[128:132])[0]   # LAST_VALID_PNT
            Vgain = struct.unpack('<f', sDesc[156:160])[0]   # VERTICAL_GAIN
            Voffs = struct.unpack('<f', sDesc[160:164])[0]   # VERTICAL_OFFSET
            dt = struct.unpack('<f', sDesc[176:180])[0]      # HORIZ_INTERVAL
            t0 = struct.unpack('<d', sDesc[180:188])[0]      # HORIZ_OFFSET
            nPts = struct.unpack('<i', sDesc[116:120])[0]    # WAVE_ARRAY_COUNT

            # Request waveform data
            self.write('C%d:WF? DAT1;' % channel, bCheckError=False)
            sData = self.read(ignore_termination=True)

            # Parse data header
            # Format: C1:WF DAT1,#9NNNNNNNNN<data>
            head_start = sData.find(b'#9')
            if head_start < 0:
                head_start = sData.find(b'#8')
                if head_start < 0:
                    self.log('Error: Could not find data header')
                    return InstrumentQuantity.getTraceDict([])
                byte_count_len = 8
                data_start = head_start + 2 + 8
            else:
                byte_count_len = 9
                data_start = head_start + 2 + 9

            # Extract the actual waveform data
            # LeCroy sends 16-bit signed integers
            data_bytes = sData[data_start:]

            # Determine actual number of points available
            nAvailable = len(data_bytes) // 2  # 2 bytes per 16-bit word

            # Use the minimum of what we expect and what we have
            nPoints = min(iLast - iFirst + 1, nAvailable)

            # Convert binary data to numpy array
            # LeCroy uses little-endian 16-bit signed integers
            vData = np.frombuffer(data_bytes[:nPoints*2], dtype='<i2')

            # Apply gain and offset to convert to voltage
            vData = vData.astype(float) * Vgain + Voffs

            # Create trace dictionary with timing information
            value = InstrumentQuantity.getTraceDict(vData, dt=dt, t0=t0)

        except Exception as e:
            self.log('Error reading trace: %s' % str(e))
            value = InstrumentQuantity.getTraceDict([])

        return value


if __name__ == '__main__':
    pass
