import numpy as np
import pytest
import tempfile
from pathlib import Path

import act


@pytest.mark.skipif(not act.io.icartt._ICARTT_AVAILABLE, reason='ICARTT is not installed.')
def test_read_icartt():
    result = act.io.icartt.read_icartt(act.tests.EXAMPLE_AAF_ICARTT)
    assert 'pitch' in result
    assert len(result['time'].values) == 14087
    assert result['true_airspeed'].units == 'm/s'
    assert 'Revision' in result.attrs
    np.testing.assert_almost_equal(result['static_pressure'].mean(), 708.75, decimal=2)

def test_write_icartt():
    with tempfile.TemporaryDirectory() as tmpdirname:
        write_file = Path(tmpdirname)
        ds = act.io.read_arm_netcdf(act.tests.EXAMPLE_MET_SAIL, cleanup_qc=True)
        ds.load()
        # Write the data to the temporary directory
        act.io.icartt.write_icartt(ds, path=str(write_file))
        # Verify output file name is correct
        nfile = str(write_file) + '/' + 'gucmetM1.b1.20230301.000000.ict'
        assert os.path.isfile(nfile) == True
        # Read the file back in
        icartt = act.io.icartt.read_icartt(nfile)
        # Check to make sure all three time variables are 
        # there since it is not 1Hz data
        assert "Time_Start" in list(icartt.keys())
        assert "Time_Stop" in list(icartt.keys())
        # Check the number of dependent variables
        assert icartt.attrs['Dependent_Var_Num'] == 30
        # Verify data handling hasn't changed
        np.testing.assert_almost_equal(icarrt.temp_mean.mean(), -11.456925, decimal=2)

        # close the data objects
        ds.close()
        icartt.close()
        del ds
        del icartt