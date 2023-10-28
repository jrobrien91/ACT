"""
Modules for Reading/Writing the International Consortium for Atmospheric
Research on Transport and Transformation (ICARTT) file format standards V2.0

References:
    ICARTT V2.0 Standards/Conventions:
    - https://www.earthdata.nasa.gov/s3fs-public/imported/ESDS-RFC-029v2.pdf

"""
import numpy as np
import xarray as xr
import datetime

try:
    import icartt
    _ICARTT_AVAILABLE = True
    _format = icartt.Formats.FFI1001
except ImportError:
    _ICARTT_AVAILABLE = False
    _format = None


def read_icartt(filename, format=_format,
                return_None=False, **kwargs):
    """

    Returns `xarray.Dataset` with stored data and metadata from a user-defined
    query of ICARTT from a single datastream. Has some procedures to ensure
    time is correctly fomatted in returned Dataset.

    Parameters
    ----------
    filename : str
        Name of file to read.
    format : str
        ICARTT Format to Read: FFI1001 or FFI2110.
    return_None : bool, optional
        Catch IOError exception when file not found and return None.
        Default is False.
    **kwargs : keywords
        keywords to pass on through to icartt.Dataset.

    Returns
    -------
    ds : xarray.Dataset (or None)
        ACT Xarray dataset (or None if no data file(s) found).

    Examples
    --------
    This example will load the example sounding data used for unit testing.

    .. code-block :: python

        import act
        ds = act.io.icartt.read_icartt(act.tests.sample_files.AAF_SAMPLE_FILE)
        print(ds.attrs['_datastream'])

    """
    if not _ICARTT_AVAILABLE:
        raise ImportError(
            "ICARTT is required to use to read ICARTT files but is not installed")

    ds = None

    # Create an exception tuple to use with try statements. Doing it this way
    # so we can add the FileNotFoundError if requested. Can add more error
    # handling in the future.
    except_tuple = (ValueError,)
    if return_None:
        except_tuple = except_tuple + (FileNotFoundError, OSError)

    try:
        # Read data file with ICARTT dataset.
        ict = icartt.Dataset(filename, format=format, **kwargs)

    except except_tuple as exception:
        # If requested return None for File not found error
        if type(exception).__name__ == 'FileNotFoundError':
            return None

        # If requested return None for File not found error
        if (type(exception).__name__ == 'OSError'
                and exception.args[0] == 'no files to open'):
            return None

    # Define the Uncertainty for each variable. Note it may not be calculated.
    # If not calculated, assign 'N/A' to the attribute
    uncertainty = ict.normalComments[6].split(':')[1].split(',')

    # Define the Upper and Lower Limit of Detection Flags
    ulod_flag = ict.normalComments[7].split(':')[1]
    ulod_value = ict.normalComments[8].split(':')[1].split(',')
    llod_flag = ict.normalComments[9].split(':')[1]
    llod_value = ict.normalComments[10].split(':')[1].split(',')

    # Convert ICARTT Object to Xarray Dataset
    ds_container = []
    # Counter for uncertainty/LOD values
    counter = 0

    # Loop over ICART variables, convert to Xarray DataArray, Append.
    for key in ict.variables:
        # Note time is the only independent variable within ICARTT
        # Short name for time must be "Start_UTC" for ICARTT files.
        if key != 'Start_UTC':
            if key == 'qc_flag':
                key2 = 'quality_flag'
            else:
                key2 = key
            da = xr.DataArray(ict.data[key],
                              coords=dict(time=ict.times),
                              name=key2, dims=['time'])
            # Assume if Uncertainity does not match the number of variables,
            # values were not set within the file. Needs to be string!
            if len(uncertainty) != len(ict.variables):
                da.attrs['uncertainty'] = 'N/A'
            else:
                da.attrs['uncertainty'] = uncertainty[counter]

            # Assume if ULOD does not match the number of variables within the
            # the file, ULOD values were not set.
            if len(ulod_value) != len(ict.variables):
                da.attrs['ULOD_Value'] = 'N/A'
            else:
                da.attrs['ULOD_Value'] = ulod_value[counter]

            # Assume if LLOD does not match the number of variables within the
            # the file, LLOD values were not set.
            if len(llod_value) != len(ict.variables):
                da.attrs['LLOD_Value'] = 'N/A'
            else:
                da.attrs['LLOD_Value'] = llod_value[counter]
            # Define the meta data:
            da.attrs['units'] = ict.variables[key].units
            da.attrs['mvc'] = ict.variables[key].miss
            da.attrs['scale_factor'] = ict.variables[key].scale
            da.attrs['ULOD_Flag'] = ulod_flag
            da.attrs['LLOD_Flag'] = llod_flag
            # Append to ds container
            ds_container.append(da.to_dataset(name=key2))
            # up the counter
            counter += 1

    # Concatenate each of the Xarray DataArrays into a single Xarray DataSet
    ds = xr.merge(ds_container)

    # Assign ICARTT Meta data to Xarray DataSet
    ds.attrs['PI'] = ict.PIName
    ds.attrs['PI_Affiliation'] = ict.PIAffiliation
    ds.attrs['Platform'] = ict.dataSourceDescription
    ds.attrs['Mission'] = ict.missionName
    ds.attrs['DateOfCollection'] = ict.dateOfCollection
    ds.attrs['DateOfRevision'] = ict.dateOfRevision
    ds.attrs['Data_Interval'] = ict.dataIntervalCode
    ds.attrs['Independent_Var'] = str(ict.independentVariable)
    ds.attrs['Dependent_Var_Num'] = len(ict.dependentVariables)
    ds.attrs['PI_Contact'] = ict.normalComments[0].split('\n')[0].split(':')[-1]
    ds.attrs['Platform'] = ict.normalComments[1].split(':')[-1]
    ds.attrs['Location'] = ict.normalComments[2].split(':')[-1]
    ds.attrs['Associated_Data'] = ict.normalComments[3].split(':')[-1]
    ds.attrs['Instrument_Info'] = ict.normalComments[4].split(':')[-1]
    ds.attrs['Data_Info'] = ict.normalComments[5][11:]
    ds.attrs['DM_Contact'] = ict.normalComments[11].split(':')[-1]
    ds.attrs['Project_Info'] = ict.normalComments[12].split(':')[-1]
    ds.attrs['Stipulations'] = ict.normalComments[13].split(':')[-1]
    ds.attrs['Comments'] = ict.normalComments[14].split(':')[-1]
    ds.attrs['Revision'] = ict.normalComments[15].split(':')[-1]
    ds.attrs['Revision_Comments'] = ict.normalComments[15 + 1].split(':')[-1]

    # Assign Additional ARM meta data to Xarray DatatSet
    ds.attrs['_datastream'] = filename.split('/')[-1].split('_')[0]

    # Return Xarray Dataset
    return ds

def write_icartt(ds, 
                 format=_format,
                 return_None=False,
                 prop_keywords=None,
                 special=None,
                 normal=None,
                 normal_keywords=None):
    """
    Returns ICARTT formatted file with stored data and metadata from a 
    single datastream. Has some procedures to allow for user-defined changes to
    ICARTT metadata.

    Parameters
    ----------
    ds : xarray.Dataset
        ACT Xarray dataset.
    format : str
        ICARTT Format to Write: FFI1001 or FFI2110.
    return_None : bool, optional
        Catch IOError exception when file not found and return None.
        Default is False.
    prop_keywords : dict (or None)
        User-defined dictionary of ICARTT file properties to include within
        initial configuration of the ICARTT file. If set to None, information
        will be taken from ACT Xarray DataSet.
        ICARTT file properties include:
            PI_NAME, PI_AFFILIATION, DATA_SOURCE_DESCRIPTION, MISSION_NAME,
            DATE_OF_COLLECTION, DATE_OF_REVISION
    special : List (or None)
        User entered special comments (strings) to be added to the metadata.
        Each indice within entered list will be it's own line. 
        Typically reserved for instructions to end-user for publications,
        or something particularly important on date.
    normal : List (or None)
        User entered normal comments, where each indice within the entered list
        will be it's own line.
        These are separate comment lines from the required keywords.
    normal_keywords : dict (or None)
        User-defined dictionary of ICARTT standard keywords to include within 
        the ICARTT standard normal comment lines. If set to None, information
        will be taken from ACT Xarray DataSet. 
        ICARTT standard Keywords include: 
            PI_CONTACT_INFO, PLATFORM, LOCATION, ASSOCIATED_DATA,
            INSTRUMENT_INFO, DATA_INFO, UNCERTAINTY, ULOD_FLAG, ULOD_VALUE,
            LLOD_FLAG, LLOD_VALUE, DM_CONTACT_INFO, PROJECT_INFO,
            STIPULATIONS_ON_USE, OTHER_COMMENTS, REVISION

    Returns
    -------
    ict : icartt.Dataset
    """
    if not _ICARTT_AVAILABLE:
        raise ImportError(
            "ICARTT is required to use to write ICARTT files but is not installed")

    # Create an exception tuple to use with try statements. Doing it this way
    # so we can add the FileNotFoundError if requested. Can add more error
    # handling in the future.
    except_tuple = (ValueError,)
    if return_None:
        except_tuple = except_tuple + (FileNotFoundError, OSError)

    try:
        # Create an empty ICARTT DataSet
        ict = icartt.Dataset(format=format)

    except except_tuple as exception:
        # If requested return None for File not found error
        if type(exception).__name__ == 'FileNotFoundError':
            return None

        # If requested return None for File not found error
        if (type(exception).__name__ == 'OSError'
                and exception.args[0] == 'no files to open'):
            return None

    # Define ICARTT file properties
    # Principle Investigator Name
    ict.PIName = ds.platform_id + " Instrument Mentors"
    # Principle Investigator Affiliation
    ict.PIAffiliation = "Atmospheric Radiation Measurement (ARM)"
    # Data Source Description
    ict.dataSourceDescription = ds.datastream
    # Mission Name
    ict.missionName = ds.location_description
     # Date of Collection
    ict.dateOfCollection = (int(ds._file_dates[0][0:4]), 
                            int(ds._file_dates[0][4:6]), 
                            int(ds._file_dates[0][6:8])
                            )
    # Date of Revision
    ict.dateOfRevision = datetime.datetime.utcnow().timetuple()[:3]
    
    # Check for user defined entries to change defaults.
    if prop_keywords:
        for key in pro_keywords:
            print('prop_keywords ', key)
            if key == "PI_NAME":
                ict.PIName = pro_keywords[key]
            if key == "PI_AFFILIATION":
                ict.PIAffiliation = pro_keywords[key]
            if key == "DATA_SOURCE_DESCRIPTION":
                ict.dataSourceDescription = pro_keywords[key]
            if key == "MISSION_NAME":
                ict.missionName = pro_keywords[key]
            if key == "DATE_OF_COLLECTION":
                ict.dateOfCollection = pro_keywords[key]
            if key == "DATE_OF_REVISION":
                ict.dateOfRevision = pro_keywords[key]
    
    # never seen this set to anything but zero    
    ict.dataIntervalCode = [0]
    
    # Define the independent variable (has to be time)
    ict.independentVariable = icartt.Variable(
        "Time",
        "Time offset from midnight",
        "Time",
        "Time",
        vartype=icartt.VariableType.IndependentVariable,
        scale=1.0,
        miss=-9999999,
    )

    # Define the dependent variables
    ignore = ["time", "base_time", "time_offset", "time_bounds"]
    for var in ds.variables:
        if var not in ignore and var[0:2] != 'qc':
            if "standard_name" in ds[var].attrs:
                ict.dependentVariables[var] = icartt.Variable(
                    var,
                    ds[var].units,
                    ds[var].standard_name,
                    ds[var].long_name
                )
            else:
                ict.dependentVariables[var] = icartt.Variable(
                    var,
                    ds[var].units,
                    ds[var].long_name,
                    ds[var].long_name
                )

    # Normal and Special Comments
    if special:
        for note in special:
            ict.specialComments.append(special[note])
    else:
        ict.specialComments.append("")

    if normal:
        for note in normal:
            ict.normalComments.freeform.append(normal[note])
    else:
        ict.normalComments.freeform.append(
            "DATA WERE OBTAINED FROM THE ATMOSPHERIC RADIATION MEASUREMENT "
            + "(ARM) USER FACILITY "
        )
        ict.normalComments.freeform.append(
            " A U.S. DEPARTMENT OF ENERGY (DOE) OFFICE OF SCIENCE USER "
            + "FACILITY MANAGED BY THE BIOLOGICAL AND "
            + "ENVIRONMENTAL RESEARCH PROGRAM. "
        )

    # Keywords Are Required to be Set - set to N/A as default
    ict.normalComments.keywords["PLATFORM"].append(ds.platform_id)
    ict.normalComments.keywords["LOCATION"].append(ds.site_id)
    ict.normalComments.keywords["ASSOCIATED_DATA"].append(
        "additional instruments at " + ds.facility_id + " site"
    )
    ict.normalComments.keywords["INSTRUMENT_INFO"].append(ds.datastream)
    ict.normalComments.keywords["DATA_INFO"].append(ds.dod_version)
    ict.normalComments.keywords["UNCERTAINTY"].append("See Mentor Handbook")
    # Upper and Lower Limit of Detection
    ict.normalComments.keywords["ULOD_FLAG"].append("-7777")
    ict.normalComments.keywords["ULOD_VALUE"].append("Not Calculated")
    ict.normalComments.keywords["LLOD_FLAG"].append("-8888")
    ict.normalComments.keywords["LLOD_VALUE"].append("Not Calculated")
    ict.normalComments.keywords["DM_CONTACT_INFO"].append(
        "DM Contact Information Found at arm.gov"
    )
    ict.normalComments.keywords["PROJECT_INFO"].append(ds.location_description)
    ict.normalComments.keywords["STIPULATIONS_ON_USE"].append(
        "The Atmospheric Radiation Measurement (ARM) user facility should "
        + "be acknowledged in publications as the origin of field studies "
        + "or data used in the research according to the following guidelines."
    )
    ict.normalComments.keywords["REVISION"].append(
        "R0: Initial ICARTT Created File."
    )
    # Check for user defined entries to change defaults.
    if normal_keywords:
        for key in normal_keywords:
                ict.normalComments.keywords[key].append(
                    normal_keywords[key]
                )
            
    # Apparently this is needed
    ict.endDefineMode()

    # Add the data, first remove unsupported qc variables and time bounds
    qc_list = ["base_time", "time_offset", "time_bounds"]
    for var in ds.variables:
        if var[0:2] == 'qc':
            qc_list.append(var)
    # Remove the unsupported variables before conversion
    ds = ds.drop_vars(qc_list)
    # Data are needed in csv format for writing, convert to Pandas Dataframe
    # then convert to CSV for output
    df = ds.to_dataframe()
    csv = df.to_csv()

    # iterate over the csv file and add to the ICARTT Dataset
    for line in csv.split('\n')[1:]:
        ict.data.add(np.array(line))
    
    # write
    ict.write(delimiter=',')
    
