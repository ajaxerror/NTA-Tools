##### Version
# 1.0
##### Created
# 5/25/2021
# Matt Boyce
##### Updated
# None

##### How to use script:
# This script is designed to be run through CLI to convert
# mzXML of MS1 data into a long-form .csv. The .csv will be
# stored in the same location as the mzXML file using the
# same file name, but appended with .csv in lieu of
# .mzXML.
#
# To run the script, call the file from command line with the
# following notation: -file 'file-path', where file-path is
# the absolute location of the mzXML file.


##### Import Dependencies
import pyopenms as pyms
import sys
import csv

def parse_mzXML(file_path):
    """
    Takes in the file path of a mzXML file, reads it and converts it to a .csv
    in the same location and name

    :param file_path: str
    :return: None
    """
    exp = pyms.MSExperiment()
    pyms.MzXMLFile().load(file_path, exp)
    mzXML_to_csv(exp, file_path)

def mzXML_to_csv(mzXML, file_path):
    """
    Reads each spectrum of a MSExperiment object, which corresponds to all
    mz, intensity pairs at a single retention time. Spectrum data are
    stored to a local .csv in a long-form format.

    :param mzXML:pyopenms.MSExperiment object
    :param file_path: str
    :return: None
    """
    with open(file_path.replace('.mzXML', '.csv'), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['rt', 'mz', 'intensity'])
        for spec in mzXML:
            mz, intensity = spec.get_peaks()
            rt = [spec.getRT()] * len(intensity)
            writer.writerows([list(row) for row in zip(rt, mz, intensity)])

if __name__ == '__main__':
    for idx, arg in enumerate(sys.argv):
        if arg == '-file':
            exp = parse_mzXML(sys.argv[idx+1])