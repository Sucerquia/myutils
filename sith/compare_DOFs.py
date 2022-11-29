from SITH.Utilities import Extractor
from pathlib import Path
import sys


def extract_DOFs(fchk_name):
    workingPath = Path.cwd()
    referencePath = workingPath / fchk_name

    with referencePath.open() as rFile:
        rData = rFile.readlines()
    rExtractor = Extractor(referencePath, rData)
    rExtractor._extract()
    # Create Geometry objects from reference and deformed data
    ref = rExtractor.getGeometry()

    return ref.dimIndices

if __name__ == '__main__':
    fchk1 = sys.argv[1]
    fchk2 = sys.argv[2]
    constrain1 = int(sys.argv[3])
    constrain2 = int(sys.argv[4])
    print(f"\n\
    ++++++++ Compare DOFs: VERBOSE - {fchk1} vs {fchk2} removing ({constrain1}, \
    {constrain2}) ++++++++")
    dofs1 = extract_DOFs(fchk1)
    dofs2 = extract_DOFs(fchk2)
    try:
        dofs1.remove((constrain1, constrain2))
    except ValueError:
        pass
    try:
        dofs2.remove((constrain1, constrain2))
    except ValueError:
        pass
    assert len(dofs1) == len(dofs2), "there are DOF in one that is not in the other"

    for i in range(len(dofs1)):
        if dofs1[i] != dofs2[i]:
            raise Exception("there are discrepancies DOFs in the files 1 and 2. ")



