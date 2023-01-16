from SITH.Utilities import Extractor
from pathlib import Path


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


def compare(fchk1, fchk2, constrain1, constrain2):
    """
    Compare the degrees of freedom of two configurations defined in fchk files.
    """
    print(f"\n\
    ++++++++ Compare DOFs: VERBOSE - {fchk1} vs {fchk2} removing ({constrain1}, \
    {constrain2}) ++++++++")
    dofs1 = extract_DOFs(fchk1)
    dofs2 = extract_DOFs(fchk2)
    constrain1 = int(constrain1)
    constrain2 = int(constrain2)
    try:
        dofs1.remove((constrain1, constrain2))
    except ValueError:
        pass
    try:
        dofs2.remove((constrain1, constrain2))

    except ValueError:
        pass
    #assert len(dofs1) == len(dofs2), "there are DOF in one that is not in the other"

    for i in range(len(dofs1)):
        if dofs1[i] != dofs2[i]:
            print(dofs1[i], dofs2[i])
            raise Exception("there are discrepancies DOFs in the files 1 and 2. ")
    print(f"\n\
    ++++++++ Compare DOFs: VERBOSE - comparison successful ++++++++")
    return 0
