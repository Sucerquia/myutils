from SITH.SITH import SITH
from myutils.sith.sith import Sith
import sys


# DEPRECTED
def save_extradofs(drelaxed, dstretched):
    """"
    This code finds the degrees of freedom in the
    optmimization calculation that don't belong to the
    forces code """
    sith = SITH(drelaxed, dstretched)

    sith.extractData()
    sith.energyAnalysis()
    print(sith._deformed[0].dims)

    a = Sith(sith)

    print(sith._deformed[0].dims)

    with open('extra_DOFs.dat', 'w') as extra_DOF:
        for i in a._deformed[0].dimIndices:
            if (i not in sith._deformed[0].dimIndices) \
                and (i[::-1] not in sith._deformed[0].dimIndices):
                for j in i:
                    extra_DOF.write(f" {j}")
                extra_DOF.write("\n")


