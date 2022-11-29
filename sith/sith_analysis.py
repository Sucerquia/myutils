import sys
sys.path.insert(0, '/hits/basement/mbm/sucerquia/utils')
from miscellaneous import all_hydrogen_atoms
from SITH.SithWriter import writeSummary
from SITH.SITH import SITH
import glob


def create_sith_analysis(drelaxed, dstreched):
    print("\
    JEDI analysis will be applied using the SITH package using the next relaxed\
    file and stretched directory:")
    print("   - ", drelaxed)
    print("   - ", dstreched)
    sith = SITH(drelaxed, dstreched)
    config = glob.glob('./optimization/*.xyz')
    Hatoms_Ala = all_hydrogen_atoms(config[0])
    sith.setKillAtoms(Hatoms_Ala)
    sith.extractData()
    sith.energyAnalysis()

    return writeSummary(sith, includeXYZ=True)

if __name__ == '__main__':
    if '-h' in sys.argv:
        print(create_sith_analysis.__doc__)
    else:
        create_sith_analysis(*sys.argv[1:])