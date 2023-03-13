from SITH.SITH import SITH

def save_dofs(fchk_file):
    sith = SITH(fchk_file, fchk_file)

    sith.extractData()
    sith.energyAnalysis()

    with open("all_dofs.dat", "w") as dofs_file:
        for dof in sith._reference.dimIndices:
            for i in dof:
                dofs_file.write(f"{i} ")
            dofs_file.write("\n")
    
