# Code taken from:
# https://github.com/chunxiangzheng/gaussian_log_file_converter/
import re


code = {"1": "H", "2": "He", "3": "Li", "4": "Be", "5": "B",
        "6": "C", "7": "N", "8": "O",  "9": "F", "10": "Ne",
        "11": "Na", "12": "Mg", "13": "Al", "14": "Si", "15": "P",
        "16": "S", "17": "Cl", "18": "Ar", "19": "K", "20": "Ca",
        "21": "Sc", "22": "Ti", "23": "V", "24": "Cr", "25": "Mn",
        "26": "Fe", "27": "Co", "28": "Ni", "29": "Cu", "30": "Zn",
        "31": "Ga", "32": "Ge", "33": "As", "34": "Se", "35": "Br",
        "36": "Kr", "37": "Rb", "38": "Sr", "39": "Y", "40": "Zr",
        "41": "Nb", "42": "Mo", "43": "Tc", "44": "Ru", "45": "Rh",
        "46": "Pd", "47": "Ag", "48": "Cd", "49": "In", "50": "Sn",
        "51": "Sb", "52": "Te", "53": "I", "54": "Xe", "55": "Cs",
        "56": "Ba", "57": "La", "58": "Ce", "59": "Pr", "60": "Nd",
        "61": "Pm", "62": "Sm", "63": "Eu", "64": "Gd", "65": "Tb",
        "66": "Dy", "67": "Ho", "68": "Er", "69": "Tm", "70": "Yb",
        "71": "Lu", "72": "Hf", "73": "Ta", "74": "W", "75": "Re",
        "76": "Os", "77": "Ir", "78": "Pt", "79": "Au", "80": "Hg",
        "81": "Tl", "82": "Pb", "83": "Bi", "84": "Po", "85": "At",
        "86": "Rn", "87": "Fr", "88": "Ra", "89": "Ac", "90": "Th",
        "91": "Pa", "92": "U", "93": "Np", "94": "Pu", "95": "Am",
        "96": "Cm", "97": "Bk", "98": "Cf", "99": "Es", "100": "Fm",
        "101": "Md", "102": "No", "103": "Lr", "104": "Rf", "105": "Db",
        "106": "Sg", "107": "Bh", "108": "Hs", "109": "Mt", "110": "Ds",
        "111": "Rg", "112": "Uub", "113": "Uut", "114": "Uuq", "115": "Uup",
        "116": "Uuh", "117": "Uus", "118": "Uuo"}


def getEnergy(structure):
    for line in structure.split("\n"):
        if line.startswith(" SCF Done:"):
            arr = line.split("=")
            return float(re.split(" +", arr[1].strip())[0])
    return 1000.0


def findInList(dataList, target):
    for i in range(0, len(dataList)):
        if dataList[i].find(target) != -1:
            return i
    return -1


def getCoordinates(dataList):
    start = findInList(dataList, "Standard orientation")
    dataList = dataList[start + 5:]
    dataList = dataList[: findInList(dataList, "-----")]
    return dataList


# add2executable
def log2xyz(finput, foutput=None):
    """
    Extract the configuration of minumum energy in optimization process in a
    xyz file from a .log gaussian file.

    Parameters
    ==========
    finput: str
        path to the log file.
    foutput: str (optional)
        name of the output file without extension.

    Note: if foutput is not given, the name output will be the same than the
    input but with xyz extension.
    """
    infoBlock = ""
    optimized = False
    optimized_structure = ""

    with open(finput, "r") as fin:
        isStructure = True
        isInfo = True
        structures = []
        currentStructure = ""
        for line in fin:
            if line.startswith(" GradGrad"):
                if isInfo:
                    isInfo = False
                if currentStructure != "":
                    structures.append((getEnergy(currentStructure),
                                       currentStructure))
                    currentStructure = ""
                isStructure = not isStructure
            elif isInfo:
                infoBlock += line
            elif isStructure:
                currentStructure += line
            else:
                if line.find("Optimized") != -1:
                    optimized = True

        if optimized:
            optimized_structure = currentStructure
        else:
            if currentStructure != "":
                structures.append((getEnergy(currentStructure),
                                   currentStructure))
            structures = sorted(structures, key=lambda item: item[0])
            optimized_structure = structures[0][1]

    if foutput:
        prefix = foutput
    else:
        prefix = finput.strip(".log")
    foutput = prefix + ".xyz"

    with open(foutput, "w") as fout:
        dataList = optimized_structure.split("\n")
        atoms = getCoordinates(dataList)
        fout.write(str(len(atoms)) + "\n\n")
        for atom in atoms:
            arr = atom.split()
            symbol = code.get(arr[1], 'X')
            fout.write("  %s %16.7f %16.7f %16.7f\n" % (symbol, float(arr[3]),
                       float(arr[4]), float(arr[5])))
