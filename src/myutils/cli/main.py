from importlib import import_module
import subprocess
from pathlib import Path
import sys
import numpy as np


def output_terminal(cmd, print_output=True, skip_error=False, print_cmd=False,
                    **kwargs):
    """
    Runs a command in a terminal and save the output in a list
    of strings

    Parameters
    ==========
    cmd: str
        bash command to be executed in the terminal.
    print_output: bool. Default=True
        True for printing the output besides of returning it. Default False.
    skip_error: bool. Default=False
        True for continuing running although the command fails. Default False.
    **kwargs:
        additional options for subprocess.Popen
    print_cmd: Default=False
        print the command that was just executed.

    Return
    ======
    (list) [#linesStr] output of the executed command, line by line.
    """
    if print_cmd:
        print(cmd)
    p = subprocess.Popen(cmd,
                         shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         text=True,
                         **kwargs)

    out = ""
    output = ""
    while not (output == '' and p.poll() is not None):
        output = p.stdout.readline()
        if output:
            out += output
            if print_output:
                print(output, end="")
    return_code = p.wait()

    if not skip_error:
        assert not return_code, f"ERROR executing the command \"{cmd}\"  " + \
            "with output_terminal with the next message:\n" + \
            p.stderr.read().strip()

    return out


pymodules = {
    'proline_state': 'myutils.sith.sith_tools',
    'gen_randpep': 'myutils.sith.sith_tools',
    'protonate': 'myutils.sith.protonate',
    'log2xyz': 'myutils.sith.g09_xyz',
    'reduce_structs': 'myutils.sith.from_extreme.info_from_opt',
    'info_from_opt': 'myutils.sith.from_extreme.info_from_opt',
    'change_prolines_dofs': 'myutils.sith.change_dofs',
    'function_doc': 'myutils.miscellaneous',
    'args_and_defaults': 'myutils.miscellaneous',
    'optimized_e': 'myutils.miscellaneous',
    'time_g09': 'myutils.miscellaneous',
    'compute_bde': 'myutils.g_values.postpro',
    'extract_system_info': 'myutils.g_values.postpro',
    'create_table_per_system': 'myutils.g_values.postpro',
    'create_table_per_method': 'myutils.g_values.postpro',
    'extract_gvals': 'myutils.g_values.postpro',
    'nrand_rad': 'myutils.g_values.g_valsetup',
    'ext_xyz_from_npz': 'myutils.g_values.g_valsetup',
    'heavya_idx_from_ref': 'myutils.g_values.g_valsetup',
    'rad_loc': 'myutils.g_values.g_valsetup',
    'iHFC_fromxyz': 'myutils.g_values.g_valsetup',
    'mat2dat': 'myutils.g_values.best_fit',
    'best_fit': 'myutils.g_values.best_fit',
    'spect_w_experiment': 'myutils.g_values.best_fit',
    'create_fit_file': 'myutils.g_values.best_fit',
    'create_amber_data': 'myutils.gromacs.ff_parameters',
    'create_grappa_data': 'myutils.gromacs.ff_parameters',
    'abstracts': 'myutils.conferences.book_creator',
    'methods_in_class': 'myutils.cli.pkg_structure.documentation_tools',
    'nh_neighs': 'myutils.ase_utils.tools',
    'h_link2': 'myutils.ase_utils.tools',
    'F_max_stretch': 'myutils.ase_utils.tools',
    'distance': 'myutils.ase_utils.tools',
    'all_xyz2pdb': 'myutils.ase_utils.tools',
    'conf2pdb': 'myutils.ase_utils.tools',
    'diff_bonds': 'myutils.ase_utils.tools',
    'extract_bonds': 'myutils.ase_utils.tools',
    'change_distance': 'myutils.ase_utils.tools',
    'shake_except': 'myutils.ase_utils.tools',
}

sh_executers = {
    'workflow': './sith/workflow.sh',
    'swap_atoms_in_com': './sith/swap_atoms_in_com.sh',
    'stretching': './sith/stretching.sh',
    'proline_mod': './sith/proline_mod.sh',
    'workflow_from_extreme': './sith/from_extreme/workflow_from_extreme.sh',
    'workflow_from_extreme2': './sith/from_extreme/second_version/workflow_from_extreme2.sh',
    'continuous_path': './sith/from_extreme/second_version/continuous_path.sh',
    'resubmit_failed': './sith/from_extreme/resubmit_failed.sh',
    'rearange_files': './sith/from_extreme/rearange_files.sh',
    'prepare_and_submit': './sith/from_extreme/prepare_and_submit.sh',
    'opt_from_xyzs': './sith/from_extreme/opt_from_xyzs.sh',
    'opt_and_forces': './sith/from_extreme/opt_and_forces.sh',
    'extr_dofs': './sith/from_extreme/extr_dofs.sh',
    'after_optimization': './sith/from_extreme/after_optimization.sh',
    'find_forces': './sith/find_forces.sh',
    'extract_forces': './sith/extract_forces.sh',
    'compute_forces': './sith/compute_forces.sh',
    'clean_ds': './sith/clean_ds.sh',
    'pkges_installer': './pkges_installer.sh',
    'spectrum_from_gval': './g_values/spectrum_from_gval.sh',
    'gval_workflow': './g_values/gval_workflow.sh',
    'gvals_template': './g_values/gvals_template.sh',
    'gval_postpro': './g_values/gval_postpro.sh',
    'gval_dataset': './g_values/gval_dataset.sh',
    'extract_EPRspec': './g_values/extract_EPRspec.sh',
    'extract_bonds_vmd': './g_values/extract_bonds_vmd.sh',
    'pulling_with_ff': './gromacs/pulling_with_ff.sh',
    'pulling': './gromacs/pulling.sh',
    'peptide_pulling': './gromacs/peptide_pulling.sh',
    'extract_distance': './gromacs/extract_distance.sh',
    'equilibrate_pdb': './gromacs/equilibrate_pdb.sh',
    'constraint_run': './gromacs/constraint_run.sh',
    'classical_minimization': './gromacs/classical_minimization.sh',
    'classical_energies': './gromacs/classical_energies.sh',
    'analysis': './gromacs/analysis.sh',
    'python_doc_fixer': './cli/pkg_structure/python_doc_fixer.sh',
    'files_tree': './cli/pkg_structure/files_tree.sh',
    'doc_pythonfile': './cli/pkg_structure/doc_pythonfile.sh',
    'doc_modules': './cli/pkg_structure/doc_modules.sh',
    'check_tests': './cli/pkg_structure/check_tests.sh',
    'check_structure': './cli/pkg_structure/check_structure.sh',
    'bash_style': './cli/pkg_structure/bash_style.sh',
    'bash_basic_structure': './cli/pkg_structure/bash_basic_structure.sh',
    'add_python_doc': './cli/pkg_structure/add_python_doc.sh',
    'generate_main': './cli/generate_main.sh',
    'basics': './basics.sh',
    'single_g09': './bash_scripts/single_g09.sh',
    'find_blocks': './bash_scripts/find_blocks.sh',
    'bash-template': './bash_scripts/bash-template.sh',
}

other_files = {
    'create_mol_png': './g_values/create_mol_png.tcl',
    'pulling_temp': './gromacs/pulling_temp.mdp',
    'nvt': './gromacs/nvt.mdp',
    'npt': './gromacs/npt.mdp',
    'minim': './gromacs/minim.mdp',
    'ions': './gromacs/ions.mdp',
    'constraint': './gromacs/constraint.mdp',
}


def _read_arguments():
    """
    Function that reads args and kwargs.

    Return
    ======
    (tuple) args
    """
    if len(sys.argv) == 1:
        return ((), {})
    argument = '_reader_args'
    values = ''
    args_dict = {}
    for entry in np.array(sys.argv[1:]):
        if '--' == entry[:2]:
            if argument == '_reader_args':
                args_dict[argument] = values.split(" ")[2:]
                argument = entry.replace('--', '')
                values = ''
            else:
                args_dict[argument] = eval(values)
                argument = entry.replace('--', '')
                values = ''
        else:
            values += ' ' + entry
    if argument == '_reader_args':
        args_dict[argument] = values.split(" ")[2:]
    else:
        args_dict[argument] = eval(values)
    args = args_dict['_reader_args']
    del args_dict['_reader_args']
    
    if '--' == sys.argv[1][:2]:
        args = ()

    return (args, args_dict)


def main():
    """
    This function run each time myutils is called from the terminal.

    Return
    ======
    (None)
    """
    # Help menu of this code
    if sys.argv[1] == '-h' or sys.argv[1] == '--help' or sys.argv[1] == 'help':
        functions = list(pymodules.keys()) + list(sh_executers.keys())
        functions.append('tests')
        functions.sort()

        print("\n"
              "'myutils' contains the next set of tools:\n\n")
        for function in functions:
            print("    -   " + function)

        print("\nTo use one of the functions from the terminal, use:\n\n"
              "  'myutils <function> <arg1> <arg2>...'\n\n"
              "For detailed information of any function, use \"-h\" as first"
              " argument (<arg1>).")

    elif sys.argv[1] == 'tests':
        testdir = Path(__file__).parent
        cmd = f"cd {str(testdir)}/../tests ; pytest -v --color=yes" + \
            ' '.join(sys.argv[2:])
        output_terminal(cmd)

    # python module from terminal
    elif sys.argv[1] in pymodules.keys():
        module = import_module(pymodules[sys.argv[1]])
        method = getattr(module, sys.argv[1])

        if '-h' in sys.argv:
            print(method.__doc__)
        elif '-path' in sys.argv:
            print(pymodules[sys.argv[1]])
        else:
            arguments = _read_arguments()
            output = method(*arguments[0], **arguments[1])
            if output is not None:
                print(output)

    # bash codes
    elif sys.argv[1] in sh_executers.keys():
        if '-path' in sys.argv[2:]:
            path = (str(Path(__file__).parent)[:-3]
                    + sh_executers[sys.argv[1]][2:])
            print(path)
        else:
            command = str(Path(__file__).parent)[:-3] + \
                sh_executers[sys.argv[1]][2:] + ' ' + \
                ' '.join(sys.argv[2:])

            output_terminal(command, print_output=True)

    # other files
    elif sys.argv[1] in other_files.keys():
        print(str(Path(__file__).parent)[:-3] + other_files[sys.argv[1]][2:])

    # own path
    elif sys.argv[1] == 'path':
        print(str(Path(__file__).parent)[:-3])
    
    # open documentation
    elif sys.argv[1] == 'doc':
        command = "xdg-open " + str(Path(__file__).parent)[:-3] + \
            "../../doc/_build/html/index.html"
        output_terminal(command)

    # Not recognized keyword
    else:
        print(f"ERROR: keyword {sys.argv[1]} not recognized as part of"
              " myutils. Use 'myutils -h' to see the options you can use.")


if __name__ == "__main__":
    main()
