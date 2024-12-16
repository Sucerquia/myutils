from importlib import import_module
from myutils.miscellaneous import output_terminal
from pathlib import Path
import sys


pymodules = {
    'F_max_stretch': 'myutils.ase_utils.tools',
    'distance': 'myutils.ase_utils.tools',
    'all_xyz2pdb': 'myutils.ase_utils.tools',
    'conf2pdb': 'myutils.ase_utils.tools',
    'diff_bonds': 'myutils.ase_utils.tools',
    'extract_bonds': 'myutils.ase_utils.tools',
    'change_distance': 'myutils.ase_utils.tools',
    'create_table_per_system': 'myutils.g_values.analysis.postpro',
    'create_table_per_method': 'myutils.g_values.analysis.postpro',
    'extract_gvals': 'myutils.g_values.analysis.postpro',
    'create_amber_data': 'myutils.gromacs.ff_parameters',
    'create_grappa_data': 'myutils.gromacs.ff_parameters',
    'methods_in_class': 'myutils.cli.pkg_structure.documentation_tools',
    'protonate': 'myutils.sith.protonate',
    'proline_state': 'myutils.sith.sith_tools',
    'gen_randpep': 'myutils.sith.sith_tools',
    'log2xyz': 'myutils.sith.g09_xyz',
    'reduce_structs': 'myutils.sith.from_extreme.info_from_opt',
    'info_from_opt': 'myutils.sith.from_extreme.info_from_opt',
    'function_doc': 'myutils.miscellaneous',
    'args_and_defaults': 'myutils.miscellaneous',
    'optimized_e': 'myutils.miscellaneous',
    'time_g09': 'myutils.miscellaneous',
}

sh_executers = {
    'single_g09': './bash_scripts/single_g09.sh',
    'bash-template': './bash_scripts/bash-template.sh',
    'find_blocks': './bash_scripts/find_blocks.sh',
    'create_g09_BDEs': './g_values/basic_scripts/create_g09_BDEs.sh',
    'clean_directories': './g_values/basic_scripts/clean_directories.sh',
    'gval_workflow': './g_values/basic_scripts/gval_workflow.sh',
    'submit': './g_values/basic_scripts/submit.sh',
    'constraint_run': './gromacs/constraint_run.sh',
    'classical_minimization': './gromacs/classical_minimization.sh',
    'pulling_with_ff': './gromacs/pulling_with_ff.sh',
    'equilibrate_pdb': './gromacs/equilibrate_pdb.sh',
    'classical_energies': './gromacs/classical_energies.sh',
    'analysis': './gromacs/analysis.sh',
    'pulling': './gromacs/pulling.sh',
    'extract_distance': './gromacs/extract_distance.sh',
    'peptide_pulling': './gromacs/peptide_pulling.sh',
    'generate_main': './cli/generate_main.sh',
    'bash_style': './cli/pkg_structure/bash_style.sh',
    'python_doc_fixer': './cli/pkg_structure/python_doc_fixer.sh',
    'doc_pythonfile': './cli/pkg_structure/doc_pythonfile.sh',
    'files_tree': './cli/pkg_structure/files_tree.sh',
    'check_structure': './cli/pkg_structure/check_structure.sh',
    'bash_basic_structure': './cli/pkg_structure/bash_basic_structure.sh',
    'add_python_doc': './cli/pkg_structure/add_python_doc.sh',
    'check_tests': './cli/pkg_structure/check_tests.sh',
    'doc_modules': './cli/pkg_structure/doc_modules.sh',
    'extract_forces': './sith/extract_forces.sh',
    'proline_mod': './sith/proline_mod.sh',
    'workflow': './sith/workflow.sh',
    'find_forces': './sith/find_forces.sh',
    'clean_ds': './sith/clean_ds.sh',
    'stretching': './sith/stretching.sh',
    'compute_forces': './sith/compute_forces.sh',
    'reduce_ds': './sith/from_extreme/reduce_ds.sh',
    'rearange_files': './sith/from_extreme/rearange_files.sh',
    'forces_from_xyzs': './sith/from_extreme/forces_from_xyzs.sh',
    'submit_forces_after_opt': './sith/from_extreme/submit_forces_after_opt.sh',
    'resubmit_failed': './sith/from_extreme/resubmit_failed.sh',
    'extr_dofs': './sith/from_extreme/extr_dofs.sh',
    'after_optimization': './sith/from_extreme/after_optimization.sh',
    'workflow_from_extreme': './sith/from_extreme/workflow_from_extreme.sh',
    'opt_and_forces': './sith/from_extreme/opt_and_forces.sh',
    'prepare_and_submit': './sith/from_extreme/prepare_and_submit.sh',
    'basics': './basics.sh',
    'pkges_installer': './pkges_installer.sh',
}

other_files = {
    'pulling_temp': './gromacs/pulling_temp.mdp',
    'nvt': './gromacs/nvt.mdp',
    'minim': './gromacs/minim.mdp',
    'ions': './gromacs/ions.mdp',
    'npt': './gromacs/npt.mdp',
    'constraint': './gromacs/constraint.mdp',
}


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
              "This package contains a set of tools you can use for different"
              "functions. \n To use any function from the terminal, use"
              "    myutils <function> <arg1> <arg2> ... "
              "where <function> can be one of the next options:")
        for function in functions:
            print("    -   " + function)

        print("\nFor detailed information of any function, use \"-h\" as first"
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

        else:
            output = method(*sys.argv[2:])
            if output is not None:
                print(output)

    # bash codes
    elif sys.argv[1] in sh_executers.keys():
        if '-path' in sys.argv[2:]:
            path = str(Path(__file__).parent)[:-3] + \
                    sh_executers[sys.argv[1]][2:]
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
