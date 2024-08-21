import os

from biosimulators_smoldyn.combine import (
    preprocess_sed_task, exec_sed_task, exec_sed_doc, exec_sedml_docs_in_combine_archive
)
import smoldyn

__all__ = [
    '__version__',
    'get_simulator_version',
    'preprocess_sed_task',
    'exec_sed_task',
    'exec_sed_doc',
    'exec_sedml_docs_in_combine_archive',
]


current_dir = os.path.dirname(__file__)
version_file_path = os.path.join(current_dir, '_VERSION')

with open(version_file_path, 'r') as f:
    __version__ = f.read().strip()


def get_simulator_version():
    """ Get the version of Smoldyn

    Returns:
        :obj:`str`: version
    """
    return smoldyn.__version__
