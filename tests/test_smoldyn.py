import os
from tempfile import mkdtemp

from biosimulators_utils.config import Config

from biosimulators_smoldyn.combine import exec_sedml_docs_in_combine_archive


TEST_CONFIG = Config(VERBOSE=False, LOG=False)


def test_exec_sedml_docs_in_combine_archive():
    archive_fp = './fixtures/Min1.omex'
    dest = mkdtemp()

    results, log = exec_sedml_docs_in_combine_archive(
        archive_filename=archive_fp,
        out_dir=dest,
        config=TEST_CONFIG
    )

    print(f'Dest files: {os.listdir(dest)}')


test_exec_sedml_docs_in_combine_archive()
