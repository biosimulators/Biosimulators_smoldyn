import os
from tempfile import mkdtemp

from biosimulators_utils.config import Config
from biosimulators_utils.log.data_model import StandardOutputErrorCapturerLevel

from biosimulators_smoldyn.combine import exec_sedml_docs_in_combine_archive, exec_sed_doc


TEST_CONFIG = Config(
    VERBOSE=True,
    LOG=True,
    COLLECT_SED_DOCUMENT_RESULTS=True,
    BUNDLE_OUTPUTS=True
)
TEST_FIXTURES_DIR = './fixtures'
TEST_DEST_DIR = './artifacts'


def show_artifacts(output_dirname: str):
    dest_files = os.listdir(output_dirname)
    print(f'Dest files in test: {dest_files}')


def test_exec_sed_doc():
    dirname = 'lotka-volterra'
    output_dirname = os.path.join(TEST_DEST_DIR, dirname)
    archive_dir = os.path.join(TEST_FIXTURES_DIR, dirname)
    sedml_fp = os.path.join(archive_dir, 'simulation.sedml')
    exec_sed_doc(
        doc=sedml_fp,
        working_dir=archive_dir,
        base_out_path=output_dirname,
        config=TEST_CONFIG,
        log_level=StandardOutputErrorCapturerLevel.c
    )

    return show_artifacts(output_dirname)


def test_exec_sedml_docs_in_combine_archive():
    dirname = 'lotka-volterra-omex'
    archive_fp = './fixtures/lotka-volterra.omex'
    output_dirname = os.path.join(TEST_DEST_DIR, dirname)
    exec_sedml_docs_in_combine_archive(
        archive_filename=archive_fp,
        out_dir=output_dirname,
        config=TEST_CONFIG
    )

    return show_artifacts(output_dirname)


if __name__ == '__main__':
    # test_exec_sed_doc()
    test_exec_sedml_docs_in_combine_archive()
