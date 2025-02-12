import contextlib, tempfile
from pathlib import Path
from functools import partial
from subprocess import check_output as _call
run_shell_cmd = partial(_call, shell=True)


@contextlib.contextmanager
def temprary_filename(mode='w+b', suffix=None):
    tmp_name = None
    try:
        tmp_file = tempfile.NamedTemporaryFile(mode=mode, suffix=suffix, delete=False)
        tmp_name = tmp_file.name
        tmp_file.close()
        yield tmp_name
    finally:
        Path(tmp_name).unlink(missing_ok=True)