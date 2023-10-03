from SOPRANO.utils.app_utils import get_annotated_input_options
from SOPRANO.utils.path_utils import Directories


def _check_detected(dir_method, options_method, ext):
    tmp_name = f"pytest.{ext}"
    other_name = "pytest.other"

    tmp_path = dir_method(tmp_name)
    other_path = dir_method(other_name)
    tmp_path.touch(exist_ok=False)
    other_path.touch(exist_ok=False)

    options = options_method()

    assert tmp_name in options.keys()
    assert other_name not in options.keys()

    tmp_path.unlink(missing_ok=False)
    other_path.unlink(missing_ok=False)


def test_get_annotated_input_options():
    _check_detected(
        Directories.app_annotated_inputs, get_annotated_input_options, "anno"
    )
