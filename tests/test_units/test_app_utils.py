from SOPRANO.utils.app_utils import get_annotated_input_options
from SOPRANO.utils.path_utils import Directories


def test_get_annotated_input_options():
    tmp_name = "pytest.anno"
    bed_name = "pytest.bed"

    tmp_path = Directories.app_annotated_inputs(tmp_name)
    bed_path = Directories.app_annotated_inputs(bed_name)
    tmp_path.touch(exist_ok=False)
    bed_path.touch(exist_ok=False)

    options = get_annotated_input_options()

    assert tmp_name in options.keys()
    assert bed_name not in options.keys()

    tmp_path.unlink(missing_ok=False)
    bed_path.unlink(missing_ok=False)
