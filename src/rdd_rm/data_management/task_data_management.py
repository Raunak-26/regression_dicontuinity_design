import pytask
from rdd_rm.config import BLD
from rdd_rm.config import SRC



@pytask.mark.r(script=SRC / "data_management" / "clean_data.r", serializer="yaml")
@pytask.mark.depends_on(
    {
        "data_info": SRC / "data_management" / "data_info.yaml",
        "data": SRC / "data" / "data.csv",
    }
)
@pytask.mark.produces(BLD / "r" / "data" / "data_clean.csv")
def task_clean_data_r():
    pass

