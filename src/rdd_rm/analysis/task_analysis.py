import pytask
from rdd_rm.config import BLD
from rdd_rm.config import GROUPS
from rdd_rm.config import SRC



@pytask.mark.r(script=SRC / "analysis" / "fit_model.r", serializer="yaml")
@pytask.mark.depends_on(
    {
        "data": BLD / "r" / "data" / "data_clean.csv",
        "data_info": SRC / "data_management" / "data_info.yaml",
    }
)
@pytask.mark.produces(BLD / "r" / "models" / "model.rds")
def task_fit_model_r():
    pass


for group in GROUPS:

    kwargs = {
        "group": group,
        "produces": BLD / "r" / "predictions" / f"{group}.csv",
    }

    @pytask.mark.depends_on(
        {
            "data": BLD / "r" / "data" / "data_clean.csv",
            "model": BLD / "r" / "models" / "model.rds",
        }
    )
    @pytask.mark.task(id=group, kwargs=kwargs)
    @pytask.mark.r(script=SRC /  "analysis" / "predict.r", serializer="yaml")
    def task_predict_r():
        pass

