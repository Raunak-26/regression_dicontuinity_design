import pytask
from rdd_rm.config import BLD
from rdd_rm.config import GROUPS
from rdd_rm.config import SRC



for group in GROUPS:

    kwargs = {
        "group": group,
        "depends_on": {"predictions": BLD / "r" / "predictions" / f"{group}.csv"},
        "produces": BLD / "r" / "figures" / f"smoking_by_{group}.png",
    }

    @pytask.mark.depends_on(
        {
            "data_info": SRC / "data_management" / "data_info.yaml",
            "data": BLD / "r" / "data" / "data_clean.csv",
        }
    )
    @pytask.mark.task(id=group, kwargs=kwargs)
    @pytask.mark.r(script=SRC / "final" / "plot_regression.r", serializer="yaml")
    def task_plot_regression_r():
        pass


@pytask.mark.depends_on({"model": BLD / "r" / "models" / "model.rds", "SRC": SRC})
@pytask.mark.produces(BLD / "r" / "tables" / "estimation_results.tex")
@pytask.mark.r(script=SRC / "final" / "create_estimation_table.r", serializer="yaml")
def task_create_estimation_table_r():
    pass

