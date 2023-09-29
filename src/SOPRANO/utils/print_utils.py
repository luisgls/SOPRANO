from datetime import datetime


def time_output():
    now = datetime.now()
    return now.strftime("%d/%m/%Y %H:%M:%S")


def task_output(msg):
    print(f"[{time_output()}] {msg}")


def title_output():
    print(
        """
███████  ██████  ██████  ██████   █████  ███    ██  ██████  
██      ██    ██ ██   ██ ██   ██ ██   ██ ████   ██ ██    ██ 
███████ ██    ██ ██████  ██████  ███████ ██ ██  ██ ██    ██ 
     ██ ██    ██ ██      ██   ██ ██   ██ ██  ██ ██ ██    ██ 
███████  ██████  ██      ██   ██ ██   ██ ██   ████  ██████  
"""
    )


def line_output(n=60):
    print("-" * n)


def startup_output(**kwargs):
    title_output()
    line_output()
    print("Selection On PRotein ANnotated regiOns")
    line_output()

    if kwargs:
        # Parameters used in pipeline
        for k, v in kwargs.items():
            print("-> {0:.<30}".format(k) + f"{v}")

        line_output()
