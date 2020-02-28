from os import system
from pynestml.frontend.pynestml_frontend import to_nest, install_nest


to_nest(input_path="modules",
        target_path="../tmp/module",
        logging_level="INFO",
        module_name="nestml_module",
        suffix="_nestml")

install_nest("../tmp/module",
             "/home/abolfazl/prog/install-dir/nest-simulator-2.18.0_build")
