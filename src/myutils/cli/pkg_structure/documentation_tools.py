import inspect
from importlib import import_module


# add2executable
def methods_in_class(module, class_name):
    """
    Shows the methods that are in a class.

    Parameters
    ==========
    module: str
        name of he module.
    class_name: str
        name of the class.

    Return
    ======
    (list) all modules in the class.
    """
    module = import_module(module)

    imported_class = getattr(module, class_name)

    # Get all methods of the class
    methods = [member[0] for member
               in inspect.getmembers(imported_class,
                                     predicate=inspect.isfunction)
               if member[0] != '__init__']
    for method in methods:
        print(method)
    return methods
